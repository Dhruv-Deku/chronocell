"""
Variant impact engine, version 2 (Phase B1): explicit joins, two chromosomes in one model, complex events
and copy number. It builds on chronocell.perturb (unchanged; the 02 · 4D dynamics → 04 Variant impact panel
still uses it as before) and adds:

1. Explicit joins. A rearrangement is described by its breakend pairs (VCF breakend ALT notation, BEDPE
   strands, or typed by hand), not only by deleted intervals: the Gate 4b post hoc showed that deleted
   intervals alone do not say how the kept pieces are joined. Each loaded window ("source") is cut at the
   breakpoints; each join links one side of one cut to one side of another:
       side "L" = the piece that ENDS at the cut (sequence to its left),
       side "R" = the piece that STARTS at the cut.
   VCF:   t[p[ -> (L, R)   t]p] -> (L, L)   ]p]t -> (R, L)   [p[t -> (R, R)
   BEDPE: strand "+" -> L,  strand "-" -> R  (LUMPY / SVTyper convention)
   The derivative chromosomes are the walks through these adjacencies from each window end; pieces reached
   from no window end (acentric fragments) are reported as lost; a breakend used twice is rejected.
2. Two chromosomes in one model. Each source has its own population model (pair-variance matrix in nm²);
   different sources are independent (no data say otherwise); a source without contact data can be a
   homogeneous chain. The derived chain is built exactly as perturb.derive does (pieces keep their joint law,
   reused pieces become independent copies, each junction is a new bond), generalised to several sources.
3. Copy number. The cell carries `reference_copies` normal homologs of each source and `copies` of each
   derivative; Hi-C reads mapped to the reference sum contacts over all of them, so the predicted change of a
   reference pair is log2( after / before ) with before = ploidy x p_reference. Copy number can come from the
   VCF (GT 0/1 -> one derivative + one normal homolog; 1/1 -> two derivatives; INFO CN), from a CNV BED file,
   or be estimated from Hi-C coverage (`estimate_copy_number`, labelled "not validated").

Outputs: contact changes in reference coordinates (within each source, and new contacts between sources);
insulation boundaries lost / gained and domains fused (reference coordinates) and domains that span a junction
on the derivative ("neo-domains"); genes affected (copy number, broken, inverted, moved next to a junction,
contacts changed); enhancer–promoter pairs gained / lost / new across a junction; a confidence from refits on
resampled counts; and a transparent ranking score. Every output is a prediction of a mechanism simulator: it is
labelled "not validated" in the app until Gate 4d (validation/RESULTS.md) passes.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np

from . import domains as DOM, ensemble as ENS

MIN_P = 0.01            # a contact counts if its probability per copy is at least this before or after
MIN_LOG2 = 1.0          # ... and changes at least two-fold
JUNCTION_NEIGHBOURHOOD_BP = 500_000     # a gene this close to a new junction is "moved next to a new partner"
BOUNDARY_TOLERANCE = 2  # beads: a boundary within this distance after the variant is the same boundary


# ======================================================================================
# Description of the rearranged genome
# ======================================================================================
@dataclass(frozen=True)
class Segment:
    source: str
    a: int                  # first bead (inclusive), source bead coordinates
    b: int                  # end bead (exclusive)
    reverse: bool = False

    @property
    def beads(self) -> np.ndarray:
        r = np.arange(self.a, self.b)
        return r[::-1] if self.reverse else r

    def label(self) -> str:
        return f"{self.source}:{self.a}-{self.b}{'(-)' if self.reverse else ''}"


@dataclass(frozen=True)
class Join:
    source1: str
    cut1: int
    side1: str              # "L" | "R"
    source2: str
    cut2: int
    side2: str
    label: str = ""


@dataclass
class Derivative:
    segments: list[Segment]
    copies: float = 1.0
    name: str = ""
    junctions: list[int] = field(default_factory=list)    # derivative bead index where each novel join starts

    @property
    def n(self) -> int:
        return int(sum(s.b - s.a for s in self.segments))

    def label(self) -> str:
        return self.name or " + ".join(s.label() for s in self.segments)


@dataclass
class Karyotype:
    derivatives: list[Derivative]
    reference_copies: dict[str, float]          # normal homologs kept, per source
    ploidy: float = 2.0
    lost: list[Segment] = field(default_factory=list)       # acentric pieces, not in any derivative
    notes: list[str] = field(default_factory=list)


_SEG = re.compile(r"^\s*([A-Za-z0-9_.]+)\s*:\s*(\d+)\s*-\s*(\d+)\s*(\(\s*-\s*\)|\(\s*\+\s*\))?\s*$")


def parse_segments(text: str) -> list[Segment]:
    """'A:0-120 + B:200-400(-) + A:150-300' -> segments (bead coordinates, end exclusive; (-) = reversed)."""
    out = []
    for part in text.split("+"):
        if not part.strip():
            continue
        m = _SEG.match(part)
        if not m:
            raise ValueError(f"Cannot read segment '{part.strip()}': expected SOURCE:START-END, optionally (-).")
        a, b = int(m.group(2)), int(m.group(3))
        if b <= a:
            raise ValueError(f"Segment '{part.strip()}' is empty.")
        out.append(Segment(m.group(1), a, b, bool(m.group(4)) and "-" in m.group(4)))
    if not out:
        raise ValueError("No segments given.")
    return out


def breakend_sides(kind: str, orientation: str) -> tuple[str, str] | None:
    """Sides joined by a breakend record: VCF bracket ALT or BEDPE strand pair. None if not stated."""
    o = (orientation or "").strip()
    if "[" in o or "]" in o:
        if o.endswith("["):
            return "L", "R"          # t[p[
        if o.endswith("]"):
            return "L", "L"          # t]p]
        if o.startswith("]"):
            return "R", "L"          # ]p]t
        if o.startswith("["):
            return "R", "R"          # [p[t
        return None
    if len(o) == 2 and set(o) <= {"+", "-"}:
        return ("L" if o[0] == "+" else "R"), ("L" if o[1] == "+" else "R")
    return None


def joins_for_variant(v, windows: dict[str, tuple[str, int, int, int]]) -> list[Join] | str:
    """Joins of one svio.Variant in the loaded windows {source: (chrom, bin0, n, resolution)}, or the reason it
    cannot be placed. DEL / DUP / INV are given as the joins they imply."""
    def place(chrom: str, pos: int):
        for name, (c, bin0, n, res) in windows.items():
            cut = int(round(pos / res)) - bin0
            if c == chrom and 0 < cut < n:
                return name, cut
        return None
    if v.kind == "BND":
        sides = breakend_sides(v.kind, v.orientation)
        if sides is None:
            return "breakend orientation not given (VCF bracket notation or BEDPE strands needed)"
        p1, p2 = place(v.chrom, v.start), place(v.mate_chrom, v.mate_pos or 0)
        if p1 is None or p2 is None:
            return "a breakend is outside the loaded windows"
        return [Join(p1[0], p1[1], sides[0], p2[0], p2[1], sides[1], v.vid or "BND")]
    pa, pb = place(v.chrom, v.start), place(v.chrom, v.end)
    if pa is None or pb is None or pa[0] != pb[0]:
        return "the variant is not inside one loaded window"
    (s, a), (_, b) = pa, pb
    if b - a < 1:
        return "shorter than one bead at this resolution"
    if v.kind == "DEL":
        return [Join(s, a, "L", s, b, "R", v.vid or "DEL")]
    if v.kind == "INV":
        return [Join(s, a, "L", s, b, "L", v.vid or "INV"), Join(s, a, "R", s, b, "R", v.vid or "INV")]
    if v.kind == "DUP":
        return "tandem duplication: use the copy-number path (a breakend would be used twice)"
    return f"unsupported kind {v.kind}"


def walk(sizes: dict[str, int], joins: list[Join]) -> tuple[list[Derivative], list[Segment], list[str]]:
    """Derivative chromosomes implied by the joins: walks from each window end through the adjacencies
    (a join, or the reference adjacency where neither side of a cut is joined). Returns (derivatives that
    contain a novel join, lost acentric pieces, notes)."""
    cuts = {s: sorted({0, n}) for s, n in sizes.items()}
    for j in joins:
        for src, c in ((j.source1, j.cut1), (j.source2, j.cut2)):
            if src not in sizes:
                raise ValueError(f"Unknown source '{src}'.")
            if not 0 < c < sizes[src]:
                raise ValueError(f"Cut {c} is not inside {src} (1..{sizes[src] - 1}).")
            cuts[src] = sorted(set(cuts[src]) | {c})
    segs = {s: list(zip(cs[:-1], cs[1:])) for s, cs in cuts.items()}

    def node(src: str, c: int, side: str) -> tuple[str, int, str]:
        starts = [a for a, _ in segs[src]]
        if side == "L":
            k = starts.index(next(a for a, b in segs[src] if b == c))
            return src, k, "E"
        return src, starts.index(c), "S"

    edge: dict[tuple, tuple] = {}
    novel: set[frozenset] = set()
    for j in joins:
        n1, n2 = node(j.source1, j.cut1, j.side1), node(j.source2, j.cut2, j.side2)
        for nd in (n1, n2):
            if nd in edge:
                raise ValueError(f"A breakend side is joined twice ({j.label or 'join'}): give the derivative "
                                 "chromosomes as explicit segments instead.")
        if n1 == n2:
            raise ValueError("A breakend side is joined to itself.")
        edge[n1], edge[n2] = n2, n1
        novel.add(frozenset((n1, n2)))
    for src, ss in segs.items():
        for k in range(len(ss) - 1):
            e, s = (src, k, "E"), (src, k + 1, "S")
            if e not in edge and s not in edge:
                edge[e], edge[s] = s, e
    seen: set[tuple[str, int]] = set()
    derivs, notes = [], []

    def run(start: tuple) -> Derivative:
        out, junc, pos, nd = [], [], 0, start
        while True:
            src, k, end = nd
            if (src, k) in seen:
                notes.append("A walk met a piece twice (a cycle); it was stopped there.")
                break
            seen.add((src, k))
            a, b = segs[src][k]
            out.append(Segment(src, a, b, reverse=(end == "E")))
            pos += b - a
            exit_node = (src, k, "S" if end == "E" else "E")
            nxt = edge.get(exit_node)
            if nxt is None:
                break
            if frozenset((exit_node, nxt)) in novel:
                junc.append(pos)
            nd = nxt
        return Derivative(_merge(out), 1.0, "", junc)

    starts = [(src, 0, "S") for src in segs] + [(src, len(ss) - 1, "E") for src, ss in segs.items()]
    for start in starts:                    # every window start first, so derivatives read pter -> qter
        if (start[0], start[1]) in seen:
            continue
        d = run(start)
        if d.junctions:
            derivs.append(d)
    lost = [Segment(src, a, b) for src, ss in segs.items() for k, (a, b) in enumerate(ss) if (src, k) not in seen]
    if lost:
        notes.append(f"{len(lost)} piece(s) reached from no window end (acentric) are treated as lost.")
    return derivs, lost, notes


def _merge(segs: list[Segment]) -> list[Segment]:
    """Join consecutive segments that continue each other in the same direction."""
    out: list[Segment] = []
    for s in segs:
        if out:
            p = out[-1]
            if p.source == s.source and p.reverse == s.reverse and (
                    (not s.reverse and p.b == s.a) or (s.reverse and s.b == p.a)):
                out[-1] = Segment(s.source, min(p.a, s.a), max(p.b, s.b), s.reverse)
                continue
        out.append(s)
    return out


def karyotype_from_joins(sizes: dict[str, int], joins: list[Join], zygosity: str = "heterozygous",
                         ploidy: float = 2.0) -> Karyotype:
    """Derivatives implied by joins. Heterozygous: one copy of each derivative and one normal homolog of each
    source; homozygous: two copies of each derivative and no normal homolog of the sources they use."""
    derivs, lost, notes = walk(sizes, joins)
    if not derivs:
        raise ValueError("These joins leave every window unchanged.")
    hom = zygosity == "homozygous"
    for d in derivs:
        d.copies = ploidy if hom else 1.0
    ref = {s: (0.0 if hom else ploidy - 1.0) for s in sizes}
    return Karyotype(derivs, ref, ploidy, lost, notes)


def karyotype_from_segments(sizes: dict[str, int], derivatives: list[list[Segment]], copies: list[float] | None = None,
                            reference_copies: dict[str, float] | None = None, ploidy: float = 2.0) -> Karyotype:
    """Explicit derivative chromosomes (e.g. a chained complex event), with their copy numbers."""
    ds = []
    for k, segs in enumerate(derivatives):
        for s in segs:
            if s.source not in sizes or s.b > sizes[s.source] or s.a < 0:
                raise ValueError(f"Segment {s.label()} is outside its source.")
        junc, pos = [], 0
        for p, s in zip(segs[:-1], segs[1:]):
            pos += p.b - p.a
            cont = p.source == s.source and p.reverse == s.reverse and ((not s.reverse and p.b == s.a) or (s.reverse and s.b == p.a))
            if not cont:
                junc.append(pos)
        ds.append(Derivative(list(segs), float(copies[k]) if copies else 1.0, "", junc))
    ref = reference_copies if reference_copies is not None else {s: ploidy - 1.0 for s in sizes}
    return Karyotype(ds, dict(ref), ploidy)


def karyotype_from_copy_number(n: int, segments_cn: list[tuple[int, int, float]], source: str = "A",
                               ploidy: float = 2.0) -> Karyotype:
    """Copy-number segments [(a, b, CN)] of one source. Losses: the segment is missing from (ploidy - CN)
    homologs (deletions joining its flanks); gains: (CN - ploidy) extra copies, each an unlinked copy (their
    position is unknown, so they add only their own internal contacts)."""
    ks: list[Derivative] = []
    whole = Segment(source, 0, n)
    losses = [(a, b, cn) for a, b, cn in segments_cn if cn < ploidy]
    gains = [(a, b, cn) for a, b, cn in segments_cn if cn > ploidy]
    ref = ploidy
    for a, b, cn in losses:
        k = ploidy - cn
        segs = ([Segment(source, 0, a)] if a > 0 else []) + ([Segment(source, b, n)] if b < n else [])
        if not segs:
            raise ValueError("A loss covering the whole window cannot be modelled.")
        ks.append(Derivative(segs, float(k), f"loss of {a}-{b} (CN {cn:g})", [a] if 0 < a and b < n else []))
        ref -= k
    for a, b, cn in gains:
        ks.append(Derivative([Segment(source, a, b)], float(cn - ploidy), f"extra copies of {a}-{b} (CN {cn:g})", []))
    if ref < 0:
        raise ValueError("Overlapping losses remove more homologs than the ploidy.")
    notes = [] if not gains else ["Extra copies are modelled unlinked: their genomic position is not known."]
    return Karyotype(ks, {source: ref}, ploidy, [], notes + ([] if losses or gains else ["No copy-number change."]))


# ======================================================================================
# Exact derived ensemble over several sources
# ======================================================================================
def _centre(S: np.ndarray) -> np.ndarray:
    S = np.asarray(S, dtype=np.float64)
    return S - S.mean(0, keepdims=True) - S.mean(1, keepdims=True) + S.mean()


def homogeneous_chain(S_like: np.ndarray, m: int) -> np.ndarray:
    """Pair variances of a homogeneous chain of m beads with S_like's mean variance per separation (nm²)."""
    n = len(S_like)
    mean = np.array([np.mean(np.diagonal(S_like, k)) if k < n else np.nan for k in range(m)])
    if m > n:
        slope = (mean[n - 1] - mean[max(1, n // 2)]) / max(1, n - 1 - max(1, n // 2))
        mean[n:] = mean[n - 1] + slope * (np.arange(n, m) - (n - 1))
    out = mean[np.abs(np.subtract.outer(np.arange(m), np.arange(m)))]
    np.fill_diagonal(out, 0.0)
    return out


def derive_segments(S: dict[str, np.ndarray], segments: list[Segment], s_bond: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pair-variance matrix (nm²) of a derivative chromosome made of `segments`, with each derived bead's
    source index (order of S) and source bead. The first use of a bead keeps its source's joint law (all
    first-use pieces of one source share it); a later use of the same beads is an independent copy."""
    names = list(S)
    blocks, offset, ext = [], {}, 0
    for nm in names:
        blocks.append(-0.5 * _centre(S[nm]))
        offset[nm] = ext
        ext += len(S[nm])
    used = {nm: np.zeros(len(S[nm]), bool) for nm in names}
    cols = []
    for seg in segments:
        b = seg.beads
        if used[seg.source][b].any():
            u = np.unique(b)
            blocks.append(-0.5 * _centre(S[seg.source][np.ix_(u, u)]))
            cols.append(ext + (b - u.min()))
            ext += len(u)
        else:
            used[seg.source][b] = True
            cols.append(offset[seg.source] + b)
    sigma = np.zeros((ext, ext))
    pos = 0
    for blk in blocks:
        k = len(blk)
        sigma[pos:pos + k, pos:pos + k] = blk
        pos += k
    total = sum(len(c) for c in cols)
    V = np.zeros((total, ext))
    junction = np.zeros(total)
    src = np.zeros(total, np.int64)
    bead = np.zeros(total, np.int64)
    base = np.zeros(ext)
    r = 0
    for k, (seg, c) in enumerate(zip(segments, cols)):
        m = len(c)
        V[r:r + m] = base[None, :]
        V[np.arange(r, r + m), c] += 1.0
        V[r:r + m, c[0]] -= 1.0
        junction[r:r + m] = k
        src[r:r + m] = names.index(seg.source)
        bead[r:r + m] = seg.beads
        base = base.copy()
        base[c[-1]] += 1.0
        base[c[0]] -= 1.0
        r += m
    G = V @ sigma @ V.T
    g = np.diag(G)
    Sd = np.clip(g[:, None] + g[None, :] - 2.0 * G, 0.0, None) + s_bond * np.abs(junction[:, None] - junction[None, :])
    np.fill_diagonal(Sd, 0.0)
    return Sd, src, bead


def _p(S_nm2: np.ndarray, r_c: float) -> np.ndarray:
    sig = np.sqrt(np.clip(S_nm2, 0.0, None)) / r_c
    p = ENS.contact_probability_from_sigma(np.where(sig > 0, sig, 1e-9))
    np.fill_diagonal(p, np.nan)
    return p


@dataclass
class Source:
    name: str
    S_nm2: np.ndarray                 # pair variances (nm²)
    r_c_nm: float
    chrom: str = ""
    bin0: int = 0
    resolution: int = 0
    signal: np.ndarray | None = None   # activity track per bead (enhancer proxy), optional
    note: str = ""

    @property
    def n(self) -> int:
        return len(self.S_nm2)


def source_from_result(name: str, res, window: tuple[int, int] | None = None, **meta) -> Source:
    """A Source from a fitted population (v3.3 or v4), optionally a window of it (bead coordinates)."""
    from .perturb import pair_variance_from_result
    S, r_c = pair_variance_from_result(res, window)
    return Source(name, np.asarray(S, float) * r_c ** 2, float(r_c), **meta)


@dataclass
class ImpactV2:
    sources: list[Source]
    karyotype: Karyotype
    p_before: np.ndarray              # (N, N) per-copy contact probability, reference coordinates (all sources)
    p_after: np.ndarray               # per copy, summed over the cell's chromosomes / ploidy
    log2_fc: np.ndarray               # within-source pairs; NaN between sources
    copy_number: np.ndarray           # (N,) copies per reference bead after the variant
    derivatives: list[dict]           # per derivative: label, copies, n, junctions, p (derivative coordinates)
    offsets: dict[str, int]
    interval: dict | None = None
    notes: list[str] = field(default_factory=list)

    def index(self, source: str, bead: int) -> int:
        return self.offsets[source] + int(bead)

    def locate(self, k: int) -> tuple[Source, int]:
        for s in self.sources:
            o = self.offsets[s.name]
            if o <= k < o + s.n:
                return s, k - o
        raise IndexError(k)

    def new_contacts(self) -> np.ndarray:
        """Per-copy contact probability between different sources after the variant (0 before)."""
        out = np.where(np.isnan(self.log2_fc) & np.isfinite(self.p_after) & ~np.isfinite(self.p_before), self.p_after, np.nan)
        return out


def impact(sources: list[Source], karyotype: Karyotype) -> ImpactV2:
    """Contact changes in reference coordinates for a karyotype of derivatives built from `sources`."""
    names = [s.name for s in sources]
    S = {s.name: s.S_nm2 for s in sources}
    r_c = sources[0].r_c_nm
    s_bond = float(np.median(np.concatenate([np.diag(s.S_nm2, 1) for s in sources])))
    offsets, N = {}, 0
    for s in sources:
        offsets[s.name] = N
        N += s.n
    ploidy = karyotype.ploidy
    before = np.full((N, N), np.nan)
    after = np.zeros((N, N))
    cn = np.zeros(N)
    for s in sources:
        o, n = offsets[s.name], s.n
        p = _p(s.S_nm2, r_c)
        before[o:o + n, o:o + n] = p
        rc = karyotype.reference_copies.get(s.name, ploidy)
        after[o:o + n, o:o + n] += rc * np.nan_to_num(p)
        cn[o:o + n] += rc
    derived = []
    for d in karyotype.derivatives:
        Sd, src, bead = derive_segments(S, d.segments, s_bond)
        Pd = _p(Sd, r_c)
        g = np.array([offsets[names[i]] for i in src]) + bead
        P0 = np.nan_to_num(Pd)
        np.add.at(after, (g[:, None], g[None, :]), d.copies * P0)
        np.add.at(cn, g, d.copies)
        derived.append({"label": d.label(), "copies": d.copies, "n": d.n, "junctions": list(d.junctions),
                        "p": Pd, "ref_index": g})
    after = after / ploidy
    present = cn > 0
    after[~present, :] = np.nan
    after[:, ~present] = np.nan
    np.fill_diagonal(after, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        fc = np.log2(after / before)
    same = np.zeros((N, N), bool)
    for s in sources:
        o = offsets[s.name]
        same[o:o + s.n, o:o + s.n] = True
    fc = np.where(same, fc, np.nan)
    after = np.where(same | (after > 0), after, np.nan)
    return ImpactV2(sources, karyotype, before, after, fc, cn, derived, offsets, None, list(karyotype.notes))


# ======================================================================================
# Read-outs
# ======================================================================================
def _sparse_upper(P: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(P)
    iu = np.triu_indices(n, 1)
    v = P[iu]
    ok = np.isfinite(v) & (v > 0)
    return iu[0][ok], iu[1][ok], v[ok]


def domain_boundaries(P: np.ndarray, resolution: int) -> tuple[list[int], np.ndarray, int]:
    n = len(P)
    w = DOM.window_for(resolution or 10_000, n)
    ci, cj, cm = _sparse_upper(P)
    ins = DOM.insulation(ci, cj, cm, n, w)
    return DOM.boundaries(ins, w), ins, w


def boundary_changes(imp: ImpactV2) -> list[dict]:
    """Per source: insulation boundaries lost / gained (reference coordinates) and domains fused; per
    derivative: domains that span a junction (neo-domains)."""
    rows = []
    for s in imp.sources:
        o, n = imp.offsets[s.name], s.n
        pb = imp.p_before[o:o + n, o:o + n]
        pa = imp.p_after[o:o + n, o:o + n]
        b0, _, _ = domain_boundaries(pb, s.resolution)
        b1, _, _ = domain_boundaries(np.where(np.isfinite(pa), pa, 0.0), s.resolution)
        cn = imp.copy_number[o:o + n]
        for b in b0:
            if cn[b] == 0:
                rows.append({"source": s.name, "bead": b, "change": "boundary deleted"})
            elif all(abs(b - x) > BOUNDARY_TOLERANCE for x in b1):
                rows.append({"source": s.name, "bead": b, "change": "boundary lost (domains fused)"})
        for b in b1:
            if all(abs(b - x) > BOUNDARY_TOLERANCE for x in b0):
                rows.append({"source": s.name, "bead": b, "change": "boundary gained"})
    for k, d in enumerate(imp.derivatives):
        if not d["junctions"]:
            continue
        res = imp.sources[0].resolution
        bd, _, _ = domain_boundaries(np.nan_to_num(d["p"]), res)
        tads = DOM.tads_from_boundaries(bd, d["n"])
        for a, b in tads:
            for j in d["junctions"]:
                if a + 2 <= j <= b - 2:
                    rows.append({"source": f"derivative {k + 1}", "bead": j, "change": "domain spans a junction (neo-domain)",
                                 "domain": (a, b)})
    return rows


def affected_genes(imp: ImpactV2, genes_by_source: dict[str, list[tuple[str, int, int, int]]]) -> list[dict]:
    """Genes per source as (name, tss_bead, first_bead, last_bead), in source bead coordinates. Effects:
    copy number, breakpoint inside, inverted, moved next to a new junction, strongest contact change."""
    rows = []
    near = {}
    for d in imp.derivatives:
        res = imp.sources[0].resolution or 10_000
        reach = max(1, JUNCTION_NEIGHBOURHOOD_BP // res)
        for j in d["junctions"]:
            for k in range(max(0, j - reach), min(d["n"], j + reach)):
                near[int(d["ref_index"][k])] = True
    cuts = {s.name: set() for s in imp.sources}
    inv = {s.name: np.zeros(s.n, bool) for s in imp.sources}
    for d in imp.karyotype.derivatives:
        for p, q in zip(d.segments[:-1], d.segments[1:]):
            cuts[p.source].add(p.a if p.reverse else p.b)
            cuts[q.source].add(q.b if q.reverse else q.a)
        for s in d.segments:
            if s.reverse:
                inv[s.source][s.a:s.b] = True
    pmax = np.fmax(np.nan_to_num(imp.p_before), np.nan_to_num(imp.p_after))
    for s in imp.sources:
        o = imp.offsets[s.name]
        for name, tss, g0, g1 in genes_by_source.get(s.name, []):
            if not 0 <= tss < s.n:
                continue
            effects = []
            cn_body = imp.copy_number[o + max(0, g0):o + min(s.n, g1 + 1)]
            ploidy = imp.karyotype.ploidy
            if cn_body.size and cn_body.min() == 0:
                effects.append("deleted" if cn_body.max() == 0 else "partly deleted")
            elif cn_body.size and cn_body.min() < ploidy:
                effects.append(f"copy loss (CN {cn_body.min():g})")
            elif cn_body.size and cn_body.max() > ploidy:
                effects.append(f"copy gain (CN {cn_body.max():g})")
            if any(g0 < c <= g1 for c in cuts[s.name]):
                effects.append("breakpoint inside the gene")
            if inv[s.name][tss]:
                effects.append("inverted")
            if near.get(o + tss):
                effects.append("next to a new junction")
            k = o + tss
            row = np.where((pmax[k] >= MIN_P), imp.log2_fc[k], np.nan)
            row[max(o, k - 1):k + 2] = np.nan
            strongest = float(row[np.nanargmax(np.abs(row))]) if np.isfinite(row).any() else float("nan")
            newc = imp.new_contacts()[k]
            new_max = float(np.nanmax(newc)) if np.isfinite(newc).any() else 0.0
            if np.isfinite(strongest) and abs(strongest) >= MIN_LOG2:
                effects.append("contacts changed")
            if new_max >= MIN_P:
                effects.append("new contacts with the partner")
            if effects:
                rows.append({"source": s.name, "gene": name, "tss_bead": int(tss), "effects": effects,
                             "strongest_log2_change": strongest, "max_new_contact_p": new_max})
    return rows


def enhancer_promoter(imp: ImpactV2, promoters: list[tuple[str, str, int]], enhancers: list[tuple[str, int]],
                      min_p: float = MIN_P, min_abs_log2: float = MIN_LOG2) -> list[dict]:
    """Promoter (source, gene, bead) x enhancer (source, bead) pairs: lost / gained (same source) and new
    across a junction (different sources, before 0)."""
    out = []
    for src_p, gene, pb in promoters:
        kp = imp.index(src_p, pb)
        for src_e, eb in enhancers:
            ke = imp.index(src_e, eb)
            if src_e == src_p and abs(eb - pb) < 2:
                continue
            b, a = imp.p_before[kp, ke], imp.p_after[kp, ke]
            row = {"gene": gene, "promoter": (src_p, int(pb)), "enhancer": (src_e, int(eb)),
                   "p_before": None if not np.isfinite(b) else float(b), "p_after": None if not np.isfinite(a) else float(a)}
            if src_e != src_p:
                if np.isfinite(a) and a >= min_p:
                    out.append(dict(row, change="new across a junction", log2_fc=float("inf")))
                continue
            if np.isfinite(b) and not np.isfinite(a):
                out.append(dict(row, change="removed", log2_fc=float("-inf")))
                continue
            fc = imp.log2_fc[kp, ke]
            if np.isfinite(fc) and abs(fc) >= min_abs_log2 and max(b, a) >= min_p:
                r = dict(row, change="gained" if fc > 0 else "lost", log2_fc=float(fc))
                if imp.interval is not None:
                    r["log2_fc_90ci"] = (float(imp.interval["lo"][kp, ke]), float(imp.interval["hi"][kp, ke]))
                out.append(r)
    return sorted(out, key=lambda r: -abs(r["log2_fc"]) if math.isfinite(r["log2_fc"]) else -1e9)


def bootstrap(refit_sources, karyotype: Karyotype, reps: int = 8, level: float = 0.9, progress=None) -> dict:
    """refit_sources(r) -> list[Source] refitted on resampled counts (replicate r). Percentile intervals of the
    within-source log2 change, and the per-copy contact probability after, over the replicates."""
    fcs, pas = [], []
    for r in range(reps):
        imp = impact(refit_sources(r), karyotype)
        fcs.append(imp.log2_fc)
        pas.append(imp.p_after)
        if progress is not None:
            progress(r + 1, reps)
    q = (1 - level) / 2
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        fc, pa = np.stack(fcs), np.stack(pas)
        return {"lo": np.nanquantile(fc, q, axis=0), "hi": np.nanquantile(fc, 1 - q, axis=0),
                "p_after_lo": np.nanquantile(pa, q, axis=0), "p_after_hi": np.nanquantile(pa, 1 - q, axis=0),
                "reps": reps, "level": level}


def ranking_score(genes: list[dict], ep: list[dict], boundaries: list[dict], interval: bool) -> dict:
    """Transparent heuristic for 'most likely to disrupt gene regulation' (NOT validated):
    3 per gene deleted / broken, 2 per enhancer-promoter pair lost, removed, gained or new (only those whose
    interval excludes no change, when refits were run), 2 per boundary lost or neo-domain, 1 per other affected gene."""
    hard = sum(1 for g in genes if any(e in ("deleted", "partly deleted", "breakpoint inside the gene") for e in g["effects"]))
    soft = len(genes) - hard

    def sure(r):
        if not interval or "log2_fc_90ci" not in r:
            return True
        lo, hi = r["log2_fc_90ci"]
        return lo > 0 or hi < 0
    ep_n = sum(1 for r in ep if sure(r))
    bnd = sum(1 for b in boundaries if b["change"] in ("boundary lost (domains fused)", "domain spans a junction (neo-domain)"))
    return {"score": 3 * hard + 2 * ep_n + 2 * bnd + soft, "genes_deleted_or_broken": hard, "ep_pairs": ep_n,
            "boundaries_lost_or_neo": bnd, "other_genes": soft}


# ======================================================================================
# Copy number from VCF genotypes, CNV BED, or Hi-C coverage
# ======================================================================================
def vcf_genotypes(data: bytes | str) -> dict[int, dict]:
    """Per VCF data line number: genotype of the first sample (FORMAT GT) and INFO/FORMAT CN when present."""
    from .svio import _text
    out = {}
    for ln, line in enumerate(_text(data).splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        info = dict(kv.partition("=")[::2] for kv in f[7].split(";"))
        rec: dict = {}
        if "CN" in info:
            try:
                rec["cn"] = float(info["CN"])
            except ValueError:
                pass
        if len(f) >= 10:
            keys, vals = f[8].split(":"), f[9].split(":")
            fmt = dict(zip(keys, vals))
            gt = fmt.get("GT", "").replace("|", "/")
            if gt:
                rec["gt"] = gt
            if "CN" in fmt:
                try:
                    rec["cn"] = float(fmt["CN"])
                except ValueError:
                    pass
        out[ln] = rec
    return out


def zygosity(gt: str | None) -> str:
    if not gt:
        return "heterozygous"
    alleles = [a for a in gt.replace("|", "/").split("/") if a not in (".", "")]
    return "homozygous" if alleles and all(a != "0" for a in alleles) else "heterozygous"


def read_cnv_bed(data: bytes | str, chrom: str, bin0: int, n: int, resolution: int) -> list[tuple[int, int, float]]:
    """chrom start end copy_number [...] -> [(a, b, CN)] in local bead coordinates of the window."""
    from .svio import _default_norm, _text
    out = []
    for line in _text(data).splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = re.split(r"\t|\s+", line.strip())
        if len(f) < 4:
            continue
        try:
            c, s, e, cn = _default_norm(f[0]), int(f[1]), int(f[2]), float(f[3])
        except (ValueError, KeyError):
            continue
        if c != chrom:
            continue
        a, b = max(0, s // resolution - bin0), min(n, -(-e // resolution) - bin0)
        if b > a:
            out.append((a, b, cn))
    return out


def estimate_copy_number(ci: np.ndarray, cj: np.ndarray, cm: np.ndarray, n: int, ploidy: float = 2.0,
                         smooth: int = 5, min_bins: int = 3) -> tuple[np.ndarray, list[tuple[int, int, float]]]:
    """Copy number per bin from Hi-C coverage: ploidy x (row sum / median row sum), median-smoothed over
    `smooth` bins, rounded; segments of >= min_bins with CN != ploidy. Coverage also follows mappability, GC
    and restriction-site density, so this estimate is labelled not validated."""
    cov = np.zeros(n)
    np.add.at(cov, np.asarray(ci, np.int64), cm)
    np.add.at(cov, np.asarray(cj, np.int64), cm)
    ok = cov > 0
    if ok.sum() < 10:
        return np.full(n, np.nan), []
    med = float(np.median(cov[ok]))
    raw = np.where(ok, ploidy * cov / med, np.nan)
    h = max(1, smooth // 2)
    sm = np.array([np.nanmedian(raw[max(0, i - h):i + h + 1]) if np.isfinite(raw[max(0, i - h):i + h + 1]).any() else np.nan
                   for i in range(n)])
    est = np.round(sm)
    segs, i = [], 0
    while i < n:
        if not np.isfinite(est[i]) or est[i] == ploidy:
            i += 1
            continue
        j = i
        while j < n and est[j] == est[i]:
            j += 1
        if j - i >= min_bins:
            segs.append((i, j, float(est[i])))
        i = j
    return est, segs
