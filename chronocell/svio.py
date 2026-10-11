"""
Structural-variant files: VCF (SVTYPE / symbolic ALT / breakend ALT) and BEDPE.

Every variant becomes a `Variant` with 0-based, half-open reference coordinates:
    DEL / DUP / INV : chrom, start, end
    BND / TRA       : chrom, start (this side's breakpoint), mate_chrom, mate_pos, orientation
Chromosome names are normalised by the caller-supplied function (chronocell.genome.normalize_chrom
by default), so 'chr9', '9', 'Chr9' and RefSeq/GenBank accessions all work. Lines that cannot be
read are returned as `skipped` with the reason; nothing is silently dropped.
"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field
from typing import Callable

KINDS = ("DEL", "DUP", "INV", "BND")


@dataclass
class Variant:
    kind: str                   # DEL | DUP | INV | BND
    chrom: str
    start: int                  # 0-based
    end: int                    # half-open (= start + 1 for a breakend)
    mate_chrom: str | None = None
    mate_pos: int | None = None
    orientation: str = ""       # BEDPE strands or VCF bracket form, informative only
    vid: str = ""
    source_line: int = 0

    @property
    def length(self) -> int:
        return self.end - self.start if self.kind != "BND" else 0


@dataclass
class VariantFile:
    variants: list[Variant] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    format: str = ""


def _text(data: bytes | str) -> str:
    if isinstance(data, str):
        return data
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8", "replace")


def _info(s: str) -> dict[str, str]:
    out = {}
    for kv in s.split(";"):
        k, _, v = kv.partition("=")
        out[k.strip()] = v.strip()
    return out


_BND = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def read_vcf(data: bytes | str, norm: Callable[[str], str] | None = None) -> VariantFile:
    norm = norm or _default_norm
    vf = VariantFile(format="VCF")
    for ln, line in enumerate(_text(data).splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            vf.skipped.append(f"line {ln}: fewer than 8 VCF columns")
            continue
        chrom, pos, vid, ref, alt, _, _, info = f[:8]
        try:
            chrom = norm(chrom)
            p0 = int(pos) - 1
        except (ValueError, KeyError) as exc:
            vf.skipped.append(f"line {ln}: {exc}")
            continue
        inf = _info(info)
        kind = inf.get("SVTYPE", "").upper()
        m = re.match(r"<([A-Z]+)", alt)
        if not kind and m:
            kind = m.group(1)
        bnd = _BND.search(alt)
        if bnd and kind in ("", "BND", "TRA"):
            kind = "BND"
        kind = {"TRA": "BND", "CNV": "", "DUP:TANDEM": "DUP"}.get(kind, kind)
        if kind.startswith("DUP"):
            kind = "DUP"
        if kind not in KINDS:
            vf.skipped.append(f"line {ln}: unsupported SV type '{kind or alt[:20]}'")
            continue
        try:
            if kind == "BND":
                if bnd:
                    mate_chrom, mate_pos = norm(bnd.group(1)), int(bnd.group(2)) - 1
                else:
                    mate_chrom, mate_pos = norm(inf["CHR2"]), int(inf["END"]) - 1
                vf.variants.append(Variant("BND", chrom, p0, p0 + 1, mate_chrom, mate_pos, alt, vid, ln))
            else:
                if "END" in inf:
                    end = int(inf["END"])
                elif "SVLEN" in inf:
                    end = p0 + 1 + abs(int(inf["SVLEN"].split(",")[0]))
                else:
                    raise ValueError("no END or SVLEN")
                # symbolic alleles: POS is the base before the event (VCF 4.x), so the event starts at POS
                start = p0 + 1 if alt.startswith("<") else p0
                if end <= start:
                    raise ValueError(f"END {end} is not after POS")
                vf.variants.append(Variant(kind, chrom, start, end, vid=vid, source_line=ln))
        except (ValueError, KeyError) as exc:
            vf.skipped.append(f"line {ln}: {exc}")
    return vf


def read_bedpe(data: bytes | str, norm: Callable[[str], str] | None = None) -> VariantFile:
    """chrom1 start1 end1 chrom2 start2 end2 [name score strand1 strand2 [type]]. Type, when absent,
    comes from the strands: same chromosome (+,-) deletion, (-,+) duplication, (+,+)/(-,-) inversion;
    different chromosomes: translocation (BND)."""
    norm = norm or _default_norm
    vf = VariantFile(format="BEDPE")
    for ln, line in enumerate(_text(data).splitlines(), 1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = re.split(r"\t|\s+", line.strip())
        if len(f) < 6:
            vf.skipped.append(f"line {ln}: fewer than 6 BEDPE columns")
            continue
        try:
            c1, s1, e1, c2, s2, e2 = norm(f[0]), int(f[1]), int(f[2]), norm(f[3]), int(f[4]), int(f[5])
        except (ValueError, KeyError) as exc:
            if ln == 1:
                continue                                    # header row
            vf.skipped.append(f"line {ln}: {exc}")
            continue
        name = f[6] if len(f) > 6 else f"sv{ln}"
        st1, st2 = (f[8], f[9]) if len(f) > 9 else ("", "")
        kind = f[10].upper() if len(f) > 10 else ""
        kind = {"TRA": "BND", "TRANSLOCATION": "BND", "DELETION": "DEL", "DUPLICATION": "DUP",
                "INVERSION": "INV"}.get(kind, kind)
        if c1 != c2:
            vf.variants.append(Variant("BND", c1, (s1 + e1) // 2, (s1 + e1) // 2 + 1, c2, (s2 + e2) // 2,
                                       f"{st1}{st2}", name, ln))
            continue
        if kind not in KINDS:
            kind = {("+", "-"): "DEL", ("-", "+"): "DUP", ("+", "+"): "INV", ("-", "-"): "INV"}.get((st1, st2), "")
        if kind not in ("DEL", "DUP", "INV"):
            vf.skipped.append(f"line {ln}: type unknown (give strands or an 11th type column)")
            continue
        a, b = sorted(((s1 + e1) // 2, (s2 + e2) // 2))
        if b <= a:
            vf.skipped.append(f"line {ln}: breakpoints coincide")
            continue
        vf.variants.append(Variant(kind, c1, a, b, orientation=f"{st1}{st2}", vid=name, source_line=ln))
    return vf


def read_any(data: bytes | str, name: str, norm: Callable[[str], str] | None = None) -> VariantFile:
    """VCF if the name or the content says so (header line, or an SVTYPE= INFO field), else BEDPE."""
    low = name.lower()
    text = _text(data)
    if (low.endswith((".vcf", ".vcf.gz")) or text.lstrip().startswith("##fileformat=VCF") or "#CHROM" in text
            or "SVTYPE=" in text):
        return read_vcf(text, norm)
    return read_bedpe(text, norm)


OPERATION = {"DEL": "deletion", "DUP": "duplication", "INV": "inversion", "BND": "translocation"}


def label(v: Variant) -> str:
    if v.kind == "BND":
        return f"BND {v.chrom}:{v.start + 1:,} - {v.mate_chrom}:{(v.mate_pos or 0) + 1:,}" + (f" ({v.vid})" if v.vid else "")
    return f"{v.kind} {v.chrom}:{v.start + 1:,}-{v.end:,} ({v.length / 1e3:,.0f} kb)" + (f" ({v.vid})" if v.vid else "")


def to_scenario(v: Variant, chrom_name: str, resolution: int, bin0: int, n: int,
                partners: tuple[str, ...] | list[str] = ()) -> tuple[str, dict] | str:
    """One variant as an operation on the loaded window, in local bead coordinates: (operation, params),
    or the reason it cannot be applied there. A breakend joins this chromosome from pter to the
    breakpoint with the partner from its breakpoint to qter (orientation is not modelled)."""
    if v.kind in ("DEL", "DUP", "INV"):
        if v.chrom != chrom_name:
            return f"on {v.chrom}, not {chrom_name}"
        a = v.start // resolution - bin0
        b = -(-v.end // resolution) - bin0
        if b <= 0 or a >= n:
            return "outside the loaded window"
        a, b = max(0, a), min(n, b)
        if v.kind == "INV" and b - a < 2:
            return "inversion shorter than two beads at this resolution"
        if v.kind == "DEL" and (a == 0 and b == n):
            return "deletes the whole loaded window"
        params: dict = {"a": int(a), "b": int(b)}
        if v.kind == "DUP":
            params["copies"] = 1
        return OPERATION[v.kind], params
    if v.chrom == chrom_name and v.mate_chrom != chrom_name:
        here, partner, ppos = v.start, v.mate_chrom, v.mate_pos
    elif v.mate_chrom == chrom_name and v.chrom != chrom_name:
        here, partner, ppos = v.mate_pos, v.chrom, v.start
    elif v.chrom == chrom_name:
        return "both breakends on this chromosome: give it as DEL, DUP or INV"
    else:
        return f"joins {v.chrom} and {v.mate_chrom}, not {chrom_name}"
    if partners and partner not in partners:
        return f"partner {partner} is not a main chromosome of this assembly"
    bp = int(here) // resolution - bin0
    if not 2 <= bp <= n - 1:
        return "breakpoint outside the loaded window"
    return "translocation", {"breakpoint": int(bp), "partner": partner, "partner_start_bp": int(ppos or 0)}


def deletion_segments(variants: list[Variant], chrom_name: str, resolution: int, bin0: int, n: int) -> list[tuple[int, int]]:
    """Every deletion on this chromosome inside the window, as local half-open bead intervals."""
    out = []
    for v in variants:
        if v.kind == "DEL":
            sc = to_scenario(v, chrom_name, resolution, bin0, n)
            if isinstance(sc, tuple):
                out.append((sc[1]["a"], sc[1]["b"]))
    return sorted(out)


def _default_norm(name: str) -> str:
    try:
        from .genome import normalize_chrom
        return normalize_chrom(name)
    except ImportError:          # pragma: no cover
        n = name.strip()
        return n if n.startswith("chr") else f"chr{n}"
