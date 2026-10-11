"""
Pipelines shared by the command line (chronocell.cli), the REST API (chronocell.api), batch mode, the local
job queue (chronocell.jobs) and the app (Phase B). Each takes contacts (chronocell.contacts_io.Contacts)
and returns plain dicts / DataFrames plus files written to an output folder; nothing here imports Streamlit.

    analyze(contacts)                      the analysis suite (chronocell.analysis)
    diff(contacts_a, contacts_b)           two-condition differential analysis (chronocell.differential)
    impact(sources, variants / joins / segments / copy number)
                                           variant impact engine v2 (chronocell.sv_engine), with genes and
                                           enhancer-promoter pairs when an assembly and a signal are given
Every output carries its standing: the analysis suite and the differential statistics are measured
procedures (their accuracy gates are Gates 6 and 7), the variant impact is a mechanism simulator.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import analysis as AN, contacts_io as IO, differential as DF, sv_engine as SV

STANDING = {
    "analyze": "Analysis suite: loops, domains, compartments; accuracy against reference calls is Gate 6 "
               "(validation/RESULTS.md).",
    "diff": "Differential statistics: false-discovery control is Gate 7 (validation/RESULTS.md); without "
            "replicates no statistics are reported.",
    "impact": "Variant impact: mechanism simulator, not validated (Gate 4d in validation/RESULTS.md).",
}
P_ADJACENT = 0.3        # adjacent contact probability for fits made here (the Gate 1 / SV-test setting)
B0_NM = 50.0


# ======================================================================================
# Analysis suite
# ======================================================================================
def analyze(c: IO.Contacts, out: Path | None = None, loops: bool = True, orient: np.ndarray | None = None,
            loop_fdr: float = 0.1) -> dict:
    t0 = time.time()
    rep = AN.run_suite(c.ci, c.cj, c.cm, c.n, c.resolution, c.weights, orient, loops, loop_fdr)
    pos = lambda k: c.start + int(k) * c.resolution                       # noqa: E731
    loops_df = pd.DataFrame([{"chrom1": c.chrom, "start1": pos(l.i), "end1": pos(l.i) + c.resolution, "chrom2": c.chrom,
                              "start2": pos(l.j), "end2": pos(l.j) + c.resolution, "observed": l.observed,
                              "expected_donut": round(l.expected_donut, 3), "oe_donut": round(l.oe_donut, 3),
                              "q_donut": l.q["donut"], "cluster_size": l.cluster_size} for l in rep.loops],
                            columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "observed",
                                     "expected_donut", "oe_donut", "q_donut", "cluster_size"])
    bnd = pd.DataFrame([{"chrom": c.chrom, "start": pos(b), "end": pos(b) + c.resolution, "method": m}
                        for m, bs in (("insulation", rep.insulation_boundaries), ("topdom", rep.topdom_boundaries))
                        for b in bs])
    dom = pd.DataFrame([{"chrom": c.chrom, "start": pos(a), "end": pos(b), "inside_over_flank": round(r, 3)}
                        for a, b, r in rep.arrowhead_domains], columns=["chrom", "start", "end", "inside_over_flank"])
    tracks = pd.DataFrame({"chrom": c.chrom, "start": [pos(k) for k in range(c.n)],
                           "end": [pos(k) + c.resolution for k in range(c.n)], "insulation": rep.insulation,
                           "compartment_eigenvector": rep.compartment})
    ps = pd.DataFrame({"separation_bp": rep.p_s["separation_bp"], "p": rep.p_s["p"]})
    summary = rep.summary() | {"region": f"{c.chrom}:{c.start:,}-{c.start + c.n * c.resolution:,}",
                               "seconds": round(time.time() - t0, 2), "standing": STANDING["analyze"]}
    result = {"summary": summary, "loops": loops_df, "boundaries": bnd, "domains": dom, "tracks": tracks, "p_s": ps,
              "notes": c.notes + rep.notes}
    if out is not None:
        write_analysis(result, Path(out))
    return result


def write_analysis(result: dict, out: Path) -> list[Path]:
    from . import exports as EX
    out.mkdir(parents=True, exist_ok=True)
    files = [EX.write_text(out / "loops.bedpe", EX.loops_bedpe(result["loops"])),
             EX.write_text(out / "loops_juicebox_2d.txt", EX.juicebox_2d(result["loops"])),
             EX.write_text(out / "boundaries.bed", EX.bed(result["boundaries"], name_col="method")),
             EX.write_text(out / "domains.bed", EX.bed(result["domains"])),
             EX.write_text(out / "insulation.bedGraph", EX.bedgraph(result["tracks"], "insulation")),
             EX.write_text(out / "compartments.bedGraph", EX.bedgraph(result["tracks"], "compartment_eigenvector"))]
    result["p_s"].to_csv(out / "p_of_s.tsv", sep="\t", index=False)
    (out / "summary.json").write_text(json.dumps({"summary": result["summary"], "notes": result["notes"]}, indent=1,
                                                 default=float), encoding="utf-8")
    return files + [out / "p_of_s.tsv", out / "summary.json"]


# ======================================================================================
# Differential analysis
# ======================================================================================
def diff(cond_a: list[IO.Contacts], cond_b: list[IO.Contacts], fdr: float = 0.05, max_sep_bp: int = 2_000_000,
         min_count: float = 5.0, out: Path | None = None, loops: bool = True) -> dict:
    ref = cond_a[0]
    for c in list(cond_a) + list(cond_b):
        if (c.n, c.resolution, c.chrom, c.start) != (ref.n, ref.resolution, ref.chrom, ref.start):
            raise IO.ContactFileError("Every replicate must cover the same region at the same resolution.")
    t0 = time.time()
    res = DF.compare([c.dense() for c in cond_a], [c.dense() for c in cond_b], ref.resolution, fdr, max_sep_bp,
                     min_count, loops=loops)
    summary = res.summary() | {"replicates_a": len(cond_a), "replicates_b": len(cond_b), "seconds": round(time.time() - t0, 2),
                               "region": f"{ref.chrom}:{ref.start:,}-{ref.start + ref.n * ref.resolution:,}",
                               "standing": STANDING["diff"]}
    result = {"summary": summary, "pixels": res.pixels, "significant": res.significant, "loops": res.loops,
              "boundaries": res.boundaries, "compartments": res.compartments, "notes": res.notes, "result": res,
              "chrom": ref.chrom, "start": ref.start, "resolution": ref.resolution}
    if out is not None:
        write_diff(result, Path(out))
    return result


def write_diff(result: dict, out: Path) -> list[Path]:
    from . import exports as EX
    out.mkdir(parents=True, exist_ok=True)
    c, s, r = result["chrom"], result["start"], result["resolution"]
    files = [EX.write_text(out / "differential_pixels.bedpe", DF.to_bedpe(result["significant"], c, s, r))]
    result["pixels"].assign(chrom=c, start1=s + result["pixels"]["i"] * r, start2=s + result["pixels"]["j"] * r) \
        .to_csv(out / "differential_pixels_all.tsv", sep="\t", index=False)
    for name in ("loops", "boundaries", "compartments"):
        df = result[name]
        if df is not None:
            df.to_csv(out / f"differential_{name}.csv", index=False)
            files.append(out / f"differential_{name}.csv")
    (out / "summary.json").write_text(json.dumps({"summary": result["summary"], "notes": result["notes"]}, indent=1,
                                                 default=float), encoding="utf-8")
    return files + [out / "differential_pixels_all.tsv", out / "summary.json"]


# ======================================================================================
# Variant impact
# ======================================================================================
def fit_source(name: str, c: IO.Contacts, b0_nm: float = B0_NM, p_adjacent: float = P_ADJACENT, progress=None,
               signal: np.ndarray | None = None) -> tuple[SV.Source, object]:
    """Population model of a contact window (v3.3 up to 400 beads, v4 above) as an engine Source."""
    from . import ensemble as ENS, population as POP
    off = c.ci != c.cj
    if c.n > POP.V33_MAX_BEADS:
        res = POP.fit_population_from_counts(c.ci[off], c.cj[off], c.cm[off], c.n, None, b0_nm=b0_nm, p_adjacent=p_adjacent,
                                             cfg=POP.config_for(c.n), progress=progress)
    else:
        res = ENS.fit_from_counts(c.ci[off], c.cj[off], c.cm[off], c.n, None, b0_nm=b0_nm, p_adjacent=p_adjacent,
                                  progress=progress)
    src = SV.source_from_result(name, res, chrom=c.chrom, bin0=c.start // c.resolution, resolution=c.resolution,
                                signal=signal)
    return src, res


def _genes_for(src: SV.Source, assembly: str | None) -> list[tuple[str, int, int, int]]:
    from . import genes as G
    if not src.chrom or not src.resolution:
        return []
    start, end = src.bin0 * src.resolution, (src.bin0 + src.n) * src.resolution
    tab = G.in_region(src.chrom, start, end, assembly)
    out = []
    for _, g in tab.iterrows():
        out.append((str(g["name"]), int(g["tss"]) // src.resolution - src.bin0, int(g["start"]) // src.resolution - src.bin0,
                    (int(g["end"]) - 1) // src.resolution - src.bin0))
    return out


def impact(sources: list[SV.Source], karyotype: SV.Karyotype, assembly: str | None = None,
           refit=None, reps: int = 0, progress=None) -> dict:
    """Run the engine and every read-out. refit(r) -> list[Source] enables the confidence (reps refits)."""
    t0 = time.time()
    imp = SV.impact(sources, karyotype)
    if refit is not None and reps > 0:
        imp.interval = SV.bootstrap(refit, karyotype, reps=reps, progress=progress)
    genes = SV.affected_genes(imp, {s.name: _genes_for(s, assembly) for s in sources})
    promoters = [(s.name, g[0], g[1]) for s in sources for g in _genes_for(s, assembly) if 0 <= g[1] < s.n]
    enhancers = []
    for s in sources:
        if s.signal is not None and np.isfinite(s.signal).any():
            thr = float(np.nanpercentile(s.signal, 90))
            enhancers += [(s.name, int(k)) for k in np.flatnonzero(np.nan_to_num(s.signal, nan=-np.inf) >= thr)]
    ep = SV.enhancer_promoter(imp, promoters, enhancers) if enhancers else []
    bnd = SV.boundary_changes(imp)
    score = SV.ranking_score(genes, ep, bnd, imp.interval is not None)
    clin = annotate_genes([g["gene"] for g in genes])
    for g in genes:
        g.update(clin.get(g["gene"], {}))
    summary = {"derivatives": [d["label"] for d in imp.derivatives], "copies": [d["copies"] for d in imp.derivatives],
               "lost_pieces": [s.label() for s in karyotype.lost], "genes_affected": len(genes), "ep_pairs": len(ep),
               "boundary_changes": len(bnd), "ranking": score, "seconds": round(time.time() - t0, 2),
               "confidence": f"{imp.interval['reps']} refits" if imp.interval else "not estimated",
               "standing": STANDING["impact"]}
    return {"summary": summary, "impact": imp, "genes": genes, "ep": ep, "boundaries": bnd, "notes": imp.notes}


def annotate_genes(names: list[str]) -> dict[str, dict]:
    """Information-only gene annotations available offline: ChronoCell's curated cancer / neuro lists, and the
    ClinVar gene summary if it was downloaded (chronocell.annotations). Never a diagnosis."""
    from . import genes as G
    out = {n: {"curated_list": G.category(n)} for n in names if G.category(n)}
    try:
        from . import annotations as AO
        cv = AO.clinvar_genes(download=False)
    except (ImportError, OSError, ValueError):
        cv = None
    if cv is not None:
        for n in names:
            if n in cv:
                out.setdefault(n, {})["clinvar_pathogenic_alleles"] = int(cv[n])
    return out


def rank_variants(results: list[tuple[str, dict]]) -> pd.DataFrame:
    """'Most likely to disrupt gene regulation' (heuristic, not validated): variants by ranking score."""
    rows = [{"variant": label, **r["summary"]["ranking"], "genes_affected": r["summary"]["genes_affected"]}
            for label, r in results]
    df = pd.DataFrame(rows)
    return df.sort_values("score", ascending=False).reset_index(drop=True) if len(df) else df
