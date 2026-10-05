"""
02 · 4D dynamics → 05 Variant impact engine v2 (Phase B1): rearrangements given by their joins (breakend pairs),
derivative chromosomes typed as segments, copy number (CNV BED, or estimated from Hi-C coverage), or a variant
file with genotypes; an optional second chromosome (partner window with its own contacts, or a homogeneous chain
when there are none). Outputs per variant: contact changes, domain boundaries lost / gained and domains spanning a
junction, genes affected (with ClinVar counts when downloaded), enhancer–promoter pairs, a confidence from refits,
and a ranking across the variants of a file. Mechanism simulator, not validated: labelled so everywhere.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from ui.common import Dataset, banner, esc, html, readout, warning_card

ss = st.session_state
MAX_BEADS = 1500
_JOIN = re.compile(r"^\s*([\w.]+):([\d,_]+):([LR])\s*-\s*([\w.]+):([\d,_]+):([LR])\s*$")


def _window_source(ds: Dataset, frame: int):
    """(Source A, lo, hi, population key) from the largest population of this dataset and frame."""
    from chronocell import sv_engine as SV
    from ui import variant_impact as VI
    pops = VI.populations(ds, frame)
    if not pops:
        return None
    lo, hi, key, res = pops[0]
    if hi - lo > MAX_BEADS:
        c = (lo + hi) // 2
        w = (max(lo, c - MAX_BEADS // 2), min(hi, c + MAX_BEADS // 2))
        src = SV.source_from_result("A", res, (w[0] - lo, w[1] - lo), chrom=ds.chrom.name, bin0=ds.bin0 + w[0],
                                    resolution=ds.chrom.resolution, signal=None if ds.signal_is_placeholder else
                                    np.asarray(ds.epi[w[0]:w[1]], float))
        return src, w[0], w[1], key
    src = SV.source_from_result("A", res, None, chrom=ds.chrom.name, bin0=ds.bin0 + lo, resolution=ds.chrom.resolution,
                                signal=None if ds.signal_is_placeholder else np.asarray(ds.epi[lo:hi], float))
    return src, lo, hi, key


def parse_joins(text: str, sources: list) -> list:
    """'A:120:L-B:40:R, ...' (bead coordinates) or 'chr9:131000000:L-chr22:23200000:R' (bp, mapped to the windows)."""
    from chronocell import sv_engine as SV
    by_name = {s.name: s for s in sources}
    joins = []
    for part in [p for p in re.split(r"[,;\n]", text or "") if p.strip()]:
        m = _JOIN.match(part)
        if not m:
            raise ValueError(f"Cannot read join '{part.strip()}': write SOURCE:POSITION:SIDE-SOURCE:POSITION:SIDE.")
        ends = []
        for name, pos, side in ((m.group(1), m.group(2), m.group(3)), (m.group(4), m.group(5), m.group(6))):
            pos = int(pos.replace(",", "").replace("_", ""))
            if name in by_name:
                ends.append((name, pos, side))
                continue
            hit = [s for s in sources if s.chrom == name and s.bin0 * s.resolution <= pos < (s.bin0 + s.n) * s.resolution]
            if not hit:
                raise ValueError(f"{name}:{pos:,} is not inside a loaded window.")
            ends.append((hit[0].name, int(round(pos / hit[0].resolution)) - hit[0].bin0, side))
        joins.append(SV.Join(*ends[0], *ends[1], part.strip()))
    return joins


def _partner(ds: Dataset, src_a):
    """Optional second chromosome: a contact file region fitted here, or a homogeneous chain."""
    from chronocell import sv_engine as SV
    kind = st.segmented_control("Second chromosome", ["None", "Contact file", "No data (homogeneous chain)"],
                                default="None", required=True, key="ve_partner_kind",
                                help="For translocations and other events joining two chromosomes.")
    if kind == "None":
        return None
    if kind == "No data (homogeneous chain)":
        c1, c2, c3 = st.columns(3)
        chrom = c1.text_input("Partner chromosome", "chr22", key="ve_pc")
        start = c2.number_input("Partner window start (bp)", 0, value=23_000_000, step=100_000, key="ve_ps")
        n = c3.number_input("Partner beads", 20, 1000, 200, 10, key="ve_pn")
        S = SV.homogeneous_chain(src_a.S_nm2, int(n))
        return SV.Source("B", S, src_a.r_c_nm, chrom=chrom, bin0=int(start) // src_a.resolution, resolution=src_a.resolution,
                         note="homogeneous chain (no partner data)")
    path = st.text_input("Partner map: local path or URL (.hic / .mcool / .cool / .pairs)", key="ve_pf")
    c1, c2, c3 = st.columns(3)
    chrom = c1.text_input("Partner chromosome", "chr22", key="ve_pc2")
    start = c2.number_input("From (bp)", 0, value=23_000_000, step=100_000, key="ve_ps2")
    end = c3.number_input("To (bp)", 0, value=25_000_000, step=100_000, key="ve_pe2")
    key = f"{path}|{chrom}|{start}|{end}|{src_a.resolution}"
    if ss.get("ve_partner", {}).get("key") == key:
        return ss.ve_partner["source"]
    if path and st.button("Fit partner window", key="ve_pfit", icon=":material/scatter_plot:"):
        from chronocell import contacts_io as IO, pipelines as PL
        try:
            with st.spinner("Reading and fitting the partner window…"):
                c = IO.read_region(path, chrom, int(start), int(end), src_a.resolution)
                src, _ = PL.fit_source("B", c)
        except (IO.ContactFileError, ValueError, OSError) as exc:
            warning_card("The partner window could not be fitted", str(exc))
            return None
        ss.ve_partner = {"key": key, "source": src}
        return src
    return None


def render(ds: Dataset, frame: int, b0: float) -> None:
    from chronocell import pipelines as PL, svio, sv_engine as SV
    from ui import analysis_view as AV
    banner("<b>Mechanism simulator, not validated.</b> Every change on this panel is a prediction of the population "
           "model rearranged exactly as the joins say; it is not a measurement.", "warn")
    AV.standing_banner("4d")
    AV.standing_banner("4c")
    if not ds.has_contacts:
        html('<p class="cc-note">Needs contacts: the engine rearranges a population model fitted to them.</p>')
        return
    got = _window_source(ds, frame)
    if got is None:
        html('<p class="cc-note">Build a population model of a window first (01 · 3D structure → 03 Model & convergence, '
             'or "Fit population around the variant" in 04 above).</p>')
        return
    src_a, lo, hi, pop_key = got
    ch = ds.chrom
    html(f'<p class="cc-note">Window A: {ch.name}:{float(ch.bin_start(ds.bin0 + lo)) / 1e6:.2f}–'
         f'{float(ch.bin_end(ds.bin0 + hi - 1)) / 1e6:.2f} Mb ({hi - lo:,} beads, beads numbered from 0).</p>')
    src_b = _partner(ds, src_a)
    sources = [src_a] + ([src_b] if src_b is not None else [])
    sizes = {s.name: s.n for s in sources}
    mode = st.segmented_control("Describe the variant by", ["Joins", "Segments", "Copy number", "Variant file"],
                                default="Joins", required=True, key="ve_mode")
    zyg = st.segmented_control("Zygosity", ["heterozygous", "homozygous"], default="heterozygous", required=True,
                               key="ve_zyg", help="Heterozygous: one derivative and one normal homolog.")
    cases: list[tuple[str, SV.Karyotype]] = []
    try:
        if mode == "Joins":
            txt = st.text_area("Joins (breakend pairs)", "A:40:L-A:60:R", key="ve_joins", height=70,
                               help="SOURCE:POSITION:SIDE - SOURCE:POSITION:SIDE. Side L = the piece ending at the cut, R = "
                                    "the piece starting there. Positions are beads of the window, or chromosome:bp "
                                    "(e.g. chr9:130700000:L-chr22:23290000:R). A deletion of beads 40-59 is A:40:L-A:60:R; "
                                    "an inversion adds A:40:R-A:60:L... as two joins.")
            joins = parse_joins(txt, sources)
            if joins:
                cases.append((txt.strip(), SV.karyotype_from_joins(sizes, joins, zyg)))
        elif mode == "Segments":
            txt = st.text_area("Derivative chromosome(s): segments joined in order; one derivative per line",
                               "A:0-40 + A:60-120(-) + A:40-60", key="ve_segs", height=70)
            derivs = [SV.parse_segments(t) for t in txt.splitlines() if t.strip()]
            if derivs:
                ref = {s: (0.0 if zyg == "homozygous" else 1.0) for s in sizes}
                cases.append((txt.strip(), SV.karyotype_from_segments(sizes, derivs, [2.0 if zyg == "homozygous" else 1.0]
                                                                      * len(derivs), ref)))
        elif mode == "Copy number":
            how = st.segmented_control("Copy number from", ["CNV BED file", "Hi-C coverage (estimate)"],
                                       default="CNV BED file", required=True, key="ve_cnsrc")
            if how == "CNV BED file":
                up = st.file_uploader("CNV BED: chrom start end copy_number", type=["bed", "txt", "tsv"], key="ve_cnv")
                segs = SV.read_cnv_bed(up.getvalue(), ch.name, src_a.bin0, src_a.n, src_a.resolution) if up else []
            else:
                m = (ds.ci >= lo) & (ds.ci < hi) & (ds.cj >= lo) & (ds.cj < hi)
                est, segs = SV.estimate_copy_number(ds.ci[m] - lo, ds.cj[m] - lo, ds.cm[m], hi - lo)
                html('<p class="cc-note">Estimated from Hi-C coverage (row sums / median, smoothed, rounded). Coverage also '
                     'follows mappability and GC: <b>not validated</b>.</p>')
            if segs:
                st.dataframe(pd.DataFrame(segs, columns=["first bead", "end bead", "copy number"]), hide_index=True, key="ve_cn_df")
                cases.append(("copy number", SV.karyotype_from_copy_number(src_a.n, segs, "A")))
            else:
                html('<p class="cc-note">No copy-number change in this window.</p>')
        else:
            up = st.file_uploader("Variant file (VCF with GT, or BEDPE)", type=["vcf", "gz", "bedpe", "txt", "tsv"], key="ve_vcf")
            if up is not None:
                data = up.getvalue()
                vf = svio.read_any(data, up.name)
                gts = SV.vcf_genotypes(data) if vf.format == "VCF" else {}
                windows = {s.name: (s.chrom, s.bin0, s.n, s.resolution) for s in sources}
                skipped = []
                for v in vf.variants:
                    js = SV.joins_for_variant(v, windows)
                    if isinstance(js, str):
                        skipped.append(f"{svio.label(v)}: {js}")
                        continue
                    z = SV.zygosity(gts.get(v.source_line, {}).get("gt")) if gts else zyg
                    try:
                        cases.append((svio.label(v), SV.karyotype_from_joins(sizes, js, z)))
                    except ValueError as exc:
                        skipped.append(f"{svio.label(v)}: {exc}")
                if skipped:
                    html('<p class="cc-note">Not applied: ' + esc("; ".join(skipped[:8])) + ("…" if len(skipped) > 8 else "")
                         + "</p>")
    except ValueError as exc:
        warning_card("The variant description could not be used", str(exc))
        return
    refits = st.toggle("Confidence from 8 refits on resampled counts (slower)", False, key="ve_boot")
    if not cases or not st.button("Run the engine", key="ve_go", type="primary", icon=":material/hub:"):
        last = ss.get("ve_last")
        if not last or last.get("pop") != pop_key:
            return
    else:
        refit = None
        if refits:
            from chronocell import contacts_io as IO
            m = (ds.ci >= lo) & (ds.ci < hi) & (ds.cj >= lo) & (ds.cj < hi)
            base = IO.Contacts(ds.ci[m] - lo, ds.cj[m] - lo, np.asarray(ds.cm[m], float), hi - lo, ch.name,
                               int(ch.bin_start(ds.bin0 + lo)), ch.resolution)
            rng = np.random.default_rng(0)

            def refit(r):
                c = IO.Contacts(base.ci, base.cj, rng.poisson(base.cm).astype(float), base.n, base.chrom, base.start, base.resolution)
                return [PL.fit_source("A", c, b0_nm=b0, signal=src_a.signal)[0]] + sources[1:]
        results = []
        bar = st.progress(0.0, text="Applying…")
        for k, (label, kt) in enumerate(cases):
            results.append((label, PL.impact(sources, kt, ch.assembly, refit, 8 if refits else 0)))
            bar.progress((k + 1) / len(cases), text=f"Variant {k + 1} / {len(cases)}")
        bar.empty()
        ss.ve_last = {"pop": pop_key, "results": results, "lo": lo}
    last = ss.ve_last
    results = last["results"]
    ranking = PL.rank_variants(results)
    if len(results) > 1:
        html('<p class="cc-eyebrow">Ranking · "most likely to disrupt gene regulation" (a transparent heuristic, not validated)</p>')
        st.dataframe(ranking, hide_index=True, width="stretch", key="ve_rank")
    pick = st.selectbox("Variant", range(len(results)), format_func=lambda i: results[i][0], key="ve_pick") \
        if len(results) > 1 else 0
    label, r = results[pick]
    s = r["summary"]
    readout([("Derivative chromosomes", f"{len(s['derivatives'])}", "; ".join(s["derivatives"])[:120]),
             ("Pieces lost (acentric)", f"{len(s['lost_pieces'])}", ", ".join(s["lost_pieces"])[:80]),
             ("Genes affected", f"{s['genes_affected']:,}", ""), ("Enhancer–promoter pairs", f"{s['ep_pairs']:,}", ""),
             ("Domain boundary changes", f"{s['boundary_changes']:,}", ""),
             ("Ranking score", f"{s['ranking']['score']}", "heuristic"), ("Confidence", s["confidence"], "")])
    for n in r["notes"]:
        html(f'<p class="cc-note">{esc(n)}</p>')
    t1, t2, t3, t4 = st.tabs(["Genes", "Boundaries & domains", "Enhancer–promoter", "Contact change map"])
    with t1:
        g = pd.DataFrame([{**x, "effects": ", ".join(x["effects"])} for x in r["genes"]])
        st.dataframe(g, hide_index=True, width="stretch", height=260, key="ve_genes")
        from chronocell import annotations as AO
        if AO.clinvar_genes(False) is None:
            if st.button("Add ClinVar counts (download the NCBI gene summary, MD5-checked; information only)", key="ve_clinvar"):
                with st.spinner("Downloading the ClinVar gene summary…"):
                    try:
                        AO.clinvar_genes.cache_clear()
                        AO.clinvar_genes(True)
                        st.rerun()
                    except OSError as exc:
                        warning_card("ClinVar could not be downloaded", str(exc))
        html('<p class="cc-note">ClinVar and the curated lists are information only, never a diagnosis.</p>')
    with t2:
        st.dataframe(pd.DataFrame(r["boundaries"]), hide_index=True, width="stretch", key="ve_bnd")
    with t3:
        if not r["ep"]:
            html('<p class="cc-note">None, or no measured activity track loaded (H3K27ac / ATAC as the enhancer proxy).</p>')
        st.dataframe(pd.DataFrame(r["ep"]), hide_index=True, width="stretch", key="ve_ep")
    with t4:
        imp = r["impact"]
        n = sources[0].n
        fc = imp.log2_fc[:n, :n]
        import plotly.graph_objects as go
        from chronocell import theme as T
        lim = float(np.nanpercentile(np.abs(fc[np.isfinite(fc)]), 99)) if np.isfinite(fc).any() else 1.0
        fig = go.Figure(go.Heatmap(z=fc, zmin=-lim, zmax=lim, colorscale="RdBu_r", colorbar=dict(title="log2 after/before")))
        fig.update_layout(height=420, margin=dict(l=10, r=10, t=10, b=10), yaxis=dict(autorange="reversed"))
        st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="ve_map")
    out = Path(tempfile.mkdtemp(prefix="chronocell_impact_"))
    ranking.to_csv(out / "variant_ranking.csv", index=False)
    pd.DataFrame([{"variant": l, **{k: (", ".join(v) if isinstance(v, list) else v) for k, v in gg.items()}}
                  for l, rr in results for gg in rr["genes"]]).to_csv(out / "genes.csv", index=False)
    pd.DataFrame([{"variant": l, **e} for l, rr in results for e in rr["ep"]]).to_csv(out / "ep_pairs.csv", index=False)
    pd.DataFrame([{"variant": l, **b} for l, rr in results for b in rr["boundaries"]]).to_csv(out / "boundaries.csv", index=False)
    import json
    (out / "summary.json").write_text(json.dumps({"variants": [{"variant": l, **rr["summary"]} for l, rr in results],
                                                  "standing": PL.STANDING["impact"]}, default=float), encoding="utf-8")
    AV.write_run(out, "impact (app)", {"window": f"{ch.name} beads {lo}-{lo + sources[0].n}", "mode": mode,
                                       "zygosity": zyg, "refits": 8 if refits else 0}, ds)
    from chronocell import report_html as RH
    c1, c2, c3 = st.columns(3)
    c1.download_button("Genes (CSV)", (out / "genes.csv").read_bytes(), "impact_genes.csv", width="stretch",
                       icon=":material/download:", key="ve_dl_g")
    c2.download_button("Ranking (CSV)", (out / "variant_ranking.csv").read_bytes(), "impact_ranking.csv", width="stretch",
                       icon=":material/download:", key="ve_dl_r")
    c3.download_button("Report (HTML)", RH.build_html(out), "chronocell_variant_report.html", "text/html", width="stretch",
                       icon=":material/description:", key="ve_dl_rep")
