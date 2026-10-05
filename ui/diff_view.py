"""
03 · Compare → "Differential analysis" (Phase B2): two conditions, each with one or more replicate maps
(.hic, .mcool, .cool, .pairs: uploaded, or local paths / URLs for large files, which are read region by region),
compared on one region at a chosen resolution. Differential contacts with a replicate-aware moderated t and
false-discovery control; loop gain / loss, boundary changes and compartment switches; tables and BEDPE / CSV
exports and an HTML report. Without at least two replicates per condition no statistics are shown, only fold
changes, and the page says so. The standing (Gate 7) is read from its result file.
"""

from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from chronocell import theme as T
from ui.common import CACHE_DIR, banner, esc, html, readout, warning_card

ss = st.session_state
UPLOADS = CACHE_DIR / "uploads"


def _sources(label: str, key: str) -> list[tuple[str, str]]:
    """(display name, local path or URL) of one condition's replicates."""
    ups = st.file_uploader(f"{label}: replicate files", type=["hic", "mcool", "cool", "pairs", "gz"],
                           accept_multiple_files=True, key=f"{key}_up",
                           help="Each file is one replicate. Files stay on this computer (.chronocell_cache/uploads).")
    paths = st.text_area(f"{label}: or local paths / URLs, one per line", key=f"{key}_paths", height=72,
                         help="For large maps: they are read region by region (a .hic URL by HTTP range requests).")
    out = []
    for f in ups or []:
        data = f.getvalue()
        UPLOADS.mkdir(parents=True, exist_ok=True)
        dest = UPLOADS / f"{hashlib.sha256(data).hexdigest()[:16]}_{Path(f.name).name}"
        if not dest.exists():
            dest.write_bytes(data)
        out.append((f.name, str(dest)))
    for line in (paths or "").splitlines():
        if line.strip():
            out.append((Path(line.strip()).name, line.strip()))
    return out


def render(default_chrom: str = "chr21", default_start: int = 28_000_000, default_end: int = 30_000_000,
           assembly: str | None = None) -> None:
    from ui import analysis_view as AV
    html('<p class="cc-note">Compare two conditions (for example untreated and treated, or healthy and tumour) on one '
         'region. Give each condition its replicate maps. With two or more replicates per condition, differential '
         'contacts get a replicate-aware test with false-discovery control; with one, only fold changes are shown.</p>')
    AV.standing_banner("7")
    ca, cb = st.columns(2)
    with ca:
        A = _sources("Condition A", "df_a")
    with cb:
        B = _sources("Condition B", "df_b")
    c1, c2, c3, c4, c5 = st.columns([1, 1.2, 1.2, 1, 0.8])
    chrom = c1.text_input("Chromosome", default_chrom, key="df_chrom")
    start = c2.number_input("From (bp)", 0, value=int(default_start), step=100_000, key="df_start")
    end = c3.number_input("To (bp)", 0, value=int(default_end), step=100_000, key="df_end")
    res = c4.selectbox("Resolution", [5_000, 10_000, 25_000, 50_000, 100_000], index=1, key="df_res",
                       format_func=lambda r: f"{r // 1000} kb")
    fdr = c5.number_input("FDR", 0.001, 0.5, 0.05, 0.01, key="df_fdr")
    if not st.button("Run differential analysis", key="df_go", type="primary", disabled=not (A and B),
                     icon=":material/compare_arrows:"):
        if not ss.get("df_last"):
            if not (A and B):
                html('<p class="cc-note">Add at least one map to each condition.</p>')
            return
    else:
        from chronocell import contacts_io as IO, pipelines as PL
        if end - start > 2_000 * res:
            warning_card("Region too large", f"At most 2,000 bins per region ({2_000 * res / 1e6:g} Mb at this resolution).")
            return
        try:
            with st.spinner("Reading the maps (region only)…"):
                CA = [IO.read_region(p, chrom, int(start), int(end), int(res), "none", assembly) for _, p in A]
                CB = [IO.read_region(p, chrom, int(start), int(end), int(res), "none", assembly) for _, p in B]
            out = Path(tempfile.mkdtemp(prefix="chronocell_diff_"))
            with st.spinner("Testing every pixel…"):
                r = PL.diff(CA, CB, float(fdr), out=out)
            AV.write_run(out, "diff (app)", {"region": f"{chrom}:{start}-{end}", "resolution": res, "fdr": fdr,
                                             "condition_a": [n for n, _ in A], "condition_b": [n for n, _ in B]},
                         extra_inputs=[{"name": n, "path": p if "://" in p else Path(p).name} for n, p in A + B])
        except (IO.ContactFileError, ValueError, OSError) as exc:
            warning_card("The differential analysis could not run", str(exc))
            return
        ss.df_last = {"result": r, "out": str(out)}
    last = ss.get("df_last")
    if not last:
        return
    r, out = last["result"], Path(last["out"])
    s = r["summary"]
    if not s["statistics"]:
        banner("<b>No statistics.</b> Each condition needs at least two replicates for a test; the tables show fold "
               "changes only.", "warn")
    readout([("Pixels tested", f"{s['tested_pixels']:,}", f"{s['region']} · {s['resolution'] // 1000} kb"),
             ("Significant", f"{s['significant']:,}" if s["statistics"] else "—", f"FDR {s['fdr']:g}"),
             ("Gained · lost in B", f"{s['gained']:,} · {s['lost']:,}" if s["statistics"] else "—", ""),
             ("Replicates", f"{s['replicates_a']} vs {s['replicates_b']}", "")])
    for n in r["notes"]:
        html(f'<p class="cc-note">{esc(n)}</p>')
    t1, t2, t3, t4, t5 = st.tabs(["Map", "Differential contacts", "Loops", "Boundaries", "Compartments"])
    with t1:
        px = r["pixels"]
        n = int(max(px["j"].max(), px["i"].max()) + 1) if len(px) else 0
        if n:
            m = np.full((n, n), np.nan)
            m[px["i"], px["j"]] = px["log2_fc"]
            m[px["j"], px["i"]] = px["log2_fc"]
            import plotly.graph_objects as go
            lim = float(np.nanpercentile(np.abs(px["log2_fc"]), 98)) or 1.0
            fig = go.Figure(go.Heatmap(z=m, zmin=-lim, zmax=lim, colorscale="RdBu_r", colorbar=dict(title="log2 B/A")))
            fig.update_layout(height=420, margin=dict(l=10, r=10, t=10, b=10), yaxis=dict(autorange="reversed"))
            st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="df_map")
            html('<p class="cc-note">log₂ (B / A) of distance-normalised counts per pixel tested; blank: too few reads.</p>')
    with t2:
        tab = r["significant"] if s["statistics"] else r["pixels"].reindex(r["pixels"]["log2_fc"].abs().sort_values(ascending=False).index)
        st.dataframe(tab.head(500), hide_index=True, width="stretch", height=280, key="df_px")
    with t3:
        st.dataframe(r["loops"] if r["loops"] is not None else pd.DataFrame(), hide_index=True, width="stretch", key="df_loops")
    with t4:
        st.dataframe(r["boundaries"] if r["boundaries"] is not None else pd.DataFrame(), hide_index=True, width="stretch",
                     key="df_bnd")
    with t5:
        st.dataframe(r["compartments"] if r["compartments"] is not None else pd.DataFrame(), hide_index=True,
                     width="stretch", key="df_cmp")
    from chronocell import report_html as RH
    cols = st.columns(4)
    for k, (label, fname) in enumerate((("Significant pixels (BEDPE)", "differential_pixels.bedpe"),
                                        ("All pixels (TSV)", "differential_pixels_all.tsv"),
                                        ("Loops (CSV)", "differential_loops.csv"))):
        p = out / fname
        cols[k].download_button(label, p.read_bytes() if p.exists() else b"", fname, width="stretch",
                                icon=":material/download:", key=f"df_dl_{k}")
    cols[3].download_button("Report (HTML)", RH.build_html(out), "chronocell_differential_report.html", "text/html",
                            width="stretch", icon=":material/description:", key="df_dl_report")
