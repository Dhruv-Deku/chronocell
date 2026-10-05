"""
01 · 3D structure → 06 Analysis suite (Phase B3): loops (HiCCUPS-like), domains (insulation, TopDom-like,
Arrowhead-like), compartments, insulation and P(s) of the loaded contacts in the current window, with exports for
IGV and Juicebox and an HTML report. Nothing is computed until "Run" is pressed. The standing (Gate 6) is read
from its result file.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from ui.common import Dataset, banner, contacts_are_synthetic, esc, html, readout

ss = st.session_state


def contacts_of(ds: Dataset, lo: int, hi: int):
    from chronocell.contacts_io import Contacts
    m = (ds.ci >= lo) & (ds.ci < hi) & (ds.cj >= lo) & (ds.cj < hi)
    ch = ds.chrom
    return Contacts(np.minimum(ds.ci[m], ds.cj[m]) - lo, np.maximum(ds.ci[m], ds.cj[m]) - lo, np.asarray(ds.cm[m], float),
                    hi - lo, ch.name, int(ch.bin_start(ds.bin0 + lo)), int(ch.resolution))


def write_run(out: Path, command: str, params: dict, ds: Dataset | None = None, extra_inputs: list | None = None) -> None:
    """Reproducibility record of an app run: inputs with SHA-256 (as the existing export records them), versions,
    parameters."""
    import datetime as dt
    from chronocell import provenance as PV
    inputs = [{"kind": i[0], "name": i[1], "sha256": i[2]} for i in (ds.inputs if ds is not None else []) if len(i) >= 3]
    rec = {"command": command, "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "parameters": params,
           "inputs": inputs + list(extra_inputs or []), "software": PV.software_versions()}
    (Path(out) / "run.json").write_text(json.dumps(rec, indent=1, default=str), encoding="utf-8")


def standing_banner(gate: str) -> None:
    from chronocell import report_html as RH
    lab, cls = RH.gate_standing(gate)
    banner(f"<b>Standing.</b> {esc(lab)}.", "info" if cls == "ok" else "warn")


def render(ds: Dataset, lo: int, hi: int) -> None:
    html('<p class="cc-note">Calls loops, domains and compartments on the <b>measured contacts</b> of this window '
         '(raw counts; KR or ICE balancing optional). Loops follow the HiCCUPS recipe (local expected from four '
         'neighbourhoods, Poisson test, false-discovery rate); domains come from insulation, a TopDom-like and an '
         'Arrowhead-like caller. Exports open in IGV and Juicebox.</p>')
    if not ds.has_contacts:
        html('<p class="cc-note">Needs contacts (Data → Graph, or a state\'s Hi-C / Micro-C file).</p>')
        return
    if contacts_are_synthetic(ds):
        banner("Input is the <b>SYNTHETIC</b> reference contact map: the calls illustrate the method only.", "warn")
    standing_banner("6")
    c1, c2 = st.columns(2)
    norm = c1.segmented_control("Balancing", ["none", "kr", "ice"], default="none", required=True, key="an_norm",
                                help="none: raw counts (as the validated settings); kr: Knight-Ruiz; ice: iterative correction.")
    loops_on = c2.toggle("Call loops", True, key="an_loops")
    if not st.button("Run analysis suite", key="an_go", type="primary", icon=":material/insights:"):
        if ss.get("an_last", {}).get("key") != f"{ds.key}:{lo}:{hi}":
            return
    else:
        from chronocell import normalize as NZ, pipelines as PL
        c = contacts_of(ds, lo, hi)
        off = c.ci != c.cj
        if norm == "kr":
            c.weights = NZ.kr_balance(c.ci[off], c.cj[off], c.cm[off], c.n).weights
        elif norm == "ice":
            b = NZ.ice_balance(c.ci[off], c.cj[off], c.cm[off], c.n).bias
            c.weights = np.where(np.isfinite(b), 1.0 / b, np.nan)
        with st.spinner("Running the analysis suite…"):
            out = Path(tempfile.mkdtemp(prefix="chronocell_analysis_"))
            r = PL.analyze(c, out, loops=loops_on, orient=np.asarray(ds.gc[lo:hi], float))
            write_run(out, "analyze (app)", {"region": r["summary"]["region"], "resolution": c.resolution,
                                             "normalization": norm, "loops": loops_on}, ds)
        ss.an_last = {"key": f"{ds.key}:{lo}:{hi}", "result": r, "out": str(out)}
    last = ss.get("an_last")
    if not last:
        return
    r, out = last["result"], Path(last["out"])
    s = r["summary"]
    readout([("Loops<small>HiCCUPS-like, FDR 0.1 per neighbourhood</small>", f"{s['loops']:,}", ""),
             ("Boundaries<small>insulation · TopDom-like</small>", f"{s['insulation_boundaries']:,} · {s['topdom_boundaries']:,}", ""),
             ("Domains<small>Arrowhead-like corner score</small>", f"{s['arrowhead_domains']:,}", ""),
             ("A compartment<small>eigenvector > 0, signed by GC content</small>",
              "—" if s["a_fraction"] is None else f"{100 * s['a_fraction']:.0f}", "% of bins"),
             ("P(s) slope", "—" if s["p_s_slope"] is None else f"{s['p_s_slope']:.2f}", "log-log"),
             ("Time", f"{s['seconds']:.1f}", "s")])
    for n in r["notes"]:
        html(f'<p class="cc-note">{esc(n)}</p>')
    t1, t2, t3 = st.tabs(["Loops", "Domains & boundaries", "P(s)"])
    with t1:
        st.dataframe(r["loops"], hide_index=True, width="stretch", height=260, key="an_loops_df")
    with t2:
        st.dataframe(pd.concat([r["boundaries"].assign(kind="boundary"), r["domains"].assign(kind="domain")], ignore_index=True),
                     hide_index=True, width="stretch", height=260, key="an_dom_df")
    with t3:
        ps = r["p_s"].dropna()
        ps = ps[ps["separation_bp"] > 0]
        st.line_chart(pd.DataFrame({"log10 separation (bp)": np.log10(ps["separation_bp"]),
                                    "log10 P(s)": np.log10(ps["p"])}).set_index("log10 separation (bp)"), height=240)
    from chronocell import report_html as RH
    cols = st.columns(5)
    for k, (label, fname) in enumerate((("Loops (BEDPE)", "loops.bedpe"), ("Juicebox 2D", "loops_juicebox_2d.txt"),
                                        ("Boundaries (BED)", "boundaries.bed"), ("Insulation (bedGraph)", "insulation.bedGraph"))):
        p = out / fname
        cols[k].download_button(label, p.read_bytes() if p.exists() else b"", fname, width="stretch",
                                icon=":material/download:", key=f"an_dl_{k}")
    cols[4].download_button("Report (HTML)", RH.build_html(out), "chronocell_analysis_report.html", "text/html",
                            width="stretch", icon=":material/description:", key="an_dl_report")
