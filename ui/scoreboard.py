"""
07 Scoreboard (October 2026): every pre-registered accuracy test on held-out real data in one place, what it asks in
plain words, what was measured, and the verdict; then charts of the retests and of the newest tests. Everything is read
from validation/results_*.json through validation/report.py (the same code that writes RESULTS.md), so nothing here is
typed in.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from chronocell import theme as T
from ui.common import banner, esc, html, readout

VALIDATION = Path(__file__).resolve().parent.parent / "validation"
PASS, FAIL, OTHER = "#2E7D4F", T.TERRACOTTA, T.MUTED
FIGS: list = []                     # charts drawn on this run, for the downloadable evidence report
AREAS = {"Structure & imaging": ("v3.3", "Gate 1", "Gate 2", "Gate 3", "Gate 5"),
         "Analysis & perturbations": ("Gate 4", "Gate 6", "Gate 7"),
         "Drug lab": ("Gate 8",),
         "Quantum lab": ("Gate Q",)}
PLAIN = {   # what each test asks, in plain words (the measured numbers come from the result files)
    "v3.3": "How much of the real, reproducible folding pattern does the model recover?",
    "Gate 1:": "Does modelling the whole chromosome at once match the windowed model?",
    "Gate 1b": "Do distances predicted from sequencing data rank like real microscope distances?",
    "Gate 1c": "Can the predicted sizes be calibrated to real nanometres?",
    "Gate 2:": "Do the stated 90 % ranges really contain 90 % of real distances?",
    "Gate 2b": "Same, when the input is sequencing Hi-C.",
    "Gate 2c": "Can the app tell which predicted distances are wrong?",
    "Gate 2d": "Do recalibrated ranges hold in every distance band?",
    "Gate 2e": "A second attempt at telling which distances are wrong.",
    "Gate 3:": "Is ChronoCell as good as PASTIS, a published method?",
    "Gate 3b": "Would a learned correction help?",
    "Gate 4:": "Does it predict what happens when cohesin is removed?",
    "Gate 4b": "Does it predict the effect of DNA deletions?",
    "Gate 4c": "Cohesin loss again, on six new regions.",
    "Gate 4d": "Rearrangements with before-and-after data (needs 3 events).",
    "Gate 5:": "Can it predict the fold from DNA sequence alone?",
    "Gate 5b": "Same, with more inputs.",
    "Gate 5m": "Does the human predictor work on mouse cells?",
    "Gate 6:": "Does it find the same DNA loops as the reference caller?",
    "Gate 6b": "Loops again: settings from three cell types, tested on two new ones.",
    "Gate 7": "Does the change-finder avoid false alarms?",
    "Gate 8": "Does the drug simulator match real drug-treated cells?",
    "Gate Q1": "Does the quantum optimiser find the best answer?",
    "Gate Q2:": "Does the quantum route find DNA 'rooms' as well as classical tools?",
    "Gate Q2b": "Rooms again, on two new cell types.",
    "Gate Q3": "Does the quantum chemistry get small molecules exactly right?",
    "Gate Q4:": "Does a quantum classifier predict which genes are on?",
    "Gate Q4b": "Gene classifier again, on a new cell type.",
    "Gate Q5": "Can a quantum classifier flag heart-risk (hERG) drugs?",
    "Gate Q6:": "Quantum chemistry on stretched molecules.",
    "Gate Q6b": "Stretched molecules with a smarter circuit (ADAPT-VQE).",
    "Gate Q6c": "Same, scored against the right spin state.",
    "Gate Q6d": "Same, with an escape from flat spots and several starting points.",
    "Gate Q7:": "Can quantum matching dock a drug into its protein?",
    "Gate Q7b": "Docking with a better score, on a new set.",
    "Gate Q7c": "Docking with a quantum-seeded hybrid search, on new complexes.",
    "Gate Q8": "A 21-property drug safety profile (ADMET) from a quantum model.",
    "Gate Q9": "Can error mitigation rescue a noisy quantum computer?",
}


def _report():
    spec = importlib.util.spec_from_file_location("chronocell_report_ui", VALIDATION / "report.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load(name: str):
    try:
        return json.loads((VALIDATION / name).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


RETESTS = ("Gate 6b", "Gate Q2b", "Gate Q4b", "Gate Q6b", "Gate Q6c", "Gate Q6d", "Gate Q7b", "Gate Q7c")


def status_of(verdict: str) -> str:
    v = verdict.lower()
    if "blocked" in v or "not run" in v or "baseline" in v or "see results" in v or "not validated" in v:
        return "other"
    if "matches within" in v or ("transfer" in v and "do not" in v):
        return "other"                              # mixed: part holds, part does not (see the verdict text)
    if "fail" in v or "too narrow" in v or "poorly" in v or "below" in v or "no usable" in v:
        return "fail" if "pass" not in v.split("fail")[0] else "pass"
    if "pass" in v or "beats" in v or "matches" in v or "transfer" in v or "usable" in v:
        return "pass"
    return "other"


def rows() -> pd.DataFrame:
    rep = _report()
    md = [rep.readme_accuracy(), rep.summary_q(), rep.summary_qd(), rep.summary_r2()]
    if hasattr(rep, "summary_new"):
        md.append(rep.summary_new())
    out, seen = [], set()
    for block in md:
        for ln in block.splitlines():
            if not ln.startswith("| ") or ln.startswith("| Test") or set(ln) <= set("|-"):
                continue
            cells = [c.strip() for c in ln.strip("|").split(" | ")]
            if len(cells) < 3 or cells[0] in seen or cells[1] == "not run":
                continue
            seen.add(cells[0])
            test, measured, verdict = cells[0], " | ".join(cells[1:-1]), cells[-1]
            area = next((a for a, keys in AREAS.items() if any(test.startswith(k) for k in keys)), "Quantum lab")
            key = next((k for k in sorted(PLAIN, key=len, reverse=True) if test.startswith(k)), None)
            out.append({"area": area, "test": test, "question": PLAIN.get(key, ""), "measured": measured, "verdict": verdict,
                        "status": status_of(verdict), "retest": test.startswith(RETESTS)})
    return pd.DataFrame(out)


def _chip(status: str) -> str:
    col = {"pass": PASS, "fail": FAIL}.get(status, OTHER)
    lab = {"pass": "PASS", "fail": "FAIL"}.get(status, "OTHER")
    return (f'<span style="display:inline-block;padding:1px 9px;border-radius:10px;font-size:11px;font-weight:600;'
            f'letter-spacing:.04em;color:#fff;background:{col}">{lab}</span>')


def overview(df: pd.DataFrame) -> None:
    n = len(df)
    p, f, o = (int((df["status"] == s).sum()) for s in ("pass", "fail", "other"))
    readout([("Tests on held-out real data", f"{n}", "each pre-registered, each run once"),
             ("Passed", f"{p}", f"{100 * p / max(n, 1):.0f} %"), ("Failed", f"{f}", "kept on the record"),
             ("Blocked / baseline / mixed", f"{o}", "see the verdict text")])
    fig = go.Figure()
    for s, col in (("pass", PASS), ("fail", FAIL), ("other", OTHER)):
        cnt = [int(((df["area"] == a) & (df["status"] == s)).sum()) for a in AREAS]
        fig.add_trace(go.Bar(y=list(AREAS), x=cnt, orientation="h", name=s, marker_color=col, text=cnt,
                             textposition="inside"))
    fig.update_layout(**T.plot_layout(230, barmode="stack", showlegend=True, legend=dict(orientation="h", y=1.18),
                                      margin=dict(l=170, r=10, t=30, b=30), xaxis=dict(title="tests")))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="sb_overview")


def table(df: pd.DataFrame) -> None:
    c1, c2 = st.columns([1.2, 1])
    area = c1.selectbox("Area", ["All"] + list(AREAS), key="sb_area")
    show = c2.segmented_control("Show", ["All", "Passed", "Failed", "Retests"], default="All", required=True, key="sb_show")
    d = df if area == "All" else df[df["area"] == area]
    if show == "Passed":
        d = d[d["status"] == "pass"]
    elif show == "Failed":
        d = d[d["status"] == "fail"]
    elif show == "Retests":
        d = d[d["retest"]]
    for _, r in d.iterrows():
        html(f'<div class="cc-callout" style="margin:6px 0"><div style="display:flex;gap:10px;align-items:center">'
             f'{_chip(r["status"])}<b>{esc(r["test"])}</b></div>'
             + (f'<div style="margin-top:4px"><i>{esc(r["question"])}</i></div>' if r["question"] else "")
             + f'<div style="margin-top:4px"><b>Measured:</b> {esc(r["measured"])}</div>'
             f'<div><b>Verdict:</b> {esc(r["verdict"])}</div></div>')


# ------------------------------------------------------------------------------------------- charts
def retest_pairs() -> list[dict]:
    """Each original test and its retests as points 'distance from the pass line' (right of 0 = passed), read from the
    result files. Units differ by row (F1, AUC, share of molecules, share of ligands); only the side of 0 compares."""
    out = []
    g = _load("results_gateq.json") or {}
    r2 = _load("results_round2.json") or {}

    def q2_margin(rows):
        vals = []
        for rr in rows.values():
            rr = rr.get("rows", rr)
            best = max(rr["insulation (classical)"]["f1"], rr["topdom (classical)"]["f1"])
            vals.append(rr["QAOA (simulated quantum)"]["f1"] - (best - 0.05))
        return min(vals)
    if "q2" in g:
        pts = [("Q2", q2_margin(g["q2"]["rows"]))] + ([("Q2b", q2_margin(r2["q2b"]["per_cell"]))] if "q2b" in r2 else [])
        out.append({"family": "Domains: F1 vs best classical − 0.05", "points": pts})
    if "q4" in g:
        pts = [("Q4", g["q4"]["auc"]["qsvm"] - g["q4"]["auc"]["rbf"] + 0.03)]
        if "q4b" in r2:
            pts.append(("Q4b", r2["q4b"]["auc"]["qsvm"] - r2["q4b"]["auc"]["rbf"] + 0.03))
        out.append({"family": "Genes: AUC vs RBF − 0.03", "points": pts})
    g6, g6b = _load("results_gate6.json"), _load("results_gate6b.json")
    if g6:
        m = lambda r: min(v["chronocell"]["f1"] - max(v["chromosight"].get("f1", 0), v["mustache"].get("f1", 0))   # noqa: E731
                          for v in r["sets"].values())
        out.append({"family": "Loops: F1 vs best tool", "points": [("6", m(g6))] + ([("6b", m(g6b))] if g6b else [])})
    qd = _load("results_qdrug.json") or {}
    if "q6" in qd:
        share = lambda errs: float(np.mean([abs(e) <= 1.6 for e in errs])) - 1.0                      # noqa: E731
        pts = [("Q6", share([x["error_mEh"] for x in qd["q6"]["vqe"]]))]
        for k, f in (("q6b", "adapt_error_mEh"), ("q6c", "adapt_error_vs_singlet_mEh"), ("q6d", "error_vs_singlet_mEh")):
            if k in r2:
                pts.append((k[:2].upper() + k[2:], share([x[f] for x in r2[k]["rows"]])))
        out.append({"family": "Molecules: share within chemical accuracy − 1", "points": pts})
    if "q7" in qd:
        pts = [("Q7", qd["q7"]["success_qaoa"] - qd["q7"]["success_random_search"])]
        if "q7b" in r2:
            pts.append(("Q7b", r2["q7b"]["success_qaoa"] - r2["q7b"]["success_random"]))
        out.append({"family": "Docking: share docked vs random search", "points": pts})
    return out


def chart_retests() -> None:
    rows_ = retest_pairs()
    if not rows_:
        return
    fig = go.Figure()
    n = len(rows_)
    for k, r in enumerate(rows_):
        y = n - k
        xs = [v for _, v in r["points"]]
        fig.add_trace(go.Scatter(x=xs, y=[y] * len(xs), mode="lines", line=dict(color=T.RULE_STRONG, width=3),
                                 showlegend=False, hoverinfo="skip"))
        for i, (lab, x) in enumerate(r["points"]):
            fig.add_trace(go.Scatter(x=[x], y=[y], mode="markers+text", text=[lab],
                                     textposition="top center" if i % 2 == 0 else "bottom center",
                                     marker=dict(size=14, color=PASS if x >= 0 else FAIL, line=dict(color=T.INK, width=0.6)),
                                     showlegend=False, hovertemplate=f"{lab}: %{{x:+.3f}}<extra></extra>"))
    fig.add_vline(x=0, line=dict(color=T.INK, dash="dot"), annotation_text="pass line")
    fig.update_layout(**T.plot_layout(110 + 80 * n, margin=dict(l=290, r=20, t=30, b=40),
                                      yaxis=dict(tickvals=[n - k for k in range(n)], ticktext=[r["family"] for r in rows_],
                                                 range=[0.4, n + 0.6]),
                                      xaxis=dict(title="distance from the pass line (right of 0 = passed)")))
    FIGS.append(fig)
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="sb_retests")
    html('<p class="cc-note">Each row: the original test and its retests (each with a new method on new data), in '
         'order; green passed, red failed. Units differ between rows: only which side of the line counts. Original '
         'verdicts always stay on the record.</p>')


def chart_molecules() -> None:
    sets = []
    q = _load("results_qdrug.json")
    if q and "q6" in q:
        sets.append(("Q6: fixed UCCSD", [abs(x["error_mEh"]) for x in q["q6"]["vqe"]]))
    r2 = _load("results_round2.json") or {}
    if "q6b" in r2:
        sets.append(("Q6b: ADAPT-VQE", [abs(x["adapt_error_mEh"]) for x in r2["q6b"]["rows"]]))
    if "q6c" in r2:
        sets.append(("Q6c: ADAPT vs singlet", [abs(x["adapt_error_vs_singlet_mEh"]) for x in r2["q6c"]["rows"]]))
    if "q6d" in r2:
        sets.append(("Q6d: ADAPT + escape", [abs(x["error_vs_singlet_mEh"]) for x in r2["q6d"]["rows"]]))
    if not sets:
        return
    fig = go.Figure()
    for k, (lab, errs) in enumerate(sets):
        e = np.maximum(np.array(errs, float), 1e-4)
        fig.add_trace(go.Box(y=e, name=lab, boxpoints="all", jitter=0.4, pointpos=0, marker=dict(size=6),
                             line=dict(color=T.ACCENT if k else T.MUTED)))
    fig.add_hline(y=1.6, line=dict(color=T.TERRACOTTA, dash="dot"), annotation_text="chemical accuracy 1.6 mHa")
    fig.update_layout(**T.plot_layout(330, showlegend=False, yaxis=dict(type="log", title="error vs exact (mHa)")))
    FIGS.append(fig)
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="sb_mol")
    html('<p class="cc-note">Each dot is one held-out molecule / bond length. Below the dotted line = chemically '
         'accurate. Each round used new molecules.</p>')


def chart_loops() -> None:
    rows_ = []
    for fname, lab in (("results_gate6.json", "Gate 6"), ("results_gate6b.json", "Gate 6b")):
        r = _load(fname)
        for cell, v in (r or {}).get("sets", {}).items():
            for tool in ("chronocell", "chromosight", "mustache"):
                if "f1" in v.get(tool, {}):
                    rows_.append({"test": f"{lab} · {cell}", "tool": tool, "f1": v[tool]["f1"]})
    if not rows_:
        return
    d = pd.DataFrame(rows_)
    fig = go.Figure()
    for tool, col in (("chronocell", T.ACCENT), ("chromosight", T.OCHRE), ("mustache", T.MUTED)):
        s = d[d["tool"] == tool]
        fig.add_trace(go.Bar(x=s["test"], y=s["f1"], name=tool, marker_color=col))
    fig.update_layout(**T.plot_layout(300, barmode="group", showlegend=True, legend=dict(orientation="h", y=1.15),
                                      yaxis=dict(title="F1 vs ENCODE reference loops", range=[0, 1])))
    FIGS.append(fig)
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="sb_loops")


def chart_admet() -> None:
    r = _load("results_admet.json")
    if not r:
        html('<p class="cc-note">Gate Q8 (ADMET profile) has not been run yet.</p>')
        return
    from chronocell.quantum import admet as A
    e = r["endpoints"]
    names = [n for g in A.GROUPS.values() for n in g if n in e]
    lab = [A.ENDPOINTS[n][1] for n in names]
    fig = go.Figure()
    for k, col in (("quantum", T.ACCENT), ("rbf8", T.OCHRE), ("rbf17", T.MUTED)):
        fig.add_trace(go.Bar(x=lab, y=[e[n]["metric"][k] for n in names], name={"quantum": "quantum kernel",
                             "rbf8": "classical RBF, same 8 inputs", "rbf17": "classical RBF, all 17"}[k], marker_color=col))
    fig.update_layout(**T.plot_layout(380, barmode="group", showlegend=True, legend=dict(orientation="h", y=1.12),
                                      margin=dict(l=50, r=10, t=30, b=140), xaxis=dict(tickangle=-40),
                                      yaxis=dict(title="test AUC (yes/no) or Spearman ρ (values)")))
    FIGS.append(fig)
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="sb_admet")
    html(f'<p class="cc-note">Official held-out scaffold splits of the TDC ADMET benchmark: {r["passed"]} of {r["of"]} '
         'endpoints met the rule (the quantum model within 0.03 of the classical model on the same inputs, its interval '
         'above chance).</p>')


def chart_mitigation() -> None:
    r = _load("results_mitigation.json")
    if not r:
        html('<p class="cc-note">Gate Q9 (error mitigation) has not been run yet.</p>')
        return
    k = r["method_key"]
    lab = [f"{x['molecule']} {x['bond_scale']}×" for x in r["rows"]]
    noisy = [max(abs(x["methods"][k]["noisy_mEh"]), 1e-4) for x in r["rows"]]
    mit = [max(abs(x["methods"][k]["mitigated_mEh"]), 1e-4) for x in r["rows"]]
    fig = go.Figure([go.Bar(x=lab, y=noisy, name="noisy (symmetry-verified)", marker_color=T.TERRACOTTA),
                     go.Bar(x=lab, y=mit, name="after zero-noise extrapolation", marker_color=T.ACCENT)])
    fig.add_hline(y=1.6, line=dict(color=T.INK, dash="dot"), annotation_text="chemical accuracy")
    fig.update_layout(**T.plot_layout(340, barmode="group", showlegend=True, legend=dict(orientation="h", y=1.12),
                                      yaxis=dict(type="log", title="energy error (mHa)"), xaxis=dict(tickangle=-40)))
    FIGS.append(fig)
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="sb_mit")


def chart_docking() -> None:
    pts = []
    q = _load("results_qdrug.json")
    if q and "q7" in q:
        pts.append(("Q7 (PoseBusters)", q["q7"]["success_qaoa"], q["q7"]["success_random_search"]))
    r2 = _load("results_round2.json") or {}
    for k, lab in (("q7b", "Q7b (Astex)"), ("q7c", "Q7c (PoseBusters new)")):
        if k in r2:
            pts.append((lab, r2[k]["success_qaoa"], r2[k]["success_random"]))
    if not pts:
        return
    fig = go.Figure([go.Bar(x=[p[0] for p in pts], y=[100 * p[1] for p in pts], name="quantum route", marker_color=T.ACCENT),
                     go.Bar(x=[p[0] for p in pts], y=[100 * p[2] for p in pts], name="random search, same score",
                            marker_color=T.MUTED)])
    fig.update_layout(**T.plot_layout(300, barmode="group", showlegend=True, legend=dict(orientation="h", y=1.15),
                                      yaxis=dict(title="ligands docked within 2 Å (%)")))
    FIGS.append(fig)
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="sb_dock")


def evidence_html(df: pd.DataFrame, figures: list) -> str:
    """A self-contained summary page (tables inline; charts load plotly.js from its CDN) for sharing the record."""
    import datetime as dt
    import html as H
    n = len(df)
    p, f, o = (int((df["status"] == s).sum()) for s in ("pass", "fail", "other"))
    css = ("body{font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;"
           "color:#1c1e1b}h1{font-size:22px}h2{font-size:16px;border-bottom:1px solid #ddd;padding-bottom:4px;margin-top:28px}"
           "table{border-collapse:collapse;width:100%;font-size:12.5px}td,th{border:1px solid #ddd;padding:4px 6px;"
           "vertical-align:top;text-align:left}th{background:#f4f4f1}.pass{color:#2e7d4f;font-weight:700}"
           ".fail{color:#c24a1e;font-weight:700}.other{color:#62645f;font-weight:700}")
    rows = "".join(f"<tr><td>{H.escape(r['area'])}</td><td><b>{H.escape(r['test'])}</b><br><i>{H.escape(r['question'])}</i></td>"
                   f"<td>{H.escape(r['measured'])}</td><td class='{r['status']}'>{r['status'].upper()}</td>"
                   f"<td>{H.escape(r['verdict'])}</td></tr>" for _, r in df.iterrows())
    figs = "".join(fig.to_html(full_html=False, include_plotlyjs="cdn" if i == 0 else False) for i, fig in enumerate(figures))
    return (f"<!doctype html><html><head><meta charset='utf-8'><title>ChronoCell-5D evidence</title><style>{css}</style>"
            f"</head><body><h1>ChronoCell-5D: accuracy on held-out real data</h1><p>Generated "
            f"{dt.datetime.now().strftime('%Y-%m-%d %H:%M')} from validation/results_*.json. Every test was written down "
            f"before its data were read and run once; failures stay on the record.</p>"
            f"<p><b>{n}</b> tests · <span class='pass'>{p} passed</span> · <span class='fail'>{f} failed</span> · "
            f"<span class='other'>{o} blocked / baseline / mixed</span></p><h2>All tests</h2><table><tr><th>Area</th>"
            f"<th>Test</th><th>Measured</th><th>Verdict</th><th>Verdict text</th></tr>{rows}</table>"
            f"<h2>Charts</h2>{figs}</body></html>")


def render() -> None:
    html('<p class="cc-note">Every accuracy claim in ChronoCell-5D is a test written down <b>before</b> its data were '
         'read, run <b>once</b> on held-out real data, and kept on the record whether it passed or failed. A failed '
         'test is never re-run; a retest (b, c, d) uses a new method and new data. This page reads the result files '
         'directly.</p>')
    try:
        df = rows()
    except Exception as exc:                      # a missing or half-written result file must not break the page
        banner(f"The result files could not be read: {esc(exc)}", "warn")
        return
    overview(df)
    FIGS.clear()
    t0, t1, t2, t3, t4, t5, t6 = st.tabs(["Before → after", "All tests", "Molecules (quantum chemistry)", "DNA loops",
                                          "Drug safety (ADMET)", "Error mitigation", "Docking"])
    with t0:
        chart_retests()
    with t1:
        table(df)
    with t2:
        chart_molecules()
    with t3:
        chart_loops()
    with t4:
        chart_admet()
    with t5:
        chart_mitigation()
    with t6:
        chart_docking()
    st.download_button("Download the evidence report (HTML)", evidence_html(df, list(FIGS)),
                       file_name="chronocell_evidence.html", mime="text/html", key="sb_report",
                       icon=":material/download:",
                       help="Every test with its question, measurement and verdict, and the charts above, in one page.")
