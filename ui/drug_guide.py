"""
04 Drug lab → Drug guide (October 2026): what each drug class is, where it acts on the fold in view, how every class
does at full dose, which pairs combine well, and molecule cards from PubChem. Everything is computed only when its
button is pressed; the information is general and educational, not medical advice.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from chronocell import drug_info as DI, therapy as TH, theme as T
from ui.common import banner, esc, html, readout, warning_card

ss = st.session_state
TARGET_PLAIN = {"compact_low_signal": "packed, quiet chromatin", "low_signal": "chromatin with few active marks",
                "hubs": "hyper-active enhancer hubs", "loops": "loop anchors", "compact": "packed chromatin (any activity)",
                "mid_signal": "moderately active stretches (poised enhancers)", "high_signal": "the most active stretches"}
DIRECTION = {1: "opens", -1: "compacts", 0: "re-draws loops"}
STATUS_SCORE = {"approved": 3, "clinical": 2, "research": 1, "hypothetical": 0}


def _status_level(evidence: str) -> int:
    e = evidence.lower()
    return 3 if "approved" in e else 2 if "clinical" in e else 1 if "research" in e else 0


def _card(d: TH.Drug, info: DI.ClassInfo | None, core: bool) -> str:
    ex = "; ".join(f"{n} ({s})" for n, s in info.examples) if info else esc(d.examples)
    return (f'<div class="cc-callout"><h4>{esc(d.name)}{"" if core else " · <small>extended set</small>"}</h4>'
            f'<b>Acts on:</b> {esc(info.target_protein) if info else "—"}. {esc(info.what_it_does) if info else ""}<br>'
            f'<b>In the simulator:</b> {DIRECTION[d.direction]} {esc(TARGET_PLAIN.get(d.target, d.target))}'
            f'{"" if d.strength == 1.0 else f" (strength {d.strength:g} × the core classes)"}.<br>'
            f'<b>Status (2025):</b> {esc(info.status) if info else esc(d.evidence)}<br>'
            f'<b>Examples:</b> {esc(ex)}<br><b>Safety themes:</b> {esc(info.safety) if info else "—"}</div>')


def classes_tab() -> None:
    html('<p class="cc-note">Twelve classes: the four core ones (the default set) and eight more (choose <b>Extended</b> '
         'above). Each is reduced to <i>where</i> it acts and <i>which way</i> it pushes; the cards say what is known about '
         'the real drugs. General information, not medical advice; check current labels.</p>')
    rows = []
    for k, d in TH.ALL_DRUGS.items():
        rows.append({"class": d.name, "key": k, "target": TARGET_PLAIN.get(d.target, d.target), "direction": d.direction,
                     "strength": d.strength, "status": _status_level(d.evidence), "core": k in TH.DRUGS})
    df = pd.DataFrame(rows)
    targets = list(dict.fromkeys(df["target"]))
    fig = go.Figure()
    for core, sym in ((True, "circle"), (False, "diamond")):
        sub = df[df["core"] == core]
        fig.add_trace(go.Scatter(
            x=[targets.index(t) + (0.12 if not core else -0.12) for t in sub["target"]], y=sub["direction"] * sub["strength"],
            mode="markers+text", text=sub["class"], textposition="top center", textfont=dict(size=9),
            marker=dict(size=10 + 6 * sub["status"], symbol=sym, color=[T.ACCENT if c else T.TERRACOTTA for c in sub["core"]],
                        line=dict(color=T.INK, width=0.6)),
            hovertext=[f"{r['class']}: {DIRECTION[r['direction']]} {r['target']}" for _, r in sub.iterrows()], hoverinfo="text",
            name="core" if core else "extended"))
    fig.add_hline(y=0, line=dict(color=T.RULE_STRONG, dash="dot"))
    fig.update_layout(**T.plot_layout(380, showlegend=True, legend=dict(orientation="h", y=1.12),
                                      xaxis=dict(tickvals=list(range(len(targets))), ticktext=targets, tickangle=-20),
                                      yaxis=dict(title="opens (+) · compacts (−) × strength")))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="dg_landscape")
    html('<p class="cc-note">Drug landscape: where each class acts (x) and which way and how hard it pushes (y). Bigger '
         'marks = further along (approved > clinical trials > research). Blue circles: core set; red diamonds: extended set.</p>')
    cols = st.columns(2)
    for i, (k, d) in enumerate(TH.ALL_DRUGS.items()):
        with cols[i % 2]:
            html(_card(d, DI.CLASSES.get(k), k in TH.DRUGS))


def where_tab(x: np.ndarray, sig: np.ndarray, valid: np.ndarray, b0: float, mb: np.ndarray) -> None:
    html('<p class="cc-note">How strongly each class reaches each stretch of the region being treated (0 = not at all, '
         '1 = fully), from the activity signal and how crowded each bead is.</p>')
    if not st.button("Map where every class acts", key="dg_where_go"):
        if "dg_where" not in ss or ss.dg_where["n"] != len(x):
            return
    else:
        W = np.stack([TH.target_weights(d, x, sig, valid, b0) for d in TH.ALL_DRUGS.values()])
        ss.dg_where = {"W": W, "n": len(x)}
    W = ss.dg_where["W"]
    names = [d.name for d in TH.ALL_DRUGS.values()]
    fig = go.Figure(go.Heatmap(z=W, x=mb, y=names, colorscale=[[0, T.PAPER_RAISED], [0.5, "#9AA2E6"], [1, "#1A2175"]],
                               zmin=0, zmax=1, colorbar=dict(title="reach", thickness=10),
                               hovertemplate="%{y}<br>%{x:.2f} Mb · %{z:.2f}<extra></extra>"))
    fig.update_layout(**T.plot_layout(420, margin=dict(l=210, r=10, t=10, b=40), xaxis=dict(title="Mb")))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="dg_where_fig")
    sig_n = (sig - np.nanmin(sig)) / ((np.nanmax(sig) - np.nanmin(sig)) or 1)
    fig2 = go.Figure(go.Scatter(x=mb, y=sig_n, mode="lines", line=dict(color=T.TERRACOTTA, width=1.2), fill="tozeroy"))
    fig2.update_layout(**T.plot_layout(110, margin=dict(l=210, r=10, t=4, b=24), yaxis=dict(title="signal", showticklabels=False)))
    st.plotly_chart(fig2, theme=None, width="stretch", config=T.PLOT_CONFIG, key="dg_where_sig")


def all_classes_tab(x, sig, valid, b0, xh, eff, sig_ref) -> None:
    if xh is None:
        html('<p class="cc-note">Needs a healthy baseline covering the same DNA (as the ranking above).</p>')
        return
    html('<p class="cc-note">Every class, core and extended, at full dose: how far each moves this fold back toward healthy.</p>')
    if st.button("Test all twelve classes", key="dg_all_go"):
        with st.spinner("Simulating twelve classes at full dose…"):
            ss.dg_all = TH.compare_drugs(x, sig, valid, b0, xh, eff, signal_ref=sig_ref, drugs=TH.ALL_DRUGS)
    df = ss.get("dg_all")
    if df is None:
        return
    colors = [T.ACCENT if k in TH.DRUGS else T.TERRACOTTA for k in df["key"]]
    fig = go.Figure(go.Bar(x=df["restoration_pct"], y=df["drug"], orientation="h", marker_color=colors,
                           hovertext=df["evidence"], hoverinfo="text+x"))
    fig.update_layout(**T.plot_layout(380, margin=dict(l=230, r=10, t=10, b=40), xaxis=dict(title="restored at full dose (%)"),
                                      yaxis=dict(autorange="reversed")))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="dg_all_fig")
    html('<p class="cc-note">Blue: core classes; red: extended. Mechanism simulation, not a prediction for patients.</p>')


def pairs_tab(x, sig, valid, b0, xh, eff, sig_ref) -> None:
    if xh is None:
        html('<p class="cc-note">Needs a healthy baseline covering the same DNA.</p>')
        return
    from chronocell.quantum import problems as PR
    keys = list(TH.DRUGS)
    html('<p class="cc-note">Pairs of core classes given one after the other at full dose: does the pair restore more than '
         'either drug alone (synergy, green) or less (red)?</p>')
    if st.button("Simulate every pair of core classes", key="dg_pairs_go"):
        bar = st.progress(0.0)
        single = {}
        for k in keys:
            single[k] = PR.combo_restoration(x, sig, valid, b0, xh, {k: 1.0}, eff, 10, sig_ref)
        M = np.full((len(keys), len(keys)), np.nan)
        n = 0
        for i, a in enumerate(keys):
            M[i, i] = single[a]
            for j in range(i + 1, len(keys)):
                M[i, j] = M[j, i] = PR.combo_restoration(x, sig, valid, b0, xh, {a: 1.0, keys[j]: 1.0}, eff, 10, sig_ref)
                n += 1
                bar.progress(n / 6)
        bar.empty()
        ss.dg_pairs = {"M": M, "single": single}
    pr = ss.get("dg_pairs")
    if pr is None:
        return
    M = pr["M"]
    best = np.array([[max(pr["single"][a], pr["single"][b]) for b in keys] for a in keys])
    gain = M - best
    names = [TH.DRUGS[k].name for k in keys]
    fig = go.Figure(go.Heatmap(z=gain, x=names, y=names, zmid=0, colorscale="RdYlGn", text=np.round(M, 0),
                               texttemplate="%{text}%", hovertemplate="%{y} + %{x}: %{text}% (vs best single %{z:+.0f})<extra></extra>"))
    fig.update_layout(**T.plot_layout(360, margin=dict(l=170, r=10, t=10, b=120), xaxis=dict(tickangle=-30)))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key="dg_pairs_fig")


def molecule_tab() -> None:
    html('<p class="cc-note">Look up a compound on <b>PubChem</b> (only the name you pick is sent; the answer is cached on this '
         'computer): its computed properties, the drug-likeness rules, a 2D drawing and a 3D model.</p>')
    names = sorted({n for c in DI.CLASSES.values() for n in c.lookup})
    c1, c2 = st.columns([1.2, 1])
    name = c1.selectbox("Compound", names, index=names.index("vorinostat") if "vorinostat" in names else 0, key="dg_mol_name")
    other = c2.text_input("…or any compound name", key="dg_mol_other", placeholder="e.g. tucidinostat")
    if st.button("Fetch from PubChem", key="dg_mol_go", icon=":material/download:"):
        q = other.strip() or name
        try:
            with st.spinner(f"Looking up {q}…"):
                p = DI.properties(q)
                try:
                    img = DI.image_png(int(p["CID"]))
                except DI.PubChemError:
                    img = None                    # the drawing is optional: a slow image service must not hide the card
                try:
                    conf = DI.conformer_3d(int(p["CID"]))
                except DI.PubChemError:
                    conf = None
            ss.dg_mol = {"name": q, "p": p, "img": img, "conf": conf}
        except DI.PubChemError as exc:
            warning_card("PubChem lookup failed", str(exc))
            return
    m = ss.get("dg_mol")
    if not m:
        return
    p = m["p"]
    dl = DI.drug_likeness(p)
    readout([("Formula", esc(p.get("MolecularFormula", "—")), f"PubChem CID {p.get('CID')}"),
             ("Molecular weight", f"{float(p.get('MolecularWeight', 0)):.1f}", "g/mol"),
             ("XLogP (lipophilicity)", f"{p.get('XLogP', '—')}", ""), ("Polar surface (TPSA)", f"{p.get('TPSA', '—')}", "Å²"),
             ("H-bond donors / acceptors", f"{p.get('HBondDonorCount')} / {p.get('HBondAcceptorCount')}", ""),
             ("Rule-of-five violations", f"{dl['ro5_violations']}", "0-1 typical of oral drugs")])
    c1, c2 = st.columns([1, 1.3])
    with c1:
        if m["img"] is not None:
            st.image(m["img"], caption=f"{m['name']} (PubChem 2D depiction)", width=300)
        else:
            html('<p class="cc-note">PubChem did not return the 2D drawing in time; the properties above are complete.</p>')
        st.dataframe(pd.DataFrame([{"rule": k, "met": "yes" if v else "no"} for k, v in dl["rules"].items()]),
                     hide_index=True, width="stretch", key="dg_mol_rules")
    with c2:
        if m["conf"] is not None:
            from ui.quantum_lab import molecule_figure
            el, xyz, bonds = m["conf"]
            st.plotly_chart(molecule_figure(el, xyz, bonds, height=380, title="PubChem 3D conformer"), theme=None,
                            width="stretch", config=T.PLOT_CONFIG, key="dg_mol_3d")
        else:
            html('<p class="cc-note">PubChem has no 3D conformer for this compound (common for large or flexible molecules).</p>')
    html(f'<p class="cc-note">IUPAC name: {esc(str(p.get("IUPACName", "—"))[:300])}</p>')


def gate8_banner() -> None:
    import json
    from pathlib import Path
    p = Path(__file__).resolve().parent.parent / "validation" / "results_gate8.json"
    try:
        r = json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
    except (OSError, ValueError):
        r = None
    if r is None:
        banner("<b>Measured standing.</b> Gate 8 (the simulator against chromatin tracing after real drug treatment) has "
               "not been run.", "info")
        return
    banner(f"<b>Measured standing (Gate 8).</b> Against chromatin tracing of IMR-90 cells after real drug treatment, "
           f"{r['passing_drugs']} of {len(r['drugs'])} drugs met the pre-registered rule: <b>{r['verdict']}</b>. "
           "Details in validation/RESULTS.md.", "info" if r["verdict"] == "pass" else "warn")


def render(x, sig, valid, b0, xh, eff, sig_ref, mb) -> None:
    gate8_banner()
    t1, t2, t3, t4, t5 = st.tabs(["Drug classes", "Where they act", "All classes at full dose", "Pairs", "Molecule (PubChem)"])
    with t1:
        classes_tab()
    with t2:
        where_tab(x, sig, valid, b0, mb)
    with t3:
        all_classes_tab(x, sig, valid, b0, xh, eff, sig_ref)
    with t4:
        pairs_tab(x, sig, valid, b0, xh, eff, sig_ref)
    with t5:
        molecule_tab()
