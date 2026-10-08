"""
Two more quantum tools (October 2026), shown in 07 Quantum lab and in the Drug lab's Quantum section:

ADMET profile      21 absorption / distribution / metabolism / excretion / toxicity properties of a drug from a
                   quantum-kernel model next to classical models (chronocell/quantum/admet.py; Gate Q8).
Noise & mitigation what a noisy quantum computer does to a molecule's energy, and how zero-noise extrapolation and
                   symmetry verification recover it (chronocell/quantum/noisy.py; Gate Q9).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from chronocell import theme as T
from ui.common import banner, esc, html, readout, warning_card

ss = st.session_state
VALIDATION = Path(__file__).resolve().parent.parent / "validation"


def _standing(fname: str, label: str) -> None:
    try:
        r = json.loads((VALIDATION / fname).read_text(encoding="utf-8")) if (VALIDATION / fname).exists() else None
    except (OSError, ValueError):
        r = None
    if not r or "pass" not in r:
        banner(f"<b>{esc(label)}</b>: its pre-registered test has not been run yet.", "info")
        return
    extra = f" ({r['passed']} of {r['of']} endpoints)" if "passed" in r else f" ({r['within']} of {r['cases']} cases)" \
        if "within" in r else ""
    banner(f"<b>{esc(label)}</b> on held-out data: <b>{'pass' if r['pass'] else 'fail'}</b>{extra}. Details in "
           "validation/RESULTS.md and on 08 Scoreboard.", "info" if r["pass"] else "warn")


# ======================================================================================
# ADMET profile
# ======================================================================================
@st.cache_resource(show_spinner="Downloading the TDC ADMET benchmark (1.5 MB, MD5-checked) and training 21 models…")
def _admet_models():
    from chronocell.quantum import admet as A
    z = A.archive()
    cfg = A.settings()
    return {n: A.train(z, n, cfg.get(n)) for n in A.ENDPOINTS}


def _smiles_for(pick: str, typed: str) -> str | None:
    from chronocell import drug_info as DI
    if typed.strip():
        return typed.strip()
    if pick and pick != "(type a SMILES)":
        p = DI.properties(pick)
        return p.get("IsomericSMILES") or p.get("SMILES") or p.get("ConnectivitySMILES")
    return None


def admet_panel(kp: str = "qadmet") -> None:
    from chronocell import drug_info as DI
    from chronocell.quantum import admet as A, molfeat as MF
    html('<p class="cc-note">Most drug candidates fail not because they miss their target but because of <b>ADMET</b>: '
         'how the body Absorbs, Distributes, Metabolises and Excretes them, and their Toxicity. This tool predicts 21 such '
         'properties from the molecule\'s structure, each learned from the Therapeutics Data Commons benchmark: a '
         '<b>quantum-kernel model</b> (the molecule\'s 8 main descriptor directions set the angles of an 8-qubit '
         'feature map) next to classical models on the same inputs.</p>')
    banner("<b>A screen for teaching and triage, not a safety assessment.</b> The models learn from a few hundred to a few "
           "thousand measured compounds each; real decisions need laboratory and clinical data.", "warn")
    _standing("results_admet.json", "Quantum ADMET profile (Gate Q8)")
    c1, c2 = st.columns([1, 1.4])
    pick = c1.selectbox("Drug", ["(type a SMILES)"] + sorted({n for c in DI.CLASSES.values() for n in c.lookup}),
                        index=0, key=f"{kp}_pick")
    typed = c2.text_input("SMILES", key=f"{kp}_smiles", placeholder="e.g. ONC(=O)CCCCCCC(=O)Nc1ccccc1 (vorinostat)",
                          help="Picked drugs are looked up on PubChem (only the name is sent).")
    if st.button("Build the ADMET profile", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        try:
            smi = _smiles_for(pick, typed)
            if not smi:
                warning_card("No molecule", "Pick a drug or type a SMILES.")
                return
            d = MF.descriptors(smi)
            X = np.array([[d[f] for f in MF.FEATURES]], float)
            models = _admet_models()
            res = {}
            for n, m in models.items():
                pr = m.predict(X)
                res[n] = {k: float(v[0]) for k, v in pr.items()}
            ss[f"{kp}_res"] = {"smiles": smi, "name": pick if not typed.strip() else "typed molecule", "pred": res,
                               "desc": d}
        except Exception as exc:
            warning_card("The profile could not be built", str(exc))
            return
    r = ss.get(f"{kp}_res")
    if not r:
        return
    models = _admet_models()
    html(f'<p class="cc-note"><b>{esc(r["name"])}</b> · <code>{esc(r["smiles"][:120])}</code></p>')
    rows = []
    for grp, names in A.GROUPS.items():
        for n in names:
            task, lab, what, unit, logt = A.ENDPOINTS[n]
            p = r["pred"][n]
            if task == "cls":
                q, c = p["quantum"], p["rbf17"]
                call = lambda v: f"{'yes' if v > 0 else 'no'} ({unit if v > 0 else 'not ' + unit})"   # noqa: E731
                rows.append({"group": grp, "property": lab, "means": what, "quantum": call(q), "classical (17 inputs)": call(c),
                             "agree": "✓" if (q > 0) == (c > 0) else "✗", "score_q": q})
            else:
                m = models[n]
                vq, vc = m.predict_value(p["quantum"]), m.predict_value(p["rbf17"])
                pct = float((m.ytr <= p["quantum"]).mean() * 100)
                rows.append({"group": grp, "property": lab, "means": what, "quantum": f"{float(vq):.3g} {unit}",
                             "classical (17 inputs)": f"{float(vc):.3g} {unit}",
                             "agree": f"percentile {pct:.0f} of training", "score_q": pct / 50 - 1})
    df = pd.DataFrame(rows)
    st.dataframe(df.drop(columns=["score_q"]), hide_index=True, width="stretch", key=f"{kp}_table")
    fig = go.Figure(go.Bar(x=df["score_q"], y=df["property"], orientation="h",
                           marker_color=[T.ACCENT if v > 0 else T.TERRACOTTA for v in df["score_q"]],
                           hovertext=df["quantum"], hoverinfo="text"))
    fig.update_layout(**T.plot_layout(560, margin=dict(l=200, r=10, t=10, b=40), yaxis=dict(autorange="reversed"),
                                      xaxis=dict(title="quantum model: decision value (yes/no) or percentile − 50 % (values)")))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_bars")
    html('<p class="cc-note">Blue: the property is present / above the training median; red: absent / below. For yes/no '
         'properties the bar is the SVM\'s decision value (further from zero = more confident). "Agree" compares the '
         'quantum model with the classical model on all 17 descriptors.</p>')


# ======================================================================================
# Noise & error mitigation
# ======================================================================================
@st.cache_data(show_spinner=False, max_entries=32)
def _mitigation_case(name: str, scale: float, p2: float, sv: bool, method: str) -> dict:
    from chronocell.quantum import molecules as M, noisy as N
    _, fn, ch = M.LIBRARY[name]
    ints = M.integrals(fn(scale), ch)
    qm = M.qubit_hamiltonian(M.active_space(ints, M.rhf(ints), 2, 2))
    r = M.vqe_uccsd(qm)
    c = N.vqe_circuit(qm, r.generators, r.thetas, f"{name} {scale}x UCCSD")
    keep = np.zeros(2 ** qm.n_qubits)
    keep[M.sector(qm)] = 1.0
    H = qm.H.toarray()
    noise = {"p1": p2 / 10, "p2": p2}
    raw = N.expectation(c, H, noise)
    z = N.zne(c, H, noise, scales=(1, 3, 5), method=method, ideal=r.energy, keep=keep if sv else None)
    cnt = c.counts()
    return {"fci": r.fci, "ideal": r.energy, "hf": r.hf, "raw": raw, "values": z.values, "mitigated": z.mitigated,
            "noisy_sv": z.noisy, "gates": cnt["gates"], "cx": cnt["cx_equivalent"], "depth": cnt["depth"],
            "qasm": c.to_qasm()}


def mitigation_panel(kp: str = "qnoise") -> None:
    from chronocell.quantum import molecules as M
    html('<p class="cc-note">Real quantum computers make mistakes: every gate has a small chance of scrambling the qubits '
         'it touches. Here a molecule\'s VQE circuit (4 qubits, about 200 gates) runs on a <b>simulated noisy chip</b> '
         'gate by gate. Two standard repairs: <b>symmetry verification</b> throws away results with the wrong number of '
         'electrons, and <b>zero-noise extrapolation</b> runs the circuit with the noise deliberately tripled and '
         'quintupled (by adding gate pairs that cancel), then extrapolates back to zero noise.</p>')
    _standing("results_mitigation.json", "Error mitigation (Gate Q9)")
    names = [n for n in M.LIBRARY if n not in ("H2",)]
    c1, c2, c3 = st.columns(3)
    name = c1.selectbox("Molecule", names, index=names.index("LiH"), key=f"{kp}_name",
                        format_func=lambda n: f"{n} · {M.LIBRARY[n][0]}")
    scale = c2.slider("Bond length (× equilibrium)", 0.8, 2.5, 1.0, 0.05, key=f"{kp}_scale")
    p2 = c3.select_slider("Two-qubit gate error", [0.0005, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02], value=0.003,
                          key=f"{kp}_p2", help="0.003 is about the best of today's superconducting chips (2024-2026); "
                                               "one-qubit gates get a tenth of it.")
    c4, c5 = st.columns(2)
    sv = c4.toggle("Symmetry verification", True, key=f"{kp}_sv")
    method = c5.segmented_control("Extrapolation", ["richardson", "linear", "exp"], default="richardson", required=True,
                                  key=f"{kp}_method")
    if not st.button("Run on the noisy simulated chip", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        if ss.get(f"{kp}_key") != (name, scale, p2, sv, method):
            return
    ss[f"{kp}_key"] = (name, scale, p2, sv, method)
    try:
        with st.spinner("Density-matrix simulation at noise ×1, ×3, ×5…"):
            r = _mitigation_case(name, float(scale), float(p2), bool(sv), method)
    except ValueError as exc:
        warning_card("This molecule does not fit", str(exc))
        return
    err = lambda e: f"{1e3 * (e - r['ideal']):+.2f}"                                          # noqa: E731
    readout([("Noise-free circuit", f"{r['ideal']:.5f}", f"Ha · exact (active space) {r['fci']:.5f}"),
             ("Noisy chip, raw", f"{r['raw']:.5f}", f"Ha · error {err(r['raw'])} mHa"),
             ("Noisy" + (" + symmetry verified" if sv else ""), f"{r['noisy_sv']:.5f}", f"Ha · error {err(r['noisy_sv'])} mHa"),
             ("After zero-noise extrapolation", f"{r['mitigated']:.5f}", f"Ha · error {err(r['mitigated'])} mHa"),
             ("Chemical accuracy", "yes" if abs(r["mitigated"] - r["ideal"]) <= 1.6e-3 else "no", "≤ 1.6 mHa"),
             ("Circuit", f"{r['gates']} gates", f"{r['cx']} CNOT-equivalents · depth {r['depth']}")])
    s = np.array([1, 3, 5], float)
    v = np.array(r["values"])
    xs = np.linspace(0, 5.2, 60)
    deg = 2 if method == "richardson" else 1
    fit = np.polyval(np.polyfit(s, v, deg), xs) if method != "exp" else None
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=s, y=v, mode="markers", marker=dict(size=11, color=T.TERRACOTTA), name="measured on the noisy chip"))
    if fit is not None:
        fig.add_trace(go.Scatter(x=xs, y=fit, mode="lines", line=dict(color=T.ACCENT, dash="dash"), name="extrapolation"))
    fig.add_trace(go.Scatter(x=[0], y=[r["mitigated"]], mode="markers", marker=dict(size=13, symbol="star", color=T.ACCENT),
                             name="mitigated (noise → 0)"))
    fig.add_hline(y=r["ideal"], line=dict(color=T.INK, dash="dot"), annotation_text="noise-free answer")
    fig.update_layout(**T.plot_layout(340, showlegend=True, legend=dict(orientation="h", y=1.15),
                                      xaxis=dict(title="noise scale (1 = the chip as it is)"), yaxis=dict(title="energy (Ha)")))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_fig")
    if st.button("Sweep the chip's error rate", key=f"{kp}_sweep"):
        rows = []
        bar = st.progress(0.0)
        grid = [0.0005, 0.001, 0.002, 0.003, 0.005, 0.01]
        for i, p in enumerate(grid):
            x = _mitigation_case(name, float(scale), p, bool(sv), method)
            rows.append({"p2": p, "raw": abs(x["raw"] - x["ideal"]) * 1e3, "mitigated": abs(x["mitigated"] - x["ideal"]) * 1e3})
            bar.progress((i + 1) / len(grid))
        bar.empty()
        ss[f"{kp}_sweep_rows"] = {"key": (name, scale, sv, method), "rows": rows}
    sw = ss.get(f"{kp}_sweep_rows")
    if sw and sw["key"] == (name, scale, sv, method):
        d = pd.DataFrame(sw["rows"])
        fig = go.Figure([go.Scatter(x=d["p2"], y=np.maximum(d["raw"], 1e-4), mode="lines+markers", name="raw",
                                    line=dict(color=T.TERRACOTTA)),
                         go.Scatter(x=d["p2"], y=np.maximum(d["mitigated"], 1e-4), mode="lines+markers", name="mitigated",
                                    line=dict(color=T.ACCENT))])
        fig.add_hline(y=1.6, line=dict(color=T.INK, dash="dot"), annotation_text="chemical accuracy")
        fig.update_layout(**T.plot_layout(300, showlegend=True, legend=dict(orientation="h", y=1.15),
                                          xaxis=dict(title="two-qubit gate error", type="log"),
                                          yaxis=dict(title="energy error (mHa)", type="log")))
        st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_sweepfig")
    st.download_button("Download the circuit (OpenQASM 2.0)", r["qasm"], file_name=f"{name}_{scale}x_uccsd.qasm",
                       key=f"{kp}_qasm")
