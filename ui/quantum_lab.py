"""
07 Quantum lab (Research mode) and the quantum panels shown inside the other workspaces.

Every panel turns a ChronoCell problem into the form a quantum computer takes (a QUBO, a qubit Hamiltonian, a
feature map, a walk), solves it on the statevector SIMULATOR of chronocell/quantum, and shows the classical answer
on the same input next to it. Nothing is computed until a Run button is pressed. Results are labelled
"simulated quantum"; the measured standing comes from Gate Q (validation/results_gateq.json).
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from chronocell import theme as T
from chronocell.quantum import LABEL, chem as CH, kernels as KQ, lattice as LT, problems as PR, qubo as QB, sim, walk as QW
from ui.common import Dataset, banner, contacts_are_synthetic, esc, html, readout, warning_card

ss = st.session_state
VALIDATION = Path(__file__).resolve().parent.parent / "validation"
PROBLEMS = {
    "tad": ("TAD boundaries", "QAOA finds where the DNA's 'rooms' (domains) begin and end."),
    "lattice": ("Lattice fold", "QAOA folds a short stretch of DNA on a grid to match its contacts."),
    "drug": ("Drug combination", "Pick the mix and doses of drug classes that restore the fold most."),
    "vqe": ("Molecule (VQE)", "Compute a small molecule's energy the way quantum chemists do on quantum computers."),
    "qsvm": ("Gene classifier", "A quantum-kernel machine predicts whether a gene is switched on."),
    "walk": ("Quantum walk", "Watch a signal spread through the contact network, quantum vs classical."),
    "swap": ("Similarity (swap test)", "A quantum circuit measures how alike two folds are."),
    "mol": ("Drug molecules (VQE)", "Energies of drug-like chemical groups from an active-space VQE."),
    "safety": ("Heart safety", "A quantum-kernel classifier screens molecules for hERG blocking."),
    "dock": ("Docking", "QAOA matches a drug's groups to its protein pocket (max clique)."),
    "admet": ("ADMET profile", "21 absorption, metabolism and toxicity properties of a drug from a quantum-kernel model."),
    "noise": ("Noise & mitigation", "Run a molecule on a simulated noisy chip and repair the answer (zero-noise extrapolation)."),
}
KIND_COLOUR = {"exact": T.INK, "classical": T.TERRACOTTA, "quantum-inspired": T.OCHRE, "simulated quantum": T.ACCENT}


# ======================================================================================
# Shared pieces
# ======================================================================================
def frozen_settings(name: str = "QUANTUM_GATEQ") -> dict | None:
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("chronocell_frozen", VALIDATION / "frozen.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return getattr(mod, name, None)
    except Exception:
        return None


def gate_q() -> dict | None:
    p = VALIDATION / "results_gateq.json"
    try:
        return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
    except (OSError, ValueError):
        return None


PART_LABEL = {"q1": "Q1 QAOA solves the domain QUBO", "q2": "Q2 quantum domain calls vs classical callers",
              "q3": "Q3 VQE within chemical accuracy", "q4": "Q4 quantum-kernel gene classifier vs classical",
              "q5": "Q5 quantum heart-safety screen vs classical", "q6": "Q6 molecule energies vs independent reference",
              "q7": "Q7 quantum docking vs random search"}


def standing(parts: tuple[str, ...] | None = None) -> None:
    r = gate_q()
    if r is None or "overall" not in r:
        banner(f"<b>{esc(LABEL.capitalize())}.</b> Gate Q has not been run, so these panels have no measured standing yet.",
               "info")
        return
    parts = parts or tuple(PART_LABEL)
    items = []
    for k in parts:
        if k in r.get("overall", {}):
            items.append(f"{PART_LABEL[k]}: <b>{'pass' if r['overall'][k] else 'fail'}</b>")
    banner(f"<b>{esc(LABEL.capitalize())}.</b> Gate Q (held-out data): " + " · ".join(items) +
           ". Details in validation/RESULTS.md.", "info" if all(r["overall"].get(k, False) for k in parts) else "warn")


R2_LABEL = {"q2b": "Q2b domain calls, new cell lines", "q4b": "Q4b gene classifier, new cell line",
            "q6b": "Q6b ADAPT-VQE, new molecules", "q6c": "Q6c ADAPT-VQE vs the lowest singlet, new molecules",
            "q6d": "Q6d ADAPT-VQE with escape and several starts, new molecules", "q7b": "Q7b docking with refinement, Astex set"}


def round2_standing(parts: tuple[str, ...]) -> None:
    """Round 2 retests (validation/results_round2.json): shown only for the parts that have been run."""
    p = VALIDATION / "results_round2.json"
    try:
        r = json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
    except (OSError, ValueError):
        r = {}
    done = [k for k in parts if k in r and "pass" in r[k]]
    if not done:
        return
    items = [f"{R2_LABEL[k]}: <b>{'pass' if r[k]['pass'] else 'fail'}</b>" for k in done]
    banner("<b>Round 2</b> (new method, new held-out data; the original verdict above stands): " + " · ".join(items) +
           ".", "info" if all(r[k]["pass"] for k in done) else "warn")


def _layout(height: int, **kw) -> dict:
    return T.plot_layout(height, **kw)


def circuit_figure(c: sim.Circuit, max_cols: int = 48, max_qubits: int = 20) -> go.Figure:
    """Gate diagram (qubit 0 on top, as Qiskit draws it)."""
    nq = min(c.n, max_qubits)
    layer = np.zeros(c.n, dtype=int)
    items = []
    truncated = False
    for op in c.ops:
        g, qs = op[0], op[1]
        lo_q, hi_q = min(qs), max(qs)
        span = list(range(lo_q, hi_q + 1))
        col = int(layer[span].max())
        if col >= max_cols:
            truncated = True
            break
        layer[span] = col + 1
        items.append((col, g, qs, op[2]))
    fig = go.Figure()
    ncol = int(layer.max()) if items else 1
    for q in range(nq):
        fig.add_trace(go.Scatter(x=[-0.6, ncol - 0.4], y=[q, q], mode="lines", line=dict(color=T.RULE_STRONG, width=1),
                                 hoverinfo="skip"))
    bx, by, bt, cx_, cy_, tx, ty = [], [], [], [], [], [], []
    for col, g, qs, ps in items:
        if g == "barrier":
            fig.add_trace(go.Scatter(x=[col, col], y=[-0.5, nq - 0.5], mode="lines",
                                     line=dict(color=T.GHOST, width=1, dash="dot"), hoverinfo="skip"))
            continue
        if any(q >= nq for q in qs):
            continue
        if len(qs) > 1:
            fig.add_trace(go.Scatter(x=[col, col], y=[min(qs), max(qs)], mode="lines",
                                     line=dict(color=T.INK_2 if hasattr(T, "INK_2") else T.INK, width=1.3), hoverinfo="skip"))
        if g == "cx":
            cx_.append(col), cy_.append(qs[0])
            tx.append(col), ty.append(qs[1])
        elif g == "cz":
            for q in qs:
                cx_.append(col), cy_.append(q)
        elif g == "cswap":
            cx_.append(col), cy_.append(qs[0])
            for q in qs[1:]:
                bx.append(col), by.append(q), bt.append("×")
        elif g == "swap":
            for q in qs:
                bx.append(col), by.append(q), bt.append("×")
        else:
            lab = g.upper() + (f"<br>{ps[0]:.2f}" if ps else "")
            if g == "rzz":
                lab = f"ZZ<br>{ps[0]:.2f}"
            for q in qs:
                bx.append(col), by.append(q), bt.append(lab)
    fig.add_trace(go.Scatter(x=bx, y=by, mode="markers+text", text=bt, textfont=dict(size=8, color=T.INK),
                             marker=dict(symbol="square", size=26, color=T.ACCENT_SOFT, line=dict(color=T.ACCENT, width=1)),
                             hoverinfo="text"))
    fig.add_trace(go.Scatter(x=cx_, y=cy_, mode="markers", marker=dict(size=8, color=T.INK), hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=tx, y=ty, mode="markers", marker=dict(symbol="circle-cross-open", size=16, color=T.INK,
                                                                      line=dict(width=1.3)), hoverinfo="skip"))
    fig.update_layout(**_layout(max(140, 34 * nq + 40), margin=dict(l=40, r=10, t=10, b=10),
                                xaxis=dict(visible=False, range=[-0.8, max(ncol, 6) - 0.2]),
                                yaxis=dict(autorange="reversed", tickvals=list(range(nq)),
                                           ticktext=[f"q{q}" for q in range(nq)], showgrid=False)))
    if truncated or c.n > nq:
        fig.add_annotation(x=1, y=1.04, xref="paper", yref="paper", showarrow=False, font=dict(size=10, color=T.MUTED),
                           text=f"first {ncol} layers{f', first {nq} qubits' if c.n > nq else ''} shown")
    return fig


def circuit_block(c: sim.Circuit, key: str, noise: dict | None = None) -> None:
    cnt = c.counts()
    st.plotly_chart(circuit_figure(c), theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{key}_circ")
    F = sim.fidelity_estimate(cnt, noise or sim.NOISE_DEFAULT)
    readout([("Qubits", f"{cnt['qubits']}", ""), ("Gates", f"{cnt['gates']:,}", f"depth {cnt['depth']}"),
             ("Two-qubit gates<small>CX equivalents after compiling</small>", f"{cnt['cx_equivalent']:,}", ""),
             ("Survival on today's hardware<small>approximate noise model</small>", f"{100 * F:.1f}", "%")])
    st.download_button("Circuit (OpenQASM 2.0)", c.to_qasm(), f"{key}.qasm", "text/plain", icon=":material/download:",
                       key=f"{key}_qasm", help="Opens in Qiskit (QuantumCircuit.from_qasm_str) and in IBM Quantum's "
                                                "composer; running it on a real device needs your own IBM Quantum account.")


def energy_bars(probs: np.ndarray, energies: np.ndarray, n: int, key: str, top: int = 12) -> None:
    order = np.argsort(-probs)[:top]
    emin = float(energies.min())
    opt = np.isclose(energies[order], emin, atol=1e-9 * max(1.0, abs(emin)))
    fig = go.Figure(go.Bar(x=[sim.bitstring(int(i), n) for i in order], y=probs[order],
                           marker_color=[T.ACCENT if o else T.GHOST for o in opt],
                           hovertext=[f"energy {energies[i]:.3f}" for i in order], hoverinfo="text+y"))
    fig.update_layout(**_layout(230, margin=dict(l=52, r=10, t=10, b=90), xaxis=dict(tickangle=-60, tickfont=dict(size=9)),
                                yaxis=dict(title="probability")))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{key}_bars")
    html('<p class="cc-note">Most likely measurement outcomes (qubit 0 on the right). Blue = a minimum-energy answer.</p>')


def convergence(history: list, key: str, label: str = "objective (energy units)") -> None:
    fig = go.Figure(go.Scatter(y=history, mode="lines", line=dict(color=T.ACCENT, width=1.4)))
    fig.update_layout(**_layout(200, xaxis=dict(title="circuit evaluations"), yaxis=dict(title=label)))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{key}_conv")


def qubo_heatmap(q: QB.QUBO, key: str) -> None:
    M = q.sym / 2 + np.diag(q.lin)
    lim = float(np.percentile(np.abs(M[M != 0]), 98)) if np.any(M != 0) else 1.0
    fig = go.Figure(go.Heatmap(z=M, zmid=0, zmin=-lim, zmax=lim, colorscale="RdBu_r", showscale=False,
                               hovertemplate="x%{y} · x%{x}: %{z:.3f}<extra></extra>"))
    fig.update_layout(**_layout(260, margin=dict(l=40, r=10, t=10, b=30), yaxis=dict(autorange="reversed")))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{key}_qubo")
    html(f'<p class="cc-note">The puzzle as a quantum computer sees it: {q.n} binary variables (one qubit each); the '
         'diagonal is what each variable costs on its own, off-diagonal cells what pairs cost together (red = '
         'penalised, blue = rewarded).</p>')


def noise_toggle(key: str) -> dict | None:
    on = st.toggle("Add hardware noise", False, key=f"{key}_noise",
                   help="Approximate model of today's superconducting devices: each gate fails with probability "
                        "0.1 % (one qubit) or 1 % (two qubits), each read-out bit flips with probability 2 %.")
    return dict(sim.NOISE_DEFAULT) if on else None


def dense_contacts(ds: Dataset, a: int, b: int) -> np.ndarray:
    m = (ds.ci >= a) & (ds.ci < b) & (ds.cj >= a) & (ds.cj < b)
    M = np.zeros((b - a, b - a))
    np.add.at(M, (ds.ci[m] - a, ds.cj[m] - a), np.asarray(ds.cm[m], float))
    return M + np.triu(M, 1).T


def _mb(ds: Dataset, local_bin: float) -> float:
    return (float(ds.chrom.bin_start(ds.bin0)) + local_bin * ds.chrom.resolution) / 1e6


# ======================================================================================
# 1 · TAD boundaries
# ======================================================================================
def tad_panel(ds: Dataset, lo: int, hi: int, kp: str = "qtad") -> None:
    html('<p class="cc-note">Each gap between neighbouring bins gets one qubit: 1 = a domain wall here. The cost '
         'rewards keeping enriched contacts inside a domain (modularity), and a penalty forbids domains smaller than '
         'the minimum size. QAOA, the gate-model quantum optimiser, searches all walls at once; exact enumeration, '
         'simulated annealing and simulated quantum annealing solve the same puzzle classically.</p>')
    if not ds.has_contacts:
        html('<p class="cc-note">Needs contacts (Data → Graph, or a state\'s Hi-C / Micro-C file).</p>')
        return
    if contacts_are_synthetic(ds):
        banner("Input is the <b>SYNTHETIC</b> reference contact map: the calls illustrate the method only.", "warn")
    standing(("q1", "q2"))
    round2_standing(("q2b",))
    fz = frozen_settings() or {}
    tad_set = fz.get("tad", {"res": 50_000, "min_size": 3, "gamma": 1.5, "boundary_cost": 0.0, "weight": "difference"})
    r2 = (frozen_settings("QUANTUM_ROUND2") or {}).get("q2b")
    choice = st.radio("Settings", ["gateq", "q2b"], horizontal=True, key=f"{kp}_settings",
                      format_func=lambda k: {"gateq": "Gate Q (original)",
                                             "q2b": "Round 2 (Q2b: chosen on three cell lines)"}[k]) if r2 else "gateq"
    if choice == "q2b":
        tad_set = r2["tad"]
    sk = "" if choice == "gateq" else "_q2b"                  # separate slider state per settings choice
    res = int(ds.chrom.resolution)
    f = max(1, int(round(tad_set["res"] / res)))
    c1, c2, c3, c4 = st.columns(4)
    nb = c1.slider("Window (bins) → qubits = bins − 1", 8, 21, 16, key=f"{kp}_nb",
                   help="Each extra bin doubles the simulator's memory; 21 bins = 20 qubits.")
    m = c2.slider("Minimum domain (bins)", 2, 5, int(tad_set["min_size"]), key=f"{kp}_m{sk}")
    gamma = c3.slider("Resolution γ", 0.5, 8.0, float(tad_set["gamma"]), 0.1, key=f"{kp}_g{sk}",
                      help="Higher: more, smaller domains.")
    p = c4.slider("QAOA depth p", 1, 8, int(fz.get("qaoa", {}).get("p", 3)), key=f"{kp}_p")
    span = 2 * m
    n_coarse = (hi - lo) // f
    if n_coarse < nb + 2 * span:
        html(f'<p class="cc-note">The region holds {n_coarse} bins of {f * res / 1000:g} kb; at least {nb + 2 * span} are '
             'needed. Widen the region.</p>')
        return
    start = st.slider(f"Window start (bin of {f * res / 1000:g} kb)", span, n_coarse - nb - span, span, key=f"{kp}_start")
    noise = noise_toggle(kp)
    run_key = f"{ds.key}:{lo}:{hi}:{f}:{nb}:{m}:{gamma}:{p}:{start}{sk}"
    if st.button("Run on the quantum simulator", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        a = lo + (start - span) * f
        b = lo + (start + nb + span) * f
        C = PR.coarsen(dense_contacts(ds, a, b), f)
        mask = PR.coverage_mask(C)
        exp = PR.expected_by_distance(C, span, mask)
        core = (span, span + nb)
        q, pos = PR.tad_qubo(C, core, gamma, m, span, exp, mask, boundary_cost=float(tad_set.get("boundary_cost", 0.0)),
                             weight=str(tad_set.get("weight", "difference")))
        with st.spinner(f"Simulating {q.n} qubits…"):
            e = q.energies()
            ex = QB.exact(q, e)
            sa = QB.anneal(q, seed=1)
            sq = QB.sqa(q, reads=16, sweeps=300, seed=1)
            r = sim.qaoa(e, p=p, shots=4096, objective=fz.get("qaoa", {}).get("objective", "cvar"), maxiter=80, seed=1)
            h, J, _ = q.ising()
            circ = sim.qaoa_circuit(h, J, r.gammas, r.betas, f"QAOA TAD boundaries p={p}")
            F = sim.fidelity_estimate(circ.counts(), noise or sim.NOISE_DEFAULT)
            noisy = sim.sample(r.probs, 4096, np.random.default_rng(7), F, (noise or sim.NOISE_DEFAULT)["readout"], q.n)
        ss[kp] = {"key": run_key, "C": C, "core": core, "q": q, "pos": pos, "e": e, "ex": ex, "sa": sa, "sq": sq, "r": r,
                  "circ": circ, "noisy": noisy, "F": F, "a": a, "f": f, "span": span, "noise": noise is not None}
    last = ss.get(kp)
    if not last or last["key"] != run_key:
        return
    q, e, ex, r = last["q"], last["e"], last["ex"], last["r"]
    emin = ex.energy
    tol = 1e-9 * max(1.0, abs(emin))
    noisy_best = int(last["noisy"][np.argmin(e[last["noisy"]])])
    pos = last["pos"]
    bp = lambda bits: [round(_mb(ds, last["a"] + k * last["f"]), 3) for k in PR.boundaries_from_bits(bits, pos)]
    rows = [("Exact (all 2^n states)", "exact", ex.energy, True, ex.seconds, bp(ex.bits)),
            ("Simulated annealing", "classical", last["sa"].energy, abs(last["sa"].energy - emin) <= tol,
             last["sa"].seconds, bp(last["sa"].bits)),
            ("Simulated quantum annealing", "quantum-inspired", last["sq"].energy, abs(last["sq"].energy - emin) <= tol,
             last["sq"].seconds, bp(last["sq"].bits)),
            (f"QAOA p={r.p} · ideal", "simulated quantum", r.best_energy, r.hit, r.seconds, bp(QB.to_bits(r.best, q.n))),
            (f"QAOA p={r.p} · with noise", "simulated quantum", float(e[noisy_best]), abs(e[noisy_best] - emin) <= tol,
             float("nan"), bp(QB.to_bits(noisy_best, q.n)))]
    tab = pd.DataFrame(rows, columns=["method", "kind", "energy", "found the optimum", "seconds", "walls (Mb)"])
    readout([("Qubits", f"{q.n}", f"{2 ** q.n:,} possible answers"),
             ("QAOA's chance of the best answer", f"{100 * r.p_optimal:.2f}", f"% (random guess {100 * ex.detail['n_optima'] / 2 ** q.n:.4f} %)"),
             ("Device", esc(r.device), f"{r.evaluations} circuit runs"),
             ("Noise survival", f"{100 * last['F']:.2g}", "% of shots unaffected")])
    st.dataframe(tab, hide_index=True, width="stretch", key=f"{kp}_race",
                 column_config={"energy": st.column_config.NumberColumn(format="%.3f"),
                                "seconds": st.column_config.NumberColumn(format="%.2f")})
    t1, t2, t3, t4 = st.tabs(["Contact map", "Puzzle (QUBO)", "Circuit", "Measurements"])
    with t1:
        C = last["C"]
        x = [round(_mb(ds, last["a"] + i * last["f"]), 3) for i in range(len(C))]
        fig = go.Figure(go.Heatmap(z=np.log1p(C), x=x, y=x, colorscale=T.CONTACT_SCALE, showscale=False,
                                   hovertemplate="%{x} · %{y} Mb<extra></extra>"))
        for name, col, bits in (("exact", T.INK, ex.bits), ("QAOA", T.ACCENT, QB.to_bits(r.best, q.n))):
            for w in bp(bits):
                fig.add_shape(type="line", x0=w, x1=w, y0=x[0], y1=x[-1], line=dict(color=col, width=2 if name == "QAOA" else 1,
                                                                                        dash="solid" if name == "QAOA" else "dot"))
        c0, c1_ = last["core"]
        fig.add_shape(type="rect", x0=x[c0], x1=x[c1_ - 1], y0=x[c0], y1=x[c1_ - 1], line=dict(color=T.MUTED, dash="dash"))
        fig.update_layout(**_layout(420, yaxis=dict(autorange="reversed", title="Mb"), xaxis=dict(title="Mb")))
        st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_map")
        html('<p class="cc-note">Blue lines: QAOA\'s walls; dotted: the exact optimum; dashed box: the window the qubits '
             'cover (the flanks only add contacts that cross its edges).</p>')
    with t2:
        qubo_heatmap(q, kp)
    with t3:
        circuit_block(last["circ"], f"{kp}_c", sim.NOISE_DEFAULT)
        convergence(r.history, kp)
    with t4:
        energy_bars(r.probs, e, q.n, kp)


# ======================================================================================
# 2 · Lattice fold
# ======================================================================================
def lattice_panel(ds: Dataset, lo: int, hi: int, coords: np.ndarray | None, kp: str = "qlat") -> None:
    html('<p class="cc-note">The textbook "folding on a lattice" problem: each bond of a short chain points along '
         'the grid (2 qubits per bond in 2D, 3 in 3D), segments that touch in the contact map are rewarded when they '
         'sit next to each other, overlaps are forbidden. The cost involves many qubits at once (a higher-order '
         'Hamiltonian), so it needs more gates than a QUBO.</p>')
    c1, c2, c3 = st.columns(3)
    dim = c1.segmented_control("Lattice", ["2D", "3D"], default="2D", required=True, key=f"{kp}_dim")
    d = 2 if dim == "2D" else 3
    beads = c2.slider("Beads (segments)", 5, 11 if d == 2 else 8, 9 if d == 2 else 7, key=f"{kp}_beads")
    p = c3.slider("QAOA depth p", 1, 6, 3, key=f"{kp}_p")
    nq = LT.n_qubits(beads, d)
    html(f'<p class="cc-note">{nq} qubits · {2 ** nq:,} possible folds.</p>')
    use_contacts = ds.has_contacts
    run_key = f"{ds.key}:{lo}:{hi}:{d}:{beads}:{p}"
    if st.button("Fold on the quantum simulator", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        if use_contacts:
            W = LT.weights_from_contacts(dense_contacts(ds, lo, hi), beads)
            src = "measured contacts" + (" (SYNTHETIC reference map)" if contacts_are_synthetic(ds) else "")
        else:
            from chronocell import domains as DOM
            x = ds.frames[0][lo:hi]
            ci, cj, cm = DOM.proximity_contacts(x, float(np.median(np.linalg.norm(np.diff(x, axis=0), axis=1))))
            M = np.zeros((hi - lo, hi - lo))
            np.add.at(M, (ci, cj), cm)
            W = LT.weights_from_contacts(M + M.T, beads)
            src = "3D proximity of the structure (no contacts loaded)"
        with st.spinner(f"Simulating {nq} qubits…"):
            t0 = time.perf_counter()
            e = LT.energies(W, d)
            t_enum = time.perf_counter() - t0
            k_ex = int(np.argmin(e))
            r = sim.qaoa(e, p=p, shots=4096, objective="cvar", maxiter=80, seed=3)
        ss[kp] = {"key": run_key, "W": W, "e": e, "exact": LT.fold_of(k_ex, W, d, e[k_ex]), "t_enum": t_enum, "r": r,
                  "qaoa": LT.fold_of(r.best, W, d, r.best_energy), "src": src, "d": d,
                  "terms": sim.diag_terms_count(e)}
    last = ss.get(kp)
    if not last or last["key"] != run_key:
        return
    r, ex, qa = last["r"], last["exact"], last["qaoa"]
    ag_q, ag_e = LT.contact_agreement(qa, last["W"]), LT.contact_agreement(ex, last["W"])
    cnt = sim.qaoa_counts(last["terms"], len(last["e"]).bit_length() - 1, r.p)
    readout([("QAOA found the best fold", "yes" if r.hit else "no", f"chance {100 * r.p_optimal:.2f} %"),
             ("Contacts realised<small>QAOA · exact</small>", f"{ag_q['contacts']} · {ag_e['contacts']}", ""),
             ("Gates it would compile to", f"{cnt['gates']:,}", f"{cnt['cx_equivalent']:,} two-qubit"),
             ("Exact enumeration", f"{last['t_enum']:.2f}", "s (classical)")])
    html(f'<p class="cc-note">Weights from {esc(last["src"])}. Pauli-Z terms by order: '
         f'{esc(", ".join(f"{k}-body {v}" for k, v in sorted(last["terms"].items())))}.</p>')
    cols = st.columns(2)
    for col, (title, fold) in zip(cols, (("QAOA (simulated quantum)", qa), ("Exact optimum (classical)", ex))):
        p3 = np.c_[fold.coords, np.zeros(len(fold.coords))] if last["d"] == 2 else fold.coords
        fig = go.Figure(go.Scatter3d(x=p3[:, 0], y=p3[:, 1], z=p3[:, 2], mode="lines+markers+text",
                                     text=[str(i) for i in range(len(p3))], textposition="top center",
                                     marker=dict(size=7, color=np.arange(len(p3)), colorscale="Viridis"),
                                     line=dict(color=T.INK, width=6)))
        for i, j in fold.contacts:
            fig.add_trace(go.Scatter3d(x=p3[[i, j], 0], y=p3[[i, j], 1], z=p3[[i, j], 2], mode="lines",
                                       line=dict(color=T.TERRACOTTA, width=3, dash="dash")))
        fig.update_layout(height=340, margin=dict(l=0, r=0, t=24, b=0), showlegend=False, paper_bgcolor="rgba(0,0,0,0)",
                          title=dict(text=title, font=dict(size=12)),
                          scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
                                     aspectmode="data"))
        col.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_{title[:4]}")
    html('<p class="cc-note">Red dashes: segments that touch on the lattice. A lattice fold is a cartoon of the '
         'contact pattern, not a structure: the 3D structure workspace builds the real model.</p>')
    energy_bars(r.probs, last["e"], len(last["e"]).bit_length() - 1, kp)


# ======================================================================================
# 3 · Drug combination
# ======================================================================================
def drug_panel(ds: Dataset, baseline: Dataset | None, b0: float, frame: int, kp: str = "qdrug") -> None:
    from chronocell import therapy as TH
    html('<p class="cc-note">Four drug classes × four dose levels (0, ⅓, ⅔, full) = 8 qubits, 256 possible '
         'treatments. The Drug lab\'s simulator is run on each drug alone and on each pair, a quadratic model of '
         'restoration is fitted to those runs, and the dose mix that maximises restoration minus a dose cost is '
         'found by QAOA and by exact search. The winning mix is then simulated again to check the model.</p>')
    banner("<b>Simulation, not a treatment recommendation.</b> The restoration numbers come from the Drug lab's "
           "mechanism simulator, which is not validated.", "warn")
    same = (baseline is not None and baseline.chrom.name == ds.chrom.name and baseline.chrom.resolution == ds.chrom.resolution
            and baseline.key != ds.key)
    if not same:
        html('<p class="cc-note">Needs a healthy baseline covering the same DNA: add a Healthy Control state (or switch on '
             'demo patients) and select the disease state.</p>')
        return
    g0, g1 = max(ds.bin0, baseline.bin0), min(ds.bin0 + ds.n, baseline.bin0 + baseline.n)
    xp = ds.frames[frame][g0 - ds.bin0:g1 - ds.bin0]
    xh = baseline.frames[0][g0 - baseline.bin0:g1 - baseline.bin0]
    sig = ds.epi[g0 - ds.bin0:g1 - ds.bin0]
    valid = ds.valid[g0 - ds.bin0:g1 - ds.bin0]
    c1, c2 = st.columns(2)
    cost = c1.slider("Dose cost (restoration % per full dose)", 0.0, 30.0, 5.0, 1.0, key=f"{kp}_cost",
                     help="How much restoration a full dose must buy to be worth giving: a stand-in for side effects.")
    eff = c2.slider("Maximum effect at full dose", 0.2, 1.0, 0.8, 0.05, key=f"{kp}_eff")
    model_key = f"{ds.key}:{baseline.key}:{frame}:{g0}:{g1}:{eff}"
    if st.button("Build the restoration model (10 simulations)", key=f"{kp}_model", icon=":material/science:"):
        lo, hi = TH.suggest_window(xp, xh, sig, 400)          # the 400 beads where the fold is most abnormal
        bar = st.progress(0.0, "Simulating…")
        mdl = PR.drug_combo_model(xp[lo:hi], sig[lo:hi], valid[lo:hi], b0, xh[lo:hi], eff, relax_iters=10,
                                  signal_ref=baseline.epi, progress=lambda f_, t: bar.progress(min(1.0, f_), t))
        bar.empty()
        ss[f"{kp}_mdl"] = {"key": model_key, "model": mdl, "lo": lo, "hi": hi}
    m = ss.get(f"{kp}_mdl")
    if not m or m["key"] != model_key:
        return
    mdl, lo, hi = m["model"], m["lo"], m["hi"]
    html(f'<p class="cc-note">Region: {ds.chrom.name}:{float(ds.chrom.bin_start(g0 + lo)) / 1e6:.2f}–'
         f'{float(ds.chrom.bin_end(g0 + hi - 1)) / 1e6:.2f} Mb ({hi - lo} beads, where the fold is most abnormal).</p>')
    q = PR.drug_combo_qubo(mdl["alpha"], mdl["beta"], mdl["eta"], cost, mdl["names"])
    e = q.energies()
    ex = QB.exact(q, e)
    run_key = f"{model_key}:{cost}"
    if ss.get(f"{kp}_run", {}).get("key") != run_key:
        r = sim.qaoa(e, p=3, shots=2048, objective="cvar", maxiter=80, seed=5)
        ss[f"{kp}_run"] = {"key": run_key, "r": r, "check": None}
    run = ss[f"{kp}_run"]
    r = run["r"]
    d_q = PR.drug_doses(QB.to_bits(r.best, q.n), len(mdl["keys"]))
    d_x = PR.drug_doses(ex.bits, len(mdl["keys"]))
    single = {k: v[-1] for k, v in mdl["singles"].items()}
    best_single = max(single, key=single.get)
    rows = [{"drug": nm, "QAOA dose": f"{dq:.0%}", "exact dose": f"{dx:.0%}"} for nm, dq, dx in zip(mdl["names"], d_q, d_x)]
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch", key=f"{kp}_doses")
    readout([("QAOA found the best mix", "yes" if r.hit else "no", f"chance {100 * r.p_optimal:.1f} %"),
             ("Model's restoration of the mix", f"{PR.predicted_restoration(mdl, d_q):.0f}", "%"),
             ("Best single drug (simulated)", f"{single[best_single]:.0f}", f"% · {esc(TH.DRUGS[best_single].name)}")])
    if st.button("Check the chosen mix in the simulator", key=f"{kp}_check"):
        with st.spinner("Simulating the combination…"):
            run["check"] = PR.combo_restoration(xp[lo:hi], sig[lo:hi], valid[lo:hi], b0, xh[lo:hi],
                                                dict(zip(mdl["keys"], d_q)), eff, 10, baseline.epi)
    if run.get("check") is not None:
        html(f'<p class="cc-meta">Simulated restoration of the chosen mix: <b>{run["check"]:.0f} %</b> (model said '
             f'{PR.predicted_restoration(mdl, d_q):.0f} %).</p>')
    with st.expander("Puzzle, circuit and measurements", expanded=False):
        qubo_heatmap(q, kp)
        h, J, _ = q.ising()
        circuit_block(sim.qaoa_circuit(h, J, r.gammas, r.betas, "QAOA drug combination"), f"{kp}_c")
        energy_bars(r.probs, e, q.n, kp)


# ======================================================================================
# 4 · Molecule (VQE)
# ======================================================================================
@st.cache_data(show_spinner=False, max_entries=64)
def _molecule_vqe(name: str, R: float, ansatz: str, noisy: bool) -> dict:
    m = CH.molecule(name, R)
    v = CH.vqe(m, ansatz, layers=2, noise=sim.NOISE_DEFAULT if noisy else None)
    return {"e_hf": m.e_hf, "e_fci": m.e_fci, "e_vqe": v.energy, "e_noisy": v.noisy_energy, "fidelity": v.fidelity,
            "electrons": v.electrons,
            "evals": v.evaluations, "history": v.history, "qasm": v.circuit.to_qasm(), "paulis": m.paulis,
            "circuit_ops": v.circuit.ops, "n": v.circuit.n}


def vqe_panel(kp: str = "qvqe") -> None:
    html('<p class="cc-note">The flagship quantum-chemistry algorithm. The molecule\'s electrons are written onto '
         '4 qubits (Jordan–Wigner), a parameterised circuit prepares a trial state, the energy is measured, and a '
         'classical optimiser turns the knobs until the energy is lowest. Everything from the integrals onwards is '
         'computed here from scratch and checked against the textbook values (Szabo & Ostlund).</p>')
    banner("These two-electron molecules are what quantum computers can handle today. A drug from the Drug lab has "
           "30–100 heavy atoms and would need thousands of error-corrected qubits.", "info")
    standing(("q3",))
    c1, c2, c3 = st.columns(3)
    name = c1.segmented_control("Molecule", list(CH.MOLECULES), default="H2", required=True, key=f"{kp}_mol")
    R = c2.slider("Bond length (Å)", 0.4, 3.0, float(CH.EQUILIBRIUM_ANGSTROM[name]), 0.01, key=f"{kp}_R")
    ansatz = c3.segmented_control("Circuit (ansatz)", ["UCCSD", "Hardware-efficient"], default="UCCSD", required=True,
                                  key=f"{kp}_ans", help="UCCSD: built from the chemistry; hardware-efficient: generic "
                                                       "layers that suit real devices.")
    noisy = noise_toggle(kp) is not None
    a = "uccsd" if ansatz == "UCCSD" else "hea"
    key = (name, float(R), a, noisy)
    if st.button("Run VQE on the quantum simulator", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        with st.spinner("Optimising the circuit…"):
            ss[f"{kp}_last"] = {"key": key, "res": _molecule_vqe(name, float(R), a, noisy)}
    last = ss.get(f"{kp}_last")
    if not last or last["key"] != key:
        return
    res = last["res"]
    err = res["e_vqe"] - res["e_fci"]
    readout([("Hartree–Fock", f"{res['e_hf']:.5f}", "Ha"), ("Exact (FCI)", f"{res['e_fci']:.5f}", "Ha"),
             ("VQE", f"{res['e_vqe']:.5f}", f"Ha · error {1e3 * err:+.3f} mHa"),
             ("Chemical accuracy", "yes" if abs(err) <= CH.CHEMICAL_ACCURACY else "no", "≤ 1.6 mHa"),
             ("Electrons in the state", f"{res['electrons']:.3f}", "the molecule has 2")])
    if a == "hea":
        html('<p class="cc-note">The hardware-efficient circuit does not conserve the number of electrons by itself; a '
             'penalty in its cost keeps it at two (without it, Gate Q’s test saw it drift to three electrons on HeH+, '
             'with energies below the true one). It can still stop in a local minimum: compare with UCCSD.</p>')
    if noisy:
        html(f'<p class="cc-meta">With hardware noise the same circuit would read <b>{res["e_noisy"]:.4f} Ha</b> '
             f'(error {1e3 * (res["e_noisy"] - res["e_fci"]):+.0f} mHa; {100 * res["fidelity"]:.0f} % of runs unaffected).</p>')
    t1, t2, t3 = st.tabs(["Energy curve", "Circuit", "Hamiltonian"])
    with t1:
        if st.button("Compute the curve (12 bond lengths)", key=f"{kp}_curve"):
            Rs = np.round(np.linspace(0.45, 2.8, 12), 3)
            ss[f"{kp}_c"] = {"key": f"{name}:{a}:{noisy}", "rows": [{"R": float(x), **{k: v for k, v in
                             _molecule_vqe(name, float(x), a, noisy).items() if k.startswith("e_")}} for x in Rs]}
        cur = ss.get(f"{kp}_c")
        if cur and cur["key"] == f"{name}:{a}:{noisy}":
            df = pd.DataFrame(cur["rows"])
            fig = go.Figure()
            for col, lab, colr, dash in (("e_hf", "Hartree–Fock", T.MUTED, "dot"), ("e_fci", "Exact (FCI)", T.INK, "solid"),
                                         ("e_vqe", "VQE (simulated quantum)", T.ACCENT, "solid")) + \
                    ((("e_noisy", "VQE with noise", T.TERRACOTTA, "dash"),) if noisy else ()):
                fig.add_trace(go.Scatter(x=df["R"], y=df[col], name=lab, mode="lines+markers",
                                         line=dict(color=colr, dash=dash, width=1.6), marker=dict(size=5)))
            fig.update_layout(**_layout(320, showlegend=True, legend=dict(orientation="h", y=1.12),
                                        xaxis=dict(title="bond length (Å)"), yaxis=dict(title="energy (hartree)")))
            st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_curvefig")
        convergence(res["history"], kp, "energy (hartree)")
    with t2:
        c = sim.Circuit(res["n"], f"VQE {name}")
        c.ops = list(res["circuit_ops"])
        circuit_block(c, f"{kp}_circ")
    with t3:
        terms = pd.DataFrame(sorted(res["paulis"].items(), key=lambda t: -abs(t[1])), columns=["Pauli string", "coefficient (Ha)"])
        st.dataframe(terms, hide_index=True, width="stretch", key=f"{kp}_paulis")
        html('<p class="cc-note">The qubit Hamiltonian: a weighted sum of Pauli strings (I, X, Y, Z on qubits 3..0). '
             'A quantum computer measures each group and adds them up.</p>')


# ======================================================================================
# 5 · Gene classifier (QSVM) and gene group
# ======================================================================================
def qsvm_panel(tab: pd.DataFrame | None, kp: str = "qsvm") -> None:
    html('<p class="cc-note">Each gene\'s features are written into the angles of a small quantum circuit (a ZZ '
         'feature map: one qubit per feature, entangled in pairs). How similar two genes look to the quantum computer '
         '(the overlap of their circuits\' states) becomes the kernel of a support-vector machine. A classical '
         'RBF-kernel SVM and logistic regression learn from the same genes.</p>')
    standing(("q4",))
    round2_standing(("q4b",))
    if tab is None or "expression" not in tab.columns or tab["expression"].notna().sum() < 30:
        html('<p class="cc-note">Needs measured expression for at least 30 genes in view (Genes → "Add measured '
             'expression", or a state\'s RNA-seq file). Gate Q4 tested the method on GM12878 → IMR-90 Hi-C features.</p>')
        return
    feats = [c for c in ("score", "crowding", "signal") if c in tab.columns]
    d = tab.dropna(subset=feats + ["expression"])
    thr = st.slider("Expressed if above (expression units)", float(d["expression"].min()), float(d["expression"].max()),
                    float(d["expression"].median()), key=f"{kp}_thr")
    y = (d["expression"].to_numpy(float) > thr).astype(int)
    if y.min() == y.max():
        html('<p class="cc-note">All genes fall on one side of the threshold.</p>')
        return
    X = d[feats].to_numpy(float)
    fz = (frozen_settings() or {}).get("qsvm", {"bandwidth": 0.05, "reps": 1, "C": 100.0})
    if st.button("Train on 70 %, test on 30 %", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        from sklearn.metrics import roc_auc_score
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVC
        from sklearn.linear_model import LogisticRegression
        tr, te = train_test_split(np.arange(len(y)), test_size=0.3, stratify=y, random_state=0)
        a, b = KQ.scale_features(X[tr], X[te])
        qs = KQ.QSVM(fz["C"], fz["bandwidth"], fz["reps"]).fit(a, y[tr])
        sc = StandardScaler().fit(X[tr])
        rbf = SVC(C=10.0, kernel="rbf", gamma="scale", class_weight="balanced").fit(sc.transform(X[tr]), y[tr])
        lr = LogisticRegression(class_weight="balanced", max_iter=2000).fit(sc.transform(X[tr]), y[tr])
        ss[kp] = {"auc": {"Quantum-kernel SVM (simulated)": roc_auc_score(y[te], qs.decision_function(b)),
                          "RBF-kernel SVM (classical)": roc_auc_score(y[te], rbf.decision_function(sc.transform(X[te]))),
                          "Logistic regression (classical)": roc_auc_score(y[te], lr.decision_function(sc.transform(X[te])))},
                  "K": KQ.kernel(a[:40], None, fz["bandwidth"], fz["reps"]), "feats": feats, "n": (len(tr), len(te)),
                  "circ": KQ.feature_map_circuit(a[0], fz["bandwidth"], fz["reps"])}
    last = ss.get(kp)
    if not last:
        return
    readout([(k, f"{v:.3f}", "AUC on the held-out 30 %") for k, v in last["auc"].items()])
    html(f'<p class="cc-note">{last["n"][0]} training and {last["n"][1]} test genes; features: {esc(", ".join(last["feats"]))}. '
         'AUC 0.5 = guessing, 1 = perfect.</p>')
    c1, c2 = st.columns(2)
    with c1:
        fig = go.Figure(go.Heatmap(z=last["K"], colorscale="Blues", showscale=False))
        fig.update_layout(**_layout(280, yaxis=dict(autorange="reversed")))
        st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_K")
        html('<p class="cc-note">Quantum kernel of the first 40 training genes (dark = similar states).</p>')
    with c2:
        circuit_block(last["circ"], f"{kp}_fm")


def gene_group_panel(ds: Dataset, tab: pd.DataFrame | None, frame: int, b0: float, offset: int = 0,
                     kp: str = "qgroup") -> None:
    html('<p class="cc-note">Choose k genes that are both active and close together in 3D: a candidate '
         'transcription hub. One qubit per candidate gene; a penalty keeps exactly k chosen.</p>')
    if tab is None or not len(tab) or "score" not in tab.columns:
        html('<p class="cc-note">No genes in view.</p>')
        return
    c1, c2, c3 = st.columns(3)
    n_c = c1.slider("Candidates (most active genes) → qubits", 6, 18, 12, key=f"{kp}_n")
    k = c2.slider("Genes to pick (k)", 2, 6, 4, key=f"{kp}_k")
    together = c3.slider("Weight on closeness", 0.0, 3.0, 1.0, 0.1, key=f"{kp}_w")
    cand = tab.dropna(subset=["score"]).sort_values("score", ascending=False)
    if "bead" not in cand.columns:
        html('<p class="cc-note">Gene positions on the fold are not available here.</p>')
        return
    cand = cand.assign(bead=cand["bead"] + int(offset))
    cand = cand[(cand["bead"] >= 0) & (cand["bead"] < ds.n)].drop_duplicates("bead").head(n_c)
    if len(cand) < k + 1:
        html('<p class="cc-note">Too few genes on the fold in view.</p>')
        return
    if st.button("Pick the group on the quantum simulator", key=f"{kp}_go", icon=":material/memory:"):
        x = ds.frames[frame][cand["bead"].to_numpy(int)]
        D = np.linalg.norm(x[:, None] - x[None], axis=-1)
        sigma = float(np.median(D[np.triu_indices(len(D), 1)])) or 1.0
        P = np.exp(-(D / sigma) ** 2)
        q = PR.gene_group_qubo(cand["score"].to_numpy(float), P, k, together, cand["gene"].tolist())
        e = q.energies()
        ex = QB.exact(q, e)
        r = sim.qaoa(e, p=3, shots=4096, objective="cvar", maxiter=80, seed=11)
        ss[kp] = {"genes": cand["gene"].tolist(), "q": q, "e": e, "ex": ex, "r": r}
    last = ss.get(kp)
    if not last:
        return
    q, r, ex = last["q"], last["r"], last["ex"]
    pick_q = [g for g, b in zip(last["genes"], QB.to_bits(r.best, q.n)) if b]
    pick_x = [g for g, b in zip(last["genes"], ex.bits) if b]
    readout([("QAOA's group", esc(", ".join(pick_q)) or "—", "simulated quantum"),
             ("Exact optimum", esc(", ".join(pick_x)), "classical"),
             ("QAOA found the optimum", "yes" if r.hit else "no", f"chance {100 * r.p_optimal:.2f} %")])
    energy_bars(r.probs, last["e"], q.n, kp)


# ======================================================================================
# 6 · Quantum walk
# ======================================================================================
def _contacts_from_coords(x: np.ndarray, b0: float) -> np.ndarray:
    from chronocell import domains as DOM
    ci, cj, cm = DOM.proximity_contacts(x, b0)
    return QW.adjacency(ci, cj, cm, len(x), observed_expected=False)


def walk_panel(ds: Dataset, before: np.ndarray, after: np.ndarray | None, b0: float, labels: tuple[str, str],
               kp: str = "qwalk") -> None:
    html('<p class="cc-note">A walker starts on one bin and moves along contacts. A <b>quantum walk</b> moves as a '
         'wave: it spreads faster (distance grows with time, not its square root) and interferes with itself, so it '
         'is sensitive to the shape of the whole network. Shown next to the classical random walk on the same '
         'contacts. Contacts here come from 3D proximity in each frame.</p>')
    n = len(before)
    step = max(1, int(math.ceil(n / 300)))
    c1, c2 = st.columns(2)
    start = c1.slider("Start bin", 0, n - 1, n // 2, key=f"{kp}_start")
    t_max = c2.slider("Duration", 2.0, 60.0, 20.0, 1.0, key=f"{kp}_t")
    if st.button("Run both walks", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        sel = np.arange(0, n, step)
        s0 = int(np.argmin(np.abs(sel - start)))
        Ab = _contacts_from_coords(before[sel], b0 * step ** (1 / 3))
        wb = QW.walks(Ab, s0, t_max, 60)
        wa = None
        if after is not None:
            wa = QW.walks(_contacts_from_coords(after[sel], b0 * step ** (1 / 3)), s0, t_max, 60)
        ss[kp] = {"wb": wb, "wa": wa, "sel": sel}
    last = ss.get(kp)
    if not last:
        return
    wb, wa, sel = last["wb"], last["wa"], last["sel"]
    on_ds = len(before) == ds.n
    mb = [round(_mb(ds, i), 2) for i in sel] if on_ds else [int(i) for i in sel]
    readout([("Position qubits", f"{wb.qubits}", f"{len(sel)} bins"),
             ("Spread at the end<small>quantum · classical</small>", f"{wb.spread()[-1]:.1f} · {wb.spread('classical')[-1]:.1f}",
              "bins")] + ([("Change in where it goes<small>quantum · classical, 0-1</small>",
                            f"{QW.compare(wb, wa)['quantum']['tv_distance']:.3f} · {QW.compare(wb, wa)['classical']['tv_distance']:.3f}", "")]
                          if wa is not None else []))
    shows = [(labels[0], wb)] + ([(labels[1], wa)] if wa is not None else [])
    cols = st.columns(len(shows))
    for col, (lab, w) in zip(cols, shows):
        for which in ("quantum", "classical"):
            P = w.quantum if which == "quantum" else w.classical
            fig = go.Figure(go.Heatmap(z=np.sqrt(P), x=mb, y=w.times, colorscale="Blues" if which == "quantum" else "Oranges",
                                       showscale=False, hovertemplate="%{x} Mb · t %{y:.1f}<extra></extra>"))
            fig.update_layout(**_layout(220, xaxis=dict(title="Mb" if on_ds else "bead"), yaxis=dict(title="time")))
            col.markdown(f'<p class="cc-eyebrow">{esc(lab)} · {which} walk</p>', unsafe_allow_html=True)
            col.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_{lab[:6]}_{which}")


def variant_set_panel(kp: str = "qvset") -> None:
    html('<p class="cc-note">From the variants ranked by the engine above, choose k to follow up: high impact, '
         'little overlap in the genes they affect. One qubit per variant.</p>')
    last = ss.get("ve_last")
    if not last or not last.get("results") or len(last["results"]) < 3:
        html('<p class="cc-note">Run "Variant impact engine v2" above on a file with at least three variants first.</p>')
        return
    from chronocell import pipelines as PL
    rank = PL.rank_variants(last["results"])
    scores = rank["score"].to_numpy(float) if "score" in rank.columns else np.arange(len(rank), 0, -1, dtype=float)
    names = rank[rank.columns[0]].astype(str).tolist()
    genes = []
    for r in last["results"]:
        g = r.get("genes") if isinstance(r, dict) else None
        genes.append({x.get("gene") for x in g} if isinstance(g, list) else set())
    n = min(len(scores), 16)
    O = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            u = genes[i] | genes[j] if i < len(genes) and j < len(genes) else set()
            O[i, j] = len(genes[i] & genes[j]) / len(u) if u else 0.0
    k = st.slider("Variants to pick", 1, max(1, n - 1), min(3, n - 1), key=f"{kp}_k")
    q = PR.variant_set_qubo(scores[:n], O, k, 1.0, names[:n])
    if st.button("Pick on the quantum simulator", key=f"{kp}_go", icon=":material/memory:"):
        e = q.energies()
        ss[kp] = {"r": sim.qaoa(e, p=3, shots=2048, objective="cvar", maxiter=80, seed=13), "ex": QB.exact(q, e), "names": names[:n]}
    res = ss.get(kp)
    if res:
        pick = lambda bits: ", ".join(nm for nm, b in zip(res["names"], bits) if b)
        readout([("QAOA's set", esc(pick(QB.to_bits(res["r"].best, q.n))), "simulated quantum"),
                 ("Exact optimum", esc(pick(res["ex"].bits)), "classical")])


# ======================================================================================
# 7 · Similarity (swap test)
# ======================================================================================
def _profile(x: np.ndarray, bins: int) -> np.ndarray:
    """Radial distance of each stretch of the fold from its centre, in `bins` equal stretches (non-negative)."""
    r = np.linalg.norm(x - x.mean(axis=0), axis=1)
    edges = np.linspace(0, len(r), bins + 1).astype(int)
    return np.array([r[a:b].mean() if b > a else 0.0 for a, b in zip(edges[:-1], edges[1:])])


def swap_panel(options: list[tuple[str, Dataset]], kp: str = "qswap") -> None:
    html('<p class="cc-note">Each fold is summarised as a profile (how far each stretch sits from the fold\'s centre), '
         'written into the amplitudes of a few qubits. The <b>swap test</b> (an extra qubit, a Hadamard, controlled '
         'swaps, a Hadamard) reads 0 with probability ½ + ½·overlap². For non-negative profiles that overlap is the '
         'squared cosine similarity, so the circuit measures a classical quantity: it is here to show how a quantum '
         'computer compares two states, with its shot noise.</p>')
    if len(options) < 2:
        html('<p class="cc-note">Needs two structures (as Compare).</p>')
        return
    names = [o[0] for o in options]
    c1, c2, c3 = st.columns(3)
    a = c1.selectbox("State A", range(len(names)), format_func=lambda i: names[i], key=f"{kp}_a")
    b = c2.selectbox("State B", range(len(names)), index=min(1, len(names) - 1), format_func=lambda i: names[i], key=f"{kp}_b")
    q = c3.slider("Qubits per state", 2, 6, 5, key=f"{kp}_q", help="Profile length 2^q.")
    shots = st.select_slider("Shots", [100, 500, 2000, 10000], 2000, key=f"{kp}_shots")
    da, db = options[a][1], options[b][1]
    g0, g1 = max(da.bin0, db.bin0), min(da.bin0 + da.n, db.bin0 + db.n)
    if da.chrom.name != db.chrom.name or g1 - g0 < 2 ** q:
        html('<p class="cc-note">The two structures must cover the same DNA.</p>')
        return
    if not st.button("Run the swap test", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        if ss.get(kp, {}).get("key") != (da.key, db.key, q, shots):
            return
    else:
        pa = _profile(da.frames[0][g0 - da.bin0:g1 - da.bin0], 2 ** q)
        pb = _profile(db.frames[0][g0 - db.bin0:g1 - db.bin0], 2 ** q)
        va, vb = KQ.amplitude_vector(pa, q), KQ.amplitude_vector(pb, q)
        pear = float(np.corrcoef(pa, pb)[0, 1]) if np.std(pa) > 0 and np.std(pb) > 0 else float("nan")
        ss[kp] = {"key": (da.key, db.key, q, shots), "res": KQ.swap_test(va, vb, int(shots), seed=0), "pear": pear,
                  "circ": KQ.swap_test_circuit(va, vb)}
    last = ss[kp]
    res = last["res"]
    readout([("Overlap from the circuit", f"{res['overlap_estimate']:.3f}", f"± {res['standard_error']:.3f} ({shots} shots)"),
             ("Exact overlap (cos²)", f"{res['overlap_exact']:.4f}", "classical"),
             ("Pearson r of the profiles", f"{last['pear']:+.3f}", "classical"), ("Qubits", f"{res['qubits']}", "")])
    circuit_block(last["circ"], f"{kp}_c")


# ======================================================================================
# Scaling and hardware notes
# ======================================================================================
def scaling_panel(kp: str = "qscale") -> None:
    df = pd.DataFrame(PR.scaling_rows())
    fig = go.Figure()
    palette = [T.ACCENT, T.TERRACOTTA, T.OCHRE, T.VIOLET, T.INK, T.MUTED]
    for colr, (name, g) in zip(palette, df.groupby("problem", sort=False)):
        fig.add_trace(go.Scatter(x=g["size"], y=g["qubits"], name=name, mode="lines+markers", line=dict(color=colr)))
    fig.add_hline(y=sim.MAX_QUBITS, line=dict(color=T.MUTED, dash="dash"),
                  annotation_text=f"this simulator: {sim.MAX_QUBITS} qubits", annotation_font_size=10)
    fig.update_layout(**_layout(320, showlegend=True, legend=dict(orientation="h", y=1.15), xaxis=dict(type="log", title="problem size"),
                                yaxis=dict(type="log", title="qubits")))
    st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_fig")
    html('<p class="cc-note">A simulator stores 2<sup>n</sup> numbers: 20 qubits take 16 MB, 30 qubits 16 GB, 40 qubits '
         '16 TB. Real quantum computers do not have that limit, but today\'s are noisy: the "survival" readouts show how '
         'quickly errors swamp circuits of the depth these problems need. A whole chromosome (thousands of bins) is far '
         'beyond any machine today; the classical methods in the other workspaces handle it in seconds.</p>')
    st.dataframe(df.assign(statevector=df["statevector_bytes"].map(_bytes)).drop(columns=["statevector_bytes"]),
                 hide_index=True, width="stretch", key=f"{kp}_tab")


def _bytes(b: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if b < 1024:
            return f"{b:.0f} {unit}"
        b /= 1024
    return f"{b:.0e} EB"


def hardware_panel() -> None:
    html('<p class="cc-note"><b>Running on a real quantum computer.</b> Every circuit here can be downloaded as '
         'OpenQASM 2.0. To run one on real hardware: create a free account on IBM Quantum Platform, open the '
         'Composer, import the .qasm file and submit it to a device (or load it in Qiskit with '
         '<code>QuantumCircuit.from_qasm_str</code>). The app never asks for or stores your account details. Expect the '
         'results to be dominated by noise for circuits with more than a few hundred two-qubit gates.</p>')
    p = VALIDATION / "results_quantum_crosscheck.json"
    if p.exists():
        try:
            r = json.loads(p.read_text(encoding="utf-8"))
            html(f'<p class="cc-note">Cross-check: the simulator\'s states were compared with Qiskit {esc(r.get("qiskit", ""))} '
                 f'on {len(r.get("circuits", []))} exported circuits; worst state fidelity '
                 f'{min(c["fidelity"] for c in r["circuits"]):.12f}.</p>')
        except (OSError, ValueError, KeyError):
            pass


# ======================================================================================
# Page 07
# ======================================================================================
def render(ds: Dataset, options: list[tuple[str, Dataset]], baseline: Dataset | None, b0: float, frame: int) -> None:
    banner(f"<b>Experimental · {esc(LABEL)}.</b> Each tool turns a ChronoCell question into a quantum algorithm and "
           "solves it on a simulator, next to the classical answer. No quantum speed-up is claimed: at these sizes the "
           "classical methods are faster.", "info")
    standing()
    prob = st.segmented_control("Problem", list(PROBLEMS), default="tad", required=True, key="qlab_problem",
                                format_func=lambda k: PROBLEMS[k][0])
    prob = prob or "tad"
    html(f'<p class="cc-purpose">{esc(PROBLEMS[prob][1])}</p>')
    if prob == "tad":
        tad_panel(ds, 0, ds.n, "qlab_tad")
    elif prob == "lattice":
        lattice_panel(ds, 0, min(ds.n, 400), ds.frames[frame][:min(ds.n, 400)], "qlab_lat")
    elif prob == "drug":
        drug_panel(ds, baseline, b0, frame, "qlab_drug")
    elif prob == "vqe":
        vqe_panel("qlab_vqe")
    elif prob == "qsvm":
        if ss.get("qlab_gene_table") is None:
            html('<p class="cc-note">Open <b>05 Genes</b> once (and add measured expression there): this panel uses the '
                 'gene table of the region shown there.</p>')
        qsvm_panel(ss.get("qlab_gene_table"), "qlab_qsvm")
    elif prob == "walk":
        after = ds.frames[-1] if ds.n_frames > 1 else None
        walk_panel(ds, ds.frames[0], after, b0, (ds.frame_labels[0] if ds.frame_labels else "frame 1",
                                                 ds.frame_labels[-1] if ds.n_frames > 1 else ""), "qlab_walk")
    elif prob == "swap":
        swap_panel(options, "qlab_swap")
    elif prob == "mol":
        molecules_panel("qlab_mol")
    elif prob == "safety":
        safety_panel("qlab_safe")
    elif prob == "admet":
        from ui import quantum_extra as QX
        QX.admet_panel("qlab_admet")
    elif prob == "noise":
        from ui import quantum_extra as QX
        QX.mitigation_panel("qlab_noise")
    else:
        docking_panel("qlab_dock")
    with st.expander("How big can these problems get?", expanded=False):
        scaling_panel()
    with st.expander("Run a circuit on a real quantum computer", expanded=False):
        hardware_panel()


# ======================================================================================
# Hooks used by the other workspaces (each inside its own expander, Research mode only)
# ======================================================================================
def structure_hooks(ds: Dataset, lo: int, hi: int, coords: np.ndarray) -> None:
    """01 3D structure → 07 Quantum."""
    t1, t2 = st.tabs(["Domain walls (QAOA)", "Lattice fold (QAOA)"])
    with t1:
        tad_panel(ds, lo, hi, "q3d_tad")
    with t2:
        lattice_panel(ds, lo, hi, coords, "q3d_lat")


def dynamics_hooks(ds: Dataset, b0: float, frame: int, traj=None) -> None:
    """02 4D dynamics → 06 Quantum: the walk compares the first and last frame of the trajectory shown."""
    t1, t2 = st.tabs(["Quantum walk", "Variant set (QAOA)"])
    with t1:
        if traj is not None and len(traj.frames) > 1:
            walk_panel(ds, traj.frames[0], traj.frames[-1], b0, (str(traj.labels[0]), str(traj.labels[-1])), "q4d_walk")
        else:
            after = ds.frames[-1] if ds.n_frames > 1 else None
            walk_panel(ds, ds.frames[0], after, b0, (ds.frame_labels[0] if ds.frame_labels else "first frame",
                                                     ds.frame_labels[-1] if ds.n_frames > 1 else ""), "q4d_walk")
    with t2:
        variant_set_panel("q4d_vset")


def genes_hooks(ds: Dataset, tab: pd.DataFrame, frame: int, b0: float, offset: int) -> None:
    """05 Genes → Quantum."""
    ss["qlab_gene_table"] = tab
    t1, t2 = st.tabs(["Gene classifier (QSVM)", "Gene group (QAOA)"])
    with t1:
        qsvm_panel(tab, "qg_qsvm")
    with t2:
        gene_group_panel(ds, tab, frame, b0, offset, "qg_group")


def drug_hooks(ds: Dataset, baseline: Dataset | None, b0: float, frame: int) -> None:
    """04 Drug lab → Quantum: two original tabs, then four drug-focused ones (molecules, heart safety, docking,
    ADMET profile)."""
    t1, t2, t3, t4, t5, t6 = st.tabs(["Drug combination (QAOA)", "Molecule energy (VQE)", "Drug molecules (active-space VQE)",
                                      "Heart safety (quantum kernel)", "Docking (QAOA max clique)",
                                      "ADMET profile (quantum kernel)"])
    with t1:
        drug_panel(ds, baseline, b0, frame, "qd_drug")
    with t2:
        vqe_panel("qd_vqe")
    with t3:
        molecules_panel("qd_mol")
    with t4:
        safety_panel("qd_safe")
    with t5:
        docking_panel("qd_dock")
    with t6:
        from ui import quantum_extra as QX
        QX.admet_panel("qd_admet")


# ======================================================================================
# Drug-focused quantum tabs (October 2026): molecules, heart safety, docking
# ======================================================================================
CPK = {"H": "#E8E8E8", "C": "#404040", "N": "#3050F8", "O": "#FF0D0D", "F": "#90E050", "Cl": "#1FF01F", "Br": "#A62929",
       "I": "#940094", "S": "#E0C030", "P": "#FF8000", "Li": "#CC80FF", "B": "#FFB5B5", "Be": "#C2FF00", "Na": "#AB5CF2"}


def molecule_figure(el: list, xyz: np.ndarray, bonds: list | None = None, height: int = 320, title: str = "",
                    extra: list | None = None) -> go.Figure:
    """Ball-and-stick 3D view; bonds from the list or from distances (< 1.25 x covalent sum)."""
    xyz = np.asarray(xyz, float)
    if bonds is None:
        rad = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39, "P": 1.07,
               "Li": 1.28, "B": 0.84, "Be": 0.96}
        bonds = [(i, j, 1) for i in range(len(el)) for j in range(i + 1, len(el))
                 if np.linalg.norm(xyz[i] - xyz[j]) < 1.25 * (rad.get(el[i], 0.8) + rad.get(el[j], 0.8))]
    fig = go.Figure()
    bx, by, bz = [], [], []
    for i, j, _ in bonds:
        bx += [xyz[i, 0], xyz[j, 0], None]
        by += [xyz[i, 1], xyz[j, 1], None]
        bz += [xyz[i, 2], xyz[j, 2], None]
    fig.add_trace(go.Scatter3d(x=bx, y=by, z=bz, mode="lines", line=dict(color="#8A8A8A", width=5), hoverinfo="skip"))
    fig.add_trace(go.Scatter3d(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], mode="markers", text=el, hoverinfo="text",
                               marker=dict(size=[5 if e == "H" else 9 for e in el], color=[CPK.get(e, "#FF1493") for e in el],
                                           line=dict(color="#202020", width=0.5))))
    for tr in extra or []:
        fig.add_trace(tr)
    fig.update_layout(height=height, margin=dict(l=0, r=0, t=24 if title else 0, b=0), showlegend=False,
                      paper_bgcolor="rgba(0,0,0,0)", title=dict(text=title, font=dict(size=12)),
                      scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False), aspectmode="data"))
    return fig


def _drug_standing(part: str, label: str) -> None:
    p = VALIDATION / "results_qdrug.json"
    try:
        r = json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
    except (OSError, ValueError):
        r = None
    if not r or part not in r.get("overall", {}):
        banner(f"<b>{esc(LABEL.capitalize())}.</b> {esc(label)}: its pre-registered test has not been run yet.", "info")
        return
    banner(f"<b>{esc(LABEL.capitalize())}.</b> {esc(label)} on held-out data: <b>{'pass' if r['overall'][part] else 'fail'}</b> "
           "(details in validation/RESULTS.md).", "info" if r["overall"][part] else "warn")


@st.cache_data(show_spinner=False, max_entries=64)
def _run_molecule(name: str, scale: float, ne: int, no: int, method: str = "uccsd") -> dict:
    from chronocell.quantum import molecules as MO
    r = MO.run(name, active=(ne, no), scale=scale, method=method)
    n_el = sum(MO.Z[a] for a, _ in r.atoms) - MO.LIBRARY[name][2]
    return {"e_hf": r.e_hf, "e_cas": r.e_cas_fci, "e_vqe": r.vqe.energy, "err": r.vqe.error, "qubits": r.n_qubits,
            "params": r.vqe.parameters, "electrons": r.vqe.electrons, "gap": r.homo_lumo_gap, "eps": r.orbital_energies.tolist(),
            "n_occ": n_el // 2, "atoms": r.atoms, "nbf": r.n_basis, "seconds": r.seconds, "history": r.vqe.history,
            "hf_converged": r.hf_converged, "method": method, "lowest_spin": r.lowest_spin, "e_singlet": r.e_singlet,
            "adapt_energies": list(getattr(r.vqe, "energies", [])), "adapt_ops": list(getattr(r.vqe, "operators", [])),
            "spin": getattr(r.vqe, "spin", None), "starts": getattr(r.vqe, "starts", None),
            "escapes": getattr(r.vqe, "escapes", None)}


def molecules_panel(kp: str = "qmol") -> None:
    from chronocell.quantum import molecules as MO
    html("<p class=\"cc-note\">The chemistry behind drug binding happens in small groups of atoms: an amine's nitrogen, a "
         "carbonyl, a nitrile, an O-H. This tool computes such fragments from first principles (STO-3G basis, "
         "Hartree–Fock), picks the most important electrons and orbitals (the <b>active space</b>), writes them onto "
         "qubits and runs <b>VQE</b> on them. The exact answer in the same space (FCI) is the reference. The integrals were "
         "checked against textbook values and against OpenFermion's independent reference data.</p>")
    _drug_standing("q6", "Molecule energies (Q6)")
    round2_standing(("q6b", "q6c", "q6d"))
    names = [n for n in MO.LIBRARY if n != "H2"]
    c1, c2, c3 = st.columns(3)
    name = c1.selectbox("Molecule", names, index=names.index("CH2O"), key=f"{kp}_name",
                        format_func=lambda n: f"{n} · {MO.LIBRARY[n][0]}")
    scale = c2.slider("Bond length (× equilibrium)", 0.7, 2.5, 1.0, 0.05, key=f"{kp}_scale",
                      help="Stretching bonds makes electrons more correlated: harder for the circuit, closer to bond breaking.")
    space = c3.selectbox("Active space (electrons, orbitals)", ["2, 2", "4, 4", "6, 6"], index=1, key=f"{kp}_space",
                         help="Qubits = 2 × orbitals. 6, 6 = 12 qubits.")
    ne, no = (int(v) for v in space.split(","))
    method = st.radio("Circuit", ["uccsd", "adapt", "adapt_multi"], horizontal=True, key=f"{kp}_method",
                      format_func=lambda m: {"uccsd": "UCCSD (fixed circuit)",
                                             "adapt": "ADAPT-VQE (grows the circuit one step at a time)",
                                             "adapt_multi": "ADAPT-VQE + escape, 4 starts (Gate Q6d)"}[m],
                      help="ADAPT-VQE adds, one at a time, the excitation that lowers the energy most, and can reuse "
                           "them; it reached chemical accuracy where the fixed UCCSD circuit did not (stretched bonds). "
                           "It takes longer: up to a few minutes at 12 qubits. The Gate Q6d version also steps off flat "
                           "spots and starts from four arrangements of the electrons, keeping the best run that ends with "
                           "paired spins: within chemical accuracy on 14 of 14 new stretched molecules, but four runs "
                           "take seconds at 8 qubits and up to about half an hour at 12.")
    key = (name, float(scale), ne, no, method)
    if st.button("Run VQE on the quantum simulator", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        try:
            with st.spinner("Integrals, Hartree–Fock, qubit Hamiltonian, VQE…"):
                ss[f"{kp}_last"] = {"key": key, "res": _run_molecule(name, float(scale), ne, no, method)}
        except ValueError as exc:
            warning_card("This active space does not fit the molecule", str(exc))
            return
    last = ss.get(f"{kp}_last")
    if not last or last["key"] != key:
        return
    r = last["res"]
    readout([("Hartree–Fock", f"{r['e_hf']:.5f}", "Ha"), ("Exact in the active space", f"{r['e_cas']:.5f}", "Ha"),
             ("VQE", f"{r['e_vqe']:.5f}", f"Ha · error {1e3 * r['err']:+.3f} mHa"),
             ("Chemical accuracy", "yes" if abs(r["err"]) <= MO.CHEMICAL_ACCURACY else "no", "≤ 1.6 mHa"),
             ("Qubits · parameters", f"{r['qubits']} · {r['params']}", f"{r['nbf']} basis functions"),
             ("HOMO–LUMO gap", f"{27.2114 * r['gap']:.1f}", "eV")])
    extra = ""
    if r.get("starts"):
        extra = (f' Starts run: {r["starts"]}; escapes from flat spots: {r["escapes"]}; spin of the result S(S+1) = '
                 f'{r["spin"]:.2f} (0 = paired; between 0 and 2 = a mix, usual where a bond is nearly broken).')
    html(f'<p class="cc-note">Correlation energy captured by the active space: {1e3 * (r["e_cas"] - r["e_hf"]):+.1f} mHa. '
         f'Electrons in the VQE state: {r["electrons"]:.3f}. {r["seconds"]:.1f} s.{extra}</p>')
    if r.get("lowest_spin", 0.0) > 1e-3:
        banner(f"<b>The exact lowest state here is not a singlet</b> (S(S+1) = {r['lowest_spin']:.1f}: unpaired "
               f"electrons, as when a bond is stretched to breaking). VQE starts from paired electrons and targets the "
               f"lowest singlet, {r['e_singlet']:.5f} Ha: its error against that state is "
               f"{1e3 * (r['e_vqe'] - r['e_singlet']):+.3f} mHa.", "info")
    if r.get("method") in ("adapt", "adapt_multi") and r.get("adapt_energies"):
        err = np.maximum(np.abs(np.array(r["adapt_energies"]) - r["e_cas"]) * 1e3, 1e-4)
        fig = go.Figure(go.Scatter(x=np.arange(1, len(err) + 1), y=err, mode="lines+markers",
                                   line=dict(color=T.ACCENT), marker=dict(size=5),
                                   hovertext=r["adapt_ops"], hoverinfo="text+y"))
        fig.add_hline(y=1e3 * MO.CHEMICAL_ACCURACY, line=dict(color=T.TERRACOTTA, dash="dot"),
                      annotation_text="chemical accuracy", annotation_position="top right")
        fig.update_layout(**_layout(260, xaxis=dict(title="operators in the circuit"),
                                    yaxis=dict(title="error vs exact (mHa)", type="log")))
        st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_adapt")
        html('<p class="cc-note">How ADAPT-VQE grows its circuit: each step adds the excitation with the steepest energy '
             'gradient (hover for which orbitals it moves electrons between) and re-optimises every angle.</p>')
    c1, c2 = st.columns(2)
    with c1:
        el = [a for a, _ in r["atoms"]]
        xyz = np.array([p for _, p in r["atoms"]], float)
        st.plotly_chart(molecule_figure(el, xyz, title=f"{name} at {scale:g} × equilibrium"), theme=None, width="stretch",
                        config=T.PLOT_CONFIG, key=f"{kp}_mol3d")
    with c2:
        eps = np.array(r["eps"]) * 27.2114
        fig = go.Figure()
        for k, e in enumerate(eps):
            fig.add_trace(go.Scatter(x=[0, 1], y=[e, e], mode="lines",
                                     line=dict(color=T.ACCENT if k < r["n_occ"] else T.TERRACOTTA, width=3),
                                     hovertext=f"orbital {k + 1}: {e:.1f} eV", hoverinfo="text"))
        fig.update_layout(**_layout(320, xaxis=dict(visible=False), yaxis=dict(title="orbital energy (eV)")))
        st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_levels")
        html('<p class="cc-note">Orbital energies: blue occupied, red empty. The gap between the highest occupied and the '
             'lowest empty orbital is a rough guide to how reactive a group is.</p>')
    slow = method == "adapt_multi" and no == 6
    if st.button("Stretch the bond: 8-point curve", key=f"{kp}_curve", disabled=slow,
                 help="Not offered for the Gate Q6d method at 12 qubits (it could take hours)." if slow else None):
        rows = []
        bar = st.progress(0.0)
        for k, s_ in enumerate(np.round(np.linspace(0.8, 2.2, 8), 2)):
            try:
                x = _run_molecule(name, float(s_), ne, no, method)
                rows.append({"scale": float(s_), "Hartree–Fock": x["e_hf"], "Exact (active space)": x["e_cas"], "VQE": x["e_vqe"]})
            except ValueError:
                pass
            bar.progress((k + 1) / 8)
        bar.empty()
        ss[f"{kp}_curve_rows"] = {"key": (name, ne, no, method), "rows": rows}
    cur = ss.get(f"{kp}_curve_rows")
    if cur and cur["key"] == (name, ne, no, method) and cur["rows"]:
        df = pd.DataFrame(cur["rows"])
        fig = go.Figure()
        for col, colr, dash in (("Hartree–Fock", T.MUTED, "dot"), ("Exact (active space)", T.INK, "solid"), ("VQE", T.ACCENT, "dash")):
            fig.add_trace(go.Scatter(x=df["scale"], y=df[col], name=col, mode="lines+markers", line=dict(color=colr, dash=dash)))
        fig.update_layout(**_layout(300, showlegend=True, legend=dict(orientation="h", y=1.15),
                                    xaxis=dict(title="bond length (× equilibrium)"), yaxis=dict(title="energy (Ha)")))
        st.plotly_chart(fig, theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_curvefig")


@st.cache_resource(show_spinner="Downloading TDC hERG (MD5-checked) and training the models…")
def _safety_model():
    from chronocell.quantum import safety as SF
    return SF.train()


def safety_panel(kp: str = "qsafe") -> None:
    from chronocell import drug_info as DI
    html("<p class=\"cc-note\">Many drugs fail because they block <b>hERG</b>, a heart potassium channel, which can disturb "
         "heart rhythm (several HDAC inhibitors carry ECG warnings). Each molecule is described by 17 numbers read from its "
         "structure (size, polarity, H-bond donors and acceptors, rings, charge, basic amines...). A <b>quantum-kernel "
         "support-vector machine</b> (the molecule's numbers set the angles of an 8-qubit feature map) learns blockers from "
         "648 measured compounds (TDC hERG), next to two classical models.</p>")
    banner("<b>A screen, not a safety assessment.</b> The models learned from a few hundred compounds; the Drug lab's "
           "epigenetic drugs are not among them. Real safety comes from patch-clamp assays and clinical ECG monitoring.", "warn")
    _drug_standing("q5", "Heart-safety classifier (Q5)")
    c1, c2 = st.columns([1, 1.4])
    pick = c1.selectbox("Drug", ["(type a SMILES)"] + sorted({n for c in DI.CLASSES.values() for n in c.lookup}),
                        key=f"{kp}_pick")
    smi = c2.text_input("SMILES", key=f"{kp}_smiles", placeholder="e.g. ONC(=O)CCCCCCC(=O)Nc1ccccc1 (vorinostat)",
                        help="Picked drugs are looked up on PubChem (only the name is sent).")
    if st.button("Screen for hERG blocking", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        try:
            if pick != "(type a SMILES)" and not smi:
                p = DI.properties(pick)
                smi = p.get("IsomericSMILES") or p.get("SMILES") or p.get("ConnectivitySMILES")
            if not smi:
                warning_card("No molecule", "Pick a drug or type a SMILES.")
                return
            model = _safety_model()
            ss[f"{kp}_last"] = {"name": pick if pick != "(type a SMILES)" else "your molecule", "smiles": smi,
                                "res": model.predict(smi)}
        except Exception as exc:
            warning_card("Could not screen this molecule", str(exc))
            return
    last = ss.get(f"{kp}_last")
    if not last:
        return
    r = last["res"]
    d = r["descriptors"]
    readout([("Quantum-kernel SVM", esc(r["qsvm_says"]), f"score {r['qsvm_score']:+.2f}"),
             ("RBF-kernel SVM (classical)", esc(r["rbf_says"]), f"score {r['rbf_score']:+.2f}"),
             ("Logistic regression (classical)", f"{100 * r['logistic_probability']:.0f}", "% blocker"),
             ("Size · polarity", f"{d['mw']:.0f} Da · {d['tpsa']:.0f} Å²", f"{d['basic_n']} basic N, {d['aromatic_rings']} aromatic rings")])
    html(f'<p class="cc-note"><b>{esc(last["name"])}</b>: <code>{esc(last["smiles"][:120])}</code></p>')
    st.dataframe(pd.DataFrame(r["nearest"]).rename(columns={"smiles": "most similar training compounds", "blocker": "measured blocker"}),
                 hide_index=True, width="stretch", key=f"{kp}_near")


def docking_panel(kp: str = "qdock") -> None:
    from chronocell.quantum import docking as DKM
    html("<p class=\"cc-note\">Docking asks where and how a drug sits in its protein pocket. Here it is a <b>matching "
         "puzzle</b>: the drug's H-bond donors, acceptors and greasy carbons are paired with pocket hot-spots where such "
         "groups are welcome; a set of pairs that one rigid placement can satisfy at once is a <b>clique</b>, and the best "
         "clique is found by QAOA (one qubit per possible pair, 20 qubits). Poses are built from the best cliques and scored. "
         "Examples: 256 drug–protein complexes from the PoseBusters benchmark (downloaded on demand, 37 MB, MD5-checked); "
         "the drug is re-docked into its own pocket.</p>")
    _drug_standing("q7", "Docking (Q7)")
    round2_standing(("q7b",))
    if f"{kp}_ids" not in ss and not st.button("Load the PoseBusters examples", key=f"{kp}_load"):
        return
    try:
        bar = st.progress(0.0, "Downloading…")
        zl = DKM.posebusters("ligands", lambda a, b: bar.progress(min(1.0, a / b) if b else 0.0))
        zp = DKM.posebusters("proteins", lambda a, b: bar.progress(min(1.0, a / b) if b else 0.0))
        bar.empty()
    except Exception as exc:
        warning_card("Could not download the examples", str(exc))
        return
    ss[f"{kp}_ids"] = DKM.list_complexes(zl)
    cid = st.selectbox("Complex (PDB id _ ligand)", ss[f"{kp}_ids"], key=f"{kp}_cid")
    p = st.slider("QAOA depth p", 1, 6, 5, key=f"{kp}_p", help="5 with the CVaR objective is the setting Gate Q7 tested.")
    if st.button("Dock on the quantum simulator", key=f"{kp}_go", type="primary", icon=":material/memory:"):
        with st.spinner("Pocket hot-spots, interaction graph, QAOA, poses…"):
            ss[f"{kp}_run"] = DKM.dock(zl, zp, cid, p=p, objective="cvar")
    run = ss.get(f"{kp}_run")
    if not run or run.cid != cid:
        return
    if not run.usable:
        banner("The supplied protein file does not contain the chain this ligand binds (it sits more than 4 Å from every "
               "protein atom); the benchmark rule skips such complexes.", "warn")

    def ok(v):
        return "—" if not np.isfinite(v) else f"{v:.1f} Å" + (" ✓" if v <= 2 else "")
    readout([("Pairs (qubits)", f"{run.qubo.n}", f"{int(run.graph.adj.sum() // 2)} compatible"),
             ("QAOA found the best clique", "yes" if run.qaoa.hit else "no", f"chance {100 * run.qaoa.p_optimal:.2f} %"),
             ("Pose error · quantum route", ok(run.rmsd_qaoa), "RMSD to the crystal; ≤ 2 Å counts as correct"),
             ("Pose error · classical clique route", ok(run.rmsd_classical), "same graph"),
             ("Pose error · random search", ok(run.rmsd_random), "1,000 placements, same score")])
    if st.button("Polish the poses: Vina-like score and refinement", key=f"{kp}_polish",
                 help="Round 2 (Gate Q7b): the quantum route's poses are scored with an AutoDock Vina-style function "
                      "(steric, hydrophobic and H-bond terms) and the best are nudged (rigid-body refinement); random "
                      "search gets the same score and the same refinement."):
        with st.spinner("Scoring and refining poses…"):
            ss[f"{kp}_pol"] = (run.cid, DKM.polish(run))
    pol = ss.get(f"{kp}_pol")
    if pol and pol[0] == run.cid:
        po = pol[1]
        readout([("Pose error · quantum route, polished", ok(po.rmsd_qaoa), f"score {po.score_qaoa:.1f}"),
                 ("Pose error · random search, polished", ok(po.rmsd_random), f"score {po.score_random:.1f}"),
                 ("Score of the crystal pose", f"{po.score_crystal:.1f}", "lower is better")])
        html('<p class="cc-note">The Vina-like score is in kcal/mol-like units (lower is better). When the crystal pose '
             'scores better than the docked one, the search missed it; when the docked pose scores better but is far '
             'from the crystal, the scoring function is at fault.</p>')
    heavy = run.lig.heavy
    centre = run.lig.xyz[heavy].mean(0)
    near = np.linalg.norm(run.prot.xyz - centre, axis=1) < 9
    pk = run.prot.xyz[near & (run.prot.el != "H")]
    extra = [go.Scatter3d(x=pk[:, 0], y=pk[:, 1], z=pk[:, 2], mode="markers", marker=dict(size=2, color="#B0B0B0"), hoverinfo="skip")]
    if len(run.sites.kind):
        extra.append(go.Scatter3d(x=run.sites.xyz[:, 0], y=run.sites.xyz[:, 1], z=run.sites.xyz[:, 2], mode="markers",
                                  marker=dict(size=6, symbol="diamond",
                                              color=[{"A": "#C24A1E", "D": "#3340D1", "H": "#A37822"}[k] for k in run.sites.kind]),
                                  text=[f"hot-spot {k}" for k in run.sites.kind], hoverinfo="text"))
    if run.pose_qaoa is not None:
        pq = run.pose_qaoa[heavy]
        extra.append(go.Scatter3d(x=pq[:, 0], y=pq[:, 1], z=pq[:, 2], mode="markers", marker=dict(size=5, color=T.ACCENT),
                                  text=["docked (quantum route)"] * len(pq), hoverinfo="text"))
    if pol and pol[0] == run.cid and pol[1].pose_qaoa is not None:
        pp = pol[1].pose_qaoa
        extra.append(go.Scatter3d(x=pp[:, 0], y=pp[:, 1], z=pp[:, 2], mode="markers",
                                  marker=dict(size=5, color=T.TERRACOTTA, symbol="square"),
                                  text=["polished (quantum route)"] * len(pp), hoverinfo="text"))
    el = [e for e, h in zip(run.lig.el, heavy) if h]
    idx = {k: i for i, k in enumerate(np.flatnonzero(heavy))}
    bonds = [(idx[i], idx[j], o) for i, j, o in run.lig.bonds if i in idx and j in idx]
    st.plotly_chart(molecule_figure(el, run.lig.xyz[heavy], bonds, height=440,
                                    title="crystal pose (atoms) · docked (blue) · hot-spots (red acceptor, blue donor, "
                                          "ochre hydrophobic) · pocket (grey)", extra=extra),
                    theme=None, width="stretch", config=T.PLOT_CONFIG, key=f"{kp}_view")
    c1, c2 = st.columns(2)
    with c1:
        qubo_heatmap(run.qubo, kp)
    with c2:
        h, J, _ = run.qubo.ising()
        circuit_block(sim.qaoa_circuit(h, J, run.qaoa.gammas, run.qaoa.betas, f"QAOA docking {run.cid}"), f"{kp}_c")
