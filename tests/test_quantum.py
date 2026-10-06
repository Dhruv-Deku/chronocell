"""Quantum lab: every simulated-quantum component is checked against exact linear algebra or published values,
and the app's existing pages are unchanged (the lab is added, in Research mode only)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from chronocell.quantum import chem as CH, kernels as KQ, lattice as LT, problems as PR, qubo as QB, sim, walk as QW

APP = str(Path(__file__).resolve().parent.parent / "app.py")


# ---------------------------------------------------------------- simulator
def test_gates_and_circuits_match_matrices():
    c = sim.Circuit(2)
    c.h(0).cx(0, 1)
    assert np.allclose(c.run(), [2 ** -0.5, 0, 0, 2 ** -0.5])
    a = sim.Circuit(3)
    b = sim.Circuit(3)
    for q in range(3):
        a.h(q)
        b.h(q)
    a.rzz(0.7, 0, 2)
    b.cx(0, 2).rz(0.7, 2).cx(0, 2)
    assert sim.state_fidelity(a.run(), b.run()) == pytest.approx(1.0)
    rng = np.random.default_rng(0)
    for lab in ("XY", "YZX", "XIZY"):
        n = len(lab)
        prep = sim.Circuit(n)
        for q in range(n):
            prep.ry(rng.random() * 3, q).rz(rng.random() * 3, q)
        psi = prep.run()
        ex = sim.Circuit(n)
        sim.pauli_exponential(ex, lab, 0.37)
        assert sim.state_fidelity(ex.run(psi), expm(-1j * 0.37 * sim.pauli_matrix(lab)) @ psi) == pytest.approx(1.0)
    qasm = c.to_qasm()
    assert qasm.startswith("OPENQASM 2.0;") and "cx q[0],q[1];" in qasm and "measure q -> c;" in qasm


def test_walsh_coefficients_reconstruct_the_cost():
    d = np.random.default_rng(1).normal(size=32)
    w = sim.walsh_coefficients(d)
    z = 1 - 2 * sim.bits_of(np.arange(32), 5)
    rec = np.array([sum(w[S] * np.prod([z[x, i] for i in range(5) if S >> i & 1]) for S in range(32)) for x in range(32)])
    assert np.allclose(rec, d)


def test_noise_model_sampling():
    rng = np.random.default_rng(0)
    p = np.zeros(16)
    p[5] = 1.0
    assert np.all(sim.sample(p, 200, rng) == 5)
    s = sim.sample(p, 4000, rng, fidelity=0.0)
    assert len(np.unique(s)) == 16                           # fully depolarised: uniform
    s = sim.sample(p, 4000, rng, readout=0.5)
    assert np.mean(s == 5) < 0.2                             # read-out flips scramble the bits


# ---------------------------------------------------------------- QUBO and solvers
def _random_qubo(n, rng, band=None):
    Q = np.triu(rng.normal(size=(n, n)), 1)
    if band is not None:
        Q[np.abs(np.subtract.outer(np.arange(n), np.arange(n))) > band] = 0
    return QB.QUBO(rng.normal(size=n), Q, rng.normal())


def test_qubo_energies_ising_and_solvers_agree():
    rng = np.random.default_rng(3)
    q = _random_qubo(12, rng)
    e = q.energies()
    x = rng.integers(0, 2, 12)
    assert e[QB.to_index(x)] == pytest.approx(float(q.energy(x)))
    h, J, c = q.ising()
    z = 1 - 2 * x
    assert h @ z + sum(v * z[i] * z[j] for (i, j), v in J.items()) + c == pytest.approx(float(q.energy(x)))
    ex = QB.exact(q, e)
    assert QB.anneal(q, seed=0).energy == pytest.approx(ex.energy)
    assert QB.sqa(q, seed=0).energy == pytest.approx(ex.energy)
    for _ in range(30):
        qb = _random_qubo(int(rng.integers(3, 14)), rng, band=int(rng.integers(1, 5)))
        assert QB.banded_exact(qb).energy == pytest.approx(QB.exact(qb).energy)


def test_qaoa_fast_path_equals_the_gate_circuit_and_finds_the_optimum():
    rng = np.random.default_rng(4)
    q = _random_qubo(8, rng)
    e = q.energies()
    h, J, _ = q.ising()
    g, b = [0.3, 0.5], [0.6, 0.2]
    fast = sim._Engine(e, prefer_gpu=False).state(g, b)
    assert sim.state_fidelity(np.asarray(fast), sim.qaoa_circuit(h, J, g, b).run()) == pytest.approx(1.0)
    r = sim.qaoa(e, p=3, shots=2048, prefer_gpu=False)
    assert r.hit and r.p_optimal > 10 / 2 ** 8


def test_tad_qubo_equals_modularity_on_feasible_walls():
    """The second-order cut expansion is exact for every solution that respects the minimum domain size."""
    rng = np.random.default_rng(5)
    n = 30
    M = rng.poisson(5, (n, n)).astype(float)
    M = np.triu(M, 1) + np.triu(M, 1).T + np.diag(rng.poisson(20, n))
    m, span, gamma = 3, 6, 1.3
    core = (span, n - span)
    mask = np.ones(n, bool)
    exp = PR.expected_by_distance(M, span, mask)
    q, pos = PR.tad_qubo(M, core, gamma, m, span, exp, mask)
    Z = float(np.mean(exp[1:span + 1]))
    checked = 0
    for _ in range(200):
        bits = np.zeros(q.n, dtype=np.int8)
        v = int(rng.integers(0, m))
        while v < q.n:                                        # random walls at least m apart (feasible)
            bits[v] = rng.random() < 0.6
            v += m + int(rng.integers(0, 3)) if bits[v] else 1
        assert PR.feasible(bits, m)
        walls = set(PR.boundaries_from_bits(bits, pos))
        direct = 0.0
        for a in range(0, n):
            for bb in range(a + 1, min(n, a + span + 1)):
                if any(k in walls for k in range(a + 1, bb + 1)):
                    direct += (M[a, bb] - gamma * exp[bb - a]) / Z
        assert float(q.energy(bits)) == pytest.approx(direct)
        checked += 1
    assert checked > 50


def test_drug_combination_qubo_equals_its_model():
    alpha, beta = np.array([10.0, 5, 2, 1]), np.array([-3.0, 1, 0, 0])
    eta = np.zeros((4, 4))
    eta[0, 1], eta[1, 3] = -4, 2
    q = PR.drug_combo_qubo(alpha, beta, eta, 3.0, list("abcd"))
    e = q.energies()
    model = {"alpha": alpha, "beta": beta, "eta": eta}
    for idx in range(256):
        d = PR.drug_doses(QB.to_bits(idx, 8), 4)
        assert e[idx] == pytest.approx(-PR.predicted_restoration(model, d) + 3.0 * d.sum())


def test_k_of_n_selection_picks_k():
    rng = np.random.default_rng(6)
    q = PR.variant_set_qubo(rng.random(7), rng.random((7, 7)) * 0.3, 3)
    assert int(QB.exact(q).bits.sum()) == 3
    q = PR.gene_group_qubo(rng.random(9), rng.random((9, 9)), 4)
    assert int(QB.exact(q).bits.sum()) == 4


# ---------------------------------------------------------------- lattice, walk
def test_lattice_fold_energies():
    W = np.zeros((6, 6))
    W[0, 5] = W[5, 0] = 1.0
    e = LT.energies(W, 2)
    k = int(np.argmin(e))
    f = LT.fold_of(k, W, 2, e[k])
    assert f.valid and (0, 5) in f.contacts and e[k] == pytest.approx(-1.0)
    straight = 0                                             # all bonds +x: no contacts, no overlap
    assert e[straight] == pytest.approx(0.0)


def test_quantum_and_classical_walks_match_matrix_exponentials():
    n = 40
    i = np.arange(n - 1)
    A = QW.adjacency(i, i + 1, np.ones(n - 1), n, observed_expected=False)
    w = QW.walks(A, 20, 10.0, 21)
    deg = A.sum(1)
    Ls = np.eye(n) - A / np.sqrt(np.outer(deg, deg))
    e = np.zeros(n)
    e[20] = 1
    assert np.allclose(np.abs(expm(-1j * Ls * w.times[7]) @ e) ** 2, w.quantum[7])
    assert np.allclose(expm(-(np.eye(n) - A / deg[None, :]) * w.times[7]) @ e, w.classical[7])
    assert np.allclose(w.quantum.sum(1), 1) and np.allclose(w.classical.sum(1), 1)
    assert w.spread()[-1] > w.spread("classical")[-1]        # ballistic vs diffusive


# ---------------------------------------------------------------- chemistry
def test_chemistry_reproduces_szabo_ostlund_and_vqe_reaches_fci():
    B = CH.BOHR_ANGSTROM
    ints = CH.integrals(("H", "H"), 1.4)
    assert ints["S"][0, 1] == pytest.approx(0.6593, abs=1e-4)
    assert ints["H"][0, 0] == pytest.approx(-1.1204, abs=1e-4)
    assert ints["eri"][0, 0, 0, 0] == pytest.approx(0.7746, abs=1e-4)
    h2 = CH.molecule("H2", 1.4 * B)
    assert h2.e_hf == pytest.approx(-1.1167, abs=2e-4)
    assert h2.e_fci == pytest.approx(-1.1373, abs=2e-4)
    assert len(h2.paulis) == 15
    assert CH.molecule("HeH+", 1.4632 * B).e_hf == pytest.approx(-2.86066, abs=2e-4)
    v = CH.vqe(h2)
    assert abs(v.error) < CH.CHEMICAL_ACCURACY
    th = np.array([0.1, -0.2, 0.3])
    U = np.eye(16)
    for (_, G), t in zip(CH.uccsd_generators(), th):
        U = expm(-1j * t * G) @ U
    psi0 = np.zeros(16)
    psi0[3] = 1
    assert sim.state_fidelity(U @ psi0, CH.uccsd_circuit(h2, th).run()) == pytest.approx(1.0)


# ---------------------------------------------------------------- kernels
def test_feature_map_state_prep_and_swap_test():
    rng = np.random.default_rng(7)
    x = rng.random(4) * 3
    assert sim.state_fidelity(KQ.feature_state(x, 1.0, 2), KQ.feature_map_circuit(x, 1.0, 2).run()) == pytest.approx(1.0)
    for q in (1, 2, 3):
        v = rng.random(2 ** q)
        v /= np.linalg.norm(v)
        c = sim.Circuit(q)
        KQ.amplitude_prep(c, v, list(range(q)))
        assert np.allclose(np.abs(c.run()), v)
    a, b = KQ.amplitude_vector(rng.random(8)), KQ.amplitude_vector(rng.random(8))
    psi = KQ.swap_test_circuit(a, b).run()
    assert float(np.sum(np.abs(psi[0::2]) ** 2)) == pytest.approx(0.5 + 0.5 * float(a @ b) ** 2)
    K = KQ.kernel(rng.random((6, 3)), None, 0.5, 2)
    assert np.allclose(np.diag(K), 1) and np.allclose(K, K.T) and np.all(np.linalg.eigvalsh(K) > -1e-9)
    X = rng.random((40, 3))
    y = (X[:, 0] > 0.5).astype(int)
    a_, = KQ.scale_features(X)
    assert set(KQ.QSVM(10.0, 0.3, 1).fit(a_, y).predict(a_)) <= {0, 1}


# ---------------------------------------------------------------- app
@pytest.fixture()
def app(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    import ui.common as C
    import ui.settings as SET
    import ui.states_panel as SP
    monkeypatch.setattr(C, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "DEMO_ROOT", tmp_path / "demo")
    monkeypatch.setattr(SET, "PATH", tmp_path / "settings.json")
    at = AppTest.from_file(APP, default_timeout=900)
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at


def test_quantum_lab_is_added_in_research_mode_and_runs_on_demand(app):
    at = app
    opts = at.segmented_control(key="workspace").options
    assert any("Quantum lab" in o for o in opts)
    assert [o for o in opts if "Quantum" not in o][:6] == [o for o in opts][:6]      # the six pages keep their order
    at.segmented_control(key="workspace").set_value("Quantum lab").run()
    assert not at.exception, [e.value for e in at.exception]
    assert "qlab_tad" not in at.session_state                                    # nothing computed before Run
    at.slider(key="qlab_tad_nb").set_value(8).run()
    at.button(key="qlab_tad_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    last = at.session_state["qlab_tad"]
    assert last["q"].n == 7 and last["ex"].energy <= last["r"].best_energy + 1e-9
    at.toggle(key="research_mode").set_value(False).run()
    assert not any("Quantum lab" in o for o in at.segmented_control(key="workspace").options)


def test_compare_keeps_its_tabs_and_adds_the_quantum_tab(app):
    at = app
    at.segmented_control(key="workspace").set_value("Compare").run()
    assert not at.exception, [e.value for e in at.exception]
    labels = [t.label for t in at.tabs]
    assert labels[:3] == ["Side by side", "Self-Math PDB State Evaluator", "Differential analysis"]
    assert "Quantum similarity (swap test)" in labels


def test_quantum_panels_run_in_the_lab_and_appear_in_the_workspaces(app):
    at = app
    assert any("07   Quantum" in e.label for e in at.expander)                     # 01 3D structure → 07
    at.segmented_control(key="workspace").set_value("Quantum lab").run()
    at.segmented_control(key="qlab_problem").set_value("lattice").run()
    at.slider(key="qlab_lat_beads").set_value(6).run()
    at.button(key="qlab_lat_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    assert at.session_state["qlab_lat"]["exact"].valid
    at.segmented_control(key="qlab_problem").set_value("vqe").run()
    at.button(key="qlab_vqe_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    assert abs(at.session_state["qlab_vqe_last"]["res"]["e_vqe"] - at.session_state["qlab_vqe_last"]["res"]["e_fci"]) < 1.6e-3
    at.segmented_control(key="qlab_problem").set_value("walk").run()
    at.button(key="qlab_walk_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    assert np.allclose(at.session_state["qlab_walk"]["wb"].quantum.sum(1), 1)
    for page, label in (("4D dynamics", "06   Quantum"), ("Drug lab", "Quantum (simulated)"), ("Genes", "Quantum (simulated)")):
        at.segmented_control(key="workspace").set_value(page).run()
        assert not at.exception, [e.value for e in at.exception]
        assert any(label in e.label for e in at.expander), page
    at.segmented_control(key="workspace").set_value("Guide").run()
    assert not at.exception and any("Quantum lab, in plain words" in m.value for m in at.markdown)
