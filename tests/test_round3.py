"""New tools (October 2026, second batch): the noisy density-matrix simulator and zero-noise extrapolation, the ADMET
models, multi-start ADAPT-VQE with the energy-scan escape, and the 08 Scoreboard page. Each piece is checked against
exact expectations; every earlier default is left as it was."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chronocell.quantum import admet as A, molecules as MO, noisy as N, sim

APP = str(Path(__file__).resolve().parent.parent / "app.py")


# ---------------------------------------------------------------- noisy simulator and mitigation
def _bell() -> sim.Circuit:
    c = sim.Circuit(2, "bell")
    c.h(0)
    c.cx(0, 1)
    c.rz(0.3, 1)
    return c


def test_density_matrix_matches_the_statevector_without_noise_and_folding_keeps_the_ideal():
    c = _bell()
    psi = c.run()
    rho = N.run(c, 0.0, 0.0)
    assert np.allclose(rho, np.outer(psi, psi.conj()), atol=1e-12)
    f = N.fold(c, 3)
    assert len(list(f.gates())) == 3 * len(list(c.gates()))
    assert np.allclose(N.run(f, 0.0, 0.0), rho, atol=1e-10)                 # C C^dagger C = C
    with pytest.raises(ValueError):
        N.fold(c, 2)


def test_depolarizing_channel_and_readout():
    c = _bell()
    rho = N.run(c, 0.0, 0.0)
    full = N.depolarize(rho, (0, 1), 2, 1.0)
    assert np.allclose(full, np.eye(4) / 4, atol=1e-12)                    # fully depolarised two qubits
    half = N.depolarize(rho, (0,), 2, 0.5)
    assert abs(np.trace(half) - 1) < 1e-12 and np.allclose(half, half.conj().T)
    p = N.readout_probabilities(np.diag([1.0, 0, 0, 0]).astype(complex), 2, 0.1)
    assert np.allclose(p, [0.81, 0.09, 0.09, 0.01])


def test_zne_and_symmetry_verification_recover_a_molecule_energy():
    ints = MO.integrals(MO.LIBRARY["LiH"][1]())
    qm = MO.qubit_hamiltonian(MO.active_space(ints, MO.rhf(ints), 2, 2))
    r = MO.vqe_uccsd(qm)
    c = N.vqe_circuit(qm, r.generators, r.thetas)
    psi = c.run()
    H = qm.H.toarray()
    assert float(np.real(psi.conj() @ H @ psi)) == pytest.approx(r.energy, abs=1e-9)   # the compiled circuit is exact
    keep = np.zeros(2 ** qm.n_qubits)
    keep[MO.sector(qm)] = 1.0
    z = N.zne(c, H, {"p1": 3e-4, "p2": 3e-3}, ideal=r.energy, keep=keep)
    raw = N.expectation(c, H, {"p1": 3e-4, "p2": 3e-3})
    assert abs(raw - r.energy) > 0.01                                       # the noise matters
    assert abs(z.error_mitigated) < MO.CHEMICAL_ACCURACY < abs(z.error_noisy)


# ---------------------------------------------------------------- ADMET
def test_admet_endpoints_kernels_and_models():
    assert len(A.ENDPOINTS) == 21 and "herg" not in A.ENDPOINTS
    assert sorted(n for g in A.GROUPS.values() for n in g) == sorted(A.ENDPOINTS)
    rng = np.random.default_rng(0)
    X = rng.normal(size=(120, 17))
    prep = A.Prep.fit(X)
    a = prep.angles(X)
    assert a.shape == (120, 8) and a.min() >= 0 and a.max() <= np.pi + 1e-12
    K = A.rbf(prep.pca(X), prep.pca(X), 0.1)
    assert np.allclose(K, K.T) and np.allclose(np.diag(K), 1.0)
    P = prep.pca(X)
    y = (P[:, 0] + 0.3 * P[:, 1] > 0).astype(int)                       # a signal inside the 8 retained components
    S = A.qstates(a, 0.05, 1)
    Kq = A.qkernel(S, S)
    assert np.allclose(np.diag(Kq), 1.0)
    s = A.fit_predict("cls", Kq[:80, :80], y[:80], Kq[80:, :80], 10.0)
    assert A.score("cls", y[80:], s) > 0.8
    yr = X[:, 0] * 2 + rng.normal(0, 0.1, 120)
    sr = A.fit_predict("reg", A.rbf(X[:80], X[:80], 0.05), yr[:80], A.rbf(X[80:], X[:80], 0.05), 0.1)
    assert A.score("reg", yr[80:], sr) > 0.8


# ---------------------------------------------------------------- ADAPT-VQE: escape and multi-start
def test_energy_scan_and_multistart_never_do_worse_than_hartree_fock_start():
    ints = MO.integrals(MO.LIBRARY["H2O"][1](1.5))
    qm = MO.qubit_hamiltonian(MO.active_space(ints, MO.rhf(ints), 4, 4))
    psi = np.zeros(2 ** qm.n_qubits)
    psi[qm.hf_index] = 1.0
    gens = MO.generalized_generators(qm)
    A_, = [G for n, G in gens if n == "2->6"]
    th, e = MO.energy_scan(qm, psi, A_, (A_ @ A_).tocsr())
    e0 = float(psi @ (qm.H @ psi))
    assert e <= e0 + 1e-12
    dets = MO.lowest_determinants(qm, 2)
    assert len(dets) == 2 and qm.hf_index not in dets
    single = MO.adapt_vqe(qm, pool="gsd", grad_tol=1e-3, max_operators=40)
    multi = MO.adapt_multistart(qm, references=2, pool="gsd", grad_tol=1e-3, max_operators=40, escape=1e-5)
    assert multi.energy <= single.energy + 1e-9 and multi.starts == 3
    assert abs(multi.energy - MO.fci_singlet(qm)) < MO.CHEMICAL_ACCURACY


# ---------------------------------------------------------------- scoreboard
def test_scoreboard_status_parsing():
    from ui import scoreboard as SB
    assert SB.status_of("pass") == "pass" and SB.status_of("fail") == "fail"
    assert SB.status_of("imaging fail, hic fail") == "fail"
    assert SB.status_of("blocked; variant engine stays a mechanism simulator") == "other"
    assert SB.status_of("pass (modest)") == "pass"
    df = SB.rows()
    assert len(df) >= 30 and set(df["status"]) <= {"pass", "fail", "other"}
    assert df["test"].str.startswith("Gate 6b").any()


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


def test_scoreboard_page_and_new_quantum_tools_render(app):
    at = app
    at.segmented_control(key="workspace").set_value("Scoreboard").run()
    assert not at.exception, [e.value for e in at.exception]
    at.segmented_control(key="workspace").set_value("Quantum lab").run()
    opts = at.segmented_control(key="qlab_problem").options
    assert any("ADMET" in o for o in opts) and any("Noise" in o for o in opts)
    at.segmented_control(key="qlab_problem").set_value("noise").run()
    at.button(key="qlab_noise_go").click().run()
    assert not at.exception, [e.value for e in at.exception]


# ---------------------------------------------------------------- the film and the Guide's pictures
def test_guide_pictures_and_film_section_render(app):
    at = app
    at.segmented_control(key="workspace").set_value("Guide").run()
    assert not at.exception, [e.value for e in at.exception]
    txt = " ".join(str(m.value) for m in at.markdown)
    assert "The tour in pictures" in txt and "Watch the film" in txt
    keys = {b.key for b in at.button}
    assert {"guide_pic_ws01", "guide_pic_ws08"} <= keys and "guide_pic_ws06" not in keys      # no button to the page itself


def test_film_numbers_come_from_the_result_files():
    import importlib.util
    from pathlib import Path
    from ui import scoreboard as SB
    spec = importlib.util.spec_from_file_location("film_data", Path(__file__).resolve().parent.parent / "motion" / "build_data.py")
    BD = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(BD)
    t = BD.tests()
    df = SB.rows()
    assert t["n"] == len(df) and t["pass"] + t["fail"] + t["other"] == t["n"]
    assert t["pass"] == int((df["status"] == "pass").sum())
    hi = BD.highlights()
    assert hi["q6d"]["within"] == len(hi["q6d"]["rows"]) == hi["q6d"]["cases"]
    assert all(r["err"] <= hi["chem_accuracy"] for r in hi["q6d"]["rows"])
    assert hi["q9"]["within"] == hi["q9"]["cases"] and hi["q7b"]["qaoa"] < hi["q7b"]["random"]
    assert all(v["chronocell"] > max(v["chromosight"], v["mustache"]) for v in hi["loops"].values())
