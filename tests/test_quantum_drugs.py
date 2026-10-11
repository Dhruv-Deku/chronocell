"""Drug-class additions (October 2026): extended drug classes, drug information, the molecule engine (active-space
VQE), the SMILES descriptors, the docking clique formulation, the heart-safety model, and the app pages. Each part is
checked against exact mathematics or published values; the original Drug lab defaults are checked unchanged."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chronocell import drug_info as DI, therapy as TH
from chronocell.quantum import docking as DK, molecules as MO, molfeat as MF, qubo as QB

APP = str(Path(__file__).resolve().parent.parent / "app.py")


# ---------------------------------------------------------------- drug lab classes
def test_core_drug_classes_are_unchanged_and_the_extended_set_is_added():
    assert list(TH.DRUGS) == ["ezh2", "hdac", "bet", "ctcf"]
    assert all(d.strength == 1.0 for d in TH.DRUGS.values())
    assert len(TH.EXTRA_DRUGS) == 8 and set(TH.ALL_DRUGS) == set(TH.DRUGS) | set(TH.EXTRA_DRUGS)
    assert set(DI.CLASSES) == set(TH.ALL_DRUGS)                 # every class has an information card
    rng = np.random.default_rng(0)
    x = np.cumsum(rng.normal(0, 50, (200, 3)), 0)
    sig = rng.random(200)
    base = TH.compare_drugs(x, sig, np.ones(200, bool), 50.0, x + rng.normal(0, 5, x.shape))
    assert len(base) == 4                                       # the default ranking still covers the core four
    for k, d in TH.EXTRA_DRUGS.items():
        r = TH.simulate_treatment(x, sig, np.ones(200, bool), 50.0, k, None, 0.8, np.array([0.0, 1.0]))
        rg0, rg1 = r.metrics.rg_nm.iloc[0], r.metrics.rg_nm.iloc[-1]
        if d.direction > 0:
            assert rg1 >= rg0 - 1e-6, k                           # opening classes do not compact the fold
        elif d.direction < 0:
            assert rg1 <= rg0 + 1e-6, k


def test_drug_likeness_rules():
    p = {"MolecularWeight": "264.32", "XLogP": 1.9, "HBondDonorCount": 3, "HBondAcceptorCount": 3, "RotatableBondCount": 8,
         "TPSA": 78.4}
    dl = DI.drug_likeness(p)
    assert dl["ro5_violations"] == 0 and all(dl["rules"].values())


# ---------------------------------------------------------------- molecules
@pytest.mark.parametrize("name", ["H2O", "NH3", "CH4", "HF", "N2", "CO"])
def test_sto3g_hartree_fock_matches_szabo_ostlund(name):
    hf = MO.rhf(MO.integrals(MO.LIBRARY[name][1]()))
    assert hf.energy == pytest.approx(MO.SZABO_OSTLUND_HF[name], abs=6e-4)   # the table is printed to 1 mEh


def test_active_space_vqe_rotations_gradient_and_accuracy():
    from scipy.sparse.linalg import expm_multiply
    ints = MO.integrals(MO.LIBRARY["H2O"][1]())
    hf = MO.rhf(ints)
    qm = MO.qubit_hamiltonian(MO.active_space(ints, hf, 4, 4))
    gens = MO.uccsd_generators(qm)
    v = np.random.default_rng(1).normal(size=2 ** qm.n_qubits)
    for _, G in gens[:6]:
        assert np.allclose(MO._rot(G, (G @ G).tocsr(), 0.41, v), expm_multiply(0.41 * G, v))
    r = MO.vqe_uccsd(qm)
    assert abs(r.error) < MO.CHEMICAL_ACCURACY and r.electrons == pytest.approx(4.0)
    assert r.energy <= r.hf + 1e-9 and r.fci <= r.energy + 1e-9


def test_adapt_vqe_pool_gradients_and_accuracy():
    """ADAPT-VQE (round 2): the generalized pool holds the UCCSD generators, its gradient formula matches a finite
    difference, and the grown circuit reaches the FCI where fixed-order UCCSD stalls (HCN at 1.3x, Q6's failure)."""
    ints = MO.integrals(MO.LIBRARY["H2O"][1]())
    qm = MO.qubit_hamiltonian(MO.active_space(ints, MO.rhf(ints), 4, 4))
    ucc = {n for n, _ in MO.uccsd_generators(qm)}
    gen = MO.generalized_generators(qm)
    assert ucc <= {n for n, _ in gen} and len(gen) > len(ucc)
    psi = np.zeros(2 ** qm.n_qubits)
    psi[qm.hf_index] = 1.0
    A = dict(gen)["2->6"]
    e = lambda t: float(MO._rot(A, (A @ A).tocsr(), t, psi) @ (qm.H @ MO._rot(A, (A @ A).tocsr(), t, psi)))
    assert 2.0 * float((qm.H @ psi) @ (A @ psi)) == pytest.approx((e(1e-5) - e(-1e-5)) / 2e-5, abs=1e-6)
    ints = MO.integrals(MO.LIBRARY["HCN"][1](1.3))
    qm = MO.qubit_hamiltonian(MO.active_space(ints, MO.rhf(ints), 4, 4))
    a = MO.adapt_vqe(qm, **MO.ADAPT_SETTINGS)
    assert a.converged and abs(a.error) < MO.CHEMICAL_ACCURACY and a.electrons == pytest.approx(4.0)
    assert np.all(np.diff(a.energies) <= 1e-9)                    # each added operator never raises the energy
    assert abs(MO.vqe_uccsd(qm).error) > MO.CHEMICAL_ACCURACY       # the fixed-order circuit misses it here


def test_sector_ground_state_spin_and_lowest_singlet():
    """Q6b post hoc: stretched LiH's lowest state in the N = 2, Sz = 0 sector is a triplet; a closed-shell VQE reaches
    the lowest singlet, which fci_singlet reports. At equilibrium the two coincide."""
    ints = MO.integrals(MO.LIBRARY["LiH"][1](3.0))
    qm = MO.qubit_hamiltonian(MO.active_space(ints, MO.rhf(ints), 4, 2))
    (e0, s0), = MO.fci_states(qm, 1)
    assert s0 == pytest.approx(2.0, abs=1e-6) and e0 == pytest.approx(MO.fci(qm))
    es = MO.fci_singlet(qm)
    assert es > e0 + 1e-3
    assert MO.vqe_uccsd(qm).energy == pytest.approx(es, abs=1e-6)
    ints = MO.integrals(MO.LIBRARY["LiH"][1]())
    qm = MO.qubit_hamiltonian(MO.active_space(ints, MO.rhf(ints), 4, 2))
    assert MO.fci_singlet(qm) == pytest.approx(MO.fci(qm)) and MO.fci_states(qm, 1)[0][1] == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------- SMILES descriptors
def test_smiles_descriptors_known_values():
    a = MF.descriptors("CC(=O)Oc1ccccc1C(=O)O")                 # aspirin
    assert a["mw"] == pytest.approx(180.16, abs=0.02) and a["hbd"] == 1 and a["hba"] == 4
    assert a["tpsa"] == pytest.approx(63.6, abs=0.1) and a["rot_bonds"] == 3 and a["aromatic_rings"] == 1
    v = MF.descriptors("ONC(=O)CCCCCCC(=O)Nc1ccccc1")           # vorinostat
    assert v["mw"] == pytest.approx(264.32, abs=0.02) and v["tpsa"] == pytest.approx(78.43, abs=0.1) and v["rot_bonds"] == 8
    assert MF.descriptors("c1ccsc1")["mw"] == pytest.approx(84.14, abs=0.02)        # thiophene: no H on S
    assert MF.descriptors("c1cc[nH]c1")["hbd"] == 1
    assert MF.descriptors("CN1CCN(CC1)c1ccccc1")["basic_n"] == 1
    with pytest.raises(MF.SmilesError):
        MF.parse("C1CC")


# ---------------------------------------------------------------- docking
def test_max_weight_clique_qubo_and_pose_recovery():
    rng = np.random.default_rng(2)
    n = 9
    adj = rng.random((n, n)) < 0.5
    adj = np.triu(adj, 1)
    adj = adj | adj.T
    w = rng.random(n) + 0.5
    g = DK.Graph([(i, i) for i in range(n)], w, adj)
    q = DK.clique_qubo(g)
    e = q.energies()
    best = QB.to_bits(int(np.argmin(e)), n)
    assert DK.is_clique(g, best)
    import itertools
    brute = max(sum(w[list(c)]) for k in range(1, n + 1) for c in itertools.combinations(range(n), k)
                if all(adj[a, b] for a, b in itertools.combinations(c, 2)))
    assert -e.min() == pytest.approx(brute)
    # exact correspondences recover a rigid transform
    P = rng.normal(size=(6, 3))
    R = DK.random_rotation(rng)
    Q = P @ R.T + np.array([3.0, -1.0, 2.0])
    lf = DK.Points(P, ["A"] * 6, np.ones(6))
    hs = DK.Points(Q, ["A"] * 6, np.ones(6))
    gg = DK.Graph([(i, i) for i in range(6)], np.ones(6), ~np.eye(6, dtype=bool))
    pose = DK.pose_from_clique(P, lf, hs, gg, np.ones(6))
    assert DK.rmsd(pose, Q) < 1e-9


# ---------------------------------------------------------------- heart-safety model (no network)
def test_safety_model_trains_and_predicts(monkeypatch):
    from chronocell.quantum import safety as SF
    smiles = ["CCO", "CCN(CC)CCOC(=O)c1ccccc1", "c1ccccc1CCN1CCCCC1", "OC(=O)CC(O)(CC(O)=O)C(O)=O", "CC(=O)Nc1ccc(O)cc1",
              "CN1CCN(CC1)c1ccccc1", "Clc1ccc(cc1)C(c1ccccc1)N1CCN(C)CC1", "OCC(O)CO", "NCC(O)=O", "c1ccc2ccccc2c1CCN(C)C"] * 4
    y = np.array([0, 1, 1, 0, 0, 1, 1, 0, 0, 1] * 4)
    X, ok = MF.matrix(smiles)
    monkeypatch.setattr(SF, "training_data", lambda progress=None: (X, y[ok], [smiles[k] for k in ok]))
    m = SF.train()
    r = m.predict("CN1CCN(CC1)c1ccc(Cl)cc1")
    assert {"qsvm_score", "rbf_score", "logistic_probability", "nearest"} <= set(r) and len(r["nearest"]) == 5


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


def test_drug_lab_keeps_its_core_default_and_adds_the_extended_set_and_guide(app):
    at = app
    at.segmented_control(key="workspace").set_value("Drug lab").run()
    assert not at.exception, [e.value for e in at.exception]
    assert at.segmented_control(key="lab_set").value == "Core (4)"
    assert len(at.radio(key="lab_drug").options) == 4
    assert any("Drug guide" in e.label for e in at.expander)
    assert any("heart safety" in e.label for e in at.expander)
    at.segmented_control(key="lab_set").set_value("Extended (12)").run()
    assert not at.exception, [e.value for e in at.exception]
    assert len(at.radio(key="lab_drug").options) == 12
    at.segmented_control(key="lab_set").set_value("Core (4)").run()
    assert not at.exception and len(at.radio(key="lab_drug").options) == 4


def test_quantum_lab_offers_the_drug_problems(app):
    at = app
    at.segmented_control(key="workspace").set_value("Quantum lab").run()
    opts = at.segmented_control(key="qlab_problem").options
    assert {"mol", "safety", "dock"} <= set(opts) or any("Docking" in o for o in opts)
    at.segmented_control(key="qlab_problem").set_value("mol").run()
    at.selectbox(key="qlab_mol_name").set_value("H2O").run()
    at.selectbox(key="qlab_mol_space").set_value("2, 2").run()
    at.button(key="qlab_mol_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    r = at.session_state["qlab_mol_last"]["res"]
    assert abs(r["err"]) < MO.CHEMICAL_ACCURACY and r["qubits"] == 4
