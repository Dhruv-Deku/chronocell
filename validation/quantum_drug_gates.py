"""
Quantum lab, drug tabs: Gates Q5 (safety classifier), Q6 (molecules) and Q7 (docking), on the statevector simulator
of chronocell/quantum, each against classical methods on the same input.

Q5  hERG cardiotoxicity (a drug blocking the hERG potassium channel can disturb heart rhythm; several HDAC
    inhibitors carry ECG warnings). Descriptors from SMILES (chronocell/quantum/molfeat.py, checked against RDKit).
    Practice: TDC hERG (Wang et al. 2016, 648 compounds; 5-fold cross-validation chooses settings). Test: trained on
    all of it, scored on the independent TDC hERG_Karim set (Karim et al. 2021, 13,445 compounds), compounds also in
    the training set (same RDKit canonical SMILES) removed. Quantum-kernel SVM vs RBF-SVM and logistic regression.
    Data: Harvard Dataverse (doi:10.7910/DVN/21LKWG), MD5 published by Dataverse.
Q6  Molecules: STO-3G Hartree-Fock / FCI from chronocell/quantum/molecules.py against OpenFermion's stored reference
    data (computed with an independent package; OpenFermion 1.7.1 installed from PyPI in the isolated environment).
    Practice: the H2 curve and Szabo & Ostlund's HF energies. Test: LiH at 1.45 A (HF and FCI). Also: active-space
    UCCSD-VQE within chemical accuracy of the active-space FCI for stretched molecules not run before.
Q7  Docking: rigid re-docking of PoseBusters ligands into their prepared proteins (Zenodo 13851241, MD5), as a
    maximum-weight clique (chronocell/quantum/docking.py). Complexes split by a hash of the PDB id (1 in 5 practice);
    complexes whose ligand does not touch the supplied protein (nearest protein atom >= 4 A) are skipped. QAOA on the
    20-vertex core of each interaction graph; the 30 heaviest distinct cliques it samples become poses, the best-scoring
    pose is kept. Classical comparators on the same graph (exact enumeration of maximal cliques, same pose step) and a
    random search with the same score (1,000 placements).

    python validation/quantum_drug_gates.py --practice q5|q6|q7
    python validation/quantum_drug_gates.py --test all          # run once (rules in frozen.QUANTUM_DRUG_GATES)
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
import time
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                                                            # noqa: E402
from chronocell.quantum import docking as DK, kernels as KQ, molecules as M, molfeat as MF, qubo as QB, sim  # noqa: E402

VENV_PY = ROOT.parent / ".chronocell_cache" / "quantum-venv" / "Scripts" / "python.exe"
CACHE = D.DATA / "cache" / "quantum_drug"
DV = "https://dataverse.harvard.edu/api/access/datafile/{}?format=original"
TDC = {"herg": (4259588, "herg.tab", "03a00602415d5fee10366eb679094305"),
       "herg_karim": (6822246, "herg_karim.tab", "cf699546d80721c9a160f4a6ac9b368b")}
PB = {"ligands": ("posebusters_ligands_256.zip", "fdeac0fa6ff4aeb46ce662f3ff8753dd"),
      "proteins": ("posebusters_spruce_structures_256.zip", "faef0d4a67d05dbfddbeba61fee5df6e")}
QSVM_GRID = [{"bandwidth": s, "reps": r, "C": c} for s in (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6) for r in (1, 2)
             for c in (0.1, 1.0, 10.0, 100.0, 1000.0)]
RBF_GRID = [{"gamma": g, "C": c} for g in (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5, 1.0, 2.0)
            for c in (0.1, 1.0, 10.0, 100.0, 1000.0)]
FEATURE_SETS = {"lipinski6": ["mw", "hbd", "hba", "rot_bonds", "tpsa", "aromatic_rings"],
                "herg8": ["mw", "tpsa", "hbd", "hba", "aromatic_rings", "basic_n", "n_halogen", "fsp3"],
                "pca8": "pca8"}
DOCK = {"tau": 1.0, "max_polar": 8, "max_hydrophobic": 4, "per_type": 8, "qubits": 20, "cliques": 30, "random_poses": 1000}
QAOA_GRID = [{"p": p, "objective": o} for p in (3, 5) for o in ("expectation", "cvar")]


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _save(name: str, obj) -> None:
    (ROOT / name).write_text(json.dumps(obj, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)),
                             encoding="utf-8")


# ======================================================================================
# Q5: hERG
# ======================================================================================
def tdc(name: str) -> list[tuple[str, int]]:
    fid, fname, md5 = TDC[name]
    p = D.fetch_url(DV.format(fid), D.DATA / "tdc" / fname, md5)
    import csv
    rows = []
    text = p.read_text(encoding="utf-8")
    delim = "\t" if "\t" in text.splitlines()[0] else ","
    for f in list(csv.reader(text.splitlines(), delimiter=delim))[1:]:
        if len(f) >= 3:
            try:
                rows.append((f[-2].strip(), int(round(float(f[-1])))))
            except ValueError:
                continue
    return rows


def rdkit_reference(smiles: list[str]) -> list[dict | None]:
    """RDKit descriptors and canonical SMILES (isolated environment), for checking molfeat and de-duplicating."""
    script = r'''
import json, sys
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Crippen
RDLogger.DisableLog("rdApp.*")
out = []
for s in json.load(open(sys.argv[1])):
    m = Chem.MolFromSmiles(s)
    if m is None:
        out.append(None); continue
    out.append({"canonical": Chem.MolToSmiles(m), "mw": Descriptors.MolWt(m), "hbd": Lipinski.NHOHCount(m),
                "hba": Lipinski.NOCount(m), "tpsa": rdMolDescriptors.CalcTPSA(m),
                "rot_bonds": rdMolDescriptors.CalcNumRotatableBonds(m), "rings": rdMolDescriptors.CalcNumRings(m),
                "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(m), "heavy_atoms": m.GetNumHeavyAtoms(),
                "fsp3": rdMolDescriptors.CalcFractionCSP3(m), "logp": Crippen.MolLogP(m)})
json.dump(out, open(sys.argv[2], "w"))
'''
    CACHE.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256("\n".join(smiles).encode()).hexdigest()[:16]
    out = CACHE / f"rdkit_{key}.json"
    if not out.exists():
        inp = CACHE / f"smiles_{key}.json"
        inp.write_text(json.dumps(smiles), encoding="utf-8")
        subprocess.run([str(VENV_PY), "-c", script, str(inp), str(out)], check=True)
    return json.loads(out.read_text(encoding="utf-8"))


def featurise(rows: list[tuple[str, int]]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    X, ok = MF.matrix([s for s, _ in rows])
    y = np.array([rows[k][1] for k in ok])
    return X, y, [rows[k][0] for k in ok]


def select(X: np.ndarray, fset: str, fit: dict | None = None) -> tuple[np.ndarray, dict]:
    """Feature set -> matrix (and the fitted transform, reused on the test set)."""
    if fset == "pca8":
        if fit is None:
            mu, sd = X.mean(0), X.std(0) + 1e-9
            Z = (X - mu) / sd
            _, _, Vt = np.linalg.svd(Z, full_matrices=False)
            fit = {"mu": mu, "sd": sd, "V": Vt[:8]}
        return ((X - fit["mu"]) / fit["sd"]) @ fit["V"].T, fit
    cols = [MF.FEATURES.index(f) for f in FEATURE_SETS[fset]]
    return X[:, cols], fit or {}


def _fit_predict(kind: str, prm: dict, Xtr, ytr, Xte) -> np.ndarray:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    if kind == "qsvm":
        a, b = KQ.scale_features(Xtr, Xte)
        return KQ.QSVM(prm["C"], prm["bandwidth"], prm["reps"]).fit(a, ytr).decision_function(b)
    sc = StandardScaler().fit(Xtr)
    if kind == "rbf":
        return SVC(C=prm["C"], kernel="rbf", gamma=prm["gamma"], class_weight="balanced").fit(sc.transform(Xtr), ytr) \
            .decision_function(sc.transform(Xte))
    return LogisticRegression(class_weight="balanced", max_iter=3000).fit(sc.transform(Xtr), ytr).decision_function(sc.transform(Xte))


def _cv(kind, prm, X, y, folds=5) -> float:
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    return float(np.mean([roc_auc_score(y[te], _fit_predict(kind, prm, X[tr], y[tr], X[te]))
                          for tr, te in StratifiedKFold(folds, shuffle=True, random_state=0).split(X, y)]))


def practice_q5() -> dict:
    rows = tdc("herg")
    X, y, smi = featurise(rows)
    ref = rdkit_reference(smi)
    agree = {}
    for f in ("mw", "hbd", "hba", "tpsa", "rot_bonds", "rings", "aromatic_rings", "heavy_atoms", "fsp3"):
        a = np.array([X[k, MF.FEATURES.index(f)] for k, r in enumerate(ref) if r])
        b = np.array([r[f] for r in ref if r])
        tol = 0.5 if f in ("mw", "tpsa") else 1e-6 if f != "fsp3" else 0.02
        agree[f] = {"exact_or_within_tol": float(np.mean(np.abs(a - b) <= tol)), "pearson": float(np.corrcoef(a, b)[0, 1])}
    out = {"made": _now(), "compounds": int(len(y)), "unreadable": int(len(rows) - len(y)), "blockers": float(y.mean()),
           "featuriser_vs_rdkit": agree, "feature_sets": {}}
    for fset in FEATURE_SETS:
        Xs, _ = select(X, fset)
        qs = [{**p, "cv_auc": _cv("qsvm", p, Xs, y)} for p in QSVM_GRID]
        rb = [{**p, "cv_auc": _cv("rbf", p, Xs, y)} for p in RBF_GRID]
        lr = _cv("logistic", {}, Xs, y)
        out["feature_sets"][fset] = {"qsvm_best": max(qs, key=lambda r: r["cv_auc"]), "rbf_best": max(rb, key=lambda r: r["cv_auc"]),
                                     "logistic": lr, "qsvm_grid": qs, "rbf_grid": rb}
        print(fset, round(out["feature_sets"][fset]["qsvm_best"]["cv_auc"], 3), round(out["feature_sets"][fset]["rbf_best"]["cv_auc"], 3),
              round(lr, 3), flush=True)
    best = max(out["feature_sets"], key=lambda k: out["feature_sets"][k]["qsvm_best"]["cv_auc"])
    fs = out["feature_sets"][best]
    out["choice"] = {"features": best, "qsvm": {k: v for k, v in fs["qsvm_best"].items() if k != "cv_auc"},
                     "rbf": {k: v for k, v in fs["rbf_best"].items() if k != "cv_auc"}}
    _save("results_qdrug_practice_q5.json", out)
    print("choice", out["choice"], json.dumps(agree), flush=True)
    return out


def test_q5(R: dict) -> dict:
    from sklearn.metrics import roc_auc_score
    tr_rows, te_rows = tdc("herg"), tdc("herg_karim")
    Xtr, ytr, str_ = featurise(tr_rows)
    Xte, yte, ste = featurise(te_rows)
    can_tr = {r["canonical"] for r in rdkit_reference(str_) if r}
    can_te = rdkit_reference(ste)
    keep = np.array([r is not None and r["canonical"] not in can_tr for r in can_te])
    Xte, yte = Xte[keep], yte[keep]
    a, fit = select(Xtr, R["features"])
    b, _ = select(Xte, R["features"], fit)
    sc = {"qsvm": _fit_predict("qsvm", R["qsvm"], a, ytr, b), "rbf": _fit_predict("rbf", R["rbf"], a, ytr, b),
          "logistic": _fit_predict("logistic", {}, a, ytr, b)}
    auc = {k: float(roc_auc_score(yte, v)) for k, v in sc.items()}
    rng = np.random.default_rng(0)
    boots = {k: [] for k in sc}
    diff = []
    for _ in range(1000):
        i = rng.integers(0, len(yte), len(yte))
        if len(set(yte[i])) < 2:
            continue
        r = {k: roc_auc_score(yte[i], v[i]) for k, v in sc.items()}
        for k in r:
            boots[k].append(r[k])
        diff.append(r["qsvm"] - r["rbf"])
    ci = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] for k, v in boots.items()}
    ci["qsvm_minus_rbf"] = [float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))]
    ok = bool(auc["qsvm"] >= auc["rbf"] - R["margin"] and ci["qsvm"][0] > 0.5)
    return {"train": int(len(ytr)), "test": int(len(yte)), "test_removed_overlap": int((~keep).sum()),
            "test_blockers": float(yte.mean()), "auc": auc, "ci95": ci, "pass": ok}


# ======================================================================================
# Q6: molecules
# ======================================================================================
def openfermion_reference(names: list[str]) -> dict:
    script = r'''
import json, os, sys, openfermion
d = os.path.join(os.path.dirname(openfermion.__file__), "testing", "data")
out = {}
for n in json.loads(sys.argv[1]):
    m = openfermion.MolecularData(filename=os.path.join(d, n))
    m.load()
    out[n] = {"geometry": m.geometry, "basis": m.basis, "charge": m.charge, "multiplicity": m.multiplicity,
              "hf_energy": float(m.hf_energy) if m.hf_energy is not None else None,
              "fci_energy": float(m.fci_energy) if m.fci_energy is not None else None}
print(json.dumps(out))
'''
    r = subprocess.run([str(VENV_PY), "-c", script, json.dumps(names)], capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def full_fci(atoms, charge=0) -> tuple[float, float]:
    ints = M.integrals(atoms, charge)
    hf = M.rhf(ints)
    n = len(ints.S)
    asp = M.active_space(ints, hf, n, ints.n_electrons)
    return hf.energy, M.fci(M.qubit_hamiltonian(asp))


def practice_q6() -> dict:
    out = {"made": _now(), "szabo_ostlund": {}, "openfermion_h2": {}}
    for name, ref in M.SZABO_OSTLUND_HF.items():
        hf = M.rhf(M.integrals(M.LIBRARY[name][1]()))
        out["szabo_ostlund"][name] = {"ours": hf.energy, "reference": ref, "diff_mEh": 1e3 * (hf.energy - ref)}
    names = [f"H2_sto-3g_singlet_{r}.hdf5" for r in ("0.5", "0.7414", "1.0", "1.5", "2.0", "2.5")]
    ref = openfermion_reference(names)
    for n, v in ref.items():
        atoms = [(a, tuple(xyz)) for a, xyz in v["geometry"]]
        e_hf, e_fci = full_fci(atoms)
        out["openfermion_h2"][n] = {"hf": e_hf, "fci": e_fci, "ref_hf": v["hf_energy"], "ref_fci": v["fci_energy"],
                                    "diff_hf_mEh": 1e3 * (e_hf - v["hf_energy"]), "diff_fci_mEh": 1e3 * (e_fci - v["fci_energy"])}
        print(n, out["openfermion_h2"][n], flush=True)
    _save("results_qdrug_practice_q6.json", out)
    return out


def test_q6(R: dict) -> dict:
    ref = openfermion_reference(R["openfermion_test"])
    rows = []
    for n, v in ref.items():
        atoms = [(a, tuple(xyz)) for a, xyz in v["geometry"]]
        e_hf, e_fci = full_fci(atoms, v["charge"])
        rows.append({"file": n, "hf": e_hf, "fci": e_fci, "ref_hf": v["hf_energy"], "ref_fci": v["fci_energy"],
                     "diff_hf_mEh": 1e3 * (e_hf - v["hf_energy"]), "diff_fci_mEh": 1e3 * (e_fci - v["fci_energy"])})
        print(rows[-1], flush=True)
    vqe_rows = []
    for name, scale, act in R["vqe_cases"]:
        r = M.run(name, active=tuple(act), scale=scale)
        vqe_rows.append({"molecule": name, "bond_scale": scale, "active": act, "qubits": r.n_qubits, "e_hf": r.e_hf,
                         "e_cas_fci": r.e_cas_fci, "e_vqe": r.vqe.energy, "error_mEh": 1e3 * r.vqe.error,
                         "parameters": r.vqe.parameters, "electrons": r.vqe.electrons, "seconds": r.seconds})
        print(vqe_rows[-1], flush=True)
    ref_ok = all(abs(x["diff_hf_mEh"]) <= R["reference_tol_mEh"] and abs(x["diff_fci_mEh"]) <= R["reference_tol_mEh"] for x in rows)
    vqe_ok = all(abs(x["error_mEh"]) <= 1e3 * M.CHEMICAL_ACCURACY for x in vqe_rows)
    return {"reference": rows, "vqe": vqe_rows, "reference_pass": ref_ok, "vqe_pass": vqe_ok, "pass": bool(ref_ok and vqe_ok)}


# ======================================================================================
# Q7: docking
# ======================================================================================
def pb_zip(kind: str) -> zipfile.ZipFile:
    name, md5 = PB[kind]
    url = f"https://zenodo.org/api/records/13851241/files/{name}/content"
    return zipfile.ZipFile(D.fetch_url(url, D.DATA / "posebusters" / name, md5))


def complexes(role: str) -> list[str]:
    zl = pb_zip("ligands")
    ids = sorted(n.split("/")[-1].replace("_ligand.sdf", "") for n in zl.namelist()
                 if n.startswith("posebusters_ligands_256/") and n.endswith("_ligand.sdf"))
    pick = (lambda i: int(hashlib.sha1(i.encode()).hexdigest(), 16) % 5 == 0)
    return [i for i in ids if pick(i) == (role == "practice")]


def maximal_cliques(g: DK.Graph, k: int) -> list:
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(len(g.vertices)))
    for a, b in zip(*np.nonzero(np.triu(g.adj, 1))):
        G.add_edge(int(a), int(b))
    out = [(float(g.weight[c].sum()), sorted(c)) for c in nx.find_cliques(G) if len(c) >= 3]
    return sorted(out, reverse=True)[:k]


def sampled_cliques(g: DK.Graph, r: sim.QAOAResult, k: int) -> list:
    seen, out = set(), []
    order = np.argsort(-r.probs)
    for idx in order[:20000]:
        bits = QB.to_bits(int(idx), len(g.vertices))
        on = tuple(np.flatnonzero(bits).tolist())
        if len(on) >= 3 and on not in seen and DK.is_clique(g, bits):
            seen.add(on)
            out.append((float(g.weight[list(on)].sum()), list(on)))
    for idx in r.samples:                              # also what was actually measured
        bits = QB.to_bits(int(idx), len(g.vertices))
        on = tuple(np.flatnonzero(bits).tolist())
        if len(on) >= 3 and on not in seen and DK.is_clique(g, bits):
            seen.add(on)
            out.append((float(g.weight[list(on)].sum()), list(on)))
    return sorted(out, reverse=True)[:k]


def dock_one(cid: str, cfg: dict, qcfg: dict | None, seed: int, tag: str) -> dict | None:
    path = CACHE / f"dock_{tag}_{cid}_{hashlib.sha1(json.dumps([cfg, qcfg]).encode()).hexdigest()[:10]}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    from scipy.spatial import cKDTree
    zl, zp = pb_zip("ligands"), pb_zip("proteins")
    lig = DK.read_sdf(zl.read(f"posebusters_ligands_256/{cid}_ligand.sdf").decode(), cid)
    prot = DK.read_pdb(zp.read(f"posebusters_spruce_structures_256/{cid}_spruced.pdb").decode(errors="replace"))
    tree = cKDTree(prot.xyz[prot.el != "H"])
    heavy = lig.heavy
    L = lig.xyz[heavy]
    if tree.query(L)[0].min() >= 4.0:
        rec = {"cid": cid, "usable": False}
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(rec), encoding="utf-8")
        return rec
    centre = L.mean(0)
    lf = DK.ligand_features(lig, cfg["max_polar"], cfg["max_hydrophobic"])
    hs = DK.hotspots_directional(prot, centre, per_type=cfg["per_type"])
    rec = {"cid": cid, "usable": True, "heavy_atoms": int(heavy.sum()), "features": len(lf.kind), "hotspots": len(hs.kind)}

    def best_pose(cliques):
        poses = []
        for w, c in cliques:
            bits = np.zeros(len(g.vertices))
            bits[c] = 1
            p = DK.pose_from_clique(lig.xyz, lf, hs, g, bits)
            if p is not None:
                poses.append((DK.pose_score(p[heavy], prot, tree), p))
        if not poses:
            return float("nan")
        return DK.rmsd(max(poses, key=lambda z: z[0])[1][heavy], L)

    if len(hs.kind) >= 3 and len(lf.kind) >= 3:
        full = DK.interaction_graph(lf, hs, tau=cfg["tau"], max_vertices=80)
        g, _ = DK.core_subgraph(full, cfg["qubits"])
        q = DK.clique_qubo(g)
        e = q.energies()
        ex = QB.exact(q, e)
        rec["vertices"] = q.n
        rec["max_clique_weight"] = float(-ex.energy)
        rec["max_clique_size"] = int(ex.bits.sum())
        rec["classical_rmsd"] = best_pose(maximal_cliques(g, cfg["cliques"]))
        if qcfg is not None:
            r = sim.qaoa(e, p=qcfg["p"], shots=4096, objective=qcfg["objective"], maxiter=80, seed=seed)
            rec["qaoa_hit"] = bool(r.hit)
            rec["qaoa_p_optimal"] = r.p_optimal
            rec["qaoa_seconds"] = r.seconds
            rec["qaoa_rmsd"] = best_pose(sampled_cliques(g, r, cfg["cliques"]))
            rec["uniform_p_optimal"] = ex.detail["n_optima"] / 2 ** q.n
    rng = np.random.default_rng(seed)
    rec["random_rmsd"] = DK.rmsd(DK.random_search(L, centre, prot, cfg["random_poses"], rng), L)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec), encoding="utf-8")
    return rec


def summarise_dock(recs: list[dict]) -> dict:
    use = [r for r in recs if r.get("usable")]
    ok = lambda k: float(np.mean([np.isfinite(r.get(k, np.nan)) and r[k] <= 2.0 for r in use])) if use else float("nan")
    out = {"complexes": len(recs), "usable": len(use), "success_classical_clique": ok("classical_rmsd"),
           "success_random_search": ok("random_rmsd")}
    if any("qaoa_hit" in r for r in use):
        qq = [r for r in use if "qaoa_hit" in r]
        out.update({"success_qaoa": ok("qaoa_rmsd"), "qaoa_hit_rate": float(np.mean([r["qaoa_hit"] for r in qq])),
                    "graphs": len(qq), "mean_qaoa_p_optimal": float(np.mean([r["qaoa_p_optimal"] for r in qq])),
                    "mean_uniform_p_optimal": float(np.mean([r["uniform_p_optimal"] for r in qq])),
                    "mean_qaoa_seconds": float(np.mean([r["qaoa_seconds"] for r in qq]))})
    return out


def practice_q7() -> dict:
    ids = complexes("practice")
    out = {"made": _now(), "settings": DOCK, "qaoa": {}}
    for qcfg in QAOA_GRID:
        recs = []
        for k, cid in enumerate(ids):
            recs.append(dock_one(cid, DOCK, qcfg, k, "practice"))
        out["qaoa"][f"p{qcfg['p']}_{qcfg['objective']}"] = summarise_dock(recs)
        print(qcfg, out["qaoa"][f"p{qcfg['p']}_{qcfg['objective']}"], flush=True)
    best = max(out["qaoa"], key=lambda k: (out["qaoa"][k]["qaoa_hit_rate"], out["qaoa"][k]["success_qaoa"], -int(k[1])))
    out["choice"] = {"p": int(best.split("_")[0][1:]), "objective": best.split("_", 1)[1]}
    _save("results_qdrug_practice_q7.json", out)
    return out


def test_q7(R: dict) -> dict:
    ids = complexes("test")
    recs = []
    t0 = time.perf_counter()
    for k, cid in enumerate(ids):
        recs.append(dock_one(cid, R["settings"], R["qaoa"], k, "test"))
        if k % 10 == 0:
            print(k, cid, round(time.perf_counter() - t0), "s", flush=True)
    s = summarise_dock(recs)
    s["solver_pass"] = bool(s["qaoa_hit_rate"] >= R["min_hit_rate"])
    s["docking_pass"] = bool(s["success_qaoa"] >= s["success_random_search"])
    s["pass"] = bool(s["solver_pass"] and s["docking_pass"])
    s["per_complex"] = recs
    return s


def run_test() -> dict:
    import frozen as F
    R = F.QUANTUM_DRUG_GATES
    out = {"made": _now(), "rule": R}
    path = ROOT / "results_qdrug.json"
    if path.exists():
        out = json.loads(path.read_text(encoding="utf-8"))
    for part, fn in (("q6", test_q6), ("q5", test_q5), ("q7", test_q7)):
        if part in out:
            continue
        out[part] = fn(R[part])
        _save("results_qdrug.json", out)
        print(part, "pass" if out[part]["pass"] else "fail", flush=True)
    out["overall"] = {k: out[k]["pass"] for k in ("q5", "q6", "q7")}
    _save("results_qdrug.json", out)
    print(json.dumps(out["overall"]), flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", choices=["q5", "q6", "q7"])
    g.add_argument("--test", choices=["all"])
    a = ap.parse_args()
    if a.practice:
        {"q5": practice_q5, "q6": practice_q6, "q7": practice_q7}[a.practice]()
    else:
        run_test()


if __name__ == "__main__":
    main()
