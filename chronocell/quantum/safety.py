"""
hERG heart-safety screen for the Drug lab (quantum-kernel SVM next to classical models), app side.

Training data: TDC hERG (Wang et al., Mol Pharm 2016; 648 compounds labelled hERG blocker / non-blocker), downloaded
on demand from Harvard Dataverse and checked against the MD5 Dataverse publishes. Descriptors: molfeat.py from
SMILES. Settings: those frozen for Gate Q5 (validation/frozen.py), read when present, else the practice defaults.

A screen, not a safety assessment: the models learn from a few hundred mostly drug-like compounds, and the
epigenetic drugs of the Drug lab are not among them. The measured standing (Gate Q5, an independent 13,000-compound
set) is shown with every prediction.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import kernels as KQ, molfeat as MF

CACHE = Path(__file__).resolve().parent.parent.parent / ".chronocell_cache" / "tdc"
HERG = ("https://dataverse.harvard.edu/api/access/datafile/4259588?format=original", "herg.tab",
        "03a00602415d5fee10366eb679094305")
DEFAULT = {"features": "herg8", "qsvm": {"bandwidth": 0.005, "reps": 2, "C": 1000.0}, "rbf": {"gamma": 0.002, "C": 100.0}}


def settings() -> dict:
    try:
        import importlib.util
        p = Path(__file__).resolve().parent.parent.parent / "validation" / "frozen.py"
        spec = importlib.util.spec_from_file_location("chronocell_frozen_q5", p)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return getattr(mod, "QUANTUM_DRUG_GATES")["q5"]
    except Exception:
        return DEFAULT


def training_data(progress=None) -> tuple[np.ndarray, np.ndarray, list[str]]:
    from ..predict import download_verified
    url, name, md5 = HERG
    p = CACHE / name
    if not p.exists():
        download_verified(url, p, md5, progress)
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
    X, ok = MF.matrix([s for s, _ in rows])
    return X, np.array([rows[k][1] for k in ok]), [rows[k][0] for k in ok]


FEATURE_SETS = {"lipinski6": ["mw", "hbd", "hba", "rot_bonds", "tpsa", "aromatic_rings"],
                "herg8": ["mw", "tpsa", "hbd", "hba", "aromatic_rings", "basic_n", "n_halogen", "fsp3"]}


@dataclass
class Model:
    fit: dict
    Xtr: np.ndarray
    ytr: np.ndarray
    smiles: list
    qsvm: object
    rbf: object
    logistic: object
    scaler: object
    cfg: dict

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self.cfg["features"] == "pca8":
            return ((X - self.fit["mu"]) / self.fit["sd"]) @ self.fit["V"].T
        cols = [MF.FEATURES.index(f) for f in FEATURE_SETS[self.cfg["features"]]]
        return X[:, cols]

    def predict(self, smiles: str) -> dict:
        d = MF.descriptors(smiles)
        x = self.transform(np.array([[d[f] for f in MF.FEATURES]], float))
        a_tr = self.transform(self.Xtr)
        _, a = KQ.scale_features(a_tr, x)
        q = float(self.qsvm.decision_function(a)[0])
        r = float(self.rbf.decision_function(self.scaler.transform(x))[0])
        lp = float(self.logistic.predict_proba(self.scaler.transform(x))[0, 1])
        # nearest training compounds in descriptor space
        z = (self.Xtr - self.Xtr.mean(0)) / (self.Xtr.std(0) + 1e-9)
        zx = (np.array([d[f] for f in MF.FEATURES]) - self.Xtr.mean(0)) / (self.Xtr.std(0) + 1e-9)
        near = np.argsort(np.linalg.norm(z - zx, axis=1))[:5]
        return {"descriptors": d, "qsvm_score": q, "qsvm_says": "likely blocker" if q > 0 else "likely non-blocker",
                "rbf_score": r, "rbf_says": "likely blocker" if r > 0 else "likely non-blocker",
                "logistic_probability": lp,
                "nearest": [{"smiles": self.smiles[k], "blocker": int(self.ytr[k])} for k in near]}


def train(progress=None) -> Model:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    cfg = settings()
    X, y, smi = training_data(progress)
    fit = {}
    if cfg["features"] == "pca8":
        mu, sd = X.mean(0), X.std(0) + 1e-9
        _, _, Vt = np.linalg.svd((X - mu) / sd, full_matrices=False)
        fit = {"mu": mu, "sd": sd, "V": Vt[:8]}
    m = Model(fit, X, y, smi, None, None, None, None, cfg)
    A = m.transform(X)
    a, = KQ.scale_features(A)
    m.qsvm = KQ.QSVM(cfg["qsvm"]["C"], cfg["qsvm"]["bandwidth"], cfg["qsvm"]["reps"]).fit(a, y)
    m.scaler = StandardScaler().fit(A)
    m.rbf = SVC(C=cfg["rbf"]["C"], kernel="rbf", gamma=cfg["rbf"]["gamma"], class_weight="balanced").fit(m.scaler.transform(A), y)
    m.logistic = LogisticRegression(class_weight="balanced", max_iter=3000).fit(m.scaler.transform(A), y)
    return m


def standing() -> dict | None:
    p = Path(__file__).resolve().parent.parent.parent / "validation" / "results_qdrug.json"
    try:
        r = json.loads(p.read_text(encoding="utf-8"))
        return r.get("q5")
    except (OSError, ValueError):
        return None
