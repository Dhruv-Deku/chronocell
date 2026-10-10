"""
ADMET profile (October 2026): 21 absorption, distribution, metabolism, excretion and toxicity endpoints of the
Therapeutics Data Commons ADMET benchmark group (Huang et al., NeurIPS Datasets 2021; Harvard Dataverse
doi:10.7910/DVN/21LKWG, admet_group.zip, MD5 published by Dataverse), each with the official scaffold split
(train_val / test). hERG is left out (it trained the heart-safety screen, Gate Q5).

Molecules are described by the 17 SMILES descriptors of molfeat.py, standardised and reduced to 8 principal
components (8 qubits). Models on the same 8 numbers: a quantum-kernel model (ZZ feature map, kernels.py: an SVM for
yes/no endpoints, kernel ridge regression for measured values) and the same model with a classical RBF kernel; a
classical RBF model on all 17 descriptors is reported next to them as the stronger classical reference. Training sets
are capped at 2,500 compounds (a fixed random subset, the same for every model) to keep kernels small.

A screen for teaching and triage, not a safety assessment: measured standing in validation/RESULTS.md (Gate Q8).
"""

from __future__ import annotations

import math
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import kernels as KQ, molfeat as MF

CACHE = Path(__file__).resolve().parent.parent.parent / ".chronocell_cache" / "tdc"
ZIP = ("https://dataverse.harvard.edu/api/access/datafile/4426004", "admet_group.zip", "e3e4dda2f7f62268448fc813ceefccd5")
CAP = 2500
# name -> (task, short label, what it means, unit / positive class, log-transform the target for fitting)
ENDPOINTS = {
    "caco2_wang": ("reg", "Caco-2 permeability", "how well it crosses a gut-cell layer", "log cm/s", False),
    "hia_hou": ("cls", "Intestinal absorption", "absorbed from the gut", "absorbed", False),
    "pgp_broccatelli": ("cls", "P-gp inhibition", "blocks the P-glycoprotein drug pump", "inhibitor", False),
    "bioavailability_ma": ("cls", "Oral bioavailability", "reaches the blood when swallowed (F >= 20 %)", "bioavailable", False),
    "lipophilicity_astrazeneca": ("reg", "Lipophilicity", "oil-versus-water preference (logD7.4)", "logD", False),
    "solubility_aqsoldb": ("reg", "Water solubility", "how much dissolves in water", "log mol/L", False),
    "bbb_martins": ("cls", "Blood-brain barrier", "crosses into the brain", "penetrates", False),
    "ppbr_az": ("reg", "Plasma protein binding", "share bound to blood proteins", "%", False),
    "vdss_lombardo": ("reg", "Volume of distribution", "how widely it spreads into tissues", "L/kg", True),
    "cyp2c9_veith": ("cls", "CYP2C9 inhibition", "blocks liver enzyme CYP2C9 (drug interactions)", "inhibitor", False),
    "cyp2d6_veith": ("cls", "CYP2D6 inhibition", "blocks liver enzyme CYP2D6", "inhibitor", False),
    "cyp3a4_veith": ("cls", "CYP3A4 inhibition", "blocks liver enzyme CYP3A4 (half of all drugs use it)", "inhibitor", False),
    "cyp2c9_substrate_carbonmangels": ("cls", "CYP2C9 substrate", "broken down by CYP2C9", "substrate", False),
    "cyp2d6_substrate_carbonmangels": ("cls", "CYP2D6 substrate", "broken down by CYP2D6", "substrate", False),
    "cyp3a4_substrate_carbonmangels": ("cls", "CYP3A4 substrate", "broken down by CYP3A4", "substrate", False),
    "half_life_obach": ("reg", "Half-life", "time for the body to clear half a dose", "h", True),
    "clearance_hepatocyte_az": ("reg", "Hepatocyte clearance", "how fast liver cells clear it", "uL/min/1e6 cells", True),
    "clearance_microsome_az": ("reg", "Microsome clearance", "how fast liver enzymes clear it", "mL/min/g", True),
    "ld50_zhu": ("reg", "Acute toxicity (LD50)", "lethal dose in rats (higher = less toxic)", "-log mol/kg", False),
    "ames": ("cls", "Mutagenicity (Ames)", "damages DNA in the Ames bacterial test", "mutagenic", False),
    "dili": ("cls", "Drug-induced liver injury", "linked to liver injury", "liver injury", False),
}
GROUPS = {"Absorption": ["caco2_wang", "hia_hou", "pgp_broccatelli", "bioavailability_ma", "lipophilicity_astrazeneca",
                         "solubility_aqsoldb"],
          "Distribution": ["bbb_martins", "ppbr_az", "vdss_lombardo"],
          "Metabolism": ["cyp2c9_veith", "cyp2d6_veith", "cyp3a4_veith", "cyp2c9_substrate_carbonmangels",
                         "cyp2d6_substrate_carbonmangels", "cyp3a4_substrate_carbonmangels"],
          "Excretion": ["half_life_obach", "clearance_hepatocyte_az", "clearance_microsome_az"],
          "Toxicity": ["ld50_zhu", "ames", "dili"]}


def archive(path: Path | None = None, progress=None) -> zipfile.ZipFile:
    from ..predict import download_verified
    url, name, md5 = ZIP
    p = path or CACHE / name
    if not p.exists():
        download_verified(url, p, md5, progress)
    return zipfile.ZipFile(p)


def read_split(z: zipfile.ZipFile, name: str, part: str) -> tuple[list[str], np.ndarray]:
    import csv
    import io
    rows = list(csv.reader(io.StringIO(z.read(f"admet_group/{name}/{part}.csv").decode("utf-8"))))
    head = rows[0]
    i_s, i_y = head.index("Drug"), head.index("Y")
    smi, y = [], []
    for r in rows[1:]:
        if len(r) > max(i_s, i_y):
            try:
                y.append(float(r[i_y]))
                smi.append(r[i_s])
            except ValueError:
                continue
    return smi, np.array(y)


def featurise(smiles: list[str], y: np.ndarray) -> tuple[np.ndarray, np.ndarray, list[int]]:
    X, ok = MF.matrix(smiles)
    return X, y[ok], ok


def cap(n: int, k: int = CAP, seed: int = 0) -> np.ndarray:
    if n <= k:
        return np.arange(n)
    return np.sort(np.random.default_rng(seed).choice(n, k, replace=False))


def target(name: str, y: np.ndarray) -> np.ndarray:
    return np.log10(np.maximum(y, 0) + 1e-3) if ENDPOINTS[name][4] else y


@dataclass
class Prep:
    """Standardise -> 8 principal components (fitted on training); angles for the quantum feature map.

    Encodings of the components into angles (Gate Q8 used "minmax"; Gate Q8b "global"):
      minmax  each component to [0, pi] by its own training range (clipped outside it)
      global  every component divided by the same number, the first component's training SD: the components keep
              their relative variance, i.e. the geometry the RBF model on the same 8 numbers sees (minmax stretches
              the minor, noisier components to the same width as the main ones)"""
    mu: np.ndarray
    sd: np.ndarray
    V: np.ndarray
    lo: np.ndarray
    hi: np.ndarray
    s1: float = 1.0

    @classmethod
    def fit(cls, X: np.ndarray, k: int = 8) -> "Prep":
        mu, sd = X.mean(0), X.std(0) + 1e-9
        _, _, Vt = np.linalg.svd((X - mu) / sd, full_matrices=False)
        P = ((X - mu) / sd) @ Vt[:k].T
        return cls(mu, sd, Vt[:k], P.min(0), P.max(0), float(P[:, 0].std()) + 1e-9)

    def z(self, X: np.ndarray) -> np.ndarray:
        return (X - self.mu) / self.sd

    def pca(self, X: np.ndarray) -> np.ndarray:
        return self.z(X) @ self.V.T

    def angles(self, X: np.ndarray, encoding: str = "minmax") -> np.ndarray:
        if encoding == "global":
            return self.pca(X) / self.s1
        span = np.where(self.hi > self.lo, self.hi - self.lo, 1.0)
        return np.clip((self.pca(X) - self.lo) / span, 0, 1) * math.pi


def rbf(A: np.ndarray, B: np.ndarray, gamma: float) -> np.ndarray:
    d = (A * A).sum(1)[:, None] + (B * B).sum(1)[None, :] - 2 * A @ B.T
    return np.exp(-gamma * np.maximum(d, 0))


def qstates(A: np.ndarray, bandwidth: float, reps: int) -> np.ndarray:
    return KQ.states(A, bandwidth, reps)


def qkernel(SA: np.ndarray, SB: np.ndarray) -> np.ndarray:
    return np.abs(SA.conj() @ SB.T) ** 2


def projected_features(S: np.ndarray) -> np.ndarray:
    """Projected quantum kernel features (Huang et al., Nat Commun 2021): each qubit's reduced state, i.e. the single-
    qubit expectations <X_k>, <Y_k>, <Z_k> of the feature-map state (3 x qubits numbers, measurable on a device)."""
    S = np.asarray(S)
    n = int(round(math.log2(S.shape[1])))
    idx = np.arange(S.shape[1])
    out = []
    for k in range(n):
        b = (idx >> k) & 1
        p = np.abs(S) ** 2
        z = (p * (1 - 2 * b)[None, :]).sum(1)
        lo = idx[b == 0]
        c = (S[:, lo].conj() * S[:, lo | (1 << k)]).sum(1)
        out += [2 * c.real, 2 * c.imag, z]
    return np.stack(out, axis=1)


def make(task: str, reg: float):
    if task == "cls":
        from sklearn.svm import SVC
        return SVC(C=reg, kernel="precomputed", class_weight="balanced")
    from sklearn.kernel_ridge import KernelRidge
    return KernelRidge(alpha=reg, kernel="precomputed")


def fit_predict(task: str, Ktr: np.ndarray, ytr: np.ndarray, Kte: np.ndarray, reg: float) -> np.ndarray:
    """Kernel model on precomputed kernels: SVM (C = reg, balanced classes) or kernel ridge (alpha = reg, centred y)."""
    m = 0.0 if task == "cls" else float(ytr.mean())
    est = make(task, reg).fit(Ktr, ytr - m)
    return (est.decision_function(Kte) if task == "cls" else est.predict(Kte) + m)


def score(task: str, y: np.ndarray, s: np.ndarray) -> float:
    if task == "cls":
        from sklearn.metrics import roc_auc_score
        return float(roc_auc_score(y, s))
    from scipy.stats import spearmanr
    return float(spearmanr(y, s).statistic)


DEFAULT = {"qk": {"bandwidth": 0.1, "reps": 1, "reg": 10.0}, "rbf8": {"gamma": 0.1, "reg": 10.0},
           "rbf17": {"gamma": 0.05, "reg": 10.0}}


@dataclass
class EndpointModel:
    """The three models of one endpoint, fitted on its (capped) training set."""
    name: str
    task: str
    prep: Prep
    Xtr: np.ndarray
    ytr: np.ndarray
    cfg: dict
    smiles: list = field(default_factory=list)
    fitted: dict = field(default_factory=dict, repr=False)

    def fit(self) -> "EndpointModel":
        m = 0.0 if self.task == "cls" else float(self.ytr.mean())
        q = self.cfg["qk"]
        Sa = qstates(self.prep.angles(self.Xtr, q.get("encoding", "minmax")), q["bandwidth"], q["reps"])
        self.fitted["quantum"] = (make(self.task, q["reg"]).fit(qkernel(Sa, Sa), self.ytr - m), Sa, m)
        for k, f in (("rbf8", self.prep.pca), ("rbf17", self.prep.z)):
            A = f(self.Xtr)
            self.fitted[k] = (make(self.task, self.cfg[k]["reg"]).fit(rbf(A, A, self.cfg[k]["gamma"]), self.ytr - m), A, m)
        return self

    def predict(self, X: np.ndarray) -> dict:
        """Scores of the three models for descriptor rows X (17 columns): SVM decision values (> 0: the positive
        class) or predicted (transformed) values."""
        if not self.fitted:
            self.fit()
        out = {}
        for k, (est, ref, m) in self.fitted.items():
            if k == "quantum":
                q = self.cfg["qk"]
                K = qkernel(qstates(self.prep.angles(X, q.get("encoding", "minmax")), q["bandwidth"], q["reps"]), ref)
            else:
                f = self.prep.pca if k == "rbf8" else self.prep.z
                K = rbf(f(X), ref, self.cfg[k]["gamma"])
            out[k] = est.decision_function(K) if self.task == "cls" else est.predict(K) + m
        return out

    def predict_value(self, raw: np.ndarray) -> np.ndarray:
        """Back-transform a regression prediction to the endpoint's unit."""
        return 10 ** np.asarray(raw) - 1e-3 if ENDPOINTS[self.name][4] else np.asarray(raw)


def train(z: zipfile.ZipFile, name: str, cfg: dict | None = None) -> EndpointModel:
    task = ENDPOINTS[name][0]
    smi, y = read_split(z, name, "train_val")
    X, y, ok = featurise(smi, y)
    idx = cap(len(y))
    X, y = X[idx], y[idx]
    if task == "reg":
        y = target(name, y)
    return EndpointModel(name, task, Prep.fit(X), X, y, cfg or DEFAULT, [smi[ok[i]] for i in idx]).fit()


def settings() -> dict:
    """Per-endpoint settings frozen for Gate Q8 (validation/frozen.py), else the defaults."""
    try:
        import importlib.util
        p = Path(__file__).resolve().parent.parent.parent / "validation" / "frozen.py"
        spec = importlib.util.spec_from_file_location("chronocell_frozen_q8", p)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return getattr(mod, "QUANTUM_ADMET")["settings"]
    except Exception:
        return {}
