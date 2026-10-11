"""
Gate Q8 (October 2026): an ADMET profile from a quantum-kernel model, 21 endpoints of the TDC ADMET benchmark group,
each tested on its official held-out scaffold split (chronocell/quantum/admet.py).

Practice: train_val sets only. Per endpoint, 5-fold cross-validation (stratified for yes/no endpoints) chooses the
quantum-kernel model's settings (bandwidth, repetitions, regularisation) and, over the same folds, the two classical
models' (RBF on the same 8 principal components; RBF on all 17 descriptors). The test sets are not read.
Test (run once, rule in frozen.QUANTUM_ADMET): train on train_val (capped at 2,500 compounds), score the test set.
Metric: ROC AUC (yes/no endpoints) or Spearman rho (measured values); 95 % intervals from 1,000 resamples.

    python validation/admet_gate.py --practice
    python validation/admet_gate.py --test
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                                           # noqa: E402
from chronocell.quantum import admet as A                      # noqa: E402

ZIP_PATH = D.DATA / "tdc" / "admet_group.zip"
# grids widened once on practice (all three models alike) when first-run best values sat at their edges
QK_GRID = [{"bandwidth": b, "reps": r} for b in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0) for r in (1, 2)]
RBF_GAMMAS = (0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0)
REGS = {"cls": (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0), "reg": (0.001, 0.01, 0.1, 1.0, 10.0, 100.0)}


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _json(path: Path, obj) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)),
                   encoding="utf-8")
    tmp.replace(path)


def data(name: str, part: str) -> tuple[np.ndarray, np.ndarray]:
    z = A.archive(ZIP_PATH)
    smi, y = A.read_split(z, name, part)
    X, y, _ = A.featurise(smi, y)
    if part == "train_val":
        idx = A.cap(len(y))
        X, y = X[idx], y[idx]
    return X, (A.target(name, y) if A.ENDPOINTS[name][0] == "reg" else y)


def practice_one(name: str) -> dict:
    from sklearn.model_selection import KFold, StratifiedKFold
    task = A.ENDPOINTS[name][0]
    X, y = data(name, "train_val")
    split = (StratifiedKFold if task == "cls" else KFold)(5, shuffle=True, random_state=0)
    regs = REGS[task]
    qk = {(g["bandwidth"], g["reps"], c): [] for g in QK_GRID for c in regs}
    r8 = {(g, c): [] for g in RBF_GAMMAS for c in regs}
    r17 = {(g, c): [] for g in RBF_GAMMAS for c in regs}
    for tr, te in split.split(X, y):
        prep = A.Prep.fit(X[tr])
        a_tr, a_te = prep.angles(X[tr]), prep.angles(X[te])
        for g in QK_GRID:
            S_tr, S_te = A.qstates(a_tr, g["bandwidth"], g["reps"]), A.qstates(a_te, g["bandwidth"], g["reps"])
            K, Kt = A.qkernel(S_tr, S_tr), A.qkernel(S_te, S_tr)
            for c in regs:
                qk[(g["bandwidth"], g["reps"], c)].append(A.score(task, y[te], A.fit_predict(task, K, y[tr], Kt, c)))
        for store, f in ((r8, prep.pca), (r17, prep.z)):
            P, Pt = f(X[tr]), f(X[te])
            for gm in RBF_GAMMAS:
                K, Kt = A.rbf(P, P, gm), A.rbf(Pt, P, gm)
                for c in regs:
                    store[(gm, c)].append(A.score(task, y[te], A.fit_predict(task, K, y[tr], Kt, c)))
    mean = lambda d: {k: float(np.mean(v)) for k, v in d.items()}           # noqa: E731
    qk, r8, r17 = mean(qk), mean(r8), mean(r17)
    bq, b8, b17 = max(qk, key=qk.get), max(r8, key=r8.get), max(r17, key=r17.get)
    out = {"endpoint": name, "task": task, "n": int(len(y)), "positives": float(y.mean()) if task == "cls" else None,
           "qk_best": {"bandwidth": bq[0], "reps": bq[1], "reg": bq[2], "cv": qk[bq]},
           "rbf8_best": {"gamma": b8[0], "reg": b8[1], "cv": r8[b8]},
           "rbf17_best": {"gamma": b17[0], "reg": b17[1], "cv": r17[b17]},
           "qk_grid": [{"bandwidth": k[0], "reps": k[1], "reg": k[2], "cv": v} for k, v in qk.items()],
           "rbf8_grid": [{"gamma": k[0], "reg": k[1], "cv": v} for k, v in r8.items()],
           "rbf17_grid": [{"gamma": k[0], "reg": k[1], "cv": v} for k, v in r17.items()]}
    print(name, task, len(y), "quantum", round(qk[bq], 3), "rbf8", round(r8[b8], 3), "rbf17", round(r17[b17], 3), flush=True)
    return out


def practice() -> dict:
    names = list(A.ENDPOINTS)
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(practice_one, names))
    out = {"made": _now(), "cap": A.CAP, "endpoints": {r["endpoint"]: r for r in res},
           "choice": {r["endpoint"]: {"qk": {k: r["qk_best"][k] for k in ("bandwidth", "reps", "reg")},
                                      "rbf8": {k: r["rbf8_best"][k] for k in ("gamma", "reg")},
                                      "rbf17": {k: r["rbf17_best"][k] for k in ("gamma", "reg")}} for r in res}}
    _json(ROOT / "results_admet_practice.json", out)
    return out


def _boot(task: str, y: np.ndarray, scores: dict, reps: int = 1000, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    n = len(y)
    b = {k: [] for k in scores}
    diff = []
    for _ in range(reps):
        i = rng.integers(0, n, n)
        if task == "cls" and len(set(y[i])) < 2:
            continue
        s = {k: A.score(task, y[i], v[i]) for k, v in scores.items()}
        for k in s:
            b[k].append(s[k])
        diff.append(s["quantum"] - s["rbf8"])
    ci = {k: [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))] for k, v in b.items()}
    ci["quantum_minus_rbf8"] = [float(np.nanpercentile(diff, 2.5)), float(np.nanpercentile(diff, 97.5))]
    return ci


def test_one(args) -> dict:
    name, cfg, margin = args
    task = A.ENDPOINTS[name][0]
    Xtr, ytr = data(name, "train_val")
    Xte, yte = data(name, "test")
    m = A.EndpointModel(name, task, A.Prep.fit(Xtr), Xtr, ytr, cfg).fit()
    sc = m.predict(Xte)
    met = {k: A.score(task, yte, v) for k, v in sc.items()}
    ci = _boot(task, yte, sc)
    floor = 0.5 if task == "cls" else 0.0
    ok = bool(met["quantum"] >= met["rbf8"] - margin and ci["quantum"][0] > floor)
    print(name, {k: round(v, 3) for k, v in met.items()}, "pass" if ok else "fail", flush=True)
    return {"endpoint": name, "task": task, "train": int(len(ytr)), "test": int(len(yte)), "metric": met, "ci95": ci,
            "pass": ok}


def test() -> dict:
    import frozen as F
    R = F.QUANTUM_ADMET
    path = ROOT / "results_admet.json"
    if path.exists():
        sys.exit("results_admet.json exists: Gate Q8 runs once.")
    jobs = [(n, R["settings"][n], R["margin"]) for n in R["endpoints"]]
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(test_one, jobs))
    passed = sum(r["pass"] for r in res)
    out = {"made": _now(), "rule": R, "endpoints": {r["endpoint"]: r for r in res}, "passed": int(passed),
           "of": len(res), "pass": bool(passed >= R["min_endpoints"])}
    _json(path, out)
    print("Gate Q8:", "pass" if out["pass"] else "fail", f"{passed} of {len(res)}", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    practice() if a.practice else test()


if __name__ == "__main__":
    main()
