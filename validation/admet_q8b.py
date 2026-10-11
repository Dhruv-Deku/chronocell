"""
Gate Q8b (October 2026): the quantum-kernel drug-property model again, with Q8's failure fixed, on 19 endpoints no
ChronoCell test has used.

Why Q8 failed, found after its test (its test splits are practice data now; Q8 stays a fail): the quantum model's
encoding stretched each of the 8 principal components to [0, pi] by its own training range, so the minor, noisier
components weighed as much as the main ones, while the RBF model it was compared with sees the components with their
real variances. Q8 missed on 5 endpoints. With the "global" encoding (every component divided by the first
component's training SD; chronocell.quantum.admet.Prep.angles), settings chosen by the same 5-fold cross-validation,
the quantum model meets Q8's per-endpoint rule on more of Q8's (now seen) test splits: --posthoc writes those numbers
to validation/results_admet_q8b_posthoc.json (the reason for Q8b, not evidence).

Data (validation/admet_q8b_split.py, split frozen before any model was fitted): TDC CYP1A2 and CYP2C19 inhibition
(Veith), PAMPA permeability (NCATS), hydration free energy (FreeSolv), skin reaction, carcinogens (Lagunin), ClinTox and
the twelve Tox21 assays; one scaffold split each (80 % train / 20 % test), Tox21 one split for all assays.

Practice (train parts only; test parts not read): per endpoint, 5-fold cross-validation (stratified for yes/no) picks
the quantum model's bandwidth, repetitions and regularisation (global encoding) and the two classical models' settings,
as Q8. Test (run once, rule in frozen.QUANTUM_ADMET_Q8B): train on the train part (capped at 2,500 compounds), score the
test part. Metric: ROC AUC (yes/no) or Spearman rho (values); 95 % intervals from 1,000 resamples.

    python validation/admet_q8b.py --posthoc       # Q8's 21 endpoints with the new encoding (seen data; the reason)
    python validation/admet_q8b.py --practice      # choose settings on the 19 new endpoints' train parts
    python validation/admet_q8b.py --test          # Gate Q8b, run once
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import admet_gate as G                                          # noqa: E402
from chronocell.quantum import admet as A, molfeat as MF        # noqa: E402

DATA = ROOT / "data" / "tdc" / "q8b"
SPLIT = ROOT / "admet_q8b_split.json"
# global encoding: angles are in units of PC1's SD. The grids were widened once, all three models alike, when the first
# practice run's best values sat at their edges (quantum bandwidth 0.02 for PAMPA and ClinTox, its regularisation 1,000
# for ClinTox, an RBF width of 3 for a Tox21 assay), as Q8's were.
BW = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0)
RBF_GAMMAS = (0.001,) + G.RBF_GAMMAS + (10.0,)
REGS = {"cls": (0.001,) + G.REGS["cls"] + (10000.0,), "reg": (0.0001,) + G.REGS["reg"] + (1000.0,)}
PARTS = ROOT / "data" / "cache" / "q8b"                          # one file per endpoint and mode: a stopped run resumes
TOX21 = ["NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma", "SR-ARE", "SR-ATAD5",
         "SR-HSE", "SR-MMP", "SR-p53"]
# endpoint -> (file, label column, task)
ENDPOINTS = {"cyp1a2_veith": ("cyp1a2_veith", "Y", "cls"), "cyp2c19_veith": ("cyp2c19_veith", "Y", "cls"),
             "pampa_ncats": ("pampa_ncats", "Y", "cls"), "freesolv": ("hydrationfreeenergy_freesolv", "Y", "reg"),
             "skin_reaction": ("skin_reaction", "Y", "cls"), "carcinogens_lagunin": ("carcinogens_lagunin", "Y", "cls"),
             "clintox": ("clintox", "Y", "cls"),
             **{f"tox21_{a}": ("tox21", a, "cls") for a in TOX21}}


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def data(name: str, part: str) -> tuple[np.ndarray, np.ndarray]:
    """Descriptor rows and labels of one endpoint's train or test part (train capped at 2,500 like Q8)."""
    fname, col, task = ENDPOINTS[name]
    split = json.loads(SPLIT.read_text(encoding="utf-8"))["splits"][fname][part]
    smi_col = "X" if fname == "tox21" else "Drug"
    with (DATA / f"{fname}.csv").open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    smi, y = [], []
    for i in split:
        v = rows[i][col]
        if v not in ("", "nan", "NaN"):
            smi.append(rows[i][smi_col])
            y.append(float(v))
    X, ok = MF.matrix(smi)
    y = np.asarray(y)[ok]
    if part == "train":
        idx = A.cap(len(y))
        X, y = X[idx], y[idx]
    return X, (y.astype(int) if task == "cls" else y)


def cv_one(args) -> dict:
    """5-fold CV of the quantum model (global encoding) and both RBF models on one endpoint's training data."""
    name, X, y, task = args
    from sklearn.model_selection import KFold, StratifiedKFold
    regs = REGS[task]
    split = list((StratifiedKFold if task == "cls" else KFold)(5, shuffle=True, random_state=0).split(X, y))
    qk, r8, r17 = {}, {}, {}
    for tr, te in split:
        prep = A.Prep.fit(X[tr])
        a_tr, a_te = prep.angles(X[tr], "global"), prep.angles(X[te], "global")
        for bw in BW:
            for reps in (1, 2):
                Str, Ste = A.qstates(a_tr, bw, reps), A.qstates(a_te, bw, reps)
                K, Kt = A.qkernel(Str, Str), A.qkernel(Ste, Str)
                for c in regs:
                    qk.setdefault((bw, reps, c), []).append(A.score(task, y[te], A.fit_predict(task, K, y[tr], Kt, c)))
        for store, f in ((r8, prep.pca), (r17, prep.z)):
            P, Pt = f(X[tr]), f(X[te])
            for gm in RBF_GAMMAS:
                K, Kt = A.rbf(P, P, gm), A.rbf(Pt, P, gm)
                for c in regs:
                    store.setdefault((gm, c), []).append(A.score(task, y[te], A.fit_predict(task, K, y[tr], Kt, c)))
    mean = lambda d: {k: float(np.mean(v)) for k, v in d.items()}           # noqa: E731
    qk, r8, r17 = mean(qk), mean(r8), mean(r17)
    bq, b8, b17 = max(qk, key=qk.get), max(r8, key=r8.get), max(r17, key=r17.get)
    out = {"endpoint": name, "task": task, "n": int(len(y)), "positives": float(np.mean(y)) if task == "cls" else None,
           "choice": {"qk": {"encoding": "global", "bandwidth": bq[0], "reps": bq[1], "reg": bq[2]},
                      "rbf8": {"gamma": b8[0], "reg": b8[1]}, "rbf17": {"gamma": b17[0], "reg": b17[1]}},
           "cv": {"quantum": qk[bq], "rbf8": r8[b8], "rbf17": r17[b17]},
           "edge": {"qk_bandwidth": bq[0] in (BW[0], BW[-1]), "qk_reg": bq[2] in (regs[0], regs[-1]),
                    "rbf8_gamma": b8[0] in (RBF_GAMMAS[0], RBF_GAMMAS[-1]), "rbf17_gamma": b17[0] in (RBF_GAMMAS[0], RBF_GAMMAS[-1])}}
    print(name, task, len(y), {k: round(v, 3) for k, v in out["cv"].items()}, out["choice"]["qk"], flush=True)
    return out


def score_test(name: str, task: str, Xtr, ytr, Xte, yte, cfg: dict, margin: float) -> dict:
    m = A.EndpointModel(name, task, A.Prep.fit(Xtr), Xtr, ytr, cfg).fit()
    sc = m.predict(Xte)
    met = {k: A.score(task, yte, v) for k, v in sc.items()}
    ci = G._boot(task, yte, sc)
    floor = 0.5 if task == "cls" else 0.0
    ok = bool(met["quantum"] >= met["rbf8"] - margin and ci["quantum"][0] > floor)
    return {"endpoint": name, "task": task, "train": int(len(ytr)), "test": int(len(yte)), "metric": met, "ci95": ci,
            "pass": ok}


def _part(mode: str, name: str, fn, *args) -> dict:
    """Run fn(*args) once per (mode, endpoint) and keep the result on disk."""
    p = PARTS / f"{mode}_{name}.json"
    if p.exists():
        return json.loads(p.read_text(encoding="utf-8"))
    r = fn(*args)
    PARTS.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(r, default=float), encoding="utf-8")
    return r


def _posthoc_part(name: str) -> dict:
    return _part("posthoc", name, _posthoc_one, name)


def _posthoc_one(name: str) -> dict:
    task = A.ENDPOINTS[name][0]
    X, y = G.data(name, "train_val")
    cv = cv_one((name, X, y, task))
    import frozen as F
    cfg = {**F.QUANTUM_ADMET["settings"][name], "qk": cv["choice"]["qk"]}
    Xte, yte = G.data(name, "test")
    r = score_test(name, task, X, y, Xte, yte, cfg, F.QUANTUM_ADMET["margin"])
    r["cv"], r["qk"] = cv["cv"], cv["choice"]["qk"]
    r["q8_quantum"] = json.loads((ROOT / "results_admet.json").read_text(encoding="utf-8"))["endpoints"][name]["metric"]["quantum"]
    print(name, "post-hoc", {k: round(v, 3) for k, v in r["metric"].items()}, "Q8 quantum", round(r["q8_quantum"], 3),
          "meets Q8's rule" if r["pass"] else "misses", flush=True)
    return r


def posthoc() -> None:
    """Q8's 21 endpoints (all seen: Q8's verdict stands) with the global encoding, settings by the same CV."""
    with ProcessPoolExecutor(max_workers=12) as ex:
        res = list(ex.map(_posthoc_part, list(A.ENDPOINTS)))
    out = {"made": _now(), "note": "Post-hoc, on Q8's seen test splits: the reason for Q8b, not evidence. Q8 stays a fail.",
           "encoding": "global", "endpoints": {r["endpoint"]: r for r in res},
           "meets_q8_rule": int(sum(r["pass"] for r in res)), "of": len(res)}
    (ROOT / "results_admet_q8b_posthoc.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("post-hoc:", out["meets_q8_rule"], "of", out["of"], "meet Q8's rule with the global encoding")


def _practice_part(name: str) -> dict:
    def run():
        X, y = data(name, "train")
        return cv_one((name, X, y, ENDPOINTS[name][2]))
    return _part("practice", name, run)


def practice() -> None:
    with ProcessPoolExecutor(max_workers=12) as ex:
        res = list(ex.map(_practice_part, list(ENDPOINTS)))
    out = {"made": _now(), "cap": A.CAP, "grid": {"qk_bandwidth": BW, "reps": [1, 2], "rbf_gamma": RBF_GAMMAS,
                                                   "regs": REGS}, "endpoints": {r["endpoint"]: r for r in res}}
    (ROOT / "results_admet_q8b_practice.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")


def _test_one(args) -> dict:
    name, cfg, margin = args
    task = ENDPOINTS[name][2]
    Xtr, ytr = data(name, "train")
    Xte, yte = data(name, "test")
    r = score_test(name, task, Xtr, ytr, Xte, yte, cfg, margin)
    print(name, {k: round(v, 3) for k, v in r["metric"].items()}, "pass" if r["pass"] else "fail", flush=True)
    return r


def test() -> None:
    import frozen as F
    R = F.QUANTUM_ADMET_Q8B
    path = ROOT / "results_admet_q8b.json"
    if path.exists():
        sys.exit("results_admet_q8b.json exists: Gate Q8b runs once.")
    jobs = [(n, R["settings"][n], R["margin"]) for n in R["endpoints"]]
    with ProcessPoolExecutor(max_workers=7) as ex:
        res = list(ex.map(_test_one, jobs))
    passed = sum(r["pass"] for r in res)
    out = {"made": _now(), "rule": R, "endpoints": {r["endpoint"]: r for r in res}, "passed": int(passed),
           "of": len(res), "pass": bool(passed >= R["min_endpoints"])}
    path.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("Gate Q8b:", "pass" if out["pass"] else "fail", f"{passed} of {len(res)}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--posthoc", action="store_true")
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    posthoc() if a.posthoc else practice() if a.practice else test()


if __name__ == "__main__":
    main()
