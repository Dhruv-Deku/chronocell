"""
Gate 5b (Phase A5): prediction without contact data, with more inputs.

Gate 5's model (CTCF peaks oriented by the JASPAR motif + GC) recovers 12-27 % of the reproducible
pattern, and the cohesin control showed its signal is compartment / insulation level. Phase A5 adds
ENCODE ATAC-seq (accessibility), H3K27ac (active chromatin) and RAD21 (cohesin) peaks of the same cell
line. Per locus and mark: s = log(1 + summed signalValue of the peaks whose summit lies in the locus).
Per pair and mark: similarity (s_i - m)(s_j - m), anchor s_i + s_j, each also times log10 separation;
for RAD21 also log(1 + RAD21 peaks between the loci) and its product with log10 separation. These come on
top of Gate 5's 12 features; the same ridge regression of the trend-removed log median distance.

Training: the untreated non-IMR-90 practice regions (K562 chr21:28-30, HCT116 chr21:28-30 and 34-37), so
every test set is a held-out cell type (Gate 5 also trained on IMR-90 chr2). Feature set and ridge
penalty by leave-one-dataset-out on those three. The frozen Gate 5 model is not changed.

Optional sequence-to-structure engines (Akita, C.Origami, Orca) are not run here; the reasons are in
NOT_RUN and in validation/RESULTS.md.

    python validation/predictor_v2.py --practice
    python validation/predictor_v2.py --test        # Gate 5b, run once (rule in frozen.PREDICTOR_V2)
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D            # noqa: E402
import predictor as G5          # noqa: E402
import protocol as PR           # noqa: E402
from chronocell import predict as PD   # noqa: E402

MODEL_PATH = ROOT / "predictor_v2_model.json"
TRAIN = ["bintu_k562_28_30", "bintu_hct116_28_30", "bintu_hct116_34_37"]
MARKS = ("atac", "h3k27ac", "rad21")
FEATURE_SETS = {"ctcf": (), "ctcf+atac": ("atac",), "ctcf+h3k27ac": ("h3k27ac",), "ctcf+rad21": ("rad21",),
                "ctcf+all": MARKS}
RIDGE = [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0]
NOT_RUN = {
    "Akita (Fudenberg et al., Nat Methods 2020)": "TensorFlow model; TensorFlow has no native GPU support on Windows, "
        "and its training data include IMR-90 and HCT116 Hi-C genome-wide, so the chr21 test regions cannot be "
        "guaranteed held out.",
    "C.Origami (Tan et al., Nat Biotechnol 2023)": "needs pyBigWig (no Windows build here) and was trained on IMR-90 "
        "Hi-C including chr21, so an IMR-90 chr21 test would be in-sample.",
    "Orca (Zhou, Nat Genet 2022)": "predicts H1 and HFF Micro-C only (not the imaged cell lines) and needs the Selene "
        "stack and multi-GB weights; not attempted here.",
}


def read_marks(path: Path, chrom: str) -> tuple[np.ndarray, np.ndarray]:
    """Summit positions and signalValue of narrowPeak peaks on one chromosome."""
    pos, sig = [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[0] != chrom:
                continue
            a, summit = int(f[1]), int(f[9])
            pos.append(a + summit if summit >= 0 else (a + int(f[2])) // 2)
            sig.append(float(f[6]))
    o = np.argsort(pos)
    return np.asarray(pos, np.int64)[o], np.asarray(sig)[o]


def mark_features(d: dict, cell: str, marks) -> np.ndarray:
    """Extra pair features for the design d of validation/predictor.py (pairs d['i'] < d['j'])."""
    L = d["loci"]
    i, j = d["i"], d["j"]
    ls = np.log10(np.maximum(d["sep"].astype(float), 1.0))
    cols = []
    for m in marks:
        pos, sig = read_marks(D.mark_peaks_path(m, cell), d["chrom"])
        csum = np.concatenate([[0.0], np.cumsum(sig)])
        a = np.searchsorted(pos, L.starts)
        b = np.searchsorted(pos, L.ends)
        s = np.log1p(csum[b] - csum[a])
        sim = (s[i] - s.mean()) * (s[j] - s.mean())
        anchor = s[i] + s[j]
        cols += [sim, anchor, sim * ls, anchor * ls]
        if m == "rad21":
            between = np.log1p(np.maximum(np.searchsorted(pos, L.starts[j]) - np.searchsorted(pos, L.ends[i]), 0))
            cols += [between, between * ls]
    return np.column_stack(cols) if cols else np.zeros((len(i), 0))


def build(key: str, marks) -> dict:
    tr = D.load(key)
    d = G5.design(key, tr)
    d["X"] = np.concatenate([d["X"], mark_features(d, d["cell"], marks)], axis=1)
    return d


def practice() -> None:
    for k in TRAIN:
        assert D.REGISTRY[k].role == "practice", k
    table = {}
    best = (None, None, -np.inf)
    for fs, marks in FEATURE_SETS.items():
        data = {}
        for k in TRAIN:
            d = build(k, marks)
            med = G5.all_copy_medians(k, D.load(k))
            log_d = np.log(np.where(med[d["i"], d["j"]] > 0, med[d["i"], d["j"]], np.nan))
            d["log_d"], d["y"] = log_d, PD.trend_residual(log_d, d["sep"])
            d["ok"] = np.isfinite(d["y"])
            data[k] = d
        for lam in RIDGE:
            rhos = {}
            for held in TRAIN:
                beta, mu, sd = PD.fit_ridge([data[k]["X"][data[k]["ok"]] for k in TRAIN if k != held],
                                            [data[k]["y"][data[k]["ok"]] for k in TRAIN if k != held], lam)
                h = data[held]
                rhos[held] = PR.spearman(((h["X"][h["ok"]] - mu) / sd) @ beta, h["y"][h["ok"]])
            mean = float(np.mean(list(rhos.values())))
            table[f"{fs}|{lam:g}"] = {"held_out_trend_removed_rho": rhos, "mean": mean}
            print(f"{fs:14s} ridge {lam:<7g} " + " ".join(f"{k} {v:+.3f}" for k, v in rhos.items()) + f" | mean {mean:+.3f}",
                  flush=True)
            if mean > best[2]:
                best = (fs, lam, mean)
        if fs == best[0]:
            best_data = data
    fs, lam, _ = best
    data = best_data
    beta, mu, sd = PD.fit_ridge([data[k]["X"][data[k]["ok"]] for k in TRAIN], [data[k]["y"][data[k]["ok"]] for k in TRAIN], lam)
    trend = PD.fit_trend([data[k]["log_d"] for k in TRAIN], [data[k]["sep"].astype(float) for k in TRAIN])
    out = {"feature_set": fs, "marks": list(FEATURE_SETS[fs]), "ridge_lambda": lam, "beta": beta.tolist(), "mu": mu.tolist(),
           "sd": sd.tolist(), "trend": list(trend), "trained_on": TRAIN,
           "mark_files": {f"{m}_{c}": D.MARK_PEAKS[(m, c)][0] for m in FEATURE_SETS[fs] for c in ("K562", "HCT116", "IMR90", "A549")},
           "fitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "note": "Practice data only (non-IMR-90 practice regions). Held-out test: Gate 5b (validation/RESULTS.md)."}
    MODEL_PATH.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    (ROOT / "results_predictor_v2_practice.json").write_text(json.dumps(
        {"mode": "practice (leave-one-dataset-out over the three training regions)", "table": table,
         "chosen": {"feature_set": fs, "ridge_lambda": lam}, "not_run": NOT_RUN}, indent=1, default=float))
    print(f"chosen: {fs}, ridge {lam:g}")


def test() -> None:
    from frozen import PREDICTOR as G5CFG, PREDICTOR_V2 as R
    from validation.benchmark import methods as M, run as BR
    model = json.loads(MODEL_PATH.read_text(encoding="utf-8"))
    beta, mu, sd = np.asarray(model["beta"]), np.asarray(model["mu"]), np.asarray(model["sd"])
    a, b = model["trend"]
    g5 = json.loads((ROOT / "results_predictor.json").read_text(encoding="utf-8"))
    rows = []
    for k in list(G5CFG["test"]) + list(G5CFG["control"]):
        e = D.REGISTRY[k]
        assert e.role == "test", k
        tr = D.load(k)
        d = build(k, model["marks"])
        n = d["loci"].n
        logd = a + b * np.log10(np.maximum(d["sep"].astype(float), 1.0)) + ((d["X"] - mu) / sd) @ beta
        pred = np.exp(G5._matrix(n, d["i"], d["j"], logd))
        np.fill_diagonal(pred, 0.0)
        sep_m = PR.separation(n, tr.starts)
        r_c = 150.0 if e.kind == "bintu_csv" else None
        masks = BR._masks(n, PR.tiles(n, 400))
        for split in range(3 if e.kind == "bintu_csv" else 1):
            st = BR._stats(k, tr.xyz, split, r_c)
            sc = BR._score_pred(M.Prediction(pred), st, sep_m, masks, tr.xyz[st["b_idx"]], BR.BOOT.get(e.kind, 30), None,
                                split == 0, split)
            rows.append({"dataset": k, "role": "control" if k in G5CFG["control"] else "test", "split": split,
                         "method": "predictor_v2", "loci": n, "status": "ok", "input": "sequence + CTCF + marks", **sc})
            ap = sc.get("all_pairs", {})
            print(f"{k:26s} s{split} all {ap.get('percent_of_ceiling', np.nan):6.1f}% raw {ap.get('raw_spearman', np.nan):.3f}",
                  flush=True)
    summ = BR.summarise(rows)
    per = {}
    for k in G5CFG["test"]:
        p = summ.get(k, {}).get("sequence + CTCF + marks", {}).get("all_pairs", {}).get("predictor_v2")
        old = g5["summary"][k]["sequence + CTCF"]["all_pairs"]["predictor"]["percent_of_ceiling"]
        ci = (p or {}).get("ci95", {}).get("percent_of_ceiling", [np.nan, np.nan])
        ok = bool(p and p["percent_of_ceiling"] >= R["min_percent"] and ci[0] > 0 and p["percent_of_ceiling"] > old)
        per[k] = {"percent_of_ceiling": p and p["percent_of_ceiling"], "ci95": ci, "gate5_percent": old, "passes": ok}
    n_pass = sum(v["passes"] for v in per.values())
    verdict = "pass" if n_pass >= R["min_sets"] else "fail"
    (ROOT / "results_predictor_v2.json").write_text(json.dumps(
        {"mode": "test (run once; model frozen)", "rule": R, "model": model, "summary": summ, "per_set": per,
         "datasets_passing": n_pass, "verdict": verdict, "rows": rows, "not_run": NOT_RUN,
         "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}, indent=1, default=float))
    print(json.dumps(per, indent=1, default=float))
    print(f"Gate 5b: {n_pass} of {len(per)} -> {verdict}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        practice() if a.practice else test()


if __name__ == "__main__":
    main()
