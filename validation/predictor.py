"""
Gate 5 (Pillar 5): can sequence and CTCF binding alone predict the pattern of 3D distances between
loci where no contact data exist?

    python validation/predictor.py --practice    # ridge penalty by leave-one-dataset-out CV, fit, freeze
    python validation/predictor.py --test        # run once -> validation/results_predictor.json

Everything is pre-registered in validation/frozen.py (PREDICTOR): model, inputs, training and test
datasets, scoring, baseline and pass rule. --practice touches practice datasets only and writes the
frozen model to validation/predictor_model.json (committed before --test runs). --test refuses to run
without that file and never changes it.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                              # noqa: E402
import protocol as PR                             # noqa: E402
from frozen import PREDICTOR as CFG               # noqa: E402
from chronocell import predict as PD              # noqa: E402

MODEL_PATH = ROOT / "predictor_model.json"
OUT = ROOT / "results_predictor.json"
TUNING_LOG = ROOT / "tuning_v4.json"
CACHE = ROOT / "data" / "cache"
_SEQ: dict = {}
_GC: dict = {}
_PEAKS: dict = {}


def _seq(chrom: str) -> bytes:
    if chrom not in _SEQ:
        _SEQ[chrom] = PD.read_fasta(D.hg38_fasta_path(chrom))
    return _SEQ[chrom]


def _gc(chrom: str) -> PD.GCIndex:
    if chrom not in _GC:
        _GC[chrom] = PD.GCIndex(_seq(chrom))
    return _GC[chrom]


def _pwm() -> np.ndarray:
    return PD.log_odds(PD.read_jaspar(D.jaspar_ctcf_path().read_text()), CFG["pseudocount"])


def peaks(cell: str, chrom: str) -> dict:
    """CTCF peaks of one cell line on one chromosome: midpoint (summit), motif strand, motif score."""
    key = (cell, chrom)
    if key in _PEAKS:
        return _PEAKS[key]
    import pandas as pd
    acc = D.CTCF_PEAKS[cell][0]
    cache = CACHE / f"ctcf_{acc}_{chrom}_oriented.npz"
    if cache.exists():
        z = np.load(cache)
        out = {k: z[k] for k in z.files}
    else:
        df = pd.read_csv(D.ctcf_peaks_path(cell), sep="\t", header=None, compression="gzip")
        df = df[df[0] == chrom]
        s, e, summit = df[1].to_numpy(np.int64), df[2].to_numpy(np.int64), df[9].to_numpy(np.int64)
        strand, score = PD.orient_peaks(_seq(chrom), s, e, summit, _pwm(), CFG["motif_min_relative_score"],
                                        CFG["summit_half_width_bp"])
        mid = np.where(summit >= 0, s + summit, (s + e) // 2)
        out = {"mid": mid, "strand": strand, "score": score, "signal": df[6].to_numpy(float)}
        CACHE.mkdir(parents=True, exist_ok=True)
        np.savez(cache, **out)
    _PEAKS[key] = out
    return out


def design(key: str, tr) -> dict:
    """Locus annotations and pair features (all pairs i < j) of one dataset."""
    e = D.REGISTRY[key]
    chrom = str(tr.chrom[0])
    starts = np.asarray(tr.starts, np.int64)
    ends = starts + int(e.locus_bp)
    pk = peaks(e.cell_line, chrom)
    loci = PD.annotate(starts, ends, pk["mid"], pk["strand"], _gc(chrom))
    i, j = np.triu_indices(loci.n, 1)
    X = PD.pair_features(loci, pk["mid"], i, j)
    sep = np.rint(np.abs((starts[j] + ends[j]) - (starts[i] + ends[i])) / 2).astype(np.int64)
    return {"loci": loci, "i": i, "j": j, "X": X, "sep": sep, "chrom": chrom, "cell": e.cell_line,
            "motif_rate": float(np.mean(pk["strand"] != 0)) if len(pk["strand"]) else float("nan"),
            "peaks_on_chrom": int(len(pk["mid"]))}


def all_copy_medians(key: str, tr) -> np.ndarray:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"pred_{key}_median_all.npz"
    if path.exists():
        return np.load(path)["median"]
    med = PR.half_stats(tr.xyz, None).median
    np.savez(path, median=med)
    return med


def _log(entry: dict) -> None:
    try:
        log = json.loads(TUNING_LOG.read_text())
    except (OSError, ValueError):
        log = []
    log.append(entry)
    TUNING_LOG.write_text(json.dumps(log, indent=1, default=float))


# ======================================================================================
# Practice: choose the ridge penalty, fit, freeze
# ======================================================================================
def practice() -> None:
    for k in CFG["train"]:
        if D.REGISTRY[k].role != "practice":
            sys.exit(f"{k} is not a practice dataset; refusing.")
    data = {}
    for k in CFG["train"]:
        t0 = time.time()
        tr = D.load(k)
        d = design(k, tr)
        med = all_copy_medians(k, tr)
        log_d = np.log(np.where(med[d["i"], d["j"]] > 0, med[d["i"], d["j"]], np.nan))
        d["log_d"] = log_d
        d["y"] = PD.trend_residual(log_d, d["sep"])
        ok = np.isfinite(d["y"])
        d["ok"] = ok
        data[k] = d
        L = d["loci"]
        print(f"{k:24s} {d['cell']:7s} loci {L.n:4d} pairs {ok.sum():7,d} CTCF peaks on {d['chrom']}: {d['peaks_on_chrom']:,} "
              f"(motif >= {CFG['motif_min_relative_score']}: {100 * d['motif_rate']:.0f} %) loci with a peak: "
              f"{np.mean((L.fwd + L.rev + L.unk) > 0) * 100:.0f} % ({time.time() - t0:.0f} s)", flush=True)
    cv = {}
    for lam in CFG["ridge_grid"]:
        rhos = {}
        for held in CFG["train"]:
            Xs = [data[k]["X"][data[k]["ok"]] for k in CFG["train"] if k != held]
            ys = [data[k]["y"][data[k]["ok"]] for k in CFG["train"] if k != held]
            beta, mu, sd = PD.fit_ridge(Xs, ys, lam)
            h = data[held]
            pred = ((h["X"][h["ok"]] - mu) / sd) @ beta
            rhos[held] = PR.spearman(pred, h["y"][h["ok"]])
        cv[lam] = rhos
        print(f"ridge {lam:g}: held-out trend-removed Spearman " +
              " ".join(f"{k} {v:+.3f}" for k, v in rhos.items()) + f" | mean {np.mean(list(rhos.values())):+.3f}",
              flush=True)
    best = max(cv, key=lambda lam: np.mean(list(cv[lam].values())))
    beta, mu, sd = PD.fit_ridge([data[k]["X"][data[k]["ok"]] for k in CFG["train"]],
                                [data[k]["y"][data[k]["ok"]] for k in CFG["train"]], best)
    trend = PD.fit_trend([data[k]["log_d"] for k in CFG["train"]], [data[k]["sep"].astype(float) for k in CFG["train"]])
    model = PD.Predictor(beta.tolist(), mu.tolist(), sd.tolist(), trend, float(best), trained_on=list(CFG["train"]),
                         motif=CFG["motif"],
                         notes=f"Fitted {dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')} on practice "
                               "data only (validation/predictor.py --practice).")
    MODEL_PATH.write_text(model.to_json() + "\n", encoding="utf-8")
    in_sample = {k: PR.spearman(model.residual(data[k]["X"][data[k]["ok"]]), data[k]["y"][data[k]["ok"]])
                 for k in CFG["train"]}
    _log({"experiment": "no_hic_predictor", "ridge_cv_lodo": {str(k): v for k, v in cv.items()}, "chosen": best,
          "in_sample_trend_removed_spearman": in_sample, "trend": trend,
          "coefficients": dict(zip(model.names, model.beta)),
          "motif_rates": {k: data[k]["motif_rate"] for k in CFG["train"]}})
    print(f"chosen ridge {best:g}; trend log d = {trend[0]:.3f} + {trend[1]:.3f} log10 s; in-sample " +
          " ".join(f"{k} {v:+.3f}" for k, v in in_sample.items()))
    for name, b in zip(model.names, model.beta):
        print(f"  {name:22s} {b:+.4f}")
    print(f"-> {MODEL_PATH}")


# ======================================================================================
# Test: run once
# ======================================================================================
def _matrix(n: int, i: np.ndarray, j: np.ndarray, v: np.ndarray) -> np.ndarray:
    m = np.zeros((n, n))
    m[i, j] = v
    m[j, i] = v
    return m


def test() -> None:
    from validation.benchmark import methods as M, run as BR
    if not MODEL_PATH.exists():
        sys.exit("validation/predictor_model.json is missing: run --practice first (and commit it).")
    model = PD.Predictor.from_json(MODEL_PATH.read_text())
    rows, summary = [], {}
    keys = list(CFG["test"]) + list(CFG["control"])
    for k in keys:
        e = D.REGISTRY[k]
        if e.role != "test":
            sys.exit(f"{k} is not a test dataset; refusing.")
        tr = D.load(k)
        d = design(k, tr)
        n = d["loci"].n
        sep_m = PR.separation(n, tr.starts)
        r_c = 150.0 if e.kind == "bintu_csv" else None
        splits = 3 if e.kind == "bintu_csv" else 1
        tiles = PR.tiles(n, 400)
        masks = BR._masks(n, tiles)
        pred = np.exp(_matrix(n, d["i"], d["j"], model.log_distance(d["X"], d["sep"].astype(float))))
        a, b = model.trend
        base = np.exp(_matrix(n, d["i"], d["j"], a + b * np.log10(np.maximum(d["sep"], 1.0))))
        np.fill_diagonal(pred, 0.0)
        np.fill_diagonal(base, 0.0)
        for split in range(splits):
            st = BR._stats(k, tr.xyz, split, r_c)
            xyz_b = tr.xyz[st["b_idx"]]
            for name, mat in (("predictor", pred), ("genomic_trend_no_data", base)):
                sc = BR._score_pred(M.Prediction(mat), st, sep_m, masks, xyz_b, BR.BOOT.get(e.kind, 30), None,
                                    split == 0, split)
                rows.append({"dataset": k, "role": "control" if k in CFG["control"] else "test", "split": split,
                             "method": name, "loci": n, "status": "ok", "input": "sequence + CTCF", **sc})
                ap = sc.get("all_pairs", {})
                print(f"{k:26s} s{split} {name:22s} all {ap.get('percent_of_ceiling', np.nan):6.1f}% raw "
                      f"{ap.get('raw_spearman', np.nan):.3f} ccc {ap.get('lin_ccc_nm', np.nan):.3f}", flush=True)
    summ = BR.summarise(rows)
    passes = {}
    for k in CFG["test"]:
        s = summ.get(k, {}).get("sequence + CTCF", {}).get("all_pairs", {})
        p, g = s.get("predictor"), s.get("genomic_trend_no_data")
        if not p or not g:
            passes[k] = {"i": False, "ii": False, "note": "not scored"}
            continue
        ci = p.get("ci95", {}).get("percent_of_ceiling", [np.nan, np.nan])
        passes[k] = {"i": bool(p["percent_of_ceiling"] > 0 and ci[0] > 0),
                     "ii": bool(p["raw_spearman"] > g["raw_spearman"]),
                     "percent_of_ceiling": p["percent_of_ceiling"], "ci95": ci,
                     "raw_spearman": p["raw_spearman"], "baseline_raw_spearman": g["raw_spearman"]}
    n_pass = sum(v["i"] and v["ii"] for v in passes.values())
    verdict = "pass" if n_pass >= CFG["pass_min_datasets"] else "fail"
    out = {"model": json.loads(MODEL_PATH.read_text()), "settings": CFG, "summary": summ, "gate5": passes,
           "datasets_passing": n_pass, "verdict": verdict, "rows": rows,
           "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
    OUT.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps(passes, indent=1, default=float))
    print(f"Gate 5: {n_pass} of {len(CFG['test'])} test datasets pass both criteria -> {verdict.upper()}  -> {OUT}")


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
