"""
Gate 4 (Pillar 4), part 1: does the cohesin-depletion perturbation predict what RAD21 degradation
does to real chromatin?

Data (Bintu et al. 2018, HCT116-RAD21-mAC, chromatin tracing at 30 kb):
  PRACTICE  chr21:28-30 Mb untreated -> 6 h auxin       (the only pair used to fit parameters)
  TEST      chr21:34-37 Mb untreated -> 6 h auxin       (run once, parameters frozen)
The untreated 34-37 Mb data (a v3.3 practice set) are the INPUT for the test; the auxin 34-37 Mb data
are the held-out truth and were never used before the test run.

Prediction: untreated half A contacts (< 150 nm, as v3.3) -> v3.3 ensemble -> median map d ->
perturb.cohesin_loss_map(d) with parameters lam, c0..c2 fitted on the practice pair by linear least
squares of log d_pred against log d_auxin (pooled over the splits).

Metrics, fixed before the test run (pass criteria in brackets):
  change agreement  Spearman between predicted log change (d_pred / d_model) and measured log change
                    (auxin median / untreated half-B median)          [full model > trend-only]
  absolute          Lin's CCC and raw Spearman of d_pred vs the measured auxin medians, and trend-removed
                    Spearman as % of the auxin ceiling (auxin half A vs half B)
                                                                      [CCC full model > no change]
Baselines: no change (d_pred = d); trend only (lam = 0, shift fitted on practice); and a data-only
reference that needs imaging (measured untreated half-B medians + the same shift), to show how much
of the error comes from the reconstruction itself.

    python validation/perturbation.py --practice     # fit parameters -> chronocell/data/perturbation_params.json
    python validation/perturbation.py --test         # run once -> validation/results_perturbation.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D            # noqa: E402
import protocol as PR           # noqa: E402
from chronocell import ensemble as E, perturb as PT   # noqa: E402

R_C = 150.0
SPLITS = 3
PAIRS = {"practice": ("bintu_hct116_28_30", "bintu_hct116_28_30_auxin"),
         "test": ("bintu_hct116_34_37", "bintu_hct116_34_37_auxin")}


def _maps(untreated: str, auxin: str, split: int) -> dict:
    u, x = D.load(untreated), D.load(auxin)
    a, b = PR.split(u.n_copies, split)
    ha = PR.half_stats(u.xyz[a], R_C)
    hb = PR.half_stats(u.xyz[b], R_C)
    xa_i, xb_i = PR.split(x.n_copies, split)
    hx = PR.half_stats(x.xyz, R_C)
    hxa, hxb = PR.half_stats(x.xyz[xa_i], R_C), PR.half_stats(x.xyz[xb_i], R_C)
    res = E.fit_ensemble(ha.freq, ha.seen, r_c_nm=R_C, cfg=E.EnsembleConfig(seed=split))
    n = u.n_loci
    sep = PR.separation(n, u.starts)
    return {"model": res.median_distance_nm, "unt_b": hb.median, "aux": hx.median, "aux_a": hxa.median,
            "aux_b": hxb.median, "sep": sep, "n": n}


def _design(m: dict, use_lam: bool) -> tuple[np.ndarray, np.ndarray]:
    iu = np.triu_indices(m["n"], 1)
    log_d = np.log(np.where(m["model"] > 0, m["model"], 1.0))      # diagonal (0 nm) is never scored
    _, resid = PT.trend_residual(log_d, m["sep"])
    L = np.log10(m["sep"][iu].astype(float))
    cols = [np.ones_like(L), L, L ** 2] + ([-resid[iu]] if use_lam else [])
    X = np.stack(cols, axis=1)
    y = np.log(m["aux"][iu]) - log_d[iu]
    return X, y


def fit_params(maps: list[dict], use_lam: bool = True) -> PT.CohesinParams:
    X = np.concatenate([_design(m, use_lam)[0] for m in maps])
    y = np.concatenate([_design(m, use_lam)[1] for m in maps])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    seps = np.concatenate([m["sep"][np.triu_indices(m["n"], 1)] for m in maps]).astype(float)
    rng = (float(np.log10(seps.min())), float(np.log10(seps.max())))
    lam = float(coef[3]) if use_lam else 0.0
    return PT.CohesinParams(lam, (float(coef[0]), float(coef[1]), float(coef[2])), rng)


def evaluate(m: dict, params: dict[str, PT.CohesinParams]) -> dict:
    n = m["n"]
    iu = np.triu_indices(n, 1)
    mask = np.triu(np.ones((n, n), bool), 1)
    meas_change = np.log(m["aux"][iu]) - np.log(m["unt_b"][iu])
    preds = {"full_model": PT.cohesin_loss_map(m["model"], m["sep"], params["full"]),
             "trend_only": PT.cohesin_loss_map(m["model"], m["sep"], params["trend"]),
             "no_change": m["model"].copy(),
             "data_only_reference": PT.cohesin_loss_map(m["unt_b"], m["sep"], params["trend"])}
    ceiling = PR.scores(m["aux_a"], m["aux_b"], m["sep"], mask)["spearman_distance_corrected"]
    out = {"auxin_ceiling_trend_removed": ceiling}
    for name, p in preds.items():
        base = m["unt_b"] if name == "data_only_reference" else m["model"]
        change = np.log(p[iu]) - np.log(base[iu])
        s = PR.scores(p, m["aux"], m["sep"], mask)
        out[name] = {"change_spearman": PR.spearman(change, meas_change) if name != "no_change" else float("nan"),
                     "change_pearson": float(np.corrcoef(change, meas_change)[0, 1]) if name != "no_change" else float("nan"),
                     "change_rmse_log": float(np.sqrt(np.mean((change - meas_change) ** 2))),
                     "auxin_lin_ccc_nm": s["lin_ccc_nm"], "auxin_raw_spearman": s["spearman"],
                     "auxin_trend_removed_pct_of_ceiling": 100 * s["spearman_distance_corrected"] / ceiling,
                     "auxin_median_scale": s["median_scale_model_over_real"]}
    out["measured_change"] = {"mean_log": float(meas_change.mean()), "sd_log": float(meas_change.std())}
    return out


def _summary(rows: list[dict]) -> dict:
    keys = [k for k in rows[0] if isinstance(rows[0][k], dict) and k != "measured_change"]
    out = {}
    for k in keys:
        out[k] = {m: float(np.mean([r[k][m] for r in rows])) for m in rows[0][k]}
    out["auxin_ceiling_trend_removed"] = float(np.mean([r["auxin_ceiling_trend_removed"] for r in rows]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    if a.practice:
        unt, aux = PAIRS["practice"]
        maps = [_maps(unt, aux, s) for s in range(SPLITS)]
        full, trend = fit_params(maps, True), fit_params(maps, False)
        stamp = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        note = ("Fitted by validation/perturbation.py --practice on Bintu et al. 2018 HCT116 chr21:28-30 Mb, "
                "untreated vs 6 h auxin (practice pair), pooled over 3 splits. Not a first-principles constant.")
        params = {"cohesin_loss": {"lam": full.lam, "c": list(full.c), "log_sep_range": list(full.log_sep_range),
                                   "fitted_on": f"{unt} -> {aux}", "fitted_utc": stamp, "note": note},
                  "cohesin_loss_trend_only": {"lam": 0.0, "c": list(trend.c), "log_sep_range": list(trend.log_sep_range),
                                              "fitted_on": f"{unt} -> {aux}", "fitted_utc": stamp,
                                              "note": "Baseline: separation-dependent shift only."}}
        PT.PARAMS_PATH.write_text(json.dumps(params, indent=1) + "\n", encoding="utf-8")
        rows = [evaluate(m, {"full": full, "trend": trend}) for m in maps]
        res = {"mode": "practice (in-sample: these numbers are fits, not tests)", "params": params, "splits": rows,
               "summary": _summary(rows)}
        (ROOT / "results_perturbation_practice.json").write_text(json.dumps(res, indent=1, default=float))
        print(json.dumps(res["summary"], indent=1, default=float))
        print(f"lam = {full.lam:.3f}, c = {full.c}")
    else:
        params = {"full": PT.cohesin_params(), "trend": None}
        raw = PT.load_params().get("cohesin_loss_trend_only")
        if params["full"] is None or raw is None:
            sys.exit("Run --practice first (parameters are fitted on practice data only).")
        params["trend"] = PT.CohesinParams(0.0, tuple(raw["c"]), tuple(raw["log_sep_range"]))
        unt, aux = PAIRS["test"]
        if D.REGISTRY[aux].role != "test":
            sys.exit("misconfigured test pair")
        maps = [_maps(unt, aux, s) for s in range(SPLITS)]
        rows = [evaluate(m, params) for m in maps]
        res = {"mode": "test (held out; run once with frozen parameters)", "dataset": f"{unt} -> {aux}",
               "params_used": PT.load_params(), "splits": rows, "summary": _summary(rows),
               "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
        (ROOT / "results_perturbation.json").write_text(json.dumps(res, indent=1, default=float))
        print(json.dumps(res["summary"], indent=1, default=float))


if __name__ == "__main__":
    main()
