"""
Gate 2 (Pillar 2): are the population model's stated intervals honest?

For every pair of loci the ensemble predicts a distribution of distances across cells (a Maxwell law
with per-axis sd sigma_ij). A stated 50 %, 80 % or 90 % central interval should contain that share of
the single-copy distances measured in held-out cells (half B). This script measures that coverage,
the PIT histogram (u = F_model(d)), and coverage by genomic separation.

Recalibration (fitted on PRACTICE data only): quantile recalibration (Kuleshov, Fenner & Ermon,
ICML 2018). With G the empirical CDF of the pooled practice PIT values, the interval at stated level
p is taken between model quantiles G^-1((1-p)/2) and G^-1((1+p)/2). The map is frozen in
chronocell/data/calibration.json and applied unchanged to the test datasets.

Per-bead reliability (population.bead_reliability, computed from the input alone) is checked on
practice data against each bead's held-out error (median |log(model / measured)| over its pairs).

    python validation/calibration.py --practice      # fit recalibration -> chronocell/data/calibration.json
    python validation/calibration.py --test          # run once -> validation/results_calibration.json

Gate 2b, sequencing Hi-C input (pre-registered in frozen.HIC_CALIBRATION): the same coverage and the
same recalibration method, with the model built from Rao et al. 2014 Hi-C on the imaged loci (the Gate 1
Hi-C settings and the app's b0 anchor) instead of half A's imaging contacts.

    python validation/calibration.py --practice --input hic   # -> chronocell/data/calibration_hic.json
    python validation/calibration.py --test --input hic       # run once -> validation/results_calibration_hic.json
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
from chronocell import ensemble as E, population as P   # noqa: E402

CAL_PATH = ROOT.parent / "chronocell" / "data" / "calibration.json"
CAL_HIC_PATH = ROOT.parent / "chronocell" / "data" / "calibration_hic.json"
PRACTICE = ["bintu_k562_28_30", "bintu_hct116_28_30", "bintu_hct116_28_30_auxin", "bintu_hct116_34_37", "su_chr2",
            "su_chr2_parm_rep"]
TEST = ["bintu_imr90_28_30", "bintu_imr90_18_20", "bintu_a549_28_30", "bintu_hct116_34_37_auxin", "su_chr21",
        "su_chr21_rep"]
SEP_BINS_BP = [0, 100_000, 300_000, 1_000_000, 3_000_000, 10_000_000, 300_000_000]


def _fit(key: str, freq, seen, r_c, n, split: int):
    """The app's model for a dataset of this size: v3.3 for <= 400 loci, v4 whole-chromosome above."""
    if n <= 400:
        return E.fit_ensemble(freq, seen, r_c_nm=r_c, cfg=E.EnsembleConfig(seed=split))
    from frozen import WHOLE_CHROMOSOME
    return P.fit_population(freq, seen, r_c_nm=r_c, cfg=P.PopulationConfig(seed=split, **WHOLE_CHROMOSOME))


def run_dataset(key: str, split: int, interval_pits=None) -> dict:
    tr = D.load(key)
    a, b = PR.split(tr.n_copies, split)
    r_c = 150.0 if D.REGISTRY[key].kind == "bintu_csv" else None
    ha = PR.half_stats(tr.xyz[a], r_c)
    r_c = ha.adjacent_median if r_c is None else r_c
    hb = PR.half_stats(tr.xyz[b], r_c)
    n = tr.n_loci
    res = _fit(key, ha.freq, ha.seen, r_c, n, split)
    iu = np.triu_indices(n, 1)
    sigma = np.zeros((n, n))
    sigma[iu] = P.pair_sigma_nm(res, iu[0], iu[1])
    sigma = sigma + sigma.T
    mask = np.triu(np.ones((n, n), bool), 1)
    sep = PR.separation(n, tr.starts)
    raw = PR.coverage(tr.xyz[b], sigma, mask)
    by_sep = []
    for lo, hi in zip(SEP_BINS_BP[:-1], SEP_BINS_BP[1:]):
        m = mask & (sep >= lo) & (sep < hi)
        if m.sum() < 20:
            continue
        c = PR.coverage(tr.xyz[b], sigma, m)
        row = {"sep_bp": [lo, hi], "pairs": int(m.sum()), "coverage": c["coverage"], "measurements": c["measurements"]}
        if interval_pits is not None:
            row["coverage_recalibrated"] = PR.coverage_with_pits(tr.xyz[b], sigma, m, interval_pits)["coverage"]
        by_sep.append(row)
    out = {"dataset": key, "split": split, "loci": n, "r_c_nm": r_c, "model": res.config.get("model", "ensemble_v3_3"),
           "coverage": raw["coverage"], "levels": raw["levels"], "measurements": raw["measurements"],
           "pit_histogram": raw["pit_histogram"], "by_separation": by_sep}
    if interval_pits is not None:
        out["coverage_recalibrated"] = PR.coverage_with_pits(tr.xyz[b], sigma, mask, interval_pits)["coverage"]
    # per-bead reliability vs held-out error
    rel = P.bead_reliability(res, ha.freq, ha.seen)
    err = np.abs(np.log(np.where(res.median_distance_nm > 0, res.median_distance_nm, np.nan)) -
                 np.log(np.where(hb.median > 0, hb.median, np.nan)))
    np.fill_diagonal(err, np.nan)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        bead_err = np.nanmedian(err, axis=1)
    out["reliability_vs_error_spearman"] = PR.spearman(rel, -bead_err)
    out["reliability_summary"] = {"min": float(np.min(rel)), "median": float(np.median(rel)), "max": float(np.max(rel))}
    out["spread_cv_model"] = float(res.spread_cv)
    return out


def _hic_counts(key: str, starts: np.ndarray) -> np.ndarray:
    """Rao et al. 2014 Hi-C of the dataset's cell line on its imaged loci (raw counts)."""
    e = D.REGISTRY[key]
    if e.kind == "bintu_csv":
        return D.bintu_hic(key)[0]
    H, _, hst = D.load_su_hic(D.REGISTRY[e.paired_hic])
    idx = np.searchsorted(hst, starts)
    if not (np.all(idx < len(hst)) and np.array_equal(hst[np.minimum(idx, len(hst) - 1)], starts)):
        raise KeyError(f"{key}: the Hi-C loci do not cover the imaged loci")
    return H[np.ix_(idx, idx)]


def interval_width(interval_pits=None, level: float = 0.9) -> float:
    """Upper / lower bound of the stated `level` interval (the same for every pair: both scale with sigma)."""
    if interval_pits is None:
        lo, hi = P.central_interval(level)
    else:
        lo, hi = (float(P.maxwell_quantile(q)) for q in interval_pits(level))
    return float(hi / lo) if lo > 0 else float("inf")


def run_dataset_hic(key: str, split: int, imaging_pits=None, hic_pits=None) -> dict:
    """Gate 2b: coverage of the stated intervals when the model is built from sequencing Hi-C."""
    from benchmark.run import _hic_input
    from chronocell import physics
    from frozen import HIC_CALIBRATION as HC
    tr = D.load(key)
    _, b = PR.split(tr.n_copies, split)
    e = D.REGISTRY[key]
    p, n_eff, _ = _hic_input(_hic_counts(key, tr.starts))
    step = e.step_bp if e.kind != "bintu_csv" else 30_000
    r_h = physics.bond_length_for(step) / float(E.gaussian_median_distance(float(HC["p_adjacent"]), 1.0))
    n = tr.n_loci
    res = _fit(key, p, n_eff, r_h, n, split)
    iu = np.triu_indices(n, 1)
    sigma = np.zeros((n, n))
    sigma[iu] = P.pair_sigma_nm(res, iu[0], iu[1])
    sigma = sigma + sigma.T
    mask = np.triu(np.ones((n, n), bool), 1)
    sep = PR.separation(n, tr.starts)
    raw = PR.coverage(tr.xyz[b], sigma, mask)
    hb_med = PR.half_stats(tr.xyz[b], None).median
    ratio = res.median_distance_nm[mask] / hb_med[mask]
    out = {"dataset": key, "split": split, "loci": n, "input": "sequencing Hi-C (Rao et al. 2014)", "r_c_nm": r_h,
           "n_eff": n_eff, "model": res.config.get("model", "ensemble_v3_3"), "coverage": raw["coverage"],
           "levels": raw["levels"], "measurements": raw["measurements"], "pit_histogram": raw["pit_histogram"],
           "median_scale_model_over_real": float(np.nanmedian(ratio[np.isfinite(ratio) & (ratio > 0)])),
           "width_90_raw": interval_width(None)}
    by_sep = []
    for lo, hi in zip(SEP_BINS_BP[:-1], SEP_BINS_BP[1:]):
        m = mask & (sep >= lo) & (sep < hi)
        if m.sum() < 20:
            continue
        row = {"sep_bp": [lo, hi], "pairs": int(m.sum()), "coverage": PR.coverage(tr.xyz[b], sigma, m)["coverage"]}
        if hic_pits is not None:
            row["coverage_recalibrated"] = PR.coverage_with_pits(tr.xyz[b], sigma, m, hic_pits)["coverage"]
        by_sep.append(row)
    out["by_separation"] = by_sep
    if imaging_pits is not None:
        out["coverage_imaging_recalibration"] = PR.coverage_with_pits(tr.xyz[b], sigma, mask, imaging_pits)["coverage"]
    if hic_pits is not None:
        out["coverage_recalibrated"] = PR.coverage_with_pits(tr.xyz[b], sigma, mask, hic_pits)["coverage"]
        out["width_90_recalibrated"] = interval_width(hic_pits)
    return out


def _print_hic(r: dict) -> None:
    print(f"{r['dataset']:28s} Hi-C input  raw {np.round(r['coverage'], 3)}  imaging recal "
          f"{np.round(r.get('coverage_imaging_recalibration', []), 3)}  Hi-C recal {np.round(r.get('coverage_recalibrated', []), 3)}"
          f"  size ratio {r['median_scale_model_over_real']:.2f}", flush=True)


def main_hic(practice: bool) -> None:
    """Gate 2b (frozen.HIC_CALIBRATION): fit the Hi-C recalibration on practice data, or run the test once."""
    from frozen import HIC_CALIBRATION as HC
    cal_img = json.loads(CAL_PATH.read_text())
    img_pits = PR.isotonic_recalibration(np.asarray(cal_img["pit_histogram"]))
    split = int(HC["split"])
    if practice:
        keys = HC["practice"]
        assert all(D.REGISTRY[k].role == "practice" for k in keys)
        rows = [run_dataset_hic(k, split, img_pits) for k in keys]
        per = [np.asarray(r["pit_histogram"]) / max(1.0, float(np.sum(r["pit_histogram"]))) for r in rows]
        balanced = np.mean(per, axis=0)
        cal = {"method": "quantile recalibration (Kuleshov, Fenner & Ermon, ICML 2018) on the PIT values of held-out "
                         "single-copy distances, model built from sequencing Hi-C",
               "pit_histogram": balanced.tolist(), "bins": len(balanced), "input": "sequencing Hi-C",
               "fitted_on": keys, "fitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
               "settings": {k: HC[k] for k in ("split", "p_adjacent", "hic_zeros", "anchor")},
               "note": "Practice datasets only; each dataset weighted equally. Pre-registered in validation/frozen.py "
                       "(HIC_CALIBRATION); its held-out test is Gate 2b in validation/RESULTS.md."}
        CAL_HIC_PATH.write_text(json.dumps(cal, indent=1) + "\n", encoding="utf-8")
        hic_pits = PR.isotonic_recalibration(balanced)
        rows = [run_dataset_hic(k, split, img_pits, hic_pits) for k in keys]          # in-sample after recalibration
        (ROOT / "results_calibration_hic_practice.json").write_text(json.dumps(
            {"mode": "practice (Hi-C recalibration fitted here; in-sample)", "rows": rows}, indent=1, default=float))
    else:
        cal = json.loads(CAL_HIC_PATH.read_text())
        hic_pits = PR.isotonic_recalibration(np.asarray(cal["pit_histogram"]))
        keys = HC["test"]
        assert all(D.REGISTRY[k].role == "test" for k in keys)
        rows = [run_dataset_hic(k, split, img_pits, hic_pits) for k in keys]
        lo90, hi90 = HC["pass_90"]
        lo50, hi50 = HC["pass_50"]
        per = {r["dataset"]: bool(lo90 <= r["coverage_recalibrated"][2] <= hi90 and lo50 <= r["coverage_recalibrated"][0] <= hi50)
               for r in rows}
        verdict = "pass" if all(per.values()) else "fail"
        (ROOT / "results_calibration_hic.json").write_text(json.dumps(
            {"mode": "test (run once; Hi-C recalibration frozen)", "calibration_used": cal, "rows": rows,
             "rule": {"pass_90": HC["pass_90"], "pass_50": HC["pass_50"], "every_dataset": True},
             "per_dataset_pass": per, "verdict": verdict,
             "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}, indent=1, default=float))
        print(f"verdict: {verdict} ({sum(per.values())} of {len(per)} datasets within the pre-registered ranges)")
    for r in rows:
        _print_hic(r)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    ap.add_argument("--splits", type=int, default=1)
    ap.add_argument("--datasets", nargs="*")
    ap.add_argument("--input", choices=("imaging", "hic"), default="imaging")
    a = ap.parse_args()
    if a.input == "hic":
        main_hic(a.practice)
        return
    if a.practice:
        keys = a.datasets or PRACTICE
        assert all(D.REGISTRY[k].role == "practice" for k in keys)
        rows = [run_dataset(k, s) for k in keys for s in range(a.splits)]
        pooled = np.sum([r["pit_histogram"] for r in rows], axis=0)
        # every practice dataset weighs the same in the recalibration, whatever its number of measurements
        per = [np.asarray(r["pit_histogram"]) / max(1.0, float(np.sum(r["pit_histogram"]))) for r in rows]
        balanced = np.mean(per, axis=0)
        cal = {"method": "quantile recalibration (Kuleshov, Fenner & Ermon, ICML 2018) on the PIT values of held-out "
                         "single-copy distances",
               "pit_histogram": balanced.tolist(), "bins": len(balanced),
               "fitted_on": keys, "fitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
               "note": "Practice datasets only; each dataset weighted equally. Applies to imaging-derived contact input "
                       "at the stated radius; not validated for sequencing Hi-C input."}
        CAL_PATH.write_text(json.dumps(cal, indent=1) + "\n", encoding="utf-8")
        pits = PR.isotonic_recalibration(np.asarray(balanced))
        rows = [run_dataset(k, s, pits) for k in keys for s in range(a.splits)]       # in-sample after recalibration
        res = {"mode": "practice (recalibration fitted here; in-sample)", "rows": rows,
               "pooled_pit_histogram": pooled.tolist()}
        (ROOT / "results_calibration_practice.json").write_text(json.dumps(res, indent=1, default=float))
        for r in rows:
            print(f"{r['dataset']:28s} raw {np.round(r['coverage'], 3)} recal {np.round(r['coverage_recalibrated'], 3)} "
                  f"rel~err rho {r['reliability_vs_error_spearman']:+.2f}", flush=True)
    else:
        cal = json.loads(CAL_PATH.read_text())
        pits = PR.isotonic_recalibration(np.asarray(cal["pit_histogram"]))
        keys = a.datasets or TEST
        assert all(D.REGISTRY[k].role == "test" for k in keys)
        rows = [run_dataset(k, s, pits) for k in keys for s in range(a.splits)]
        res = {"mode": "test (run once; recalibration frozen)", "calibration_used": cal, "rows": rows,
               "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
        out = ROOT / "results_calibration.json"
        if a.datasets:
            out = ROOT / f"results_calibration_{'_'.join(a.datasets)}.json"
        out.write_text(json.dumps(res, indent=1, default=float))
        for r in rows:
            print(f"{r['dataset']:28s} raw {np.round(r['coverage'], 3)} recal {np.round(r['coverage_recalibrated'], 3)} "
                  f"rel~err rho {r['reliability_vs_error_spearman']:+.2f}", flush=True)


if __name__ == "__main__":
    main()
