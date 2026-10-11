"""
Gate 4c (Phase B1): does the cohesin-depletion perturbation (chronocell/perturb.py, parameters fitted for
Gate 4 on imaging data and unchanged here) predict what RAD21 degradation does to sequencing Hi-C, on
regions it has never seen?

Data: HCT116 RAD21-mAC cells, untreated vs 6 h auxin (Rao et al., Cell 171:305, 2017), ENCODE GRCh38
maps read remotely by region (only the needed blocks are downloaded):
    in situ Hi-C   ENCSR123UVP untreated (ENCFF750AOC)  ->  ENCSR637QCS 6 h auxin (ENCFF301BWY)   [main]
    intact Hi-C    ENCSR958BEA untreated (ENCFF528XGK)  ->  ENCSR087JOM 6 h 5-Ph-IAA (ENCFF317OIA) [secondary]

Prediction (untreated counts only): Gate 1 Hi-C settings -> v3.3 population of the 2 Mb window at 10 kb ->
median distance map d -> perturb.cohesin_loss_map(d) = d' -> predicted change of contact probability,
log2 p(d') / p(d) with p the model's Gaussian contact law. Baseline "trend only": the same with lam = 0
(the fitted separation shift alone).
Measured change: log2 of (auxin count + 1) / auxin total over (untreated count + 1) / untreated total, on
pairs at least 2 bins apart with at least MIN_READS reads in the two maps together.
Score: Spearman between predicted and measured change, for the full model and for trend only; 95 %
interval of the difference by resampling loci (200 resamples).

    python validation/cohesin_hic.py --practice     # the Gate 4 practice region, chr21:28-30 Mb
    python validation/cohesin_hic.py --test         # run once (rule in frozen.COHESIN_HIC)
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                                  # noqa: E402
import phase_a                                        # noqa: E402,F401  (registers the ENCODE Hi-C sources)
from chronocell import ensemble as E, perturb as PT   # noqa: E402

ENC = "https://www.encodeproject.org/files/{0}/@@download/{0}.hic"
PAIRS = {"in situ": ("encode_hct116_insitu", "encode_hct116_rad21_auxin6h_insitu"),
         "intact": ("encode_hct116_rad21_untreated_intact", "encode_hct116_rad21_auxin6h_intact")}
for _k, _acc in (("encode_hct116_rad21_auxin6h_insitu", "ENCFF301BWY"),
                 ("encode_hct116_rad21_untreated_intact", "ENCFF528XGK"),
                 ("encode_hct116_rad21_auxin6h_intact", "ENCFF317OIA")):
    D.HIC_SOURCES.setdefault(_k, (ENC.format(_acc), "GRCh38", "HCT116"))
RES = 10_000
WINDOW = 2_000_000
MIN_READS = 10
MIN_SEP_BINS = 2
BOOT = 200
R_C = 150.0
PRACTICE = [("chr21", 28_000_000)]
OUT_PRACTICE = ROOT / "results_cohesin_hic_practice.json"
OUT = ROOT / "results_cohesin_hic.json"


def maps(pair: str, chrom: str, start: int) -> tuple[np.ndarray, np.ndarray]:
    u, a = PAIRS[pair]
    return (D.hic_region(u, chrom, start, start + WINDOW, RES).astype(float),
            D.hic_region(a, chrom, start, start + WINDOW, RES).astype(float))


def predicted_changes(untreated: np.ndarray, seed: int = 0) -> dict:
    """Fit the untreated window and return the predicted log2 contact change, full model and trend only."""
    from benchmark.run import _hic_input
    params = PT.cohesin_params()
    if params is None:
        raise SystemExit("chronocell/data/perturbation_params.json has no cohesin parameters.")
    n = len(untreated)
    freq, n_eff, _ = _hic_input(untreated)
    res = E.fit_ensemble(freq, n_eff, r_c_nm=R_C, cfg=E.EnsembleConfig(seed=seed))
    d = np.asarray(res.median_distance_nm, dtype=np.float64)
    sep = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) * float(RES)
    trend_only = PT.CohesinParams(0.0, params.c, params.log_sep_range)

    def p_of(dist: np.ndarray) -> np.ndarray:
        sigma = np.where(dist > 0, dist, 1e-9) / (E.MAXWELL_MEDIAN * R_C)
        return E.contact_probability_from_sigma(sigma)

    p0 = p_of(d)
    out = {}
    for name, prm in (("model", params), ("trend_only", trend_only)):
        with np.errstate(divide="ignore", invalid="ignore"):
            out[name] = np.log2(p_of(PT.cohesin_loss_map(d, sep, prm)) / p0)
    out["contact_fit"] = float(res.contact_fit) if np.isfinite(res.contact_fit) else None
    out["device"] = res.config.get("device_used", "cpu")
    return out


def measured_change(u: np.ndarray, a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = len(u)
    sep = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    upper = np.triu(np.ones((n, n), bool), MIN_SEP_BINS)
    tu, ta = u[upper].sum(), a[upper].sum()
    mask = upper & (u + a >= MIN_READS) & (sep >= MIN_SEP_BINS)
    return np.log2((a + 1.0) / ta) - np.log2((u + 1.0) / tu), mask


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    from scipy.stats import spearmanr
    ok = np.isfinite(x) & np.isfinite(y)
    return float(spearmanr(x[ok], y[ok]).statistic) if ok.sum() > 10 else float("nan")


def score_region(pred: dict, meas: np.ndarray, mask: np.ndarray, seed: int = 0) -> dict:
    n = len(meas)
    ii, jj = np.nonzero(mask)
    m = meas[ii, jj]
    full, trend = pred["model"][ii, jj], pred["trend_only"][ii, jj]
    out = {"pairs": int(len(ii)), "model": _spearman(full, m), "trend_only": _spearman(trend, m)}
    out["difference"] = out["model"] - out["trend_only"]
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(BOOT):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)   # loci resampled with replacement
        pw = w[ii] * w[jj]
        keep = pw > 0
        rep = np.repeat(np.flatnonzero(keep), pw[keep].astype(int))
        if len(rep) < 20:
            continue
        diffs.append(_spearman(full[rep], m[rep]) - _spearman(trend[rep], m[rep]))
    out["difference_ci95"] = [float(v) for v in np.nanpercentile(diffs, [2.5, 97.5])]
    return out


def run_region(pair: str, chrom: str, start: int) -> dict:
    u, a = maps(pair, chrom, start)
    pred = predicted_changes(u)
    meas, mask = measured_change(u, a)
    s = score_region(pred, meas, mask)
    s.update(region=f"{chrom}:{start / 1e6:g}-{(start + WINDOW) / 1e6:g} Mb", pair=pair, contact_fit=pred["contact_fit"],
             reads_untreated=float(u.sum()), reads_auxin=float(a.sum()), device=pred["device"])
    return s


def practice() -> None:
    rows = []
    for pair in PAIRS:
        for chrom, start in PRACTICE:
            r = run_region(pair, chrom, start)
            rows.append(r)
            print(f"{pair:8s} {r['region']:22s} model {r['model']:+.3f}  trend only {r['trend_only']:+.3f}  difference "
                  f"{r['difference']:+.3f} [{r['difference_ci95'][0]:+.3f}, {r['difference_ci95'][1]:+.3f}]  pairs {r['pairs']:,}",
                  flush=True)
    OUT_PRACTICE.write_text(json.dumps({"mode": "practice (the Gate 4 practice region; nothing fitted)", "rows": rows,
                                        "settings": {"res": RES, "window": WINDOW, "min_reads": MIN_READS,
                                                     "min_sep_bins": MIN_SEP_BINS, "boot": BOOT, "r_c_nm": R_C},
                                        "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")},
                                       indent=1, default=float), encoding="utf-8")


def test() -> None:
    from frozen import COHESIN_HIC as R
    if OUT.exists():
        sys.exit(f"{OUT.name} exists: Gate 4c runs once.")
    assert (R["res"], R["window"], R["min_reads"], R["min_sep_bins"], R["boot"], R["r_c_nm"]) == \
        (RES, WINDOW, MIN_READS, MIN_SEP_BINS, BOOT, R_C), "settings differ from the pre-registration"
    rows = []
    for pair in R["pairs"]:
        for chrom, start in R["regions"]:
            part = D.DATA / "cache" / f"gate4c_{pair.replace(' ', '')}_{chrom}_{start}.json"
            if part.exists():                       # an interrupted run resumes region by region
                r = json.loads(part.read_text())
            else:
                r = run_region(pair, chrom, start)
                part.parent.mkdir(parents=True, exist_ok=True)
                part.write_text(json.dumps(r, default=float))
            r["beats_trend"] = bool(r["difference"] > 0 and r["difference_ci95"][0] > 0)
            r["tier"] = "main" if pair == R["main_pair"] else "secondary"
            rows.append(r)
            print(f"{pair:8s} {r['region']:22s} model {r['model']:+.3f}  trend only {r['trend_only']:+.3f}  difference "
                  f"{r['difference']:+.3f} [{r['difference_ci95'][0]:+.3f}, {r['difference_ci95'][1]:+.3f}]  "
                  f"{'beats' if r['beats_trend'] else 'does not beat'} trend only", flush=True)
    n_main = sum(r["beats_trend"] for r in rows if r["tier"] == "main")
    verdict = "pass" if n_main >= R["min_regions"] else "fail"
    OUT.write_text(json.dumps({"mode": "test (run once; Gate 4 parameters unchanged)", "rule": R, "rows": rows,
                               "regions_beating_trend": n_main, "verdict": verdict,
                               "params": PT.load_params().get("cohesin_loss"),
                               "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")},
                              indent=1, default=float), encoding="utf-8")
    print(f"Gate 4c: {n_main} of {len(R['regions'])} held-out regions beat trend only (needed {R['min_regions']}) -> {verdict}")


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
