"""
Gate 2c (Pillar 2): can the input alone say WHICH distances will be wrong?

Gate 2 found that a per-bead score from the fit's misfit does not predict per-bead error (TUNING.md
section 9; RESULTS.md, Gate 2). This script tests per-PAIR scores, computed from the input and the
fitted model only, against each pair's held-out error, and their per-bead summaries.

Error of a pair (scale-free): e_ij = |r_ij - median of r over its separation stratum|, with
r_ij = log(model median / half-B measured median). Removing each stratum's median keeps systematic
size and trend biases (reported elsewhere: CCC, size ratio) out of it; what is left is "this pair is
off relative to pairs at the same separation", which is what a per-pair reliability can rank.
Strata: 10 equal-count bins of genomic separation.

Candidate scores (higher = more reliable), all from the input + model, never from half B:
  input_se    minus the delta-method standard error of the pair's log target distance from its contact
              count: |d log sigma / d log f| * sqrt((1 - f) / c); c = copies in contact (imaging) or
              reads (Hi-C); unobserved pairs get c = 0.5
  misfit      minus |log sigma_model^2 - log sigma_target^2| (how far the model is from this pair's own
              input); unobserved pairs rank lowest
  combined    minus sqrt(input_se^2 + misfit^2)
Pair-level score: Spearman of score vs -e within each stratum, averaged over strata. Bead-level: each
bead's median pair score vs minus its median pair error, Spearman over beads. 95 % intervals: 200
bootstrap resamples of the loci (pairs follow their loci).

    python validation/reliability.py --practice     # all candidates, practice data -> results_reliability_practice.json
    python validation/reliability.py --test         # the frozen candidate and rule (frozen.RELIABILITY), run once
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

import datasets as D            # noqa: E402
import protocol as PR           # noqa: E402
from calibration import PRACTICE as IMG_PRACTICE, TEST as IMG_TEST, _fit, _hic_counts   # noqa: E402
from chronocell import ensemble as E, population as P                                  # noqa: E402

STRATA = 10
BOOT = 200
CANDIDATES = ("input_se", "misfit", "combined")


def _log_sigma_slope(f: np.ndarray) -> np.ndarray:
    """|d log sigma / d log f| of the Maxwell inversion (1/3 for rare contacts, larger near f = 1)."""
    f = np.clip(f, 1e-6, 0.999)
    h = 0.01
    up = np.log(E.sigma_from_contact_probability(np.minimum(f * np.exp(h), 0.9995)))
    dn = np.log(E.sigma_from_contact_probability(f * np.exp(-h)))
    return np.abs(up - dn) / (np.log(np.minimum(f * np.exp(h), 0.9995)) - np.log(f * np.exp(-h)))


def scores(res, freq: np.ndarray, counts: np.ndarray, n_obs) -> dict[str, np.ndarray]:
    """Candidate reliability scores per pair (n x n, higher = more reliable), from the input and model only."""
    n = len(freq)
    f = np.asarray(freq, dtype=np.float64)
    c = np.nan_to_num(np.asarray(counts, dtype=np.float64), nan=0.0)
    observed = np.isfinite(f) & (c > 0)
    f_use = np.where(observed, f, np.nanmedian(f[np.isfinite(f)]) if np.isfinite(f).any() else 0.1)
    se = _log_sigma_slope(f_use) * np.sqrt(np.clip(1 - f_use, 1e-6, 1) / np.maximum(c, 0.5))
    log_t, _, obs_t = P.targets_from_frequency(f, n_obs)
    iu = np.triu_indices(n, 1)
    sig = P.pair_sigma_nm(res, iu[0], iu[1]) / float(res.config.get("r_c_nm", 150.0))
    log_m = np.zeros((n, n))
    log_m[iu] = 2.0 * np.log(np.maximum(sig, 1e-12))
    log_m = log_m + log_m.T
    mis = np.where(obs_t, np.abs(log_m - log_t.astype(np.float64)), np.inf)
    out = {"input_se": -se, "misfit": -mis, "combined": -np.sqrt(se ** 2 + mis ** 2)}
    for v in out.values():
        np.fill_diagonal(v, np.nan)
    return out


def _strata(sep: np.ndarray, iu) -> np.ndarray:
    s = sep[iu].astype(np.float64)
    edges = np.unique(np.quantile(s, np.linspace(0, 1, STRATA + 1)))
    return np.clip(np.searchsorted(edges, s, side="right") - 1, 0, len(edges) - 2)


def pair_error(model_med: np.ndarray, truth_med: np.ndarray, strata: np.ndarray, iu) -> np.ndarray:
    """Scale-free error per pair i < j: |log(model / truth) - its stratum median| (NaN if unmeasured)."""
    r = np.log(np.where(model_med[iu] > 0, model_med[iu], np.nan)) - np.log(np.where(truth_med[iu] > 0, truth_med[iu], np.nan))
    e = np.full_like(r, np.nan)
    for k in np.unique(strata):
        m = (strata == k) & np.isfinite(r)
        if m.any():
            e[m] = np.abs(r[m] - np.median(r[m]))
    return e


def stratified_spearman(score: np.ndarray, err: np.ndarray, strata: np.ndarray) -> float:
    """Mean over strata of Spearman(score, -error) within the stratum (pairs with an infinite score rank last)."""
    vals, wts = [], []
    for k in np.unique(strata):
        m = (strata == k) & np.isfinite(err) & ~np.isnan(score)
        if m.sum() < 10:
            continue
        s = np.where(np.isfinite(score[m]), score[m], np.nanmin(score[m][np.isfinite(score[m])]) - 1.0
                     if np.isfinite(score[m]).any() else 0.0)
        r = PR.spearman(s, -err[m])
        if np.isfinite(r):
            vals.append(r)
            wts.append(m.sum())
    return float(np.average(vals, weights=wts)) if vals else float("nan")


def bead_level(score_p: np.ndarray, err_p: np.ndarray, iu, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Per-bead median pair score and median pair error over the bead's partners."""
    S = np.full((n, n), np.nan)
    Er = np.full((n, n), np.nan)
    S[iu] = score_p
    Er[iu] = err_p
    S = np.where(np.isnan(S), S.T, S)
    Er = np.where(np.isnan(Er), Er.T, Er)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        bs = np.nanmedian(np.where(np.isinf(S), -1e6, S), axis=1)
        be = np.nanmedian(Er, axis=1)
    return bs, be


def evaluate(score_mat: np.ndarray, model_med: np.ndarray, truth_med: np.ndarray, sep: np.ndarray, seed: int = 0) -> dict:
    """Pair- and bead-level scores with 95 % locus-bootstrap intervals."""
    n = len(model_med)
    iu = np.triu_indices(n, 1)
    st = _strata(sep, iu)
    err = pair_error(model_med, truth_med, st, iu)
    sc = score_mat[iu]
    pair = stratified_spearman(sc, err, st)
    bs, be = bead_level(sc, err, iu, n)
    bead = PR.spearman(bs, -be)
    rng = np.random.default_rng(seed)
    bp, bb = [], []
    for _ in range(BOOT):
        loci = rng.integers(0, n, n)
        a, b = np.triu_indices(n, 1)
        ia, ib = loci[a], loci[b]
        keep = ia != ib
        ia, ib = np.minimum(ia[keep], ib[keep]), np.maximum(ia[keep], ib[keep])
        idx = (ia, ib)
        st_b = _strata(sep, idx)
        err_b = pair_error(model_med, truth_med, st_b, idx)
        bp.append(stratified_spearman(score_mat[idx], err_b, st_b))
        bb.append(PR.spearman(bs[loci], -be[loci]))
    q = lambda v: [round(float(x), 4) for x in np.nanpercentile(v, [2.5, 97.5])]  # noqa: E731
    return {"pair_stratified_spearman": pair, "pair_ci95": q(bp), "bead_spearman": bead, "bead_ci95": q(bb),
            "pairs": int(np.isfinite(err).sum()), "beads": int(n)}


def run_dataset(key: str, input_kind: str, candidates=CANDIDATES, split: int = 0) -> dict:
    tr = D.load(key)
    a, b = PR.split(tr.n_copies, split)
    n = tr.n_loci
    sep = PR.separation(n, tr.starts)
    hb = PR.half_stats(tr.xyz[b], None if D.REGISTRY[key].kind != "bintu_csv" else 150.0)
    if input_kind == "imaging":
        r_c = 150.0 if D.REGISTRY[key].kind == "bintu_csv" else None
        ha = PR.half_stats(tr.xyz[a], r_c)
        r_c = ha.adjacent_median if r_c is None else r_c
        freq, n_obs, counts = ha.freq, ha.seen, ha.freq * ha.seen
    else:
        from benchmark.run import _hic_input
        from chronocell import physics
        from frozen import HIC_CALIBRATION as HC
        freq, n_obs, counts = _hic_input(_hic_counts(key, tr.starts))
        e = D.REGISTRY[key]
        r_c = physics.bond_length_for(e.step_bp if e.kind != "bintu_csv" else 30_000) / float(
            E.gaussian_median_distance(float(HC["p_adjacent"]), 1.0))
    res = _fit(key, freq, n_obs, r_c, n, split)
    sc = scores(res, freq, counts, n_obs)
    model_med = np.asarray(res.median_distance_nm, dtype=np.float64)
    out = {"dataset": key, "input": input_kind, "loci": n, "split": split, "model": res.config.get("model", "ensemble_v3_3")}
    for name in candidates:
        out[name] = evaluate(sc[name], model_med, hb.median, sep)
        r = out[name]
        print(f"{key:26s} {input_kind:7s} {name:9s} pair {r['pair_stratified_spearman']:+.3f} {r['pair_ci95']}  "
              f"bead {r['bead_spearman']:+.3f} {r['bead_ci95']}", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    from frozen import HIC_CALIBRATION as HC
    if a.practice:
        plan = [(k, "imaging") for k in IMG_PRACTICE] + [(k, "hic") for k in HC["practice"]]
        assert all(D.REGISTRY[k].role == "practice" for k, _ in plan)
        rows = [run_dataset(k, inp) for k, inp in plan]
        (ROOT / "results_reliability_practice.json").write_text(json.dumps(
            {"mode": "practice (candidates compared; in-sample)", "rows": rows,
             "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}, indent=1, default=float))
        return
    from frozen import RELIABILITY as R
    plan = [(k, "imaging") for k in IMG_TEST] + [(k, "hic") for k in HC["test"]]
    assert all(D.REGISTRY[k].role == "test" for k, _ in plan)
    rows = [run_dataset(k, inp, (R["score"],)) for k, inp in plan]
    verdict = {}
    for inp in ("imaging", "hic"):
        for level, key, ci in (("pair", "pair_stratified_spearman", "pair_ci95"), ("bead", "bead_spearman", "bead_ci95")):
            rs = [r[R["score"]] for r in rows if r["input"] == inp]
            ok = [x[key] >= R["min_spearman"] and x[ci][0] > 0 for x in rs]
            verdict[f"{inp}_{level}"] = {"passing": int(sum(ok)), "datasets": len(ok),
                                         "verdict": "pass" if all(ok) else "fail"}
    (ROOT / "results_reliability.json").write_text(json.dumps(
        {"mode": "test (run once; score and rule frozen in validation/frozen.py RELIABILITY)", "rule": R, "rows": rows,
         "verdict": verdict, "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")},
        indent=1, default=float))
    for k, v in verdict.items():
        print(f"{k:14s} {v['verdict']} ({v['passing']} of {v['datasets']})")


if __name__ == "__main__":
    main()
