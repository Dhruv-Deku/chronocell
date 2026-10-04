"""
Practice-only tuning for v4 (Pillars 1, 2, 4): every experiment here uses PRACTICE datasets only, and
every number it prints is written to validation/tuning_v4.json, then summarised in
validation/TUNING.md (abandoned settings included). Test datasets are refused.

    python validation/tune_practice.py hic        # Hi-C -> imaging settings, Su chr2
    python validation/tune_practice.py whole      # whole-chromosome fit settings, Su chr2
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D            # noqa: E402
import protocol as PR           # noqa: E402
from chronocell import ensemble as E, normalize as NZ, population as P   # noqa: E402

CACHE = ROOT / "data" / "cache"
LOG = ROOT / "tuning_v4.json"


def _practice(key: str):
    if D.REGISTRY[key].role != "practice":
        raise SystemExit(f"{key} is not a practice dataset; refusing.")
    return D.load(key)


def truth_cache(key: str, split: int, r_c_rule: str = "adjacent") -> dict:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{key}_split{split}_{r_c_rule}.npz"
    if path.exists():
        z = np.load(path)
        return {k: z[k] for k in z.files}
    tr = _practice(key)
    a, b = PR.split(tr.n_copies, split)
    ha = PR.half_stats(tr.xyz[a], None)
    hb = PR.half_stats(tr.xyz[b], ha.adjacent_median)
    out = {"freq_a": ha.freq, "seen_a": ha.seen, "median_a": ha.median, "median_b": hb.median,
           "r_c": np.array(ha.adjacent_median), "starts": tr.starts}
    np.savez(path, **out)
    return out


def _log(entry: dict) -> None:
    try:
        log = json.loads(LOG.read_text())
    except (OSError, ValueError):
        log = []
    log.append(entry)
    LOG.write_text(json.dumps(log, indent=1, default=float))


def _masks(n: int):
    tiles = PR.tiles(n, 400)
    within = np.zeros((n, n), bool)
    for a, b in tiles:
        within[a:b, a:b] = True
    upper = np.triu(np.ones((n, n), bool), 1)
    return tiles, {"within_tiles": upper & within, "cross_tiles": upper & ~within, "all_pairs": upper}


def _score(pred, t, sep, masks):
    out = {}
    for name, m in masks.items():
        if np.isfinite(pred[m]).mean() < 0.99:
            continue
        c = PR.scores(t["median_a"], t["median_b"], sep, m)["spearman_distance_corrected"]
        s = PR.scores(pred, t["median_b"], sep, m)
        out[name] = {"pct": 100 * s["spearman_distance_corrected"] / c, "raw": s["spearman"], "ccc": s["lin_ccc_nm"],
                     "scale": s["median_scale_model_over_real"]}
    return out


def hic_variants(key: str = "su_chr2", split: int = 0) -> None:
    t = truth_cache(key, split)
    H, _, hst = D.load_su_hic(D.REGISTRY[D.REGISTRY[key].paired_hic])
    n = len(H)
    sep = PR.separation(n, t["starts"])
    tiles, masks = _masks(n)
    iu = np.triu_indices(n, 1)
    counts = H.copy()
    np.fill_diagonal(counts, 0.0)
    ci, cj = iu
    cm = counts[iu]
    keep = cm > 0
    bci, bcj, bcm, notes = NZ.balanced_contacts(ci[keep], cj[keep], cm[keep], n)
    bal = np.zeros((n, n))
    bal[bci, bcj] = bcm
    bal = bal + bal.T
    masked = np.ones(n, bool)
    masked[np.unique(np.r_[bci, bcj])] = False
    variants = []
    for norm in ("raw", "ice"):
        for zeros in ("clip", "unobserved"):
            for p_adj in (0.3, 0.5, 0.7, 0.9):
                variants.append((norm, zeros, p_adj))
    b0 = float(t["r_c"]) * E.gaussian_median_distance(0.5, 1.0) / 1.0     # imaging-scale anchor, only affects nm
    for norm, zeros, p_adj in variants:
        m = (bal if norm == "ice" else counts).copy()
        p = E.counts_to_probability(m, p_adj)
        adj = np.diag(m, 1)
        n_eff = float(np.median(adj[adj > 0])) / p_adj
        if zeros == "unobserved":
            p = np.where(m > 0, p, np.nan)
        if norm == "ice":
            p[masked, :] = np.nan
            p[:, masked] = np.nan
        r_c = b0 / float(E.gaussian_median_distance(p_adj, 1.0))
        t0 = time.time()
        pred = {}
        for name, fn in (("windowed_v3_3", "w"), ("whole_v4", "v")):
            if fn == "w":
                out = np.full((n, n), np.nan)
                for a, b in tiles:
                    r = E.fit_ensemble(p[a:b, a:b], n_eff, r_c_nm=r_c, cfg=E.EnsembleConfig(seed=split))
                    out[a:b, a:b] = r.median_distance_nm
                pred[name] = out
            else:
                r = P.fit_population(p, n_eff, r_c_nm=r_c, cfg=P.PopulationConfig(seed=split, replicas=20, frames=2))
                pred[name] = r.median_distance_nm.astype(np.float64)
        res = {k: _score(v, t, sep, masks) for k, v in pred.items()}
        entry = {"experiment": "hic_to_imaging", "dataset": key, "split": split, "normalisation": norm,
                 "zero_counts": zeros, "p_adjacent": p_adj, "ice_note": notes[0] if norm == "ice" else "",
                 "seconds": time.time() - t0, "scores": res}
        _log(entry)
        w = res["windowed_v3_3"]["within_tiles"]["pct"]
        v = res["whole_v4"]
        print(f"{norm:4s} zeros={zeros:10s} p_adj={p_adj:.1f}: within windowed {w:5.1f}% whole {v['within_tiles']['pct']:5.1f}% "
              f"| cross whole {v['cross_tiles']['pct']:6.1f}% | all whole {v['all_pairs']['pct']:5.1f}% raw {v['all_pairs']['raw']:.3f}",
              flush=True)


def whole_variants(key: str = "su_chr2", split: int = 0) -> None:
    t = truth_cache(key, split)
    n = len(t["freq_a"])
    sep = PR.separation(n, t["starts"])
    tiles, masks = _masks(n)
    r_c = float(t["r_c"])
    grid = [
        ("default full rank fp32", {}),
        ("full rank fp64", {"dtype": "float64"}),
        ("3000 iterations", {"iterations": 3000}),
        ("lr 0.02", {"learning_rate": 0.02}),
        ("rank 256 + random walk", {"rank_cap": 256}),
        ("rank 128 + random walk", {"rank_cap": 128}),
        ("coarse start (k=3)", {"init": "coarse", "coarse_beads": 312}),
        ("uniform weights", {"weighting": "uniform"}),
    ]
    for name, kw in grid:
        t0 = time.time()
        r = P.fit_population(t["freq_a"], t["seen_a"], r_c_nm=r_c,
                             cfg=P.PopulationConfig(seed=split, replicas=20, frames=2, **kw))
        res = _score(r.median_distance_nm.astype(np.float64), t, sep, masks)
        entry = {"experiment": "whole_chromosome_fit", "dataset": key, "split": split, "variant": name, "settings": kw,
                 "seconds": time.time() - t0, "fit_seconds": r.config["fit_seconds"],
                 "best_misfit": r.history["best_loss"][0], "contact_fit": r.contact_fit, "scores": res}
        _log(entry)
        print(f"{name:26s} within {res['within_tiles']['pct']:5.1f}% cross {res['cross_tiles']['pct']:5.1f}% "
              f"all {res['all_pairs']['pct']:5.1f}% ccc {res['all_pairs']['ccc']:.3f} misfit {r.history['best_loss'][0]:.4f} "
              f"fit {r.config['fit_seconds']:.0f}s", flush=True)


def reliability_variants(keys=("bintu_k562_28_30", "bintu_hct116_28_30", "bintu_hct116_28_30_auxin",
                                "bintu_hct116_34_37", "su_chr2", "su_chr2_parm_rep"), split: int = 0) -> None:
    """Which per-bead reliability score, computed from the input alone, ranks held-out per-bead error?
    Candidates: (a) misfit-based (population.bead_reliability); (b) input sampling noise: the mean binomial
    standard error of log f over the bead's pairs, sqrt((1 - f) / (N f)); (c) both multiplied."""
    for key in keys:
        tr = _practice(key)
        a, b = PR.split(tr.n_copies, split)
        r_c = 150.0 if D.REGISTRY[key].kind == "bintu_csv" else None
        ha = PR.half_stats(tr.xyz[a], r_c)
        r_c = ha.adjacent_median if r_c is None else r_c
        hb = PR.half_stats(tr.xyz[b], r_c)
        n = tr.n_loci
        if n <= 400:
            res = E.fit_ensemble(ha.freq, ha.seen, r_c_nm=r_c, cfg=E.EnsembleConfig(seed=split))
        else:
            from frozen import WHOLE_CHROMOSOME
            res = P.fit_population(ha.freq, ha.seen, r_c_nm=r_c, cfg=P.PopulationConfig(seed=split, replicas=20, frames=2,
                                                                                            **WHOLE_CHROMOSOME))
        err = np.abs(np.log(np.where(res.median_distance_nm > 0, res.median_distance_nm, np.nan)) -
                     np.log(np.where(hb.median > 0, hb.median, np.nan)))
        np.fill_diagonal(err, np.nan)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            bead_err = np.nanmedian(err, axis=1)
            f = np.clip(ha.freq, 0.5 / np.maximum(ha.seen, 1), 1 - 0.5 / np.maximum(ha.seen, 1))
            se = np.sqrt((1 - f) / (np.maximum(ha.seen, 1) * f))
            np.fill_diagonal(se, np.nan)
            noise = np.nanmean(se, axis=1)
        misfit_rel = P.bead_reliability(res, ha.freq, ha.seen)
        noise_rel = 1.0 / (1.0 + noise)
        scores = {"a_misfit": PR.spearman(misfit_rel, -bead_err), "b_sampling_noise": PR.spearman(noise_rel, -bead_err),
                  "c_product": PR.spearman(misfit_rel * noise_rel, -bead_err)}
        _log({"experiment": "bead_reliability", "dataset": key, "split": split, "spearman_vs_neg_error": scores})
        print(key, {k: round(v, 3) for k, v in scores.items()}, flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["hic", "whole", "reliability"])
    ap.add_argument("--dataset", default="su_chr2")
    ap.add_argument("--split", type=int, default=0)
    a = ap.parse_args()
    if a.what == "reliability":
        reliability_variants(split=a.split)
        return
    {"hic": hic_variants, "whole": whole_variants}[a.what](a.dataset, a.split)


if __name__ == "__main__":
    main()
