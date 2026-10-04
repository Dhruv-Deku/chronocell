"""
Gate 1 (Pillar 1): does the whole-chromosome population model at least match the windowed v3.3
model on the same regions?

Protocol (validation/protocol.py), per random split of the imaged chromosome copies:
  input  = half A's contact frequencies at radius r_c (r_c = half A's median adjacent-locus
           distance: the imaging analogue of the app's "adjacent contact probability 0.5" anchor;
           rule fixed before any run, see validation/TUNING.md section 5);
  truth  = half B's median distances.
Models
  windowed v3.3   chronocell.ensemble.fit_ensemble on the fewest equal tiles of <= 400 loci (the
                  app's maximum window); it cannot predict pairs that cross a tile edge.
  whole v4        chronocell.population.fit_population on the whole chromosome at once.
  no 3D           each frequency inverted on its own (Maxwell).
  genomic only    power law in genomic separation, fitted on half A.
  ceiling         half A's own medians.
Scores are reported on (i) within-tile pairs, where both models predict ("the same regions"), and
(ii) all pairs and cross-tile pairs, which only the whole-chromosome model can predict.

A second input, sequencing Hi-C (Rao et al. 2014, IMR-90, binned on the imaged loci by Su et al.),
gives the direct Hi-C -> imaging test with the same truth.

    python validation/gate1.py --practice          # Su chr2 (+ p-arm replicate): tuning allowed
    python validation/gate1.py --test              # Su chr21 (+ replicate): run once, settings frozen
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
from chronocell import ensemble as E, population as P   # noqa: E402

SPLITS = 3


def _windowed(freq: np.ndarray, seen: np.ndarray, r_c: float, tiles: list[tuple[int, int]], seed: int,
              iterations: int | None = None) -> tuple[np.ndarray, float]:
    n = len(freq)
    out = np.full((n, n), np.nan)
    t0 = time.time()
    for a, b in tiles:
        cfg = E.EnsembleConfig(seed=seed, **({"iterations": iterations} if iterations else {}))
        res = E.fit_ensemble(freq[a:b, a:b], seen[a:b, a:b], r_c_nm=r_c, cfg=cfg)
        out[a:b, a:b] = res.median_distance_nm
    return out, time.time() - t0


def _counts_freq(counts: np.ndarray, p_adjacent: float, zeros: str = "clip") -> tuple[np.ndarray, float]:
    """Hi-C counts -> contact probabilities (the app's adjacent-pair anchor). zeros = "unobserved"
    treats zero-count pairs as missing (NaN) instead of "rarer than one in N_eff"."""
    c = np.asarray(counts, dtype=np.float64).copy()
    np.fill_diagonal(c, 0.0)
    p = E.counts_to_probability(c, p_adjacent)
    adj = np.diag(c, 1)
    n_eff = float(np.median(adj[adj > 0])) / p_adjacent
    if zeros == "unobserved":
        p = np.where(c > 0, p, np.nan)
    return p, n_eff


def score_block(pred: dict[str, np.ndarray], truth: np.ndarray, sep: np.ndarray, masks: dict[str, np.ndarray],
                ceiling: np.ndarray) -> dict:
    out = {}
    for mname, m in masks.items():
        if not m.any():
            continue
        ceil = PR.scores(ceiling, truth, sep, m)
        row = {"ceiling": ceil}
        for name, p in pred.items():
            if np.isfinite(p[m]).mean() < 0.99:          # a model that cannot predict these pairs is skipped
                continue
            s = PR.scores(p, truth, sep, m)
            s["percent_of_ceiling"] = 100 * s["spearman_distance_corrected"] / ceil["spearman_distance_corrected"]
            row[name] = s
        out[mname] = row
    return out


def run_dataset(key: str, settings: dict, splits: int, log=print) -> dict:
    tr = D.load(key)
    entry = D.REGISTRY[key]
    n = tr.n_loci
    sep = PR.separation(n, tr.starts)
    tiles = PR.tiles(n, 400)
    within = np.zeros((n, n), bool)
    for a, b in tiles:
        within[a:b, a:b] = True
    upper = np.triu(np.ones((n, n), bool), 1)
    masks = {"within_tiles": upper & within, "cross_tiles": upper & ~within, "all_pairs": upper}
    hic = None
    if entry.paired_hic:
        H, _, hst = D.load_su_hic(D.REGISTRY[entry.paired_hic])
        if not np.array_equal(hst, tr.starts):
            raise ValueError("Hi-C loci do not match the imaged loci.")
        hic = H
    runs = []
    for split in range(splits):
        a_idx, b_idx = PR.split(tr.n_copies, split)
        ha = PR.half_stats(tr.xyz[a_idx], None)
        r_c = ha.adjacent_median                       # pre-registered rule (TUNING.md section 5)
        if settings.get("radius_nm"):
            r_c = float(settings["radius_nm"])
            ha = PR.half_stats(tr.xyz[a_idx], r_c)
        hb = PR.half_stats(tr.xyz[b_idx], r_c)
        truth = hb.median
        pred: dict[str, np.ndarray] = {}
        timing: dict[str, float] = {}
        pred["windowed_v3_3"], timing["windowed_v3_3"] = _windowed(ha.freq, ha.seen, r_c, tiles, split)
        cfg = P.PopulationConfig(seed=split, **settings.get("population", {}))
        t0 = time.time()
        whole = P.fit_population(ha.freq, ha.seen, r_c_nm=r_c, cfg=cfg)
        timing["whole_v4"] = time.time() - t0
        pred["whole_v4"] = whole.median_distance_nm.astype(np.float64)
        fa = np.clip(ha.freq, 0.5 / np.maximum(ha.seen, 1), 1 - 0.5 / np.maximum(ha.seen, 1))
        pred["no_3d_direct_inversion"] = E.gaussian_median_distance(fa, r_c)
        pred["genomic_distance_only"] = PR.genomic_baseline(ha.median, sep, upper)
        run = {"split": split, "r_c_nm": r_c, "copies": [len(a_idx), len(b_idx)], "seconds": timing,
               "whole_v4_contact_fit": whole.contact_fit, "whole_v4_spread_cv": whole.spread_cv,
               "imaging_input": score_block(pred, truth, sep, masks, ha.median)}
        if hic is not None:                            # direct Hi-C -> imaging test (same truth)
            hp: dict[str, np.ndarray] = {}
            p_adj = float(settings.get("p_adjacent", 0.5))
            for anchor_name, factor in settings.get("hic_anchors", {"literature_b0": 1.0}).items():
                b0_nm = _literature_b0(entry) * float(factor)
                f_h, n_eff = _counts_freq(hic, p_adj, settings.get("hic_zeros", "clip"))
                r_h = b0_nm / float(E.gaussian_median_distance(p_adj, 1.0))
                hp_w, _ = _windowed(f_h, np.full_like(f_h, n_eff), r_h, tiles, split)
                hw = P.fit_population(f_h, n_eff, r_c_nm=r_h, cfg=cfg)
                hp[f"windowed_v3_3[{anchor_name}]"] = hp_w
                hp[f"whole_v4[{anchor_name}]"] = hw.median_distance_nm.astype(np.float64)
            run["hic_input"] = score_block(hp, truth, sep, masks, ha.median)
            run["hic_p_adjacent"] = p_adj
        runs.append(run)
        log(f"{key} split {split}: " + _line(run))
    return {"dataset": key, "role": entry.role, "loci": n, "copies": tr.n_copies, "tiles": tiles, "runs": runs,
            "summary": _summary(runs)}


def _literature_b0(entry: D.Entry) -> float:
    from chronocell import physics
    return physics.bond_length_for(entry.step_bp)


def _line(run: dict) -> str:
    w = run["imaging_input"]["within_tiles"]
    parts = [f"{k} {w[k]['percent_of_ceiling']:.1f}%" for k in ("windowed_v3_3", "whole_v4", "no_3d_direct_inversion")
             if k in w]
    if "all_pairs" in run["imaging_input"] and "whole_v4" in run["imaging_input"]["all_pairs"]:
        parts.append(f"whole all-pairs {run['imaging_input']['all_pairs']['whole_v4']['percent_of_ceiling']:.1f}%")
    return "within tiles: " + ", ".join(parts) + f" | whole {run['seconds']['whole_v4']:.0f}s"


def _summary(runs: list[dict]) -> dict:
    out: dict = {}
    for inp in ("imaging_input", "hic_input"):
        if inp not in runs[0]:
            continue
        for mname in runs[0][inp]:
            for model in runs[0][inp][mname]:
                if model == "ceiling":
                    continue
                vals = [r[inp][mname][model] for r in runs if model in r[inp][mname]]
                ceil = [r[inp][mname]["ceiling"]["spearman_distance_corrected"] for r in runs]
                m = float(np.mean([v["spearman_distance_corrected"] for v in vals]))
                out.setdefault(inp, {}).setdefault(mname, {})[model] = {
                    "trend_removed_rho": m, "ceiling": float(np.mean(ceil)), "percent_of_ceiling": 100 * m / float(np.mean(ceil)),
                    "raw_spearman": float(np.mean([v["spearman"] for v in vals])),
                    "lin_ccc_nm": float(np.mean([v["lin_ccc_nm"] for v in vals])),
                    "median_scale": float(np.mean([v["median_scale_model_over_real"] for v in vals])),
                    "splits": len(vals), "pairs": vals[0]["pairs"]}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    ap.add_argument("--splits", type=int, default=SPLITS)
    ap.add_argument("--datasets", nargs="*")
    ap.add_argument("--settings", default=None, help="JSON settings (practice only); test always uses validation/frozen.py")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.test:
        from frozen import GATE1
        settings = GATE1
        keys = a.datasets or ["su_chr21", "su_chr21_rep"]
        if any(D.REGISTRY[k].role != "test" for k in keys):
            sys.exit("--test runs test datasets only")
        out = Path(a.out or ROOT / "results_gate1.json")
    else:
        settings = json.loads(a.settings) if a.settings else {}
        keys = a.datasets or ["su_chr2"]
        if any(D.REGISTRY[k].role != "practice" for k in keys):
            sys.exit("--practice runs practice datasets only")
        out = Path(a.out or ROOT / "results_gate1_practice.json")
    t0 = time.time()
    res = {"settings": settings, "mode": "test" if a.test else "practice",
           "datasets": [run_dataset(k, settings, a.splits) for k in keys]}
    res["seconds"] = time.time() - t0
    out.write_text(json.dumps(res, indent=1, default=float))
    print(f"-> {out} ({res['seconds']:.0f} s)")


if __name__ == "__main__":
    main()
