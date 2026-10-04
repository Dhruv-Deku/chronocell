"""
Gate 4, structural variants: does the SV simulator predict the contacts of a real rearranged genome
better than a genomic-distance shift?

    python validation/sv_validation.py                   # -> validation/results_sv.json

Everything below is fixed in validation/frozen.py (SV_VALIDATION), committed before this script read
any K562 counts on chromosome 9:
  event      the two chr9 regions lost entirely in K562 (Zhou et al., Genome Res 2019; hg19)
  reference  GM12878, variant K562: Rao et al. 2014 in situ Hi-C (GEO GSE63525), raw counts, 25 kb
  model      v4 population fitted on the GM12878 window (Gate 1 Hi-C settings); the deletions are then
             applied by chronocell.perturb (deletion_pieces + derive). Nothing is fitted on K562.
  pairs      "spanning pairs": beads on different kept runs whose separation in the rearranged chain
             is 2-80 bins (50 kb - 2 Mb)
  metrics    Spearman of each prediction vs K562 counts; trend-removed Spearman for the model (each map
             divided by its own mean over within-run pairs at the same separation); 95 % intervals by a
             block bootstrap over 8 x 8-bin blocks of the pair grid
  baselines  distance shift (GM12878's mean count at the new separation) and no change (GM12878's
             count for the pair at its original separation)
  pass       (i) model Spearman > distance shift with the 95 % interval of the difference above 0, and
             (ii) the model's trend-removed Spearman interval above 0
A check that is not part of the rule: the published deleted bins should have ~no K562 reads.
"""

from __future__ import annotations

import datetime as dt
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                       # noqa: E402
import protocol as PR                      # noqa: E402
from frozen import SV_VALIDATION as CFG    # noqa: E402
from chronocell import ensemble as E, perturb as PT, physics, population as P   # noqa: E402

OUT = ROOT / "results_sv.json"


def deleted_bins(n: int, start: int, binsize: int, intervals) -> np.ndarray:
    """Bins whose midpoint lies inside a deleted interval."""
    mids = start + (np.arange(n) + 0.5) * binsize
    drop = np.zeros(n, bool)
    for lo, hi in intervals:
        drop |= (mids >= lo) & (mids < hi)
    return drop


def segments(drop: np.ndarray) -> list[tuple[int, int]]:
    edges = np.flatnonzero(np.diff(np.r_[0, drop.astype(np.int8), 0]))
    return [(int(a), int(b)) for a, b in zip(edges[::2], edges[1::2])]


def fit_reference(G: np.ndarray) -> E.EnsembleResult:
    """The v4 population on the reference counts, exactly as the Gate 1 / benchmark Hi-C input."""
    c = np.asarray(G, dtype=float).copy()
    np.fill_diagonal(c, 0.0)
    p_adj = float(CFG["p_adjacent"])
    p = E.counts_to_probability(c, p_adj)
    adj = np.diag(c, 1)
    n_eff = float(np.median(adj[adj > 0])) / p_adj
    if CFG["hic_zeros"] == "unobserved":
        p = np.where(c > 0, p, np.nan)
    b0 = physics.bond_length_for(CFG["binsize"])
    r_c = b0 / float(E.gaussian_median_distance(p_adj, 1.0))
    cfg = P.PopulationConfig(seed=CFG["seed"], replicas=20, frames=2, **CFG["population"])
    return P.fit_population(p, n_eff, r_c_nm=r_c, cfg=cfg)


def mean_by_separation(M: np.ndarray, sep: np.ndarray, mask: np.ndarray, smax: int) -> np.ndarray:
    s = sep[mask]
    v = M[mask]
    ok = np.isfinite(v) & (s <= smax)
    tot = np.bincount(s[ok], weights=v[ok], minlength=smax + 1)
    cnt = np.bincount(s[ok], minlength=smax + 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)


def block_bootstrap(blocks: np.ndarray, fns: dict, reps: int, seed: int) -> dict:
    """95 % intervals of each statistic fns[name](index array), resampling whole blocks of pairs."""
    rng = np.random.default_rng(seed)
    ids, inv = np.unique(blocks, return_inverse=True)
    members = [np.flatnonzero(inv == k) for k in range(len(ids))]
    draws = {k: [] for k in fns}
    for _ in range(reps):
        pick = rng.integers(0, len(ids), len(ids))
        idx = np.concatenate([members[k] for k in pick])
        for k, f in fns.items():
            draws[k].append(f(idx))
    return {k: [float(x) for x in np.nanpercentile(v, [2.5, 97.5])] for k, v in draws.items()}


def run() -> dict:
    t0 = time.time()
    chrom, start, end, binsize = CFG["chrom"], CFG["start"], CFG["end"], CFG["binsize"]
    G = D.hic_region(CFG["reference"], chrom, start, end, binsize)
    K = D.hic_region(CFG["variant"], chrom, start, end, binsize)
    n = len(G)
    assert K.shape == G.shape
    drop = deleted_bins(n, start, binsize, CFG["deleted"])
    segs = segments(drop)

    # check (not part of the rule): reads on the published deleted bins
    near = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 1
    cov_g = (G * near).sum(1)
    cov_k = (K * near).sum(1)
    check = {"k562_deleted_over_kept_coverage": float(np.median(cov_k[drop]) / np.median(cov_k[~drop])),
             "gm12878_deleted_over_kept_coverage": float(np.median(cov_g[drop]) / np.median(cov_g[~drop])),
             "bins_without_gm12878_reads": int((cov_g == 0).sum()),
             "deleted_bins": int(drop.sum()), "kept_bins": int((~drop).sum())}

    res = fit_reference(G)
    S, r_c = PT.pair_variance_from_result(res)
    imp = PT.variant_impact_from_variance(S, r_c, "deletion", {"segments": segs})
    pieces, _ = PT.deletion_pieces(n, segs)
    run_of = np.full(n, -1)
    rank = np.full(n, -1)
    k = 0
    for r, p in enumerate(pieces):
        run_of[p.beads] = r
        rank[p.beads] = np.arange(k, k + len(p.beads))
        k += len(p.beads)

    kept = ~drop
    upper = np.triu(np.ones((n, n), bool), 1)
    both = upper & kept[:, None] & kept[None, :]
    same_run = run_of[:, None] == run_of[None, :]
    new_sep = np.abs(np.subtract.outer(rank, rank))
    old_sep = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    smin, smax = int(CFG["min_new_separation_bins"]), int(CFG["max_new_separation_bins"])
    spanning = both & ~same_run & (new_sep >= smin) & (new_sep <= smax)
    within = both & same_run

    # predictions of K562 counts on spanning pairs
    exp_g = mean_by_separation(G, old_sep, upper, n)                    # GM12878 expected, all pairs
    pred = {"model": imp.p_after, "distance_shift": exp_g[np.minimum(new_sep, n)], "no_change": G}
    i, j = np.nonzero(spanning)
    obs = K[i, j]
    vals = {name: m[i, j] for name, m in pred.items()}
    # trend-removed: each map over its own within-run mean at the same separation
    e_model = mean_by_separation(imp.p_after, new_sep, within, smax)
    e_k = mean_by_separation(K, new_sep, within, smax)
    s_ij = new_sep[i, j]
    tr_model = vals["model"] / e_model[s_ij]
    tr_obs = obs / e_k[s_ij]

    def rho(name):
        return lambda idx: PR.spearman(vals[name][idx], obs[idx])

    stats = {"model": rho("model"), "distance_shift": rho("distance_shift"), "no_change": rho("no_change"),
             "model_minus_distance_shift": lambda idx: rho("model")(idx) - rho("distance_shift")(idx),
             "model_trend_removed": lambda idx: PR.spearman(tr_model[idx], tr_obs[idx])}
    every = np.arange(len(i))
    point = {k: float(f(every)) for k, f in stats.items()}
    B = int(CFG["block_bins"])
    blocks = (i // B) * (n // B + 1) + (j // B)
    ci = block_bootstrap(blocks, stats, int(CFG["bootstrap"]), int(CFG["seed"]))

    # all spanning pairs at any separation (secondary, no rule attached)
    span_all = both & ~same_run & (new_sep >= smin)
    ia, ja = np.nonzero(span_all)
    secondary = {name: PR.spearman(m[ia, ja], K[ia, ja]) for name, m in
                 {"model": imp.p_after, "distance_shift": exp_g[np.minimum(new_sep, n)], "no_change": G}.items()}

    crit_i = point["model_minus_distance_shift"] > 0 and ci["model_minus_distance_shift"][0] > 0
    crit_ii = ci["model_trend_removed"][0] > 0
    by_junction = {}
    for a in range(len(pieces)):
        for b in range(a + 1, len(pieces)):
            m = spanning & (run_of[:, None] == a) & (run_of[None, :] == b)
            if m.sum() >= 20:
                ii, jj = np.nonzero(m)
                by_junction[f"run{a}-run{b}"] = {
                    "pairs": int(m.sum()),
                    **{name: PR.spearman(mat[ii, jj], K[ii, jj]) for name, mat in
                       {"model": imp.p_after, "distance_shift": exp_g[np.minimum(new_sep, n)], "no_change": G}.items()}}
    url_g, asm, _ = D.HIC_SOURCES[CFG["reference"]]
    url_k = D.HIC_SOURCES[CFG["variant"]][0]
    out = {
        "settings": {k: v for k, v in CFG.items() if k != "population"} | {"population": dict(CFG["population"])},
        "data": {"reference_url": url_g, "variant_url": url_k, "assembly": asm, "licence": D.HIC_LICENSE,
                 "citation": "Rao SSP et al. (2014) Cell 159:1665-1680; deletions: Zhou B et al. (2019) Genome Res "
                             "29:472-484"},
        "window": f"{chrom}:{start:,}-{end:,} ({n} bins of {binsize // 1000} kb)",
        "deleted_runs_bins": segs, "kept_runs": [[int(p.beads[0]), int(p.beads[-1])] for p in pieces],
        "check_deleted_bins": check,
        "reference_fit": {"contact_fit_spearman": float(res.contact_fit), "best_misfit": float(res.history["best_loss"][0]),
                          "fit_seconds": float(res.config["fit_seconds"])},
        "spanning_pairs": int(len(i)),
        "spearman": point, "ci95_block_bootstrap": ci,
        "secondary_all_spanning_pairs": {"pairs": int(len(ia)), **secondary},
        "by_junction": by_junction,
        "criteria": {"i_beats_distance_shift": bool(crit_i), "ii_pattern_beyond_separation": bool(crit_ii)},
        "verdict": "validated on this event" if (crit_i and crit_ii) else "not validated (mechanism simulator)",
        "seconds": round(time.time() - t0, 1),
        "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    }
    return out


def main() -> None:
    out = run()
    OUT.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("check_deleted_bins", "spanning_pairs", "spearman", "ci95_block_bootstrap",
                                          "secondary_all_spanning_pairs", "by_junction", "criteria", "verdict")},
                     indent=1, default=float))
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
