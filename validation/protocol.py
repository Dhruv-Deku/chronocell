"""
The held-out protocol and the metrics, shared by every validation script, so all methods are scored
the same way. It is the v3.3 protocol (validation/validate_tracing.py), generalised to whole
chromosomes:

  1. Split the imaged chromosome copies at random into halves A and B (seeded).
  2. Input  = half A's contact frequencies only (pairs closer than the contact radius).
  3. Truth  = half B's measured median distance per pair; the model never sees half B.
  4. Scores (microscopy accuracy, never mixed with the contact-map fit):
     - trend-removed Spearman rho: each distance divided by the mean at its genomic separation,
       then ranked (the "distance-pattern recovery beyond the genomic-separation trend");
     - the same as a % of the ceiling = how well half A's own medians predict half B;
     - raw Spearman rho, Pearson r, Lin's concordance correlation (absolute size, nm);
     - median model / measured size;
     - calibration: share of half B's single-copy distances inside the model's stated intervals.

Everything is computed in blocks of rows, so a 935-locus chromosome with thousands of copies stays
within a few hundred MB.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

LEVELS = (0.5, 0.8, 0.9)


def split(n_copies: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    idx = np.random.default_rng(seed).permutation(n_copies)
    return idx[: n_copies // 2], idx[n_copies // 2:]


@dataclass
class HalfStats:
    freq: np.ndarray          # (n, n) fraction of observed copies with d < r_c (NaN if never observed)
    seen: np.ndarray          # (n, n) copies with both loci detected
    median: np.ndarray        # (n, n) median distance (nm), NaN if unobserved
    adjacent_median: float    # median distance between consecutive loci (nm)


def half_stats(xyz: np.ndarray, r_c_nm: float | None, rows: int = 24, max_bytes: float = 4e8) -> HalfStats:
    """Contact frequency, observation counts and median distance for every locus pair, in row blocks.
    r_c_nm None -> the contact radius is the median adjacent-locus distance of these copies."""
    x = np.asarray(xyz, dtype=np.float32)
    c, n, _ = x.shape
    rows = max(1, min(rows, int(max_bytes // max(1, c * n * 4))))
    adj = np.linalg.norm(x[:, 1:] - x[:, :-1], axis=-1)
    adj_med = float(np.nanmedian(adj))
    r = adj_med if r_c_nm is None else float(r_c_nm)
    freq = np.full((n, n), np.nan)
    seen = np.zeros((n, n))
    med = np.full((n, n), np.nan)
    for a in range(0, n, rows):
        b = min(n, a + rows)
        d = np.linalg.norm(x[:, a:b, None, :] - x[:, None, :, :], axis=-1)       # (c, block, n)
        ok = np.isfinite(d)
        k = ok.sum(0)
        seen[a:b] = k
        close = (d < r).sum(0)
        freq[a:b] = np.where(k > 0, close / np.maximum(k, 1), np.nan)
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                med[a:b] = np.nanmedian(d, axis=0)
    np.fill_diagonal(med, 0.0)
    return HalfStats(freq, seen, med, adj_med)


def separation(n: int, starts: np.ndarray | None = None) -> np.ndarray:
    """Genomic separation of every locus pair: |i - j| in loci, or |start_i - start_j| in bp."""
    if starts is None:
        return np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    s = np.asarray(starts, dtype=np.int64)
    return np.abs(np.subtract.outer(s, s))


def _rank(v: np.ndarray) -> np.ndarray:
    order = np.argsort(v, kind="mergesort")
    r = np.empty(len(v), dtype=np.float64)
    r[order] = np.arange(len(v))
    # average ties
    sv = v[order]
    edges = np.flatnonzero(np.diff(sv)) + 1
    starts = np.r_[0, edges]
    ends = np.r_[edges, len(v)]
    tie = ends - starts > 1
    for s0, e0 in zip(starts[tie], ends[tie]):
        r[order[s0:e0]] = (s0 + e0 - 1) / 2.0
    return r


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, np.float64), np.asarray(b, np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    ra, rb = _rank(a[ok]), _rank(b[ok])
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def lin_ccc(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, np.float64), np.asarray(b, np.float64)
    cov = np.mean((a - a.mean()) * (b - b.mean()))
    return float(2 * cov / (a.var() + b.var() + (a.mean() - b.mean()) ** 2))


def pair_mask(n: int, sep: np.ndarray, same_chrom: np.ndarray | None = None) -> np.ndarray:
    """Upper-triangle pairs to score (i < j, and on the same chromosome for genome-scale data)."""
    m = np.triu(np.ones((n, n), bool), 1)
    if same_chrom is not None:
        m &= same_chrom
    return m


def scores(pred: np.ndarray, truth: np.ndarray, sep: np.ndarray, mask: np.ndarray | None = None) -> dict:
    """The v3.3 microscopy metrics on the pairs in `mask` (default: all i < j)."""
    n = len(truth)
    mask = np.triu(np.ones((n, n), bool), 1) if mask is None else mask
    p, t, s = pred[mask], truth[mask], sep[mask]
    ok = np.isfinite(p) & np.isfinite(t) & (p > 0) & (t > 0)
    p, t, s = p[ok].astype(np.float64), t[ok].astype(np.float64), s[ok]
    if p.size < 3:
        return {"pairs": int(p.size)}
    _, inv = np.unique(s, return_inverse=True)
    pe = (np.bincount(inv, weights=p) / np.bincount(inv))[inv]      # observed / expected at each separation
    te = (np.bincount(inv, weights=t) / np.bincount(inv))[inv]
    return {"spearman": spearman(p, t), "pearson": float(np.corrcoef(p, t)[0, 1]),
            "spearman_distance_corrected": spearman(p / pe, t / te), "lin_ccc_nm": lin_ccc(p, t),
            "median_scale_model_over_real": float(np.median(p / t)), "pairs": int(p.size)}


def genomic_baseline(median_a: np.ndarray, sep: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Distance from genomic separation alone: a power law fitted on half A's medians."""
    ok = mask & np.isfinite(median_a) & (median_a > 0) & (sep > 0)
    k, logc = np.polyfit(np.log(sep[ok]), np.log(median_a[ok]), 1)
    out = np.exp(logc) * np.maximum(sep, 1).astype(np.float64) ** k
    np.fill_diagonal(out, 0.0)
    return out


# ======================================================================================
# Calibration of the population's predicted spread (Pillar 2)
# ======================================================================================
def coverage(xyz_b: np.ndarray, sigma_nm: np.ndarray, mask: np.ndarray, levels=LEVELS, rows: int = 24,
             recal=None, pit_bins: int = 200, max_bytes: float = 4e8) -> dict:
    """Share of half-B single-copy distances inside the model's central intervals, per stated level.

    sigma_nm (n, n): per-axis sd of each pair vector (Maxwell law). `recal` (optional) maps a stated
    level to the level actually used for the interval (a recalibration fitted on practice data).
    Also returns the histogram of PIT values u = F_model(d), the input to recalibration.
    """
    from chronocell import population as POP
    x = np.asarray(xyz_b, dtype=np.float32)
    c, n, _ = x.shape
    rows = max(1, min(rows, int(max_bytes // max(1, c * n * 4))))
    inside = {lv: 0 for lv in levels}
    total = 0
    pit = np.zeros(pit_bins)
    by_sep: dict[int, list] = {}
    bounds = {lv: POP.central_interval(recal(lv) if recal else lv) for lv in levels}
    for a in range(0, n, rows):
        b = min(n, a + rows)
        d = np.linalg.norm(x[:, a:b, None, :] - x[:, None, :, :], axis=-1)       # (c, block, n)
        m = mask[a:b] & np.isfinite(sigma_nm[a:b]) & (sigma_nm[a:b] > 0)
        sig = np.where(m, sigma_nm[a:b], np.nan)
        z = d / sig[None]                                                          # in units of sigma
        ok = np.isfinite(z)
        total += int(ok.sum())
        for lv, (lo, hi) in bounds.items():
            inside[lv] += int(((z > lo) & (z < hi) & ok).sum())
        u = POP.maxwell_cdf(np.where(ok, z, 0.0))[ok]
        pit += np.histogram(u, bins=pit_bins, range=(0.0, 1.0))[0]
    out = {"levels": list(levels), "coverage": [inside[lv] / max(total, 1) for lv in levels],
           "measurements": total, "pit_histogram": pit.tolist()}
    return out


def isotonic_recalibration(pit_histogram: np.ndarray):
    """Quantile recalibration (Kuleshov, Fenner & Ermon, ICML 2018): with G the empirical CDF of the
    PIT values u = F_model(d) on practice data, the recalibrated CDF is G(F_model(d)). A central
    interval at stated level p is therefore taken between the model quantiles G^-1((1-p)/2) and
    G^-1((1+p)/2). Returned as level -> (lower PIT, upper PIT) through `interval_pits`, and as an
    equivalent symmetric level for display."""
    h = np.asarray(pit_histogram, dtype=np.float64)
    edges = np.linspace(0.0, 1.0, len(h) + 1)
    cdf = np.concatenate([[0.0], np.cumsum(h) / h.sum()])

    def g_inv(q: float) -> float:                     # G^-1, monotone by construction
        return float(np.interp(q, cdf, edges))

    def interval_pits(level: float) -> tuple[float, float]:
        return g_inv((1 - level) / 2), g_inv((1 + level) / 2)

    return interval_pits


def coverage_with_pits(xyz_b: np.ndarray, sigma_nm: np.ndarray, mask: np.ndarray, interval_pits, levels=LEVELS,
                       rows: int = 24, max_bytes: float = 4e8) -> dict:
    """Coverage of recalibrated intervals: [F^-1(lower PIT), F^-1(upper PIT)] x sigma per pair."""
    from chronocell import population as POP
    x = np.asarray(xyz_b, dtype=np.float32)
    c, n, _ = x.shape
    rows = max(1, min(rows, int(max_bytes // max(1, c * n * 4))))
    bounds = {lv: tuple(float(POP.maxwell_quantile(q)) for q in interval_pits(lv)) for lv in levels}
    inside = {lv: 0 for lv in levels}
    total = 0
    for a in range(0, n, rows):
        b = min(n, a + rows)
        d = np.linalg.norm(x[:, a:b, None, :] - x[:, None, :, :], axis=-1)
        m = mask[a:b] & np.isfinite(sigma_nm[a:b]) & (sigma_nm[a:b] > 0)
        z = d / np.where(m, sigma_nm[a:b], np.nan)[None]
        ok = np.isfinite(z)
        total += int(ok.sum())
        for lv, (lo, hi) in bounds.items():
            inside[lv] += int(((z > lo) & (z < hi) & ok).sum())
    return {"levels": list(levels), "coverage": [inside[lv] / max(total, 1) for lv in levels], "measurements": total,
            "bounds_sigma": {str(lv): bounds[lv] for lv in levels}}


def bootstrap_truth_ci(xyz_b: np.ndarray, pred: np.ndarray, sep: np.ndarray, mask: np.ndarray, reps: int = 100,
                       seed: int = 0, max_pairs: int = 20_000, key: str = "spearman_distance_corrected") -> tuple[float, float]:
    """95 % interval of a score from resampling half B's copies (truth noise), on a pair subsample."""
    rng = np.random.default_rng(seed)
    iu = np.argwhere(mask)
    if len(iu) > max_pairs:
        iu = iu[rng.choice(len(iu), max_pairs, replace=False)]
    i, j = iu[:, 0], iu[:, 1]
    x = np.asarray(xyz_b, dtype=np.float32)
    d = np.linalg.norm(x[:, i] - x[:, j], axis=-1)                                # (copies, pairs)
    sub_sep = sep[i, j]
    p = pred[i, j]
    vals = []
    import warnings
    for _ in range(reps):
        pick = rng.integers(0, len(x), len(x))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            t = np.nanmedian(d[pick], axis=0)
        vals.append(_flat_scores(p, t, sub_sep)[key])
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def _flat_scores(p: np.ndarray, t: np.ndarray, s: np.ndarray) -> dict:
    ok = np.isfinite(p) & np.isfinite(t) & (p > 0) & (t > 0)
    p, t, s = p[ok].astype(np.float64), t[ok].astype(np.float64), s[ok]
    _, inv = np.unique(s, return_inverse=True)
    pe = (np.bincount(inv, weights=p) / np.bincount(inv))[inv]
    te = (np.bincount(inv, weights=t) / np.bincount(inv))[inv]
    return {"spearman": spearman(p, t), "spearman_distance_corrected": spearman(p / pe, t / te),
            "lin_ccc_nm": lin_ccc(p, t)}


def tiles(n: int, max_size: int = 400) -> list[tuple[int, int]]:
    """The fewest equal, non-overlapping windows of at most max_size loci covering [0, n)."""
    k = int(math.ceil(n / max_size))
    edges = np.round(np.linspace(0, n, k + 1)).astype(int)
    return [(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:])]
