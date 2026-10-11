"""
Analysis suite for contact maps (Phase B3): loops, domains, compartments, insulation and P(s), at a chosen
resolution, on one chromosome or region. chronocell/domains.py is reused unchanged (its insulation score,
boundaries, compartments and the app's quick loop list keep their current outputs); this module adds:

loops        HiCCUPS-like (Rao et al., Cell 2014): for every pixel, local expected counts from four
             neighbourhoods (donut, lower-left, horizontal, vertical) scaled by the distance-expected map,
             E_K = sum_K(balanced obs) x E_d(i, j) / sum_K(E_d), returned to the raw scale with the balancing
             weights; Poisson test of the raw count; Benjamini-Hochberg within each neighbourhood; Rao's
             enrichment thresholds (1.75 donut and lower-left, 1.5 horizontal and vertical, 2 for one of donut /
             lower-left); enriched pixels within `cluster_bp` are clustered and the most enriched kept;
             single-pixel clusters need summed q <= 0.02. Computed in diagonal tiles (no whole matrix in memory).
topdom       TopDom-like (Shin et al., NAR 2016): mean contact signal across each bin in a w-bin diamond,
             local minima, kept where the diamond is significantly lower (one-sided rank-sum p < 0.05) than the
             contacts inside the two flanking domains.
arrowhead    Arrowhead-like corner score (simplified from Rao 2014): for candidate domains between boundaries,
             mean O/E inside the domain triangle over the mean O/E of the two flanking regions of equal size;
             domains with a ratio >= min_ratio are kept (nesting allowed).
compartments chronocell.domains.compartments (O/E correlation, first eigenvector), signed by GC or gene density.
insulation   chronocell.domains.insulation / boundaries.
p_of_s       contact probability against separation (balanced or raw), with the log-log slope.

Input everywhere: sparse upper-triangle counts (ci < cj, cm) of n bins; `weights` are balancing weights
(chronocell.normalize ICE or KR), NaN for masked bins.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import domains as DOM

LOOP_PARAMS = {5_000: (4, 7), 10_000: (2, 5), 25_000: (1, 3)}      # resolution -> (peak width p, window w), Rao 2014


def _params(resolution: int) -> tuple[int, int]:
    if resolution in LOOP_PARAMS:
        return LOOP_PARAMS[resolution]
    return (2, 5) if resolution < 25_000 else (1, 3)


def kernels(p: int, w: int) -> dict[str, np.ndarray]:
    """HiCCUPS neighbourhoods on a (2w+1)^2 window centred on the pixel (rows = i, columns = j)."""
    size = 2 * w + 1
    c = w
    yy, xx = np.mgrid[0:size, 0:size]
    dy, dx = yy - c, xx - c
    inner = (np.abs(dy) <= p) & (np.abs(dx) <= p)
    donut = ~inner & (dy != 0) & (dx != 0)
    lower_left = (dy > 0) & (dx < 0) & ~inner                       # i larger (further from the diagonal side), j smaller
    horizontal = (np.abs(dy) <= 1) & ~inner
    vertical = (np.abs(dx) <= 1) & ~inner
    return {k: v.astype(float) for k, v in (("donut", donut), ("lower_left", lower_left),
                                            ("horizontal", horizontal), ("vertical", vertical))}


def expected_balanced(ci, cj, cm, n: int, weights: np.ndarray | None, max_d: int) -> np.ndarray:
    """Mean balanced count per diagonal d < max_d over bins valid on both ends."""
    ci, cj, cm = np.asarray(ci, np.int64), np.asarray(cj, np.int64), np.asarray(cm, float)
    wv = np.ones(n) if weights is None else np.asarray(weights, float)
    valid = np.isfinite(wv) & (wv > 0)
    d = cj - ci
    keep = (d < max_d) & valid[ci] & valid[cj]
    bal = cm[keep] * wv[ci[keep]] * wv[cj[keep]]
    sums = np.bincount(d[keep], weights=bal, minlength=max_d)[:max_d]
    pairs = np.array([np.sum(valid[:n - k] & valid[k:]) if k < n else 0 for k in range(max_d)], float)
    return sums / np.maximum(pairs, 1)


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    q = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


@dataclass
class LoopCall:
    i: int
    j: int
    observed: float
    expected_donut: float
    oe_donut: float
    q: dict[str, float]
    cluster_size: int


def call_loops(ci, cj, cm, n: int, resolution: int, weights: np.ndarray | None = None, max_sep_bp: int = 2_000_000,
               min_sep_bp: int = 30_000, fdr: float = 0.1, cluster_bp: int = 20_000, tile: int = 400,
               p: int | None = None, w: int | None = None) -> list[LoopCall]:
    from scipy.ndimage import correlate
    from scipy.stats import poisson
    p0, w0 = _params(resolution)
    p, w = p or p0, w or w0
    D = max(w + 2, int(max_sep_bp // resolution) + 1)
    dmin = max(w + 1, int(np.ceil(min_sep_bp / resolution)))
    ci, cj, cm = np.asarray(ci, np.int64), np.asarray(cj, np.int64), np.asarray(cm, float)
    lo_, hi_ = np.minimum(ci, cj), np.maximum(ci, cj)
    ci, cj = lo_, hi_
    wv = np.ones(n) if weights is None else np.asarray(weights, float)
    valid = np.isfinite(wv) & (wv > 0)
    wz = np.where(valid, wv, 0.0)
    Ed = expected_balanced(ci, cj, cm, n, weights, D + w + 1)
    K = kernels(p, w)
    order = np.argsort(ci, kind="stable")
    ci, cj, cm = ci[order], cj[order], cm[order]
    parts: list[dict] = []
    for s in range(0, n, tile):
        r0, r1 = max(0, s - w), min(n, s + tile + w)
        c0, c1 = r0, min(n, s + tile + D + w)
        a, b = np.searchsorted(ci, r0), np.searchsorted(ci, r1)
        ii, jj, vv = ci[a:b], cj[a:b], cm[a:b]
        sel = (jj >= c0) & (jj < c1)
        raw = np.zeros((r1 - r0, c1 - c0))
        raw[ii[sel] - r0, jj[sel] - c0] = vv[sel]
        # the lower triangle inside the tile mirrors the upper one (kernels near the diagonal read it)
        rr, cc = np.meshgrid(np.arange(r0, r1), np.arange(c0, c1), indexing="ij")
        low = cc < rr
        if low.any():
            raw[low] = raw[cc[low] - r0, rr[low] - c0]
        bal = raw * wz[r0:r1, None] * wz[None, c0:c1]
        dd = np.abs(cc - rr)
        E = Ed[np.minimum(dd, len(Ed) - 1)] * (dd < len(Ed))
        E = E * (valid[r0:r1, None] & valid[None, c0:c1])
        rows = np.arange(s, min(n, s + tile))
        target = np.zeros_like(raw, bool)
        for i in rows:
            jlo, jhi = i + dmin, min(n, i + D)
            if jhi > jlo:
                target[i - r0, jlo - c0:jhi - c0] = True
        target &= (raw > 0) & valid[r0:r1, None] & valid[None, c0:c1]
        if not target.any():
            continue
        ti, tj = np.nonzero(target)
        obs = raw[ti, tj]
        part = {"i": ti + r0, "j": tj + c0, "obs": obs}
        for name, k in K.items():
            so = correlate(bal, k, mode="constant")[ti, tj]
            se = correlate(E, k, mode="constant")[ti, tj]
            ek = np.where(se > 0, so * E[ti, tj] / np.maximum(se, 1e-300), np.nan)
            lam = ek / np.maximum(wz[ti + r0] * wz[tj + c0], 1e-300)
            part["e_" + name] = lam
            part["p_" + name] = np.where(np.isfinite(lam) & (lam > 0), poisson.sf(obs - 1, np.maximum(lam, 1e-12)), 1.0)
        parts.append(part)
    if not parts:
        return []
    allp = {k: np.concatenate([pt[k] for pt in parts]) for k in parts[0]}
    q = {k: _bh(allp["p_" + k]) for k in K}
    oe = {k: np.where(np.isfinite(allp["e_" + k]) & (allp["e_" + k] > 0), allp["obs"] / np.maximum(allp["e_" + k], 1e-300), 0.0)
          for k in K}
    ok = np.all([q[k] <= fdr for k in K], axis=0)
    ok &= (oe["donut"] >= 1.75) & (oe["lower_left"] >= 1.75) & (oe["horizontal"] >= 1.5) & (oe["vertical"] >= 1.5)
    ok &= np.maximum(oe["donut"], oe["lower_left"]) >= 2.0
    idx = np.flatnonzero(ok)
    if not len(idx):
        return []
    calls = [(int(allp["i"][t]), int(allp["j"][t]), float(allp["obs"][t]), float(allp["e_donut"][t]), float(oe["donut"][t]),
              {k: float(q[k][t]) for k in K}) for t in idx]
    from scipy.spatial import cKDTree
    pts = np.array([[c[0], c[1]] for c in calls], float)
    radius = max(1.0, cluster_bp / resolution)
    parent = list(range(len(calls)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in cKDTree(pts).query_pairs(radius + 1e-9):
        parent[find(a)] = find(b)
    groups: dict[int, list[int]] = {}
    for k in range(len(calls)):
        groups.setdefault(find(k), []).append(k)
    out = []
    for members in groups.values():
        best = max(members, key=lambda k: calls[k][4])
        i, j, o, e, oe, q = calls[best]
        if len(members) == 1 and sum(q.values()) > 0.02:
            continue
        out.append(LoopCall(i, j, o, e, oe, q, len(members)))
    return sorted(out, key=lambda c: (c.i, c.j))


def topdom(ci, cj, cm, n: int, window: int = 5, alpha: float = 0.05) -> tuple[list[int], np.ndarray]:
    """TopDom-like boundaries and the per-bin diamond signal."""
    from scipy.stats import mannwhitneyu
    ci, cj, cm = np.asarray(ci, np.int64), np.asarray(cj, np.int64), np.asarray(cm, float)
    a, b = np.minimum(ci, cj), np.maximum(ci, cj)
    near = (b - a) <= 2 * window
    a, b, v = a[near], b[near], cm[near]
    look = {}
    for x, y, z in zip(a, b, v):
        look[(int(x), int(y))] = float(z)
    sig = np.full(n, np.nan)
    for i in range(window, n - window):
        vals = [look.get((x, y), 0.0) for x in range(i - window + 1, i + 1) for y in range(i + 1, i + window + 1)]
        sig[i] = float(np.mean(vals))
    cand = []
    for i in range(window, n - window):
        if not np.isfinite(sig[i]):
            continue
        seg = sig[max(0, i - window):i + window + 1]
        if sig[i] <= np.nanmin(seg):
            cand.append(i)
    kept = []
    for i in cand:
        diamond = [look.get((x, y), 0.0) for x in range(i - window + 1, i + 1) for y in range(i + 1, i + window + 1)]
        inside = [look.get((x, y), 0.0) for x in range(i - window + 1, i + 1) for y in range(x + 1, i + 1)] + \
                 [look.get((x, y), 0.0) for x in range(i + 1, i + window + 1) for y in range(x + 1, i + window + 1)]
        if len(inside) < 3 or len(diamond) < 3:
            continue
        if np.allclose(diamond, inside[0]) and np.allclose(inside, inside[0]):
            continue
        pv = mannwhitneyu(diamond, inside, alternative="less").pvalue
        if pv < alpha:
            kept.append(i)
    return kept, sig


def _dense_oe(ci, cj, cm, n: int, lo: int, hi: int) -> np.ndarray:
    m = hi - lo
    mat = np.zeros((m, m))
    ci, cj, cm = np.asarray(ci, np.int64), np.asarray(cj, np.int64), np.asarray(cm, float)
    sel = (ci >= lo) & (ci < hi) & (cj >= lo) & (cj < hi)
    mat[ci[sel] - lo, cj[sel] - lo] = cm[sel]
    mat = mat + np.triu(mat, 1).T
    oe = np.zeros_like(mat)
    for d in range(m):
        diag = np.diagonal(mat, d)
        mu = diag.mean()
        if mu > 0:
            idx = np.arange(m - d)
            oe[idx, idx + d] = diag / mu
            oe[idx + d, idx] = diag / mu
    return oe


def arrowhead(ci, cj, cm, n: int, bounds: list[int], min_ratio: float = 1.25, max_bins: int = 400,
              min_size: int = 3) -> list[tuple[int, int, float]]:
    """Domains (a, b, ratio) between candidate boundaries whose inside O/E exceeds the flanks by min_ratio."""
    edges = sorted(set([0] + [int(x) for x in bounds if 0 < x < n] + [n]))
    out = []
    for x in range(len(edges)):
        for y in range(x + 1, len(edges)):
            a, b = edges[x], edges[y]
            size = b - a
            if size < min_size or size > max_bins:
                continue
            lo, hi = max(0, a - size), min(n, b + size)
            oe = _dense_oe(ci, cj, cm, n, lo, hi)
            A, B = a - lo, b - lo
            inside = oe[A:B, A:B][np.triu_indices(size, 1)]
            left = oe[max(0, A - size):A, A:B].ravel() if A > 0 else np.array([])
            right = oe[A:B, B:min(len(oe), B + size)].ravel() if B < len(oe) else np.array([])
            flank = np.concatenate([left, right])
            if inside.size and flank.size and flank.mean() > 0:
                ratio = float(inside.mean() / flank.mean())
                if ratio >= min_ratio:
                    out.append((a, b, ratio))
    return out


def p_of_s(ci, cj, cm, n: int, resolution: int, weights: np.ndarray | None = None,
           max_sep_bp: int | None = None) -> dict:
    """Contact probability against separation (mean per diagonal, normalised to sum 1) and its log-log slope."""
    max_d = n if max_sep_bp is None else min(n, int(max_sep_bp // resolution) + 1)
    e = expected_balanced(ci, cj, cm, n, weights, max_d)
    s = np.arange(max_d) * resolution
    ok = (np.arange(max_d) >= 1) & (e > 0)
    pr = np.where(ok, e / max(e[ok].sum(), 1e-300), np.nan)
    slope = float("nan")
    fit = ok & (np.arange(max_d) >= 2) & (np.arange(max_d) <= max(3, max_d // 10))
    if fit.sum() >= 3:
        slope = float(np.polyfit(np.log10(s[fit]), np.log10(pr[fit]), 1)[0])
    return {"separation_bp": s, "p": pr, "slope": slope}


@dataclass
class SuiteReport:
    resolution: int
    n: int
    loops: list[LoopCall]
    insulation: np.ndarray
    insulation_boundaries: list[int]
    topdom_boundaries: list[int]
    arrowhead_domains: list[tuple[int, int, float]]
    compartment: np.ndarray
    p_s: dict
    notes: list[str] = field(default_factory=list)

    def summary(self) -> dict:
        return {"resolution": self.resolution, "bins": self.n, "loops": len(self.loops),
                "insulation_boundaries": len(self.insulation_boundaries), "topdom_boundaries": len(self.topdom_boundaries),
                "arrowhead_domains": len(self.arrowhead_domains),
                "a_fraction": None if not np.isfinite(self.compartment).any() else
                round(float(np.mean(self.compartment[np.isfinite(self.compartment)] > 0)), 3),
                "p_s_slope": None if not np.isfinite(self.p_s["slope"]) else round(self.p_s["slope"], 3)}


def run_suite(ci, cj, cm, n: int, resolution: int, weights: np.ndarray | None = None, orient: np.ndarray | None = None,
              loops: bool = True, loop_fdr: float = 0.1) -> SuiteReport:
    """Every analysis on one map (sparse counts of n bins at `resolution`)."""
    notes = []
    w = DOM.window_for(resolution, n)
    ins = DOM.insulation(ci, cj, cm, n, w)
    bnd = DOM.boundaries(ins, w)
    td, _ = topdom(ci, cj, cm, n, window=max(3, w // 2))
    arr = arrowhead(ci, cj, cm, n, sorted(set(bnd) | set(td)))
    ev, _ = DOM.compartments(ci, cj, cm, n, orient)
    if orient is None:
        notes.append("Compartment sign is arbitrary (no GC or gene-density track given).")
    lp = call_loops(ci, cj, cm, n, resolution, weights, fdr=loop_fdr) if loops else []
    if weights is None and loops:
        notes.append("Loops called on raw counts (no balancing weights given).")
    return SuiteReport(resolution, n, lp, ins, bnd, td, arr, ev, p_of_s(ci, cj, cm, n, resolution, weights), notes)
