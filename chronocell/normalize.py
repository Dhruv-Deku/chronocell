"""
Contact-map normalisation: ICE matrix balancing.

ICE = "Iterative Correction and Eigenvector decomposition" (Imakaev et al., Nat Methods 2012). The
iterative-correction half assumes every bin is equally "visible" to the experiment. Biases from
mappability, GC content and restriction-site density multiply each contact: M_ij = b_i b_j T_ij.
Rows of the corrected map T are made to sum to the same total by repeated scaling:

    s_i = sum_j M_ij,  s <- s / mean(s),  M_ij <- M_ij / (s_i s_j),  b_i <- b_i s_i

This runs until the variance of the correction falls below `tol` (cooler's convention). As in
cooler's defaults:
- the first `ignore_diags` diagonals (self and nearest-neighbour ligation products) are left out of
  the sums;
- bins with too few contacts are masked (bias NaN) instead of being inflated.

Filtering is repeated until every kept bin still has `min_nnz` contacts with other kept bins.
Otherwise a bin whose partners were masked keeps a near-empty row, and the balancing oscillates
without converging on sparse maps.

(The eigenvector half, A/B compartments, lives in chronocell.domains.)

Works on sparse (i, j, count) lists in O(nnz) memory per iteration, so whole chromosomes are fine.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class BalanceResult:
    values: np.ndarray        # balanced value for every input pixel (NaN where a bin is masked)
    bias: np.ndarray          # (n,) multiplicative bias per bin, NaN = masked; balanced = raw / (bias_i bias_j)
    masked: np.ndarray        # (n,) bool, bins excluded (low coverage)
    iterations: int
    converged: bool
    row_sum_cv: float         # coefficient of variation of row sums after balancing (0 = perfect)


def ice_balance(ci: np.ndarray, cj: np.ndarray, cm: np.ndarray, n: int, ignore_diags: int = 2,
                min_nnz: int = 10, mad_max: float = 5.0, max_iter: int = 5000, tol: float = 1e-5) -> BalanceResult:
    """ICE-balance a sparse, upper- or lower-triangular contact list (each pair listed once).

    min_nnz  bins with fewer non-zero pixels to other kept bins (outside the ignored diagonals) are
             masked, repeatedly until stable.
    mad_max  bins whose log coverage lies more than mad_max median absolute deviations below the
             median are masked too (cooler's filter), so near-empty rows are not blown up.
    """
    ci = np.asarray(ci, dtype=np.int64)
    cj = np.asarray(cj, dtype=np.int64)
    cm = np.asarray(cm, dtype=np.float64)
    if not (len(ci) == len(cj) == len(cm)):
        raise ValueError("ci, cj and cm must have the same length.")
    if len(cm) and (np.any(cm < 0) or not np.all(np.isfinite(cm))):
        raise ValueError("Contact counts must be finite and non-negative.")
    use = np.abs(cj - ci) >= ignore_diags
    i, j, v = ci[use], cj[use], cm[use]

    def rowsum(w: np.ndarray) -> np.ndarray:
        return np.bincount(i, w, minlength=n) + np.bincount(j, w, minlength=n)

    nnz = np.bincount(i, (v > 0).astype(float), minlength=n) + np.bincount(j, (v > 0).astype(float), minlength=n)
    cov = rowsum(v)
    good = (nnz >= min_nnz) & (cov > 0)
    if good.sum() >= 3 and mad_max > 0:
        lc = np.log(cov[good])
        med = np.median(lc)
        mad = np.median(np.abs(lc - med))
        if mad > 0:
            ok = lc >= med - mad_max * mad
            idx = np.flatnonzero(good)
            good[idx[~ok]] = False
    for _ in range(n):                                        # re-filter on contacts among kept bins
        live_nz = good[i] & good[j] & (v > 0)
        new = good & (np.bincount(i, live_nz.astype(float), minlength=n)
                      + np.bincount(j, live_nz.astype(float), minlength=n) >= min_nnz)
        if (new == good).all():
            break
        good = new
    if good.sum() < 2:
        raise ValueError("Too few bins with enough contacts to balance.")

    bias = np.ones(n)
    live = good[i] & good[j]
    w = np.where(live, v, 0.0)
    converged, it = False, 0
    for it in range(1, max_iter + 1):
        s = rowsum(w)
        s = s / s[good].mean()
        s[~good] = 1.0
        s[s == 0] = 1.0
        w = w / (s[i] * s[j])
        bias *= s
        if float(np.var(s[good])) < tol:
            converged = True
            break
    final = rowsum(w)[good]
    row_cv = float(final.std() / final.mean()) if final.mean() > 0 else float("nan")

    bias_out = np.where(good, bias, np.nan)
    b_all = np.where(good, bias, np.nan)
    values = cm / (b_all[ci] * b_all[cj])                     # every input pixel, including ignored diagonals
    return BalanceResult(values, bias_out, ~good, it, converged, row_cv)


def balanced_contacts(ci: np.ndarray, cj: np.ndarray, cm: np.ndarray, n: int,
                      **kwargs) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """ICE-balance a contact list and return it in the same (ci, cj, cm) form the rest of ChronoCell
    uses. Pixels touching masked bins are dropped. Values are rescaled so their median matches the
    raw median: downstream steps use ratios, but feature scales (log1p counts) stay familiar."""
    res = ice_balance(ci, cj, cm, n, **kwargs)
    ok = np.isfinite(res.values) & (res.values > 0)
    raw = np.asarray(cm, dtype=np.float64)[ok]
    vals = res.values[ok]
    vals = vals * (np.median(raw) / np.median(vals)) if vals.size else vals
    notes = [f"ICE balanced in {res.iterations} iterations ({'converged' if res.converged else 'NOT converged'}; "
             f"row-sum CV {res.row_sum_cv:.1e}); {int(res.masked.sum()):,} low-coverage bins masked"]
    return np.asarray(ci)[ok], np.asarray(cj)[ok], vals, notes


# ======================================================================================
# KR balancing (Knight & Ruiz, IMA J Numer Anal 2013), as Juicer's "KR" normalisation (Phase B4)
# ======================================================================================
@dataclass
class KRResult:
    weights: np.ndarray       # (n,) balanced = raw * w_i * w_j; NaN = masked
    masked: np.ndarray
    iterations: int
    converged: bool
    row_sum_cv: float


def kr_balance(ci: np.ndarray, cj: np.ndarray, cm: np.ndarray, n: int, ignore_diags: int = 1,
               min_coverage_quantile: float = 0.02, tol: float = 1e-6, max_iter: int = 300) -> KRResult:
    """Symmetric matrix balancing by the Knight-Ruiz inner-outer Newton iteration (their bnewt). Bins with
    no contacts, or in the lowest `min_coverage_quantile` of coverage, are masked first (the method needs a
    matrix with total support). Each pair is listed once (i != j); the first `ignore_diags` diagonals are left
    out, as for ICE."""
    from scipy import sparse
    ci, cj, cm = np.asarray(ci, np.int64), np.asarray(cj, np.int64), np.asarray(cm, float)
    use = (np.abs(cj - ci) >= ignore_diags) & (cm > 0)
    i, j, v = ci[use], cj[use], cm[use]
    cov = np.bincount(i, weights=v, minlength=n) + np.bincount(j, weights=v, minlength=n)
    masked = cov <= 0
    if (~masked).sum() > 10:
        thr = np.quantile(cov[~masked], min_coverage_quantile)
        masked |= cov < thr
    keep = ~masked
    idx = -np.ones(n, np.int64)
    idx[keep] = np.arange(keep.sum())
    sel = keep[i] & keep[j]
    m = int(keep.sum())
    if m < 3:
        return KRResult(np.full(n, np.nan), masked, 0, False, float("nan"))
    a, b, w = idx[i[sel]], idx[j[sel]], v[sel]
    A = sparse.coo_matrix((np.r_[w, w], (np.r_[a, b], np.r_[b, a])), shape=(m, m)).tocsr()
    A = A / float(A.sum() / m)                    # scale so the solution is O(1)
    x, it, ok = _bnewt(A, tol, max_iter)
    rs = x * (A @ x)
    weights = np.full(n, np.nan)
    # back to the raw scale: balanced = raw * w_i * w_j with A = raw / s, so w = x / sqrt(s)
    s = float(np.r_[w, w].sum() / m)
    weights[keep] = x / np.sqrt(s)
    return KRResult(weights, masked, it, ok, float(np.std(rs) / max(np.mean(rs), 1e-300)))


def _bnewt(A, tol: float, max_iter: int, delta: float = 0.1, Delta: float = 3.0):
    n = A.shape[0]
    e = np.ones(n)
    x = e.copy()
    g, etamax = 0.9, 0.1
    eta, stop_tol, rt = etamax, tol * 0.5, tol ** 2
    v = x * (A @ x)
    rk = 1 - v
    rho_km1 = float(rk @ rk)
    rout = rold = rho_km1
    it = 0
    while rout > rt and it < max_iter:
        it += 1
        k = 0
        y = e.copy()
        innertol = max(eta ** 2 * rout, rt)
        rho_km2 = rho_km1
        p = np.zeros(n)
        Z = rk / v
        while rho_km1 > innertol:
            k += 1
            if k == 1:
                Z = rk / v
                p = Z.copy()
                rho_km1 = float(rk @ Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p
            wv = x * (A @ (x * p)) + v * p
            alpha = rho_km1 / float(p @ wv)
            ap = alpha * p
            ynew = y + ap
            if ynew.min() <= delta:
                ind = ap < 0
                gamma = float(((delta - y[ind]) / ap[ind]).min()) if ind.any() else 0.0
                y = y + gamma * ap
                break
            if ynew.max() >= Delta:
                ind = ynew > Delta
                gamma = float(((Delta - y[ind]) / ap[ind]).min())
                y = y + gamma * ap
                break
            y = ynew
            rk = rk - alpha * wv
            rho_km2 = rho_km1
            Z = rk / v
            rho_km1 = float(rk @ Z)
            if k > 10 * n:
                break
        x = x * y
        v = x * (A @ x)
        rk = 1 - v
        rho_km1 = float(rk @ rk)
        rout = rho_km1
        rat = rout / rold
        rold = rout
        res_norm = np.sqrt(rout)
        eta_o = eta
        eta = g * rat
        if g * eta_o ** 2 > 0.1:
            eta = max(eta, g * eta_o ** 2)
        eta = max(min(eta, etamax), stop_tol / max(res_norm, 1e-300))
    return x, it, bool(rout <= rt)
