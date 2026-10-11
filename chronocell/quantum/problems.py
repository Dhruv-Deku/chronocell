"""
ChronoCell problems written as QUBOs (see qubo.py for the form).

TAD boundaries (tad_qubo)
    One bit per gap between neighbouring bins of a window: x_k = 1 puts a domain boundary between bins k-1 and k.
    A contact (a, b) is "cut" when a boundary lies between its two bins. The objective is modularity restricted to
    contiguous domains: cutting a pair costs its observed-minus-expected contacts,
        E(x) = sum_{a<b, b-a<=S} w_ab cut_ab(x),   w_ab = (M_ab - gamma E[M | b-a]) / Z,
    so domains keep enriched contacts inside and boundaries fall where crossing contacts are depleted; gamma is the
    resolution (larger: more boundaries). cut_ab = 1 - prod_k (1 - x_k) is replaced by its second-order expansion
    sum_k x_k - sum_{k<l} x_k x_l, which is EXACT whenever at most two boundaries lie between a and b. A penalty
    lambda x_k x_l for boundaries closer than m bins enforces a minimum domain size; with the span S <= 2m every
    feasible solution has at most two boundaries inside any counted pair, so the QUBO equals the modularity
    objective on all feasible solutions.
Variant set (variant_set_qubo)     pick k variants with the highest impact scores and least overlap.
Drug combination (drug_combo_*)   doses of each drug class (2 bits: 0, 1/3, 2/3, 1) maximising simulated restoration
                                  minus a dose cost; restoration modelled as a quadratic in the doses fitted to the
                                  Drug lab's own simulations (singles across doses, pairs at full dose) and the chosen
                                  mix is re-simulated to check it.
Gene group (gene_group_qubo)      pick k genes that are active and close together in 3D (a candidate hub).
"""

from __future__ import annotations

import math

import numpy as np

from .qubo import QUBO


# ======================================================================================
# TAD boundaries
# ======================================================================================
def coarsen(M: np.ndarray, factor: int) -> np.ndarray:
    """Sum f x f blocks of a dense contact matrix (the trailing remainder is dropped)."""
    M = np.asarray(M, float)
    if factor <= 1:
        return M
    n = (M.shape[0] // factor) * factor
    m = M[:n, :n].reshape(n // factor, factor, n // factor, factor).sum(axis=(1, 3))
    return m


def expected_by_distance(M: np.ndarray, max_d: int | None = None, mask: np.ndarray | None = None) -> np.ndarray:
    """Mean contacts at each separation d (bins with no coverage excluded)."""
    M = np.asarray(M, float)
    n = len(M)
    ok = np.ones(n, bool) if mask is None else np.asarray(mask, bool)
    max_d = n - 1 if max_d is None else min(max_d, n - 1)
    out = np.zeros(max_d + 1)
    for d in range(max_d + 1):
        i = np.arange(n - d)
        good = ok[i] & ok[i + d]
        out[d] = float(M[i[good], i[good] + d].mean()) if good.any() else 0.0
    return out


def coverage_mask(M: np.ndarray, min_frac: float = 0.2) -> np.ndarray:
    """Bins whose total contacts reach min_frac of the median bin (gaps and unmappable bins are False)."""
    cov = np.asarray(M, float).sum(axis=1)
    med = float(np.median(cov[cov > 0])) if np.any(cov > 0) else 0.0
    return cov >= min_frac * med if med > 0 else cov > 0


def tad_qubo(M: np.ndarray, core: tuple[int, int], gamma: float = 1.0, min_size: int = 3, span: int | None = None,
             expected: np.ndarray | None = None, mask: np.ndarray | None = None, c_lambda: float = 2.0,
             boundary_cost: float = 0.0, weight: str = "difference") -> tuple[QUBO, np.ndarray]:
    """QUBO for domain boundaries inside bins [core[0], core[1]) of the dense matrix M (which may extend beyond the
    core: flanking bins add the contacts that cross the core's edges). Returns (QUBO, positions) where positions[v]
    is the bin index k of variable v (a boundary between bins k-1 and k).
    weight: "difference" (M - gamma E) / Z, or "log" log2((M + 1) / (gamma E + 1)) (robust to depth);
    boundary_cost: a cost added for every boundary, in units of the mean |linear coefficient| (fewer, stronger
    boundaries)."""
    M = np.asarray(M, float)
    N = len(M)
    c0, c1 = int(core[0]), int(core[1])
    if not (0 <= c0 < c1 <= N) or c1 - c0 < 3:
        raise ValueError("the core window must hold at least 3 bins inside the matrix")
    span = 2 * min_size if span is None else int(span)
    mask = coverage_mask(M) if mask is None else np.asarray(mask, bool)
    exp = expected_by_distance(M, span, mask) if expected is None else np.asarray(expected, float)
    Z = float(np.mean(exp[1:span + 1])) or 1.0
    pos = np.arange(c0 + 1, c1)
    nv = len(pos)
    var_of = {int(k): v for v, k in enumerate(pos)}
    lin = np.zeros(nv)
    quad = np.zeros((nv, nv))
    for a in range(max(0, c0 - span), min(N, c1 + span)):
        if not mask[a]:
            continue
        for b in range(a + 1, min(N, a + span + 1)):
            if not mask[b]:
                continue
            vs = [var_of[k] for k in range(a + 1, b + 1) if k in var_of]
            if not vs:
                continue
            if weight == "log":
                w = math.log2((M[a, b] + 1.0) / (gamma * exp[b - a] + 1.0))
            else:
                w = (M[a, b] - gamma * exp[b - a]) / Z
            for i, v in enumerate(vs):
                lin[v] += w
                for u in vs[i + 1:]:
                    quad[v, u] -= w
    lin += boundary_cost * (float(np.mean(np.abs(lin))) if nv else 0.0)
    lam = c_lambda * (float(np.abs(lin).max(initial=0.0)) + 1e-9)
    for v in range(nv):
        for u in range(v + 1, min(nv, v + min_size)):
            quad[v, u] += lam
    q = QUBO(lin, quad, 0.0, [f"gap {int(k)}" for k in pos], "TAD boundaries")
    return q, pos


def boundaries_from_bits(bits: np.ndarray, positions: np.ndarray) -> list[int]:
    return [int(positions[v]) for v in np.flatnonzero(np.asarray(bits) > 0)]


def feasible(bits: np.ndarray, min_size: int) -> bool:
    on = np.flatnonzero(np.asarray(bits) > 0)
    return bool(np.all(np.diff(on) >= min_size)) if len(on) > 1 else True


def match(pred: list[float], ref: list[float], tol: float) -> tuple[int, int, int]:
    """One-to-one matching within tol (greedy by distance): (true positives, false positives, false negatives)."""
    pairs = sorted((abs(p - r), i, j) for i, p in enumerate(pred) for j, r in enumerate(ref) if abs(p - r) <= tol)
    used_p, used_r = set(), set()
    for _, i, j in pairs:
        if i not in used_p and j not in used_r:
            used_p.add(i)
            used_r.add(j)
    tp = len(used_p)
    return tp, len(pred) - tp, len(ref) - tp


def f1(tp: int, fp: int, fn: int) -> dict:
    p = tp / (tp + fp) if tp + fp else float("nan")
    r = tp / (tp + fn) if tp + fn else float("nan")
    f = 2 * tp / (2 * tp + fp + fn) if tp + fp + fn else float("nan")
    return {"precision": p, "recall": r, "f1": f, "tp": tp, "fp": fp, "fn": fn}


# ======================================================================================
# k-of-n selection helpers
# ======================================================================================
def _k_hot(n: int, k: int, A: float) -> tuple[np.ndarray, np.ndarray, float]:
    """A (sum x - k)^2 = A[sum x_i (1 - 2k) + 2 sum_{i<j} x_i x_j + k^2]."""
    lin = np.full(n, A * (1 - 2 * k))
    quad = np.triu(np.full((n, n), 2 * A), 1)
    return lin, quad, A * k * k


def variant_set_qubo(scores: np.ndarray, overlap: np.ndarray, k: int, redundancy: float = 1.0,
                     labels: list | None = None) -> QUBO:
    """Choose k variants: maximise sum of scores, penalise pairs with overlapping effects (overlap in [0, 1])."""
    s = np.asarray(scores, float)
    n = len(s)
    sc = float(np.abs(s).max(initial=0.0)) or 1.0
    s = s / sc
    O = np.triu(np.asarray(overlap, float), 1)
    A = 1.0 + float(np.abs(s).sum()) + redundancy * float(O.sum())
    kl, kq, kc = _k_hot(n, k, A)
    return QUBO(-s + kl, redundancy * O + kq, kc, labels or [f"variant {i + 1}" for i in range(n)], "Variant set")


def gene_group_qubo(activity: np.ndarray, proximity: np.ndarray, k: int, together: float = 1.0,
                    labels: list | None = None) -> QUBO:
    """Choose k genes: active (activity, any scale) and close together in 3D (proximity in [0, 1])."""
    a = np.asarray(activity, float)
    n = len(a)
    a = (a - a.min()) / ((a.max() - a.min()) or 1.0)
    P = np.triu(np.asarray(proximity, float), 1)
    A = 1.0 + float(a.sum()) + together * float(P.sum())
    kl, kq, kc = _k_hot(n, k, A)
    return QUBO(-a + kl, -together * P + kq, kc, labels or [f"gene {i + 1}" for i in range(n)], "Gene group")


# ======================================================================================
# Drug combination
# ======================================================================================
LEVELS = np.array([0.0, 1 / 3, 2 / 3, 1.0])


def drug_doses(bits: np.ndarray, n_drugs: int) -> np.ndarray:
    b = np.asarray(bits, float).reshape(n_drugs, 2)
    return (b[:, 0] + 2 * b[:, 1]) / 3


def drug_combo_qubo(alpha: np.ndarray, beta: np.ndarray, eta: np.ndarray, dose_cost: float,
                    names: list[str]) -> QUBO:
    """Minimise -R(d) + dose_cost * sum d with R(d) = sum (alpha_k d_k + beta_k d_k^2) + sum_{k<l} eta_kl d_k d_l
    (restoration in % of the way back to healthy), d_k = (b_k0 + 2 b_k1) / 3 on two bits per drug."""
    D = len(alpha)
    n = 2 * D
    lin = np.zeros(n)
    quad = np.zeros((n, n))
    w = np.array([1 / 3, 2 / 3])
    for k in range(D):
        i0, i1 = 2 * k, 2 * k + 1
        # alpha d - cost d:  linear in bits
        lin[i0] += -(alpha[k] - dose_cost) * w[0]
        lin[i1] += -(alpha[k] - dose_cost) * w[1]
        # beta d^2 = beta (b0 + 4 b1 + 4 b0 b1) / 9
        lin[i0] += -beta[k] / 9
        lin[i1] += -4 * beta[k] / 9
        quad[i0, i1] += -4 * beta[k] / 9
    for k in range(D):
        for l in range(k + 1, D):
            e = float(eta[k, l])
            if e == 0:
                continue
            for a in range(2):
                for b in range(2):
                    quad[2 * k + a, 2 * l + b] += -e * w[a] * w[b]
    labels = [f"{nm} · {lv}" for nm in names for lv in ("+1/3", "+2/3")]
    return QUBO(lin, quad, 0.0, labels, "Drug combination")


def fit_dose_response(doses: np.ndarray, restoration: np.ndarray) -> tuple[float, float]:
    """Least-squares R(d) = alpha d + beta d^2 through the origin."""
    d = np.asarray(doses, float)
    r = np.asarray(restoration, float)
    ok = np.isfinite(r)
    A = np.stack([d[ok], d[ok] ** 2], axis=1)
    if len(A) < 2:
        return float(r[ok][-1]) if ok.any() else 0.0, 0.0
    coef, *_ = np.linalg.lstsq(A, r[ok], rcond=None)
    return float(coef[0]), float(coef[1])


def combo_restoration(x: np.ndarray, signal: np.ndarray, valid: np.ndarray, b0: float, x_healthy: np.ndarray,
                      doses: dict[str, float], efficacy: float = 0.8, relax_iters: int = 15,
                      signal_ref: np.ndarray | None = None) -> float:
    """Restoration (%) of a combination, simulated by applying the drugs one after another (Drug lab model, in
    the order of therapy.DRUGS) and measuring the RMSD to healthy against the untreated fold's."""
    from .. import physics, therapy as TH
    y = np.asarray(x, float)
    _, xh, _ = physics.kabsch_rmsd(y, np.asarray(x_healthy, float))
    r0 = physics.kabsch_rmsd(xh, y)[0]
    for key in TH.DRUGS:
        d = float(doses.get(key, 0.0))
        if d <= 0:
            continue
        res = TH.simulate_treatment(y, signal, valid, b0, key, x_healthy, efficacy, np.array([d]), relax_iters,
                                    signal_ref=signal_ref)
        y = res.frames[-1]
    r = physics.kabsch_rmsd(xh, y)[0]
    return float(np.clip(100 * (1 - r / r0), -100, 100)) if r0 > 0 else float("nan")


def drug_combo_model(x: np.ndarray, signal: np.ndarray, valid: np.ndarray, b0: float, x_healthy: np.ndarray,
                     efficacy: float = 0.8, relax_iters: int = 15, signal_ref: np.ndarray | None = None,
                     progress=None) -> dict:
    """Fit the quadratic restoration model from Drug lab simulations: each drug at 4 doses, each pair at full dose."""
    from .. import therapy as TH
    keys = list(TH.DRUGS)
    D = len(keys)
    alpha = np.zeros(D)
    beta = np.zeros(D)
    full = np.zeros(D)
    singles = {}
    for k, key in enumerate(keys):
        res = TH.simulate_treatment(x, signal, valid, b0, key, x_healthy, efficacy, LEVELS, relax_iters,
                                    signal_ref=signal_ref)
        r = res.metrics["restoration_pct"].to_numpy(float)
        singles[key] = r.tolist()
        alpha[k], beta[k] = fit_dose_response(LEVELS, r)
        full[k] = float(r[-1])
        if progress:
            progress((k + 1) / (D + D * (D - 1) / 2), f"single drug {TH.DRUGS[key].name}")
    eta = np.zeros((D, D))
    pairs = {}
    step = D
    for k in range(D):
        for l in range(k + 1, D):
            rkl = combo_restoration(x, signal, valid, b0, x_healthy, {keys[k]: 1.0, keys[l]: 1.0}, efficacy,
                                    relax_iters, signal_ref)
            pairs[f"{keys[k]}+{keys[l]}"] = rkl
            eta[k, l] = rkl - full[k] - full[l]
            step += 1
            if progress:
                progress(step / (D + D * (D - 1) / 2), f"pair {keys[k]} + {keys[l]}")
    return {"keys": keys, "names": [TH.DRUGS[k].name for k in keys], "alpha": alpha, "beta": beta, "eta": eta,
            "singles": singles, "pairs": pairs, "levels": LEVELS.tolist()}


def predicted_restoration(model: dict, doses: np.ndarray) -> float:
    d = np.asarray(doses, float)
    r = float(np.sum(model["alpha"] * d + model["beta"] * d ** 2))
    D = len(d)
    for k in range(D):
        for l in range(k + 1, D):
            r += float(model["eta"][k, l] * d[k] * d[l])
    return r


# ======================================================================================
# Sizes
# ======================================================================================
def scaling_rows() -> list[dict]:
    """How many qubits each problem needs as it grows, and the memory one statevector of that size takes."""
    rows = []
    for name, size_label, sizes, qubits in (
            ("TAD boundaries", "bins in the window", [10, 20, 30, 100, 1000], lambda s: s - 1),
            ("Lattice fold (2D)", "beads", [6, 8, 10, 12, 30], lambda s: 2 * (s - 2)),
            ("Lattice fold (3D)", "beads", [5, 6, 7, 10, 30], lambda s: 3 * (s - 2)),
            ("Drug combination", "drug classes", [2, 4, 8, 16, 50], lambda s: 2 * s),
            ("Gene group", "candidate genes", [8, 12, 20, 50, 500], lambda s: s),
            ("Quantum walk", "bins (position register)", [64, 256, 1024, 10_000, 250_000],
             lambda s: int(math.ceil(math.log2(s)))),
    ):
        for s in sizes:
            q = qubits(s)
            rows.append({"problem": name, "size": s, "size_unit": size_label, "qubits": q,
                         "statevector_bytes": 16.0 * 2.0 ** q})
    return rows
