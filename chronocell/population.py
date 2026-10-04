"""
v4 population model: the v3.3 maximum-entropy Gaussian ensemble (chronocell.ensemble), scaled from
400-bead windows to whole chromosomes, with per-pair uncertainty.

What stays the same as v3.3 (so results remain comparable)
- Contact frequency f_ij -> per-axis pair variance s_ij through the Maxwell distribution
  (ensemble.sigma_from_contact_probability), with half-count clipping and binomial reliability
  weights w_ij = N f / (1 - f).
- The objective: weighted least squares on log pair variances, minimised with Adam, keeping the
  best covariance seen.
- Exact statistics: every per-pair number is a closed-form property of the fitted Gaussian.
- Exact Langevin (Ornstein-Uhlenbeck) trajectories by mode decomposition of the covariance.

What changes, and why
1. Parameterisation. Bead positions are x = A z + w with A an N x r matrix (r <= rank cap) and w an
   independent-bond random walk with per-bond variances v_k. Pair variance:
       s_ij = |a_i - a_j|^2 + sum_{k=i}^{j-1} v_k .
   With r = N this is the v3.3 family (the random walk is then redundant); with r < N the random
   walk carries the short-range variance the truncated rank would otherwise lose. Both parts are
   covariance matrices, so every parameter value is a valid Gaussian chain ensemble.
2. Exact gradients in row blocks. dL/dA = 2 (diag(E 1) - E) A with e_ij = dL/ds_ij, computed one
   block of rows at a time, so memory is O(block x N) beyond the N x N targets, never the
   O(N^2 x r) of automatic differentiation. Cost per iteration: two GEMMs of N^2 r flops.
3. Coarse-to-fine start. For large N the start comes from a coarse fit on block centroids (the
   centroid pair variance is exactly mean_{i in B1, j in B2} s_ij - D_B1 - D_B2, with D_B the mean
   squared deviation inside a block), prolongated to the beads; otherwise from the PSD projection
   used by v3.3.
4. float32 by default, CUDA when present, CPU otherwise.

Physics, stated plainly
- Members are Gaussian chains: bond lengths fluctuate (the adjacent-pair *median* is anchored, not
  each bond) and there is no excluded volume. `excluded_volume_members` relaxes sampled members to
  fixed bond length and no overlaps for display; it changes pair distances, and the change is
  measured and reported, so population statistics always come from the exact Gaussian ensemble.
- Uncertainty here is the cell-to-cell spread the ensemble predicts (a Maxwell distribution per
  pair). It is calibrated against held-out imaging in validation/ (Pillar 2), not assumed correct.
"""

from __future__ import annotations

import json
import math
import time
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Callable

import numpy as np
import torch

from . import ensemble as ENS

MAXWELL_MEDIAN = ENS.MAXWELL_MEDIAN                  # median of |N(0, I_3)| / sigma
MAXWELL_MEAN = 2.0 * math.sqrt(2.0 / math.pi)        # mean of |N(0, I_3)| / sigma
MAXWELL_SD = math.sqrt(3.0 - 8.0 / math.pi)          # sd of |N(0, I_3)| / sigma


# ======================================================================================
# Maxwell distribution of a pair distance |r|, r ~ N(0, sigma^2 I_3)
# ======================================================================================
def maxwell_cdf(x: np.ndarray) -> np.ndarray:
    """P(|r| < x sigma) for a 3D Gaussian pair vector with per-axis sd sigma (the same 20,000-point
    table as ensemble.contact_probability_from_sigma; interpolation error < 1e-6)."""
    x = np.asarray(x, dtype=np.float64)
    return np.where(x <= 0, 0.0, np.interp(x, ENS._X_GRID, ENS._P_GRID))


def maxwell_quantile(q: np.ndarray | float) -> np.ndarray:
    """x such that P(|r| < x sigma) = q (tabulated inversion, as ensemble.sigma_from_contact_probability)."""
    q = np.clip(np.asarray(q, dtype=np.float64), ENS._P_GRID[0], ENS._P_GRID[-1])
    return np.interp(q, ENS._P_GRID, ENS._X_GRID)


def central_interval(level: float) -> tuple[float, float]:
    """(lower, upper) multipliers of sigma for the central `level` interval of the Maxwell law."""
    if not 0.0 < level < 1.0:
        raise ValueError("level must be in (0, 1)")
    return float(maxwell_quantile((1 - level) / 2)), float(maxwell_quantile((1 + level) / 2))


@dataclass(frozen=True)
class PairSummary:
    """Per-pair distance distribution predicted by the ensemble (nm)."""
    mean: float
    sd: float
    median: float
    lower: float
    upper: float
    level: float
    contact_probability: float

    def as_dict(self) -> dict:
        return asdict(self)


CALIBRATION_PATH = Path(__file__).with_name("data") / "calibration.json"
CALIBRATION_HIC_PATH = Path(__file__).with_name("data") / "calibration_hic.json"


@lru_cache(maxsize=2)
def _recalibration_cdf(kind: str = "imaging") -> tuple[np.ndarray, np.ndarray] | None:
    path = CALIBRATION_HIC_PATH if kind == "hic" else CALIBRATION_PATH
    try:
        h = np.asarray(json.loads(path.read_text(encoding="utf-8"))["pit_histogram"], dtype=np.float64)
    except (OSError, ValueError, KeyError):
        return None
    return np.linspace(0.0, 1.0, len(h) + 1), np.concatenate([[0.0], np.cumsum(h) / h.sum()])


def recalibrated_interval(level: float, kind: str = "imaging") -> tuple[float, float] | None:
    """(lower, upper) multipliers of sigma for a stated central `level`, after the quantile recalibration
    (Kuleshov et al., ICML 2018) fitted on practice data and frozen in data/calibration.json (imaging-
    derived contacts, kind "imaging") or data/calibration_hic.json (sequencing Hi-C, kind "hic"): the
    interval runs between model quantiles G^-1((1 - level)/2) and G^-1((1 + level)/2), with G the
    empirical CDF of practice PIT values. None if the calibration file is missing. Held-out tests:
    Gate 2 (imaging) and Gate 2b (Hi-C) in validation/RESULTS.md."""
    if not 0.0 < level < 1.0:
        raise ValueError("level must be in (0, 1)")
    r = _recalibration_cdf(kind)
    if r is None:
        return None
    edges, cdf = r
    u_lo, u_hi = np.interp((1 - level) / 2, cdf, edges), np.interp((1 + level) / 2, cdf, edges)
    return float(maxwell_quantile(u_lo)), float(maxwell_quantile(u_hi))


def summarise_sigma(sigma_nm: np.ndarray, r_c_nm: float, level: float = 0.8,
                    recalibrated: bool = False, kind: str = "imaging") -> dict[str, np.ndarray]:
    """Mean, sd, median and central `level` interval of the pair distance for per-axis sd sigma (nm);
    with `recalibrated`, the interval bounds use the practice-fitted recalibration for input `kind`
    ("imaging" or "hic"), if available."""
    s = np.asarray(sigma_nm, dtype=np.float64)
    lo, hi = (recalibrated and recalibrated_interval(level, kind)) or central_interval(level)
    return {"mean": MAXWELL_MEAN * s, "sd": MAXWELL_SD * s, "median": MAXWELL_MEDIAN * s, "lower": lo * s,
            "upper": hi * s, "contact_probability": ENS.contact_probability_from_sigma(s / r_c_nm)}


# ======================================================================================
# The model
# ======================================================================================
@dataclass
class GaussianChain:
    """x = A z + w: A (N x r) in units of r_c, w a random walk with per-bond variances v (N - 1,).
    Lengths are reported in nm through r_c_nm."""
    A: np.ndarray
    v: np.ndarray
    r_c_nm: float

    @property
    def n(self) -> int:
        return self.A.shape[0]

    @property
    def rank(self) -> int:
        return self.A.shape[1]

    @property
    def _cum(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.v, dtype=np.float64)])

    def pair_variance(self, i: np.ndarray, j: np.ndarray) -> np.ndarray:
        """Per-axis pair variance s_ij (units r_c^2) for arrays of index pairs."""
        i, j = np.asarray(i, np.int64), np.asarray(j, np.int64)
        d = self.A[i].astype(np.float64) - self.A[j].astype(np.float64)
        c = self._cum
        return np.einsum("pk,pk->p", d, d) + np.abs(c[j] - c[i])

    def sigma_nm(self, i: np.ndarray, j: np.ndarray) -> np.ndarray:
        return np.sqrt(np.maximum(self.pair_variance(i, j), 0.0)) * self.r_c_nm

    def pair_summary(self, i: int, j: int, level: float = 0.8) -> PairSummary:
        s = summarise_sigma(self.sigma_nm(np.array([i]), np.array([j])), self.r_c_nm, level)
        return PairSummary(*(float(s[k][0]) for k in ("mean", "sd", "median", "lower", "upper")), level,
                           float(s["contact_probability"][0]))

    def dense_variance(self, block: int = 1024, dtype=np.float32) -> np.ndarray:
        """Full N x N pair-variance matrix (units r_c^2), built block by block."""
        n = self.n
        A = self.A.astype(np.float64)
        g = np.einsum("ik,ik->i", A, A)
        c = self._cum
        out = np.empty((n, n), dtype=dtype)
        for a in range(0, n, block):
            b = min(n, a + block)
            s = g[a:b, None] + g[None, :] - 2.0 * (A[a:b] @ A.T)
            s += np.abs(c[None, :] - c[a:b, None])
            np.maximum(s, 0.0, out=s)
            out[a:b] = s
        np.fill_diagonal(out, 0.0)
        return out

    def covariance(self) -> np.ndarray:
        """Centred per-axis covariance of bead positions (units r_c^2): J (A A^T + C_rw) J."""
        n = self.n
        A = self.A.astype(np.float64)
        c = self._cum
        cov = A @ A.T + np.minimum.outer(c, c)          # Brownian covariance Cov(w_i, w_j) = C_min(i,j)
        cov -= cov.mean(axis=0, keepdims=True)
        cov -= cov.mean(axis=1, keepdims=True)
        return cov

    def sample(self, k: int, rng: np.random.Generator) -> np.ndarray:
        """k independent members (k, N, 3) in nm, drawn exactly from the ensemble."""
        z = rng.normal(size=(k, self.rank, 3))
        x = np.einsum("nr,krc->knc", self.A.astype(np.float64), z)
        steps = rng.normal(size=(k, self.n - 1, 3)) * np.sqrt(self.v)[None, :, None]
        x[:, 1:] += np.cumsum(steps, axis=1)
        x -= x.mean(axis=1, keepdims=True)
        return x * self.r_c_nm


# ======================================================================================
# Fitting
# ======================================================================================
@dataclass
class PopulationConfig:
    iterations: int = 1500            # Adam steps (v3.3 default; see validation/TUNING.md for v4 choices)
    learning_rate: float = 0.01
    weighting: str = "binomial"       # "binomial": N f / (1 - f) per pair; "uniform"
    rank_cap: int = 1024              # r = min(N, rank_cap); the random walk carries the rest
    local_term: bool = True
    init: str = "auto"                # "auto" | "pca" | "coarse"
    coarse_beads: int = 400           # beads at the coarse level of the coarse-to-fine start
    coarse_iterations: int = 600
    dense_init_max: int = 3000        # full eigendecomposition for the PSD start up to this N
    block_rows: int = 512
    dtype: str = "float32"
    device: str = "auto"
    replicas: int = 100
    frames: int = 50
    frame_interval: float = 0.05
    max_trajectory_floats: int = 60_000_000   # frames x replicas x N x 3 kept in memory (float32)
    max_eig_beads: int = 8000         # exact OU modes up to this N (one eigendecomposition)
    seed: int = 0


# Whole-chromosome defaults, chosen on practice data only and frozen in validation/frozen.py
# (WHOLE_CHROMOSOME); tests/test_v4_population.py checks the two stay identical.
WHOLE_CHROMOSOME_DEFAULTS = {"rank_cap": 256, "local_term": True, "iterations": 1500, "learning_rate": 0.01,
                             "weighting": "binomial", "dtype": "float32", "init": "auto"}
V33_MAX_BEADS = 400        # up to this size the app keeps the validated v3.3 fit (chronocell.ensemble)
MAX_BEADS = 6000           # largest population the app builds (the default resolution keeps chromosomes below it)


def config_for(n: int, **overrides) -> PopulationConfig:
    """The app's population-model settings for an n-bead window (frozen whole-chromosome defaults)."""
    return PopulationConfig(**(WHOLE_CHROMOSOME_DEFAULTS | overrides))


def _device(name: str) -> torch.device:
    if name == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(name)


def targets_from_frequency(freq: np.ndarray, n_observed: np.ndarray | float | None,
                           weighting: str = "binomial") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """v3.3 preprocessing on a dense matrix: symmetrise ignoring NaN, clip to half a count, invert to
    per-axis variance targets (units r_c^2). Returns (log_target, weight, observed) as N x N arrays,
    float32, symmetric, zero on the diagonal; weights are normalised to mean 1 over observed pairs."""
    f_in = np.asarray(freq, dtype=np.float64)
    n = len(f_in)
    if f_in.ndim != 2 or f_in.shape != (n, n) or n < 4:
        raise ValueError("freq must be a square matrix with at least 4 beads.")
    both = np.stack([f_in, f_in.T])
    k_obs = np.isfinite(both).sum(0)
    f = np.where(k_obs > 0, np.nansum(both, 0) / np.maximum(k_obs, 1), np.nan)
    observed = np.isfinite(f)
    np.fill_diagonal(observed, False)
    if observed.sum() < 2 * (n - 1):
        raise ValueError("Too few observed pairs to define an ensemble.")
    n_obs = np.broadcast_to(np.asarray(1000.0 if n_observed is None else n_observed, dtype=np.float64), (n, n))
    n_obs = np.maximum(n_obs, 1.0)
    fc = np.clip(np.where(observed, f, 0.5), 0.5 / n_obs, 1 - 0.5 / n_obs)
    log_t = 2.0 * np.log(ENS.sigma_from_contact_probability(fc))
    w = fc * n_obs / (1 - fc) if weighting == "binomial" else np.ones_like(fc)
    w = np.where(observed, w, 0.0)
    w /= w[observed].mean()
    np.fill_diagonal(log_t, 0.0)
    return log_t.astype(np.float32), w.astype(np.float32), observed


def _psd_start(log_t: np.ndarray, observed: np.ndarray, rank: int, dense_max: int, seed: int) -> np.ndarray:
    """v3.3 start: PSD projection of -1/2 J S J (S = target variances, gaps filled per separation),
    top `rank` components."""
    n = len(log_t)
    S = np.where(observed, np.exp(log_t.astype(np.float64)), np.nan)
    np.fill_diagonal(S, 0.0)
    if not observed[np.triu_indices(n, 1)].all():
        for s in range(1, n):
            d = np.diagonal(S, s)
            if np.isnan(d).any():
                fill = np.nanmedian(d) if np.isfinite(d).any() else np.nanmedian(S[np.isfinite(S) & (S > 0)])
                idx = np.arange(n - s)
                bad = np.isnan(d)
                S[idx[bad], idx[bad] + s] = fill
                S[idx[bad] + s, idx[bad]] = fill
    B = -0.5 * S
    B -= B.mean(axis=0, keepdims=True)
    B -= B.mean(axis=1, keepdims=True)
    if n <= dense_max:
        w, v = np.linalg.eigh(B)
        w, v = w[::-1][:rank], v[:, ::-1][:, :rank]
    else:                                              # randomised subspace iteration for the top components
        rng = np.random.default_rng(seed)
        Bt = torch.as_tensor(B, dtype=torch.float32)
        Q = torch.as_tensor(rng.normal(size=(n, rank + 16)), dtype=torch.float32)
        for _ in range(8):
            Q, _ = torch.linalg.qr(Bt @ Q)
        T = Q.T @ Bt @ Q
        w_, u = torch.linalg.eigh(T)
        w, v = w_.numpy()[::-1][:rank].astype(np.float64), (Q @ u).numpy()[:, ::-1][:, :rank].astype(np.float64)
    return v * np.sqrt(np.clip(w, 0.0, None))


def _coarse_targets(log_t: np.ndarray, w: np.ndarray, observed: np.ndarray, k: int):
    """Block-centroid targets for blocks of k consecutive beads: S_c(B1,B2) = mean s_ij - D_B1 - D_B2."""
    n = len(log_t)
    m = int(math.ceil(n / k))
    S = np.where(observed, np.exp(log_t.astype(np.float64)), np.nan)
    blocks = [slice(b * k, min(n, (b + 1) * k)) for b in range(m)]
    D = np.zeros(m)
    for b, sl in enumerate(blocks):
        sub = S[sl, sl]
        off = ~np.eye(sub.shape[0], dtype=bool)
        D[b] = 0.5 * np.nanmean(sub[off]) if off.any() and np.isfinite(sub[off]).any() else 0.0
    Sc = np.full((m, m), np.nan)
    Wc = np.zeros((m, m))
    for a in range(m):
        for b in range(a + 1, m):
            blk = S[blocks[a], blocks[b]]
            if np.isfinite(blk).any():
                Sc[a, b] = Sc[b, a] = max(np.nanmean(blk) - D[a] - D[b], 1e-6)
                Wc[a, b] = Wc[b, a] = w[blocks[a], blocks[b]].sum()
    obs_c = np.isfinite(Sc)
    np.fill_diagonal(obs_c, False)
    log_c = np.where(obs_c, np.log(np.where(obs_c, Sc, 1.0)), 0.0)
    if obs_c.any():
        Wc = np.where(obs_c, Wc, 0.0)
        Wc /= Wc[obs_c].mean()
    return log_c.astype(np.float32), Wc.astype(np.float32), obs_c


def _loss_grad(A: torch.Tensor, phi: torch.Tensor | None, T: torch.Tensor, W: torch.Tensor, P: float,
               block: int) -> tuple[float, torch.Tensor, torch.Tensor | None]:
    """L = (1/P) sum_{i<j} w_ij (log s_ij - t_ij)^2 and its exact gradients, one block of rows at a time.

    s_ij = |a_i - a_j|^2 + |C_j - C_i|, C = cumsum(exp(phi)). With e_ij = dL/ds_ij (symmetric):
        dL/dA     = 2 (diag(E 1) - E) A
        dL/dv_k   = sum_{i <= k < j} e_ij          (pairs whose path crosses bond k)
                  = sum_{i <= k} rowsum_i - 2 sum_{j <= k} c_j,   c_j = sum_{i < j} e_ij
        dL/dphi_k = v_k dL/dv_k
    (the second form of dL/dv: rows <= k, all columns, minus the k x k corner, which is twice the
    upper-triangle column sums). Memory beyond A, T and W is one block x N slab.
    """
    n = A.shape[0]
    dev, dtype = A.device, A.dtype
    use_v = phi is not None
    g = (A * A).sum(1)
    cum = torch.cat([torch.zeros(1, dtype=dtype, device=dev), torch.cumsum(torch.exp(phi), 0)]) if use_v else None
    gradA = torch.empty_like(A)
    rowsum = torch.empty(n, dtype=dtype, device=dev)
    colup = torch.zeros(n, dtype=dtype, device=dev) if use_v else None
    loss = torch.zeros((), dtype=torch.float64, device=dev)
    for a in range(0, n, block):
        b = min(n, a + block)
        Ab = A[a:b]
        s = (Ab @ A.T).mul_(-2.0).add_(g[a:b, None]).add_(g[None, :])
        if use_v:
            s.add_((cum[None, :] - cum[a:b, None]).abs_())
        s.clamp_(min=1e-9)
        res = torch.log(s).sub_(T[a:b])
        E = res * W[a:b]
        loss += (E * res).sum(dtype=torch.float64)
        E.div_(s).mul_(2.0 / P)                                # dL/ds_ij for this row block (both triangles)
        E.diagonal(offset=a).zero_()
        rowsum[a:b] = E.sum(1)
        gradA[a:b] = 2.0 * (rowsum[a:b, None] * Ab - E @ A)    # dL/da_i = sum_j e_ij 2 (a_i - a_j)
        if use_v:
            colup += torch.triu(E, diagonal=1 + a).sum(0)       # local row r is global row a + r
    val = float(loss.item()) / (2.0 * P)                        # each pair appears in both triangles
    if not use_v:
        return val, gradA, None
    R = (torch.cumsum(rowsum, 0) - 2.0 * torch.cumsum(colup, 0))[:-1]
    return val, gradA, R * torch.exp(phi)


def _adam_fit(log_t: np.ndarray, w: np.ndarray, A0: np.ndarray, v0: np.ndarray | None, iterations: int, lr: float,
              block: int, dev: torch.device, dtype: torch.dtype,
              progress: Callable[[int, int, dict[str, float]], None] | None = None,
              stage: str = "fit") -> tuple[np.ndarray, np.ndarray | None, dict]:
    """Adam on L = (1/P) sum_{i<j} w_ij (log s_ij - t_ij)^2 (see _loss_grad), keeping the best iterate."""
    T = torch.as_tensor(log_t, dtype=dtype, device=dev)
    W = torch.as_tensor(w, dtype=dtype, device=dev)
    P = float((w > 0).sum()) / 2.0
    A = torch.as_tensor(A0, dtype=dtype, device=dev).clone()
    use_v = v0 is not None
    phi = torch.as_tensor(np.log(np.maximum(v0, 1e-8)), dtype=dtype, device=dev).clone() if use_v else None
    mA, sA = torch.zeros_like(A), torch.zeros_like(A)
    if use_v:
        mP, sP = torch.zeros_like(phi), torch.zeros_like(phi)
    b1, b2, eps = 0.9, 0.999, 1e-8
    hist: dict[str, list[float]] = {"iteration": [], "loss": []}
    best, best_A, best_phi = float("inf"), A.clone(), (phi.clone() if use_v else None)
    for it in range(iterations + 1):
        val, gradA, gphi = _loss_grad(A, phi, T, W, P, block)
        if not math.isfinite(val):
            raise FloatingPointError(f"Population fit diverged at iteration {it} ({stage}).")
        if val < best:
            best, best_A = val, A.clone()
            if use_v:
                best_phi = phi.clone()
        if it % 10 == 0 or it == iterations:
            hist["iteration"].append(it)
            hist["loss"].append(val)
            if progress:
                progress(it, iterations, {"loss": val, "stage": stage})
        if it == iterations:
            break
        t = it + 1
        mA.mul_(b1).add_(gradA, alpha=1 - b1)
        sA.mul_(b2).addcmul_(gradA, gradA, value=1 - b2)
        A -= lr * (mA / (1 - b1 ** t)) / (torch.sqrt(sA / (1 - b2 ** t)) + eps)
        if use_v:
            mP.mul_(b1).add_(gphi, alpha=1 - b1)
            sP.mul_(b2).addcmul_(gphi, gphi, value=1 - b2)
            phi -= lr * (mP / (1 - b1 ** t)) / (torch.sqrt(sP / (1 - b2 ** t)) + eps)
    hist["best_loss"] = [best]
    v_out = torch.exp(best_phi).cpu().numpy().astype(np.float64) if use_v else None
    return best_A.cpu().numpy().astype(np.float64), v_out, hist


def _prolong(Ac: np.ndarray, n: int, k: int) -> np.ndarray:
    """Coarse rows (block centroids) -> bead rows, linearly interpolated between block centres."""
    m = Ac.shape[0]
    centres = np.arange(m) * k + (np.minimum(k, n - np.arange(m) * k) - 1) / 2.0
    beads = np.arange(n)
    return np.stack([np.interp(beads, centres, Ac[:, c]) for c in range(Ac.shape[1])], axis=1)


def fit_population(freq: np.ndarray, n_observed: np.ndarray | float | None = None, r_c_nm: float = 150.0,
                   cfg: PopulationConfig | None = None,
                   progress: Callable[[int, int, dict[str, float]], None] | None = None) -> ENS.EnsembleResult:
    """Fit a population of chains to a dense contact-frequency matrix of any size.

    Same inputs and outputs as ensemble.fit_ensemble (an EnsembleResult), plus `result.model`, the
    fitted GaussianChain, from which any per-pair statistic is exact.
    """
    cfg = cfg or PopulationConfig()
    t0 = time.time()
    dev = _device(cfg.device)
    dtype = torch.float32 if cfg.dtype == "float32" else torch.float64
    torch.manual_seed(cfg.seed)
    log_t, w, observed = targets_from_frequency(freq, n_observed, cfg.weighting)
    n = len(log_t)
    rank = int(min(n, cfg.rank_cap))
    use_local = bool(cfg.local_term and rank < n)      # with full rank the random walk would be redundant
    init = cfg.init if cfg.init != "auto" else ("coarse" if n > cfg.dense_init_max else "pca")
    history: dict[str, list] = {}
    t_init = time.time()
    if init == "coarse" and n > cfg.coarse_beads:
        k = int(math.ceil(n / cfg.coarse_beads))
        log_c, w_c, obs_c = _coarse_targets(log_t, w, observed, k)
        A_c0 = _psd_start(log_c, obs_c, min(len(log_c), rank), cfg.dense_init_max, cfg.seed)
        A_c, _, h_c = _adam_fit(log_c, w_c, A_c0, None, cfg.coarse_iterations, cfg.learning_rate, cfg.block_rows,
                                dev, dtype, progress, stage="coarse")
        history["coarse_loss"] = h_c["loss"]
        A0 = _prolong(A_c, n, k)
        if A0.shape[1] < rank:
            A0 = np.concatenate([A0, np.zeros((n, rank - A0.shape[1]))], axis=1)
        A0 += np.random.default_rng(cfg.seed).normal(scale=1e-3, size=A0.shape)
    else:
        A0 = _psd_start(log_t, observed, rank, cfg.dense_init_max, cfg.seed)
    v0 = None
    if use_local:
        adj = np.exp(np.diagonal(log_t, 1).astype(np.float64))
        d = A0[1:] - A0[:-1]
        v0 = np.maximum(adj - np.einsum("ik,ik->i", d, d), 0.05 * adj)
    t_init = time.time() - t_init
    A, v, h = _adam_fit(log_t, w, A0, v0, cfg.iterations, cfg.learning_rate, cfg.block_rows, dev, dtype, progress)
    history.update(h)
    t_fit = time.time() - t0
    model = GaussianChain(A.astype(np.float32), (v if v is not None else np.zeros(n - 1)), float(r_c_nm))
    return _finish(model, freq, observed, cfg, history, t0, t_fit, t_init, dev, init)


def _finish(model: GaussianChain, freq: np.ndarray, observed: np.ndarray, cfg: PopulationConfig, history: dict,
            t0: float, t_fit: float, t_init: float, dev: torch.device, init: str) -> ENS.EnsembleResult:
    n = model.n
    r_c = model.r_c_nm
    S = model.dense_variance()                                       # units r_c^2, float32
    sigma = np.sqrt(S, dtype=np.float32)
    median = (ENS.MAXWELL_MEDIAN * r_c) * sigma
    p_model = ENS.contact_probability_from_sigma(sigma).astype(np.float32)
    np.fill_diagonal(p_model, 1.0)
    iu = np.triu_indices(n, 1)
    f_in = np.asarray(freq, dtype=np.float64)
    f_sym = np.where(np.isfinite(f_in), f_in, f_in.T)
    obs_iu = observed[iu]
    rng = np.random.default_rng(cfg.seed)
    pick = np.flatnonzero(obs_iu)
    if pick.size > 2_000_000:
        pick = rng.choice(pick, 2_000_000, replace=False)
    contact_fit = ENS._spearman(p_model[iu][pick], f_sym[iu][pick])

    t1 = time.time()
    traj, sampled, spread, rep, eig_exact = _langevin(model, cfg, rng, median)
    return ENS.EnsembleResult(
        median_distance_nm=median, contact_probability=p_model,
        covariance_nm2=None if n > 3000 else (model.covariance() * r_c ** 2),
        couplings=None, trajectories_nm=traj, sampled_median_distance_nm=sampled, representative_nm=rep,
        contact_fit=contact_fit, spread_cv=spread, history=history, seconds=time.time() - t0,
        config=asdict(cfg) | {"r_c_nm": r_c, "device_used": str(dev), "fit_seconds": t_fit, "init_seconds": t_init,
                              "init": init, "rank": model.rank, "local_term": bool(model.v.any()),
                              "sampling_seconds": time.time() - t1, "exact_ou_modes": eig_exact, "n_beads": n,
                              "model": "population_v4"},
        model=model)


def _langevin(model: GaussianChain, cfg: PopulationConfig, rng: np.random.Generator, median: np.ndarray):
    """Exact Ornstein-Uhlenbeck trajectories of the fitted ensemble, chunked by replica.

    Up to cfg.max_eig_beads the physical process dx = -P x dt + sqrt(2) dW (P = covariance^+) is
    propagated exactly mode by mode, as in v3.3. Frames kept are capped by max_trajectory_floats.
    Above the cap, members are exact independent draws and the trajectory has one frame."""
    n = model.n
    r_c = model.r_c_nm
    frames = int(max(1, min(cfg.frames, cfg.max_trajectory_floats // max(1, cfg.replicas * n * 3))))
    eig_exact = n <= cfg.max_eig_beads
    if eig_exact:
        cov = torch.as_tensor(model.covariance(), dtype=torch.float64 if n <= 3000 else torch.float32)
        lam_t, modes_t = torch.linalg.eigh(cov)
        lam, modes = lam_t.double().numpy(), modes_t.double().numpy()
        keep = lam > lam.max() * 1e-10
        lam, modes = lam[keep], modes[:, keep]
        m = len(lam)
        # the same simulated time as v3.3 (cfg.frames x frame_interval slowest relaxation times); when fewer
        # frames are kept, each kept frame is one exact OU step over the longer interval
        dt = cfg.frame_interval * lam.max() * max(1, cfg.frames // frames)
        decay = np.exp(-dt / lam)
        kick = np.sqrt(lam * (1 - decay ** 2))
        out = np.empty((frames, cfg.replicas, n, 3), dtype=np.float32)
        chunk = max(1, int(4_000_000 // max(1, m * 3)))
        for a in range(0, cfg.replicas, chunk):
            b = min(cfg.replicas, a + chunk)
            y = rng.normal(size=(m, b - a, 3)) * np.sqrt(lam)[:, None, None]       # equilibrium start
            for f in range(frames):
                y = y * decay[:, None, None] + rng.normal(size=y.shape) * kick[:, None, None]
                x = (modes @ y.reshape(m, -1)).reshape(n, b - a, 3)                  # one GEMM per frame
                out[f, a:b] = (x.transpose(1, 0, 2) * r_c).astype(np.float32)
    else:
        out = model.sample(cfg.replicas, rng).astype(np.float32)[None]
    # finite-sample check and spread on a pair subset (all pairs when small); <= 5e7 floats (~200 MB, as v3.3)
    iu = np.triu_indices(n, 1)
    sel = np.arange(len(iu[0]))
    if sel.size > 200_000:
        sel = rng.choice(sel, 200_000, replace=False)
    ii, jj = iu[0][sel], iu[1][sel]
    flat = out.reshape(-1, n, 3)
    budget = max(cfg.replicas, int(5e7 // max(len(ii), 1)))
    take = flat[:: max(1, int(np.ceil(len(flat) / budget)))]
    d = np.linalg.norm(take[:, ii] - take[:, jj], axis=-1)
    sampled = np.zeros((n, n), dtype=np.float32)
    sampled[ii, jj] = np.median(d, axis=0)
    sampled = sampled + sampled.T
    spread = float(np.mean(d.std(0) / np.maximum(d.mean(0), 1e-12)))
    last = np.linalg.norm(out[-1][:, ii] - out[-1][:, jj], axis=-1)
    rep_i = int(np.argmin(np.abs(np.log(np.maximum(last, 1e-9)) - np.log(np.maximum(median[ii, jj], 1e-9))).mean(1)))
    return out, sampled, spread, out[-1][rep_i].astype(np.float64), eig_exact


def fit_population_from_counts(ci: np.ndarray, cj: np.ndarray, cm: np.ndarray, n: int, valid: np.ndarray | None = None,
                               b0_nm: float = 50.0, p_adjacent: float = 0.5, cfg: PopulationConfig | None = None,
                               progress: Callable[[int, int, dict[str, float]], None] | None = None) -> ENS.EnsembleResult:
    """Population model from a sequencing contact list, any size: the same two explicit assumptions as
    ensemble.fit_from_counts (adjacent contact probability p_adjacent; adjacent median distance b0)."""
    ci = np.asarray(ci, dtype=np.int64)
    cj = np.asarray(cj, dtype=np.int64)
    counts = np.zeros((n, n), dtype=np.float64)
    np.add.at(counts, (ci, cj), np.asarray(cm, dtype=np.float64))
    counts = counts + counts.T
    np.fill_diagonal(counts, 0.0)
    adj = np.diag(counts, 1)
    adj = adj[adj > 0]
    if adj.size < 3:
        raise ValueError("Too few adjacent-bead contacts to set the probability scale for a population model.")
    p = ENS.counts_to_probability(counts, p_adjacent)
    n_eff = float(np.median(adj)) / p_adjacent
    if valid is not None:
        bad = ~np.asarray(valid, bool)
        p[bad, :] = np.nan
        p[:, bad] = np.nan
    r_c = b0_nm / float(ENS.gaussian_median_distance(p_adjacent, 1.0))
    res = fit_population(p, max(n_eff, 1.0), r_c_nm=r_c, cfg=cfg, progress=progress)
    res.config.update({"input": "sequencing counts", "p_adjacent_assumed": p_adjacent, "n_effective": n_eff,
                       "anchor_b0_nm": b0_nm})
    return res


def fit_population_from_medians(median_nm: np.ndarray, p_adjacent: float = 0.5, n_observed: float = 1000.0,
                                cfg: PopulationConfig | None = None,
                                progress: Callable[[int, int, dict[str, float]], None] | None = None) -> ENS.EnsembleResult:
    """Population model whose pair medians follow a given distance map, e.g. a prediction made without
    contact data (chronocell.predict). Each median d_ij gives a per-axis spread sigma_ij = d_ij / MAXWELL_MEDIAN
    and a contact frequency at the radius r_c where adjacent beads touch with probability p_adjacent; the
    frequencies are then fitted like measured ones (v3.3 up to V33_MAX_BEADS beads, v4 above), which projects
    the map onto a valid 3D Gaussian ensemble. n_observed only sets the fit weights (no sampling noise here)."""
    d = np.asarray(median_nm, dtype=np.float64)
    n = len(d)
    if n < 3 or d.shape != (n, n):
        raise ValueError("Need a square distance map of at least 3 beads.")
    adj = np.diag(d, 1)
    adj = adj[np.isfinite(adj) & (adj > 0)]
    if adj.size == 0:
        raise ValueError("No adjacent distances to anchor the contact radius.")
    r_c = float(np.median(adj)) / float(ENS.gaussian_median_distance(p_adjacent, 1.0))
    sig = np.where(np.isfinite(d) & (d > 0), d / MAXWELL_MEDIAN, np.nan)
    freq = np.where(np.isfinite(sig), ENS.contact_probability_from_sigma(np.nan_to_num(sig, nan=1.0) / r_c), np.nan)
    np.fill_diagonal(freq, np.nan)
    if n > V33_MAX_BEADS:
        res = fit_population(freq, n_observed, r_c_nm=r_c, cfg=cfg or config_for(n), progress=progress)
    else:
        res = ENS.fit_ensemble(freq, n_observed, r_c_nm=r_c, progress=progress)
    res.config.update({"input": "distance map (no contact data)", "p_adjacent_assumed": p_adjacent,
                       "n_observed_weight": n_observed})
    return res


# ======================================================================================
# Uncertainty read-outs (any EnsembleResult: v3.3 dense covariance or v4 model)
# ======================================================================================
def pair_sigma_nm(res: ENS.EnsembleResult, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Per-axis sd (nm) of the pair vector, exact, for v3.3 and v4 results alike."""
    i, j = np.asarray(i, np.int64), np.asarray(j, np.int64)
    model = getattr(res, "model", None)
    if model is not None:
        return model.sigma_nm(i, j)
    return np.asarray(res.median_distance_nm, dtype=np.float64)[i, j] / MAXWELL_MEDIAN


def pair_summary(res: ENS.EnsembleResult, i: int, j: int, level: float = 0.8, recalibrated: bool = False,
                 kind: str = "imaging") -> PairSummary:
    """Mean, sd, median and central `level` interval (nm) of one pair's distance across the population
    (interval recalibrated on practice data for input `kind` if `recalibrated`)."""
    sig = pair_sigma_nm(res, np.array([i]), np.array([j]))
    s = summarise_sigma(sig, float(res.config.get("r_c_nm", 150.0)), level, recalibrated, kind)
    return PairSummary(*(float(s[k][0]) for k in ("mean", "sd", "median", "lower", "upper")), level,
                       float(s["contact_probability"][0]))


def rmsf_nm(members_nm: np.ndarray, reference_nm: np.ndarray | None = None) -> np.ndarray:
    """Per-bead root-mean-square fluctuation (nm) of members after optimal superposition (proper
    rotations) onto a reference (default: the first member). This is ensemble spread, i.e. how
    consistent the population is about a bead's position, not accuracy."""
    x = np.asarray(members_nm, dtype=np.float64)
    ref = x[0] if reference_nm is None else np.asarray(reference_nm, dtype=np.float64)
    ref_c = ref - ref.mean(0)
    al = np.empty_like(x)
    for k in range(len(x)):
        q = x[k] - x[k].mean(0)
        u, _, vt = np.linalg.svd(q.T @ ref_c)
        d = np.sign(np.linalg.det(u @ vt))
        rot = u @ np.diag([1.0, 1.0, d]) @ vt
        al[k] = q @ rot
    mean = al.mean(0)
    return np.sqrt(((al - mean) ** 2).sum(-1).mean(0))


def bead_reliability(res: ENS.EnsembleResult, freq: np.ndarray, n_observed: np.ndarray | float | None,
                     tau: float = 0.25) -> np.ndarray:
    """Per-bead reliability in [0, 1] from the input alone (no held-out data):
        coverage_i  share of the bead's pairs that were observed (not NaN / unassembled)
        misfit_i    weighted mean |log s_model - log s_target| over its observed pairs
        reliability_i = coverage_i * exp(-misfit_i / tau)
    tau = 0.25 is a fixed choice, not tuned: a weighted mean misfit of 0.25 in log variance lowers the
    score by a factor e. Whether the score predicts held-out error is measured in validation/ (Pillar 2)."""
    log_t, w, observed = targets_from_frequency(freq, n_observed)
    n = len(log_t)
    iu = np.triu_indices(n, 1)
    sig = pair_sigma_nm(res, iu[0], iu[1]) / float(res.config.get("r_c_nm", 150.0))
    log_m = np.zeros((n, n))
    log_m[iu] = 2.0 * np.log(np.maximum(sig, 1e-12))
    log_m = log_m + log_m.T
    err = np.abs(log_m - log_t.astype(np.float64)) * w
    cover = observed.sum(1) / max(n - 1, 1)
    wsum = w.sum(1)
    misfit = np.where(wsum > 0, err.sum(1) / np.maximum(wsum, 1e-12), np.inf)
    return np.clip(cover * np.exp(-misfit / tau), 0.0, 1.0)


def excluded_volume_members(members_nm: np.ndarray, b0_nm: float, d_min_factor: float = 0.8,
                            iters: int = 40) -> tuple[np.ndarray, dict]:
    """Relax sampled members to bond length b0 and no overlaps (synthetic.relax), and measure how much
    that moves them away from the Gaussian ensemble (median relative change of pair distances)."""
    from . import physics, synthetic
    out = np.empty_like(members_nm, dtype=np.float64)
    before = after = 0
    for k, x in enumerate(np.asarray(members_nm, dtype=np.float64)):
        before += physics.loss_steric(x, b0_nm, d_min_factor * b0_nm).overlaps
        out[k] = synthetic.relax(x, b0_nm, d_min_factor * b0_nm, iters=iters)
        after += physics.loss_steric(out[k], b0_nm, d_min_factor * b0_nm).overlaps
    n = members_nm.shape[1]
    rng = np.random.default_rng(0)
    ii = rng.integers(0, n, 20_000)
    jj = rng.integers(0, n, 20_000)
    ok = ii != jj
    d0 = np.linalg.norm(members_nm[:, ii[ok]] - members_nm[:, jj[ok]], axis=-1)
    d1 = np.linalg.norm(out[:, ii[ok]] - out[:, jj[ok]], axis=-1)
    change = float(np.median(np.abs(d1 - d0) / np.maximum(d0, 1e-9)))
    return out, {"overlaps_before": int(before), "overlaps_after": int(after),
                 "median_relative_pair_distance_change": change,
                 "median_bond_over_b0_after": float(np.median(physics.bond_lengths(out[0])) / b0_nm)}
