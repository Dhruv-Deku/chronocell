"""
QUBO / Ising problems and their classical solvers.

A QUBO over bits x_i in {0, 1}:  E(x) = sum_i a_i x_i + sum_{i<j} Q_ij x_i x_j + c.
Qubit i holds x_i; |0> is x = 0. With x_i = (1 - z_i) / 2 and z_i the eigenvalue of Z_i (+1 for |0>), the same
energy is the Ising form  E = sum_i h_i z_i + sum_{i<j} J_ij z_i z_j + const, which QAOA and annealers use.

Solvers (all on the same QUBO):
exact        enumerates all 2^n energies (n <= 26), returns every minimum.
anneal       simulated annealing (Metropolis single-bit flips, geometric temperature schedule), many reads at once.
sqa          simulated quantum annealing: path-integral Monte Carlo of the transverse-field Ising model with P
             Trotter replicas (Santoro et al., Science 2002; Martonak et al., PRB 2002). A CLASSICAL algorithm that
             imitates a quantum annealer, labelled so.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np

MAX_EXACT = 26


@dataclass
class QUBO:
    lin: np.ndarray                       # a_i
    quad: np.ndarray                      # Q_ij, upper triangle used (i < j), symmetric copies ignored
    offset: float = 0.0
    labels: list = field(default_factory=list)
    name: str = "QUBO"

    def __post_init__(self):
        self.lin = np.asarray(self.lin, float)
        q = np.triu(np.asarray(self.quad, float), 1)
        self.quad = q
        if not self.labels:
            self.labels = [f"x{i}" for i in range(self.n)]

    @property
    def n(self) -> int:
        return len(self.lin)

    @property
    def sym(self) -> np.ndarray:
        """Symmetric coupling matrix with zero diagonal (Q_ij on both sides)."""
        return self.quad + self.quad.T

    def energy(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, float)
        return x @ self.lin + np.einsum("...i,ij,...j->...", x, self.quad, x) + self.offset

    def energies(self) -> np.ndarray:
        """All 2^n energies; index k <-> bits x_i = (k >> i) & 1 (qubit i)."""
        n = self.n
        if n > MAX_EXACT:
            raise ValueError(f"{n} variables: too many to enumerate (limit {MAX_EXACT}).")
        idx = np.arange(2 ** n, dtype=np.int64)
        e = np.full(2 ** n, float(self.offset))
        bits = [((idx >> i) & 1).astype(np.float64) for i in range(n)]
        for i in range(n):
            if self.lin[i] != 0:
                e += self.lin[i] * bits[i]
        ii, jj = np.nonzero(self.quad)
        for i, j in zip(ii, jj):
            e += self.quad[i, j] * (bits[i] * bits[j])
        return e

    def ising(self) -> tuple[np.ndarray, dict, float]:
        """(h, J {(i, j): J_ij}, const) with E = sum h_i z_i + sum J_ij z_i z_j + const."""
        Qs = self.sym
        h = -self.lin / 2 - Qs.sum(axis=1) / 4
        J = {(int(i), int(j)): float(self.quad[i, j] / 4) for i, j in zip(*np.nonzero(self.quad))}
        const = float(self.offset + self.lin.sum() / 2 + self.quad.sum() / 4)
        return h, J, const

    def density(self) -> float:
        n = self.n
        return float(np.count_nonzero(self.quad) / max(1, n * (n - 1) / 2))


def to_bits(index: int, n: int) -> np.ndarray:
    return ((int(index) >> np.arange(n)) & 1).astype(np.int8)


def to_index(bits: np.ndarray) -> int:
    b = np.asarray(bits, np.int64)
    return int((b << np.arange(len(b))).sum())


@dataclass
class Solution:
    method: str
    bits: np.ndarray
    energy: float
    seconds: float
    detail: dict = field(default_factory=dict)

    @property
    def index(self) -> int:
        return to_index(self.bits)


def exact(q: QUBO, energies: np.ndarray | None = None) -> Solution:
    t0 = time.perf_counter()
    e = q.energies() if energies is None else energies
    emin = float(e.min())
    opt = np.flatnonzero(np.isclose(e, emin, rtol=0, atol=1e-9 * max(1.0, abs(emin))))
    return Solution("exact (all 2^n states)", to_bits(int(opt[0]), q.n), emin, time.perf_counter() - t0,
                    {"optima": opt.tolist()[:64], "n_optima": int(len(opt)), "states": int(len(e))})


def anneal(q: QUBO, reads: int = 64, sweeps: int = 600, t_hot: float | None = None, t_cold: float | None = None,
           seed: int = 0) -> Solution:
    """Simulated annealing (vectorised over reads)."""
    t0 = time.perf_counter()
    rng = np.random.default_rng(seed)
    n = q.n
    Qs = q.sym
    scale = float(np.abs(q.lin).max(initial=0) + np.abs(Qs).sum(axis=1).max(initial=0)) or 1.0
    t_hot = t_hot or 0.5 * scale
    t_cold = t_cold or 1e-3 * scale
    temps = t_hot * (t_cold / t_hot) ** (np.arange(sweeps) / max(1, sweeps - 1))
    x = rng.integers(0, 2, (reads, n)).astype(float)
    field_ = x @ Qs + q.lin                    # local field: dE for 0->1 is field_i
    for T in temps:
        for i in rng.permutation(n):
            dE = (1 - 2 * x[:, i]) * field_[:, i]
            acc = (dE <= 0) | (rng.random(reads) < np.exp(-np.maximum(dE, 0) / T))
            if acc.any():
                delta = (1 - 2 * x[acc, i])
                x[acc, i] += delta
                field_[acc] += delta[:, None] * Qs[i][None, :]
    e = q.energy(x)
    k = int(np.argmin(e))
    return Solution("simulated annealing (classical)", x[k].astype(np.int8), float(e[k]), time.perf_counter() - t0,
                    {"reads": reads, "sweeps": sweeps, "energies": e})


def sqa(q: QUBO, reads: int = 16, sweeps: int = 400, trotter: int = 16, temperature: float | None = None,
        gamma0: float | None = None, gamma1: float = 1e-3, seed: int = 0) -> Solution:
    """Simulated quantum annealing (path-integral Monte Carlo) on the Ising form; transverse field Gamma ramps from
    gamma0 to gamma1. Classical algorithm."""
    t0 = time.perf_counter()
    rng = np.random.default_rng(seed)
    h, J, const = q.ising()
    n, P = q.n, int(trotter)
    Jm = np.zeros((n, n))
    for (i, j), v in J.items():
        Jm[i, j] = Jm[j, i] = v
    scale = float(np.abs(h).max(initial=0) + np.abs(Jm).sum(axis=1).max(initial=0)) or 1.0
    T = temperature or 0.05 * scale
    gamma0 = gamma0 or 3.0 * scale
    gammas = gamma0 * (gamma1 / gamma0) ** (np.arange(sweeps) / max(1, sweeps - 1))
    s = rng.choice([-1.0, 1.0], size=(reads, P, n))
    PT = P * T
    for G in gammas:
        jperp = -0.5 * T * math.log(math.tanh(G / PT))           # ferromagnetic coupling between replicas (>0)
        for parity in (0, 1):
            sl = np.arange(parity, P, 2)
            for i in rng.permutation(n):
                loc = h[i] + s[:, sl, :] @ Jm[i]                  # classical local field in those replicas
                up = s[:, (sl + 1) % P, i]
                dn = s[:, (sl - 1) % P, i]
                si = s[:, sl, i]
                # H_eff = sum_k E_cl(s_k) / P - jperp sum_k s_k s_k+1; cost of flipping s_ki:
                dE = -2 * si * loc / P + 2 * jperp * si * (up + dn)
                acc = (dE <= 0) | (rng.random(dE.shape) < np.exp(-np.maximum(dE, 0) / T))
                s[:, sl, i] = np.where(acc, -si, si)
    x = ((1 - s) / 2).reshape(-1, n)
    e = q.energy(x)
    k = int(np.argmin(e))
    return Solution("simulated quantum annealing (path-integral Monte Carlo, classical)", x[k].astype(np.int8),
                    float(e[k]), time.perf_counter() - t0, {"reads": reads, "sweeps": sweeps, "trotter": P})


def bandwidth(q: QUBO) -> int:
    """Largest |i - j| with a non-zero coupling (0 for no couplings)."""
    ii, jj = np.nonzero(q.quad)
    return int((jj - ii).max()) if len(ii) else 0


def banded_exact(q: QUBO, max_band: int = 20) -> Solution:
    """Exact minimum of a QUBO whose couplings only join variables at most w apart, by dynamic programming over the
    last w bits: O(n 2^w) time instead of 2^n. Many chain-like ChronoCell problems (domain walls along DNA) are
    banded, so this classical algorithm solves them at any length."""
    t0 = time.perf_counter()
    n, w = q.n, bandwidth(q)
    if w > max_band:
        raise ValueError(f"bandwidth {w} is too wide for the dynamic programme (limit {max_band})")
    if w == 0:
        bits = (q.lin < 0).astype(np.int8)
        return Solution("exact (dynamic programming)", bits, float(q.energy(bits)), time.perf_counter() - t0,
                        {"bandwidth": 0})
    S = 1 << w
    mask = S - 1
    states = np.arange(S, dtype=np.int64)
    sbits = ((states[:, None] >> np.arange(w)) & 1).astype(float)           # bit d-1 = x_{i+1-d}
    cost = np.full(S, np.inf)
    cost[0] = 0.0
    back = np.zeros((n, S), dtype=np.int8)
    for i in range(n):
        c = np.array([q.quad[i - d, i] if i - d >= 0 else 0.0 for d in range(1, w + 1)])
        add1 = q.lin[i] + sbits @ c                                           # cost of x_i = 1 given previous bits
        new = np.full(S, np.inf)
        for x in (0, 1):
            nxt = states[(states & 1) == x]                                   # s' with newest bit x
            pa = nxt >> 1
            pb = pa | (1 << (w - 1))
            ca = cost[pa] + (add1[pa] if x else 0.0)
            cb = cost[pb] + (add1[pb] if x else 0.0)
            pick = cb < ca
            new[nxt] = np.where(pick, cb, ca)
            back[i, nxt] = pick.astype(np.int8)
        cost = new
    s = int(np.argmin(cost))
    best = float(cost[s]) + q.offset
    bits = np.zeros(n, dtype=np.int8)
    for i in range(n - 1, -1, -1):
        bits[i] = s & 1
        s = (s >> 1) | (int(back[i, s]) << (w - 1)) if back[i, s] else (s >> 1)
        s &= mask
    return Solution("exact (dynamic programming over the band)", bits, best, time.perf_counter() - t0,
                    {"bandwidth": w, "states": int(S)})
