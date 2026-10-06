"""
Continuous-time quantum walk on a contact graph (Farhi & Gutmann, PRA 1998), against the classical continuous-time
random walk on the same graph.

Graph: bins are nodes; edge weights are contacts (raw or observed/expected), diagonal removed. With D the weighted
degrees, the quantum walk evolves |psi(t)> = exp(-i L t)|start> with the symmetric normalised Laplacian
L = I - D^-1/2 A D^-1/2, and the classical walk p(t) = exp(-L_rw t) p(0) with L_rw = I - A D^-1 (probability
conserving). Both start on one bin. The quantum walk spreads ballistically (distance ~ t) where the classical one
diffuses (distance ~ sqrt t), and interference makes it sensitive to the graph's structure: after a change to the
contacts (a variant, cohesin loss) the two walks show where the signal now travels.

On a quantum computer the position needs ceil(log2 n) qubits and exp(-i L t) a Hamiltonian-simulation circuit;
here the evolution is exact (eigendecomposition).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def adjacency(ci: np.ndarray, cj: np.ndarray, cm: np.ndarray, n: int, observed_expected: bool = True) -> np.ndarray:
    A = np.zeros((n, n))
    ci, cj, cm = np.asarray(ci, int), np.asarray(cj, int), np.asarray(cm, float)
    m = (ci != cj) & (ci >= 0) & (cj >= 0) & (ci < n) & (cj < n)
    np.add.at(A, (ci[m], cj[m]), cm[m])
    A = np.maximum(A, A.T)
    if observed_expected:
        s = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
        exp = np.array([A[s == d].mean() if np.any(s == d) else 0 for d in range(n)])
        with np.errstate(divide="ignore", invalid="ignore"):
            A = np.where(exp[s] > 0, A / exp[s], 0.0)
    np.fill_diagonal(A, 0.0)
    return A


@dataclass
class Walks:
    times: np.ndarray
    quantum: np.ndarray        # (T, n) probabilities
    classical: np.ndarray      # (T, n)
    start: int
    qubits: int

    def spread(self, which: str = "quantum") -> np.ndarray:
        P = self.quantum if which == "quantum" else self.classical
        x = np.arange(P.shape[1])
        mu = P @ x
        return np.sqrt(np.maximum(P @ (x ** 2) - mu ** 2, 0))

    def time_average(self, which: str = "quantum") -> np.ndarray:
        P = self.quantum if which == "quantum" else self.classical
        return P.mean(axis=0)

    def arrival(self, threshold: float | None = None, which: str = "quantum") -> np.ndarray:
        """First time each bin's probability reaches `threshold` (default 1/n); NaN if never."""
        P = self.quantum if which == "quantum" else self.classical
        thr = 1.0 / P.shape[1] if threshold is None else threshold
        hit = P >= thr
        first = np.where(hit.any(axis=0), hit.argmax(axis=0), -1)
        return np.where(first >= 0, self.times[np.maximum(first, 0)], np.nan)


def walks(A: np.ndarray, start: int, t_max: float = 20.0, steps: int = 80) -> Walks:
    A = np.asarray(A, float)
    n = len(A)
    deg = A.sum(axis=1)
    deg = np.where(deg > 0, deg, 1.0)
    dm = 1 / np.sqrt(deg)
    Ls = np.eye(n) - dm[:, None] * A * dm[None, :]
    w, V = np.linalg.eigh(Ls)
    times = np.linspace(0, t_max, steps)
    e0 = np.zeros(n)
    e0[int(start)] = 1.0
    c0 = V.T @ e0
    Q = np.abs((V * np.exp(-1j * np.outer(times, w))[:, None, :]) @ c0) ** 2      # (T, n)
    # classical: exp(-L_rw t) = D^1/2 exp(-Ls t) D^-1/2
    Pc = (V * np.exp(-np.outer(times, w))[:, None, :]) @ (V.T @ (dm * e0))
    Pc = np.maximum(Pc * np.sqrt(deg)[None, :], 0)
    Pc /= Pc.sum(axis=1, keepdims=True)
    return Walks(times, Q, Pc, int(start), int(math.ceil(math.log2(max(2, n)))))


def compare(before: Walks, after: Walks) -> dict:
    """How the change moves where each walk spends its time (total-variation distance of time averages)."""
    out = {}
    for which in ("quantum", "classical"):
        a, b = before.time_average(which), after.time_average(which)
        out[which] = {"tv_distance": float(0.5 * np.abs(a - b).sum()),
                      "most_changed_bins": np.argsort(-np.abs(a - b))[:10].tolist()}
    return out
