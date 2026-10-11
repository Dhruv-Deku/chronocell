"""
Lattice folding of a short chain to match its contacts (the textbook "protein folding on a lattice" problem
applied to chromatin segments; Perdomo-Ortiz et al., Sci Rep 2012; Robert et al., npj Quantum Inf 2021).

Each bond after the first is a direction on the square (2D) or cubic (3D) lattice: 2 bits per bond in 2D
(+x, +y, -x, -y) and 3 bits in 3D (+-x, +-y, +-z; the two unused codes are penalised). The first bond is fixed
along +x (removes rotations). Energy of a fold:
    E = - sum_{|i-j|>=3} w_ij [beads i and j are lattice neighbours] + lambda (overlapping beads) + lambda (unused codes)
with w_ij the positive part of the observed/expected contact enrichment of segments i and j, scaled to max 1.
This cost has terms in up to all bits of a pair of beads, so it is a higher-order (HUBO) diagonal Hamiltonian,
not a QUBO: QAOA handles it as the 2^n energies (and qaoa_counts gives the multi-qubit gates it would compile to).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

DIRS = {2: np.array([[1, 0], [0, 1], [-1, 0], [0, -1]], dtype=np.int16),
        3: np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1], [0, 0, 0], [0, 0, 0]],
                    dtype=np.int16)}
BITS = {2: 2, 3: 3}
VALID = {2: 4, 3: 6}


def n_qubits(beads: int, dim: int) -> int:
    return BITS[dim] * (beads - 2)


def weights_from_contacts(M: np.ndarray, beads: int) -> np.ndarray:
    """Coarse-grain a dense contact matrix to `beads` segments and return w_ij = max(0, O/E - 1) scaled to max 1
    (expected = mean contacts at the same segment separation)."""
    M = np.asarray(M, float)
    n = len(M)
    edges = np.linspace(0, n, beads + 1).astype(int)
    C = np.zeros((beads, beads))
    for i in range(beads):
        for j in range(beads):
            C[i, j] = M[edges[i]:edges[i + 1], edges[j]:edges[j + 1]].mean()
    E = np.zeros(beads)
    for d in range(beads):
        E[d] = np.mean([C[i, i + d] for i in range(beads - d)])
    W = np.zeros((beads, beads))
    for i in range(beads):
        for j in range(i + 3, beads):
            if E[j - i] > 0:
                W[i, j] = W[j, i] = max(0.0, C[i, j] / E[j - i] - 1.0)
    mx = W.max()
    return W / mx if mx > 0 else W


def positions(index: np.ndarray, beads: int, dim: int) -> tuple[np.ndarray, np.ndarray]:
    """Bead coordinates (m, beads, dim) of fold indices, and the number of unused codes per fold."""
    idx = np.asarray(index, np.int64).reshape(-1)
    b = BITS[dim]
    m = len(idx)
    pos = np.zeros((m, beads, dim), dtype=np.int16)
    pos[:, 1, 0] = 1
    bad = np.zeros(m, dtype=np.int16)
    for t in range(beads - 2):
        code = (idx >> (b * t)) & ((1 << b) - 1)
        pos[:, t + 2] = pos[:, t + 1] + DIRS[dim][code]
        bad += (code >= VALID[dim]).astype(np.int16)
    return pos, bad


def energies(W: np.ndarray, dim: int = 2, overlap_penalty: float | None = None, chunk: int = 1 << 18) -> np.ndarray:
    """Energy of every fold (2^n_qubits)."""
    W = np.asarray(W, float)
    beads = len(W)
    q = n_qubits(beads, dim)
    lam = overlap_penalty if overlap_penalty is not None else 1.0 + 2.0 * float(W.sum(axis=1).max(initial=0.0))
    out = np.empty(2 ** q)
    pairs = [(i, j) for i in range(beads) for j in range(i + 2, beads)]
    for s in range(0, 2 ** q, chunk):
        idx = np.arange(s, min(2 ** q, s + chunk), dtype=np.int64)
        pos, bad = positions(idx, beads, dim)
        e = lam * bad.astype(float)
        for i, j in pairs:
            d = np.abs(pos[:, i].astype(np.int32) - pos[:, j]).sum(axis=1)
            e += lam * (d == 0)
            if j - i >= 3 and W[i, j] > 0:
                e -= W[i, j] * (d == 1)
        out[s:s + len(idx)] = e
    return out


@dataclass
class Fold:
    index: int
    coords: np.ndarray            # (beads, dim) lattice coordinates
    energy: float
    contacts: list                # (i, j) lattice neighbours with |i-j| >= 3
    valid: bool                   # no overlap, no unused code


def fold_of(index: int, W: np.ndarray, dim: int, energy: float) -> Fold:
    beads = len(W)
    pos, bad = positions(np.array([index]), beads, dim)
    p = pos[0].astype(int)
    cont = [(i, j) for i in range(beads) for j in range(i + 3, beads) if np.abs(p[i] - p[j]).sum() == 1]
    overl = any(np.abs(p[i] - p[j]).sum() == 0 for i in range(beads) for j in range(i + 1, beads))
    return Fold(int(index), p, float(energy), cont, bool(bad[0] == 0 and not overl))


def contact_agreement(fold: Fold, W: np.ndarray) -> dict:
    """How much of the enrichment the fold realises: sum of w over its contacts / best achievable by any k pairs."""
    w = np.asarray(W, float)
    got = sum(w[i, j] for i, j in fold.contacts)
    k = len(fold.contacts)
    top = np.sort(w[np.triu_indices(len(w), 3)])[::-1]
    return {"contacts": k, "weight_captured": float(got),
            "fraction_of_top_k": float(got / top[:k].sum()) if k and top[:k].sum() > 0 else float("nan")}
