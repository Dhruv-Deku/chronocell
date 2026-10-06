"""
Quantum chemistry of small, drug-relevant molecules (H to Ne), written from scratch, for active-space VQE.

Basis: STO-3G (Hehre, Stewart & Pople, J Chem Phys 1969): every 1s and 2sp shell is a contraction of 3 Gaussians
fitted to a Slater orbital, scaled by the standard molecular exponents zeta. Cartesian s and p functions.
Integrals: McMurchie-Davidson (Hermite expansion; Boys function), vectorised over the primitives of each
contracted pair / quartet; 8-fold symmetry for the two-electron integrals.
Restricted Hartree-Fock with DIIS. Active space: frozen doubly occupied core folded into a constant and an
effective one-electron operator; active orbitals around the HOMO / LUMO. Qubits: Jordan-Wigner, spin orbital
2 x orbital + spin = qubit (as chem.py), sparse matrices. FCI in the active space = exact diagonalisation in the
right electron-number and spin sector (the reference VQE is scored against). VQE: disentangled UCCSD (spin-
conserving singles and doubles, each a unitary exp(theta (T - T^dagger))), BFGS.

Everything here runs classically; the "quantum" part is that the VQE state is what a quantum computer would
prepare (circuit export via chem.uccsd-style Pauli exponentials is available for small active spaces).
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.special import gamma as _gamma, gammainc

BOHR = 0.529177210903                                   # angstrom per bohr
Z = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10}
ZETA_1S = {"H": 1.24, "He": 1.69, "Li": 2.69, "Be": 3.68, "B": 4.68, "C": 5.67, "N": 6.67, "O": 7.66, "F": 8.65,
           "Ne": 9.64}
ZETA_2SP = {"Li": 0.80, "Be": 1.15, "B": 1.50, "C": 1.72, "N": 1.95, "O": 2.25, "F": 2.55, "Ne": 2.88}
A1S = np.array([2.22766, 0.405771, 0.109818])
D1S = np.array([0.154329, 0.535328, 0.444635])
A2SP = np.array([0.994203, 0.231031, 0.0751386])
D2S = np.array([-0.0999672, 0.399513, 0.700115])
D2P = np.array([0.155916, 0.607684, 0.391957])
CHEMICAL_ACCURACY = 1.6e-3


# ======================================================================================
# Geometries (angstrom). Experimental geometries as used by Szabo & Ostlund, Table 3.13, where noted.
# ======================================================================================
def _bent(o: str, h: str, r: float, angle_deg: float) -> list:
    t = math.radians(angle_deg) / 2
    return [(o, (0.0, 0.0, 0.0)), (h, (r * math.sin(t), 0.0, r * math.cos(t))), (h, (-r * math.sin(t), 0.0, r * math.cos(t)))]


def _pyramid(x: str, h: str, r: float, angle_deg: float) -> list:
    """XH3 with three equal X-H bonds and H-X-H angles."""
    a = math.radians(angle_deg)
    # H atoms on a cone: angle between bonds a  ->  cone half-angle b with cos(a) = 1 - 1.5 sin^2(b)
    sb2 = (1 - math.cos(a)) / 1.5
    b = math.asin(math.sqrt(sb2))
    out = [(x, (0.0, 0.0, 0.0))]
    for k in range(3):
        ph = 2 * math.pi * k / 3
        out.append((h, (r * math.sin(b) * math.cos(ph), r * math.sin(b) * math.sin(ph), -r * math.cos(b))))
    return out


def _tetra(x: str, h: str, r: float) -> list:
    s = r / math.sqrt(3)
    return [(x, (0.0, 0.0, 0.0)), (h, (s, s, s)), (h, (-s, -s, s)), (h, (-s, s, -s)), (h, (s, -s, -s))]


def _linear(*atoms_and_positions) -> list:
    return [(a, (0.0, 0.0, z)) for a, z in atoms_and_positions]


B = BOHR
LIBRARY = {   # name -> (description, function(scale) -> atoms, net charge)
    "H2": ("hydrogen", lambda s=1.0: _linear(("H", 0.0), ("H", 1.4 * B * s)), 0),
    "LiH": ("lithium hydride", lambda s=1.0: _linear(("Li", 0.0), ("H", 3.015 * B * s)), 0),
    "HF": ("hydrogen fluoride", lambda s=1.0: _linear(("F", 0.0), ("H", 1.733 * B * s)), 0),
    "H2O": ("water (H-bond donor/acceptor)", lambda s=1.0: _bent("O", "H", 1.809 * B * s, 104.52), 0),
    "NH3": ("ammonia (a basic amine's nitrogen)", lambda s=1.0: _pyramid("N", "H", 1.913 * B * s, 106.7), 0),
    "CH4": ("methane (a hydrophobic carbon)", lambda s=1.0: _tetra("C", "H", 2.050 * B * s), 0),
    "N2": ("nitrogen", lambda s=1.0: _linear(("N", 0.0), ("N", 2.074 * B * s)), 0),
    "CO": ("carbon monoxide", lambda s=1.0: _linear(("C", 0.0), ("O", 2.132 * B * s)), 0),
    "HCN": ("hydrogen cyanide (nitrile group)", lambda s=1.0: _linear(("H", -1.066 * s), ("C", 0.0), ("N", 1.153 * s)), 0),
    "CH2O": ("formaldehyde (carbonyl group)", lambda s=1.0: [("C", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 1.208 * s)),
                                                            ("H", (0.943, 0.0, -0.587)), ("H", (-0.943, 0.0, -0.587))], 0),
}
# Szabo & Ostlund Table 3.13 (STO-3G, experimental geometries): practice anchors for the integrals.
SZABO_OSTLUND_HF = {"H2": -1.117, "CO": -111.225, "N2": -107.496, "CH4": -39.727, "NH3": -55.454, "H2O": -74.963,
                    "HF": -98.571}


# ======================================================================================
# Basis
# ======================================================================================
@dataclass
class Shell:
    centre: np.ndarray        # bohr
    l: tuple                  # (lx, ly, lz)
    alpha: np.ndarray         # exponents
    coef: np.ndarray          # contraction coefficients x primitive normalisation (and contracted renormalised)
    atom: int
    label: str


def _dfact(n: int) -> int:
    return 1 if n <= 0 else n * _dfact(n - 2)


def _prim_norm(a: np.ndarray, l: tuple) -> np.ndarray:
    L = sum(l)
    return (2 * a / math.pi) ** 0.75 * (4 * a) ** (L / 2) / math.sqrt(_dfact(2 * l[0] - 1) * _dfact(2 * l[1] - 1) * _dfact(2 * l[2] - 1))


def basis(atoms: list) -> list[Shell]:
    """Contracted Cartesian basis functions, positions given in angstrom."""
    out = []
    for k, (el, pos) in enumerate(atoms):
        c = np.array(pos, float) / BOHR
        a1 = A1S * ZETA_1S[el] ** 2
        out.append(Shell(c, (0, 0, 0), a1, D1S * _prim_norm(a1, (0, 0, 0)), k, f"{el}{k + 1} 1s"))
        if el in ZETA_2SP:
            a2 = A2SP * ZETA_2SP[el] ** 2
            out.append(Shell(c, (0, 0, 0), a2, D2S * _prim_norm(a2, (0, 0, 0)), k, f"{el}{k + 1} 2s"))
            for l, nm in (((1, 0, 0), "2px"), ((0, 1, 0), "2py"), ((0, 0, 1), "2pz")):
                out.append(Shell(c, l, a2, D2P * _prim_norm(a2, l), k, f"{el}{k + 1} {nm}"))
    for s in out:                                     # renormalise each contraction exactly
        nrm = _overlap(s, s)
        s.coef = s.coef / math.sqrt(nrm)
    return out


# ======================================================================================
# McMurchie-Davidson
# ======================================================================================
def _E(i: int, j: int, t: int, Q: float, a, b):
    """Hermite expansion coefficient, vectorised over arrays a, b (same shape)."""
    p = a + b
    q = a * b / p
    if t < 0 or t > i + j:
        return np.zeros_like(p)
    if i == j == t == 0:
        return np.exp(-q * Q * Q)
    if j == 0:
        return (1 / (2 * p)) * _E(i - 1, j, t - 1, Q, a, b) - (q * Q / a) * _E(i - 1, j, t, Q, a, b) + \
            (t + 1) * _E(i - 1, j, t + 1, Q, a, b)
    return (1 / (2 * p)) * _E(i, j - 1, t - 1, Q, a, b) + (q * Q / b) * _E(i, j - 1, t, Q, a, b) + \
        (t + 1) * _E(i, j - 1, t + 1, Q, a, b)


def boys(n: int, T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, float)
    small = T < 1e-10
    Ts = np.where(small, 1.0, T)
    val = _gamma(n + 0.5) * gammainc(n + 0.5, Ts) / (2 * Ts ** (n + 0.5))
    return np.where(small, 1.0 / (2 * n + 1) - T / (2 * n + 3), val)


def _R_table(L: int, p, PC: np.ndarray) -> dict:
    """Hermite Coulomb integrals R^{0}_{tuv} for t + u + v <= L (arrays over primitives).
    PC has shape (..., 3)."""
    X, Y, Zc = PC[..., 0], PC[..., 1], PC[..., 2]
    T = p * (X * X + Y * Y + Zc * Zc)
    memo: dict = {}
    base = {n: (-2 * p) ** n * boys(n, T) for n in range(L + 1)}

    def R(t, u, v, n):
        key = (t, u, v, n)
        if key in memo:
            return memo[key]
        if t < 0 or u < 0 or v < 0:
            return 0.0
        if t == u == v == 0:
            val = base[n]
        elif t > 0:
            val = (t - 1) * R(t - 2, u, v, n + 1) + X * R(t - 1, u, v, n + 1)
        elif u > 0:
            val = (u - 1) * R(t, u - 2, v, n + 1) + Y * R(t, u - 1, v, n + 1)
        else:
            val = (v - 1) * R(t, u, v - 2, n + 1) + Zc * R(t, u, v - 1, n + 1)
        memo[key] = val
        return val
    return {(t, u, v): R(t, u, v, 0) for t in range(L + 1) for u in range(L + 1 - t) for v in range(L + 1 - t - u)}


@dataclass
class _Pair:
    p: np.ndarray            # (P,) exponent sums
    P: np.ndarray            # (P, 3) centres
    cc: np.ndarray           # (P,) coefficient products
    E: list                  # E[d][t] arrays (P,), d = 0, 1, 2
    lsum: tuple              # la + lb per direction


def _pair(A: Shell, Bs: Shell) -> _Pair:
    a = np.repeat(A.alpha, len(Bs.alpha))
    b = np.tile(Bs.alpha, len(A.alpha))
    cc = np.repeat(A.coef, len(Bs.coef)) * np.tile(Bs.coef, len(A.coef))
    p = a + b
    P = (a[:, None] * A.centre + b[:, None] * Bs.centre) / p[:, None]
    E = []
    for d in range(3):
        Q = A.centre[d] - Bs.centre[d]
        E.append([_E(A.l[d], Bs.l[d], t, Q, a, b) for t in range(A.l[d] + Bs.l[d] + 1)])
    return _Pair(p, P, cc, E, tuple(A.l[d] + Bs.l[d] for d in range(3)))


def _overlap(A: Shell, Bs: Shell) -> float:
    pr = _pair(A, Bs)
    return float(np.sum(pr.cc * pr.E[0][0] * pr.E[1][0] * pr.E[2][0] * (math.pi / pr.p) ** 1.5))


def _overlap_prims(A: Shell, Bs: Shell, lb_shift: tuple) -> np.ndarray:
    """Primitive overlaps (unweighted) with B's angular momentum shifted (for the kinetic integral)."""
    lb = tuple(Bs.l[d] + lb_shift[d] for d in range(3))
    if min(lb) < 0:
        return 0.0
    a = np.repeat(A.alpha, len(Bs.alpha))
    b = np.tile(Bs.alpha, len(A.alpha))
    p = a + b
    out = (math.pi / p) ** 1.5
    for d in range(3):
        out = out * _E(A.l[d], lb[d], 0, A.centre[d] - Bs.centre[d], a, b)
    return out


def _kinetic(A: Shell, Bs: Shell) -> float:
    b = np.tile(Bs.alpha, len(A.alpha))
    cc = np.repeat(A.coef, len(Bs.coef)) * np.tile(Bs.coef, len(A.coef))
    l, m, n = Bs.l
    t0 = b * (2 * (l + m + n) + 3) * _overlap_prims(A, Bs, (0, 0, 0))
    t1 = -2 * b ** 2 * (_overlap_prims(A, Bs, (2, 0, 0)) + _overlap_prims(A, Bs, (0, 2, 0)) + _overlap_prims(A, Bs, (0, 0, 2)))
    t2 = -0.5 * (l * (l - 1) * _overlap_prims(A, Bs, (-2, 0, 0)) + m * (m - 1) * _overlap_prims(A, Bs, (0, -2, 0))
                 + n * (n - 1) * _overlap_prims(A, Bs, (0, 0, -2)))
    return float(np.sum(cc * (t0 + t1 + t2)))


def _nuclear(pr: _Pair, C: np.ndarray, Zc: float) -> float:
    L = sum(pr.lsum)
    R = _R_table(L, pr.p, pr.P - C[None, :])
    tot = np.zeros_like(pr.p)
    for t, et in enumerate(pr.E[0]):
        for u, eu in enumerate(pr.E[1]):
            for v, ev in enumerate(pr.E[2]):
                tot = tot + et * eu * ev * R[(t, u, v)]
    return float(-Zc * np.sum(pr.cc * (2 * math.pi / pr.p) * tot))


def _eri(p1: _Pair, p2: _Pair) -> float:
    p = p1.p[:, None]
    q = p2.p[None, :]
    alpha = p * q / (p + q)
    PQ = p1.P[:, None, :] - p2.P[None, :, :]
    L = sum(p1.lsum) + sum(p2.lsum)
    R = _R_table(L, alpha, PQ)
    E1 = [(t, u, v, (p1.E[0][t] * p1.E[1][u] * p1.E[2][v])[:, None])
          for t in range(len(p1.E[0])) for u in range(len(p1.E[1])) for v in range(len(p1.E[2]))]
    E2 = [(t, u, v, ((-1) ** (t + u + v)) * (p2.E[0][t] * p2.E[1][u] * p2.E[2][v])[None, :])
          for t in range(len(p2.E[0])) for u in range(len(p2.E[1])) for v in range(len(p2.E[2]))]
    tot = 0.0
    for t, u, v, e1 in E1:
        for tau, nu, phi, e2 in E2:
            tot = tot + e1 * e2 * R[(t + tau, u + nu, v + phi)]
    pref = 2 * math.pi ** 2.5 / (p * q * np.sqrt(p + q))
    return float(np.sum(p1.cc[:, None] * p2.cc[None, :] * pref * tot))


@dataclass
class Integrals:
    S: np.ndarray
    T: np.ndarray
    V: np.ndarray
    eri: np.ndarray           # (pq|rs), chemists' notation
    enuc: float
    labels: list
    n_electrons: int


def integrals(atoms: list, charge: int = 0) -> Integrals:
    bas = basis(atoms)
    n = len(bas)
    pairs = {(i, j): _pair(bas[i], bas[j]) for i in range(n) for j in range(i + 1)}
    S, T, V = np.zeros((n, n)), np.zeros((n, n)), np.zeros((n, n))
    centres = [(np.array(pos, float) / BOHR, Z[el]) for el, pos in atoms]
    for (i, j), pr in pairs.items():
        S[i, j] = S[j, i] = float(np.sum(pr.cc * pr.E[0][0] * pr.E[1][0] * pr.E[2][0] * (math.pi / pr.p) ** 1.5))
        T[i, j] = T[j, i] = _kinetic(bas[i], bas[j])
        V[i, j] = V[j, i] = sum(_nuclear(pr, C, Zc) for C, Zc in centres)
    eri = np.zeros((n, n, n, n))
    keys = list(pairs)
    for a, (i, j) in enumerate(keys):
        for (k, l) in keys[:a + 1]:
            v = _eri(pairs[(i, j)], pairs[(k, l)])
            for (p_, q_, r_, s_) in ((i, j, k, l), (j, i, k, l), (i, j, l, k), (j, i, l, k),
                                     (k, l, i, j), (l, k, i, j), (k, l, j, i), (l, k, j, i)):
                eri[p_, q_, r_, s_] = v
    enuc = 0.0
    for a in range(len(centres)):
        for b in range(a):
            enuc += centres[a][1] * centres[b][1] / float(np.linalg.norm(centres[a][0] - centres[b][0]))
    ne = sum(Z[el] for el, _ in atoms) - charge
    return Integrals(S, T, V, eri, enuc, [s.label for s in bas], ne)


# ======================================================================================
# Hartree-Fock
# ======================================================================================
@dataclass
class HF:
    energy: float
    C: np.ndarray
    eps: np.ndarray
    iterations: int
    converged: bool
    n_occ: int


def rhf(ints: Integrals, max_iter: int = 200, tol: float = 1e-10) -> HF:
    """Closed-shell Hartree-Fock. Two routes from the core-Hamiltonian guess, the lower energy kept: a long damped
    iteration polished by DIIS, and DIIS after a short warm-up (DIIS alone can lock onto an excited solution;
    it does for N2)."""
    a = _rhf(ints, max_iter, tol, warmup=400, warm_tol=1e-8)
    b = _rhf(ints, max_iter, tol, warmup=10, warm_tol=1e-4)
    return a if a.energy <= b.energy + 1e-9 else b


def _rhf(ints: Integrals, max_iter: int, tol: float, warmup: int, warm_tol: float) -> HF:
    H = ints.T + ints.V
    S, eri = ints.S, ints.eri
    w, U = np.linalg.eigh(S)
    X = U @ np.diag(w ** -0.5) @ U.T
    nocc = ints.n_electrons // 2
    e, Cp = np.linalg.eigh(X.T @ H @ X)
    C = X @ Cp
    D = 2 * C[:, :nocc] @ C[:, :nocc].T
    # damped warm-up: DIIS from the core guess can lock onto an excited solution (it does for N2)
    for _ in range(warmup):
        J = np.einsum("pqrs,rs->pq", eri, D)
        K = np.einsum("prqs,rs->pq", eri, D)
        e, Cp = np.linalg.eigh(X.T @ (H + J - 0.5 * K) @ X)
        C = X @ Cp
        Dn = 2 * C[:, :nocc] @ C[:, :nocc].T
        step = float(np.max(np.abs(Dn - D)))
        D = 0.5 * D + 0.5 * Dn
        if step < warm_tol:
            break
    focks, errs = [], []
    e_old = 0.0
    conv = False
    for it in range(1, max_iter + 1):
        J = np.einsum("pqrs,rs->pq", eri, D)
        K = np.einsum("prqs,rs->pq", eri, D)
        F = H + J - 0.5 * K
        err = X.T @ (F @ D @ S - S @ D @ F) @ X
        focks.append(F)
        errs.append(err)
        if len(focks) > 8:
            focks.pop(0)
            errs.pop(0)
        if len(focks) >= 2:                       # DIIS extrapolation
            m = len(focks)
            Bm = -np.ones((m + 1, m + 1))
            Bm[m, m] = 0
            for a in range(m):
                for b in range(m):
                    Bm[a, b] = float(np.sum(errs[a] * errs[b]))
            rhs = np.zeros(m + 1)
            rhs[m] = -1
            try:
                c = np.linalg.solve(Bm, rhs)[:m]
                F = sum(ci * Fi for ci, Fi in zip(c, focks))
            except np.linalg.LinAlgError:
                pass
        e, Cp = np.linalg.eigh(X.T @ F @ X)
        C = X @ Cp
        D_new = 2 * C[:, :nocc] @ C[:, :nocc].T
        J = np.einsum("pqrs,rs->pq", eri, D_new)
        K = np.einsum("prqs,rs->pq", eri, D_new)
        e_el = 0.5 * float(np.sum(D_new * (2 * H + J - 0.5 * K)))
        if abs(e_el - e_old) < tol and np.max(np.abs(D_new - D)) < 1e-7:
            D = D_new
            conv = True
            break
        e_old = e_el
        D = D_new
    F = H + np.einsum("pqrs,rs->pq", eri, D) - 0.5 * np.einsum("prqs,rs->pq", eri, D)
    e, Cp = np.linalg.eigh(X.T @ F @ X)
    C = X @ Cp
    return HF(e_el + ints.enuc, C, e, it, conv, nocc)


# ======================================================================================
# Active space, qubits, FCI, VQE
# ======================================================================================
@dataclass
class ActiveSpace:
    n_orbitals: int
    n_electrons: int
    core: list
    active: list
    e_core: float             # nuclear repulsion + frozen-core energy
    h: np.ndarray             # effective one-electron integrals (active MOs)
    g: np.ndarray             # (pq|rs) over active MOs


def active_space(ints: Integrals, hf: HF, n_orbitals: int, n_electrons: int) -> ActiveSpace:
    nocc = hf.n_occ
    n_core = nocc - n_electrons // 2
    if n_core < 0 or n_core + n_orbitals > len(hf.eps) or n_electrons % 2:
        raise ValueError("active space does not fit this molecule")
    core = list(range(n_core))
    act = list(range(n_core, n_core + n_orbitals))
    C = hf.C
    hmo = C.T @ (ints.T + ints.V) @ C
    g_all = np.einsum("pi,qj,rk,sl,pqrs->ijkl", C, C, C, C, ints.eri, optimize=True)
    e_core = ints.enuc + sum(2 * hmo[c, c] for c in core) + \
        sum(2 * g_all[c, c, d, d] - g_all[c, d, d, c] for c in core for d in core)
    h = hmo[np.ix_(act, act)].copy()
    for c in core:
        h += 2 * g_all[np.ix_(act, act, [c], [c])][:, :, 0, 0] - g_all[np.ix_(act, [c], [c], act)][:, 0, 0, :]
    g = g_all[np.ix_(act, act, act, act)]
    return ActiveSpace(n_orbitals, n_electrons, core, act, float(e_core), h, g)


def _ops(n_modes: int):
    import scipy.sparse as sp
    dim = 2 ** n_modes
    idx = np.arange(dim)
    a = []
    for m in range(n_modes):
        has = (idx >> m) & 1
        src = idx[has == 1]
        sign = (-1.0) ** np.array([bin(k & ((1 << m) - 1)).count("1") for k in src])
        a.append(sp.csr_matrix((sign, (src ^ (1 << m), src)), shape=(dim, dim)))
    return a


@dataclass
class QubitModel:
    H: object                 # scipy sparse (2^n x 2^n)
    N: object                 # electron number
    Sz: object
    n_qubits: int
    n_electrons: int
    hf_index: int
    ops: list = field(repr=False, default_factory=list)


def qubit_hamiltonian(asp: ActiveSpace) -> QubitModel:
    import scipy.sparse as sp
    nm = 2 * asp.n_orbitals
    a = _ops(nm)
    ad = [x.T.tocsr() for x in a]
    dim = 2 ** nm
    H = sp.identity(dim, format="csr") * asp.e_core
    for p in range(nm):
        for q in range(nm):
            if p % 2 == q % 2 and abs(asp.h[p // 2, q // 2]) > 1e-12:
                H = H + asp.h[p // 2, q // 2] * (ad[p] @ a[q])
    pairs = [(p, r) for p in range(nm) for r in range(nm) if p % 2 == r % 2]
    for p, r in pairs:
        for q, s in pairs:
            v = asp.g[p // 2, r // 2, q // 2, s // 2]
            if abs(v) > 1e-12 and p != q and r != s:
                H = H + 0.5 * v * (ad[p] @ ad[q] @ a[s] @ a[r])
    N = sum(ad[m] @ a[m] for m in range(nm))
    Sz = sum((0.5 if m % 2 == 0 else -0.5) * (ad[m] @ a[m]) for m in range(nm))
    hf_index = (1 << asp.n_electrons) - 1
    return QubitModel(H.tocsr(), N.tocsr(), Sz.tocsr(), nm, asp.n_electrons, hf_index, a)


def sector(qm: QubitModel) -> np.ndarray:
    idx = np.arange(2 ** qm.n_qubits)
    ne = np.array([bin(k).count("1") for k in idx])
    sz = np.array([sum(((k >> m) & 1) * (0.5 if m % 2 == 0 else -0.5) for m in range(qm.n_qubits)) for k in idx])
    return idx[(ne == qm.n_electrons) & (np.abs(sz) < 1e-9)]


def fci(qm: QubitModel) -> float:
    sec = sector(qm)
    Hs = qm.H[sec][:, sec].toarray()
    return float(np.linalg.eigvalsh(Hs)[0])


def uccsd_generators(qm: QubitModel) -> list:
    """Anti-Hermitian generators T - T^dagger (sparse) of spin-conserving singles and doubles from the HF state."""
    a = qm.ops
    ad = [x.T.tocsr() for x in a]
    n, ne = qm.n_qubits, qm.n_electrons
    occ, vir = list(range(ne)), list(range(ne, n))
    gens = []
    for i in occ:
        for v in vir:
            if i % 2 == v % 2:
                T = ad[v] @ a[i]
                gens.append((f"{i}->{v}", (T - T.T).tocsr()))
    for i, j in itertools.combinations(occ, 2):
        for v, w in itertools.combinations(vir, 2):
            if (i % 2 + j % 2) == (v % 2 + w % 2) and sorted([i % 2, j % 2]) == sorted([v % 2, w % 2]):
                T = ad[v] @ ad[w] @ a[j] @ a[i]
                if T.nnz:
                    gens.append((f"{i}{j}->{v}{w}", (T - T.T).tocsr()))
    return gens


@dataclass
class VQEResult:
    energy: float
    fci: float
    hf: float
    parameters: int
    evaluations: int
    history: list
    electrons: float

    @property
    def error(self) -> float:
        return self.energy - self.fci


def _rot(A, A2, theta: float, v: np.ndarray) -> np.ndarray:
    """exp(theta A) v for a fermionic excitation generator A = T - T^dagger (A^3 = -A on its support):
    exp(theta A) = 1 + sin(theta) A + (1 - cos(theta)) A^2."""
    return v + math.sin(theta) * (A @ v) + (1 - math.cos(theta)) * (A2 @ v)


def vqe_uccsd(qm: QubitModel, maxiter: int = 300, gens: list | None = None) -> VQEResult:
    """Disentangled UCCSD, exact closed-form rotations and an adjoint (backward-pass) gradient, BFGS."""
    from scipy.optimize import minimize
    gens = uccsd_generators(qm) if gens is None else gens
    As = [G for _, G in gens]
    A2s = [(G @ G).tocsr() for G in As]
    psi0 = np.zeros(2 ** qm.n_qubits)
    psi0[qm.hf_index] = 1.0
    hist: list[float] = []

    def state(th):
        psi = psi0.copy()
        for A, A2, t in zip(As, A2s, th):
            psi = _rot(A, A2, t, psi)
        return psi

    def f_and_grad(th):
        psi = state(th)
        lam = qm.H @ psi
        e = float(psi @ lam)
        hist.append(e)
        g = np.zeros(len(th))
        cur = psi
        for k in range(len(th) - 1, -1, -1):
            g[k] = 2.0 * float(lam @ (As[k] @ cur))
            cur = _rot(As[k], A2s[k], -th[k], cur)
            lam = _rot(As[k], A2s[k], -th[k], lam)
        return e, g

    x0 = np.zeros(len(As))
    e_hf = float(psi0 @ (qm.H @ psi0))
    if len(As):
        r = minimize(f_and_grad, x0, jac=True, method="BFGS", options={"maxiter": maxiter, "gtol": 1e-8})
        th = r.x
    else:
        th = x0
    psi = state(th)
    e = float(psi @ (qm.H @ psi))
    return VQEResult(e, fci(qm), e_hf, len(As), len(hist), hist, float(psi @ (qm.N @ psi)))


@dataclass
class MoleculeResult:
    name: str
    atoms: list
    n_basis: int
    e_hf: float
    hf_converged: bool
    homo_lumo_gap: float      # hartree
    orbital_energies: np.ndarray
    active: tuple             # (electrons, orbitals)
    n_qubits: int
    e_cas_fci: float
    vqe: VQEResult | None
    seconds: float
    labels: list


def run(name: str, atoms: list | None = None, charge: int = 0, active: tuple = (2, 2), do_vqe: bool = True,
        scale: float = 1.0) -> MoleculeResult:
    import time
    t0 = time.perf_counter()
    if atoms is None:
        _, fn, charge = LIBRARY[name]
        atoms = fn(scale)
    ints = integrals(atoms, charge)
    hf = rhf(ints)
    ne_act, no_act = active
    asp = active_space(ints, hf, no_act, ne_act)
    qm = qubit_hamiltonian(asp)
    v = vqe_uccsd(qm) if do_vqe else None
    e_cas = v.fci if v is not None else fci(qm)
    gap = float(hf.eps[hf.n_occ] - hf.eps[hf.n_occ - 1]) if hf.n_occ < len(hf.eps) else float("nan")
    return MoleculeResult(name, atoms, len(ints.S), hf.energy, hf.converged, gap, hf.eps, active, qm.n_qubits, e_cas, v,
                          time.perf_counter() - t0, ints.labels)
