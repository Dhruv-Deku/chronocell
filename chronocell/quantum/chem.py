"""
Minimal-basis quantum chemistry from scratch, mapped to qubits and solved with VQE (the quantum algorithm for
molecular energies; Peruzzo et al., Nat Commun 2014).

Molecules: H2 and HeH+ (two electrons each, one 1s function per atom: 4 spin orbitals = 4 qubits), the standard
first molecules run on quantum computers (O'Malley et al., PRX 2016; Kandala et al., Nature 2017).

Pipeline (all here, no chemistry package): STO-3G contracted Gaussians (Szabo & Ostlund, Appendix B; zeta 1.24 for
H, 2.0925 for He as in their HeH+ example) -> closed-form one- and two-electron integrals over s Gaussians ->
restricted Hartree-Fock -> molecular-orbital integrals -> second-quantised Hamiltonian -> Jordan-Wigner qubits
(spin orbital m = 2 x orbital + spin is qubit m) -> exact diagonalisation in the two-electron sector (FCI,
the reference) -> VQE with the UCCSD ansatz (or a hardware-efficient one) on the statevector simulator.

The Drug lab's drugs have 30-100 heavy atoms; molecules of that size need thousands of error-corrected qubits,
far beyond any quantum computer today. This is a demonstrator of the method, not a drug calculation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.special import erf

from . import sim

BOHR_ANGSTROM = 0.529177210903
HARTREE_KCAL = 627.509474
CHEMICAL_ACCURACY = 1.6e-3                     # hartree (1 kcal/mol)
STO3G_ALPHA = np.array([0.109818, 0.405771, 2.22766])          # zeta = 1
STO3G_COEF = np.array([0.444635, 0.535328, 0.154329])
ZETA = {"H": 1.24, "He": 2.0925}
CHARGE = {"H": 1.0, "He": 2.0}
MOLECULES = {"H2": (("H", "H"), 0), "HeH+": (("He", "H"), 1)}       # atoms, net charge
EQUILIBRIUM_ANGSTROM = {"H2": 0.7414, "HeH+": 0.7743}


def _f0(t: float) -> float:
    return 1.0 if t < 1e-12 else 0.5 * math.sqrt(math.pi / t) * erf(math.sqrt(t))


@dataclass
class Basis:
    centres: np.ndarray            # (K, 3) bohr
    alphas: list                   # per function: 3 exponents
    coefs: list                    # per function: 3 normalised contraction coefficients


def basis_for(atoms: tuple[str, ...], R_bohr: float) -> Basis:
    centres = np.array([[0.0, 0.0, 0.0], [R_bohr, 0.0, 0.0]])[:len(atoms)]
    alphas, coefs = [], []
    for a in atoms:
        al = STO3G_ALPHA * ZETA[a] ** 2
        alphas.append(al)
        coefs.append(STO3G_COEF * (2 * al / math.pi) ** 0.75)
    return Basis(centres, alphas, coefs)


def integrals(atoms: tuple[str, ...], R_bohr: float) -> dict:
    """Overlap S, core Hamiltonian H, two-electron (pq|rs) in chemists' notation, nuclear repulsion."""
    B = basis_for(atoms, R_bohr)
    K = len(atoms)
    Zs = [CHARGE[a] for a in atoms]
    S = np.zeros((K, K))
    T = np.zeros((K, K))
    V = np.zeros((K, K))
    for m in range(K):
        for n in range(K):
            A, Bc = B.centres[m], B.centres[n]
            rab2 = float(np.sum((A - Bc) ** 2))
            for a, ca in zip(B.alphas[m], B.coefs[m]):
                for b, cb in zip(B.alphas[n], B.coefs[n]):
                    p = a + b
                    mu = a * b / p
                    pre = (math.pi / p) ** 1.5 * math.exp(-mu * rab2)
                    S[m, n] += ca * cb * pre
                    T[m, n] += ca * cb * mu * (3 - 2 * mu * rab2) * pre
                    P = (a * A + b * Bc) / p
                    for C, Zc in zip(B.centres, Zs):
                        V[m, n] += ca * cb * (-2 * math.pi / p) * Zc * math.exp(-mu * rab2) * _f0(p * float(np.sum((P - C) ** 2)))
    eri = np.zeros((K, K, K, K))
    for p_ in range(K):
        for q in range(K):
            for r in range(K):
                for s in range(K):
                    eri[p_, q, r, s] = _eri(B, p_, q, r, s)
    enuc = Zs[0] * Zs[1] / R_bohr if K == 2 else 0.0
    return {"S": S, "H": T + V, "eri": eri, "enuc": enuc}


def _eri(B: Basis, p: int, q: int, r: int, s: int) -> float:
    A, Bc, C, D = B.centres[p], B.centres[q], B.centres[r], B.centres[s]
    rab2 = float(np.sum((A - Bc) ** 2))
    rcd2 = float(np.sum((C - D) ** 2))
    tot = 0.0
    for a, ca in zip(B.alphas[p], B.coefs[p]):
        for b, cb in zip(B.alphas[q], B.coefs[q]):
            pp = a + b
            P = (a * A + b * Bc) / pp
            for c, cc in zip(B.alphas[r], B.coefs[r]):
                for d, cd in zip(B.alphas[s], B.coefs[s]):
                    qq = c + d
                    Q = (c * C + d * D) / qq
                    tot += ca * cb * cc * cd * 2 * math.pi ** 2.5 / (pp * qq * math.sqrt(pp + qq)) * \
                        math.exp(-a * b / pp * rab2 - c * d / qq * rcd2) * \
                        _f0(pp * qq / (pp + qq) * float(np.sum((P - Q) ** 2)))
    return tot


def rhf(ints: dict, n_electrons: int = 2, max_iter: int = 100, tol: float = 1e-10) -> dict:
    S, Hc, eri = ints["S"], ints["H"], ints["eri"]
    w, U = np.linalg.eigh(S)
    X = U @ np.diag(w ** -0.5) @ U.T
    K = len(S)
    P = np.zeros((K, K))
    nocc = n_electrons // 2
    e_old = 0.0
    for it in range(max_iter):
        G = np.einsum("ls,mnsl->mn", P, eri) - 0.5 * np.einsum("ls,mlsn->mn", P, eri)
        F = Hc + G
        e, Cp = np.linalg.eigh(X.T @ F @ X)
        C = X @ Cp
        P = 2 * C[:, :nocc] @ C[:, :nocc].T
        e_el = 0.5 * float(np.sum(P * (Hc + F)))
        if abs(e_el - e_old) < tol and it > 0:
            break
        e_old = e_el
    return {"C": C, "orbital_energies": e, "e_hf": e_el + ints["enuc"], "iterations": it + 1}


def _annihilators(n_modes: int) -> list[np.ndarray]:
    """Jordan-Wigner annihilation operators; basis index bit m = occupation of mode m (= qubit m)."""
    dim = 2 ** n_modes
    ops = []
    for m in range(n_modes):
        a = np.zeros((dim, dim))
        for k in range(dim):
            if (k >> m) & 1:
                sign = (-1) ** bin(k & ((1 << m) - 1)).count("1")
                a[k ^ (1 << m), k] = sign
        ops.append(a)
    return ops


@dataclass
class Molecule:
    name: str
    R_angstrom: float
    e_hf: float
    e_fci: float
    hamiltonian: np.ndarray            # 16 x 16 qubit Hamiltonian (real symmetric)
    number: np.ndarray                 # electron-number operator
    hf_index: int
    n_qubits: int
    paulis: dict = field(default_factory=dict)
    enuc: float = 0.0


def molecule(name: str, R_angstrom: float) -> Molecule:
    atoms, charge = MOLECULES[name]
    ne = sum(int(CHARGE[a]) for a in atoms) - charge
    R = R_angstrom / BOHR_ANGSTROM
    ints = integrals(atoms, R)
    hf = rhf(ints, ne)
    C = hf["C"]
    h = C.T @ ints["H"] @ C
    g = np.einsum("pi,qj,rk,sl,pqrs->ijkl", C, C, C, C, ints["eri"])          # (ij|kl) in MOs
    K = len(h)
    nm = 2 * K
    a = _annihilators(nm)
    ad = [x.T for x in a]
    Hq = ints["enuc"] * np.eye(2 ** nm)
    for p in range(nm):
        for q in range(nm):
            if p % 2 == q % 2:
                Hq += h[p // 2, q // 2] * ad[p] @ a[q]
    for p in range(nm):
        for q in range(nm):
            for r in range(nm):
                for s in range(nm):
                    if p % 2 == r % 2 and q % 2 == s % 2:
                        v = g[p // 2, r // 2, q // 2, s // 2]                  # <pq|rs> = (pr|qs)
                        if abs(v) > 1e-14:
                            Hq += 0.5 * v * ad[p] @ ad[q] @ a[s] @ a[r]
    N = sum(ad[m] @ a[m] for m in range(nm))
    sector = np.flatnonzero(np.isclose(np.diag(N), ne))
    e_fci = float(np.linalg.eigvalsh(Hq[np.ix_(sector, sector)])[0])
    hf_index = (1 << ne) - 1
    return Molecule(name, float(R_angstrom), float(hf["e_hf"]), e_fci, Hq, N, hf_index, nm,
                    sim.pauli_decompose(Hq), float(ints["enuc"]))


# ======================================================================================
# VQE
# ======================================================================================
def uccsd_generators(n_qubits: int = 4, n_electrons: int = 2) -> list[tuple[str, np.ndarray]]:
    """Hermitian generators G_k = i (T_k - T_k^dagger) of spin-conserving singles and the double excitation."""
    a = _annihilators(n_qubits)
    ad = [x.T for x in a]
    occ = list(range(n_electrons))
    vir = list(range(n_electrons, n_qubits))
    gens = []
    for i in occ:
        for v in vir:
            if i % 2 == v % 2:
                T = ad[v] @ a[i]
                gens.append((f"single {i}->{v}", 1j * (T - T.T)))
    if n_electrons == 2 and len(vir) >= 2:
        T = ad[vir[0]] @ ad[vir[1]] @ a[occ[1]] @ a[occ[0]]
        gens.append((f"double {occ[0]}{occ[1]}->{vir[0]}{vir[1]}", 1j * (T - T.T)))
    return gens


def uccsd_circuit(mol: Molecule, thetas: np.ndarray) -> sim.Circuit:
    """Hartree-Fock state, then exp(-i theta_k G_k) for each generator as Pauli exponentials (the Pauli terms of one
    excitation commute, so each factor is exact)."""
    c = sim.Circuit(mol.n_qubits, f"UCCSD-VQE {mol.name} R={mol.R_angstrom:g} A")
    for q in range(mol.n_qubits):
        if (mol.hf_index >> q) & 1:
            c.x(q)
    for (name, G), th in zip(uccsd_generators(mol.n_qubits), thetas):
        for lab, coef in sim.pauli_decompose(G).items():
            sim.pauli_exponential(c, lab, th * coef)
    return c


def hea_circuit(n: int, thetas: np.ndarray, layers: int, hf_index: int = 0) -> sim.Circuit:
    """Hardware-efficient ansatz: RY on every qubit, then a CX chain, repeated; starts from the HF bitstring."""
    c = sim.Circuit(n, f"hardware-efficient ansatz ({layers} layers)")
    for q in range(n):
        if (hf_index >> q) & 1:
            c.x(q)
    t = iter(thetas)
    for _ in range(layers):
        for q in range(n):
            c.ry(next(t), q)
        for q in range(n - 1):
            c.cx(q, q + 1)
    for q in range(n):
        c.ry(next(t), q)
    return c


@dataclass
class VQEResult:
    ansatz: str
    energy: float
    thetas: np.ndarray
    evaluations: int
    history: list
    error: float                         # energy - FCI
    circuit: sim.Circuit
    noisy_energy: float = float("nan")
    fidelity: float = 1.0
    electrons: float = float("nan")       # <N> of the final state (the molecule has n_electrons)


def vqe(mol: Molecule, ansatz: str = "uccsd", layers: int = 2, seed: int = 0, noise: dict | None = None,
        maxiter: int = 400, number_penalty: float | None = None) -> VQEResult:
    """VQE. UCCSD conserves the electron count by construction. The hardware-efficient ansatz does not, and on the
    4-qubit Hamiltonian (which also holds 0-4 electron states) it can drift to a state with the wrong number of
    electrons and an energy below the molecule's: the standard remedy, a penalty number_penalty * <(N - n)^2>
    (hartree), is added to its cost by default (2.0). The reported energy never includes the penalty."""
    from scipy.optimize import minimize
    H = mol.hamiltonian
    ne = bin(mol.hf_index).count("1")
    pen = (0.0 if ansatz == "uccsd" else 2.0) if number_penalty is None else float(number_penalty)
    dN = mol.number - ne * np.eye(len(H))
    D2 = dN @ dN
    hist: list[float] = []
    if ansatz == "uccsd":
        k = len(uccsd_generators(mol.n_qubits))

        def build(th):
            return uccsd_circuit(mol, th)
        x0 = np.zeros(k)
    else:
        k = mol.n_qubits * (layers + 1)

        def build(th):
            return hea_circuit(mol.n_qubits, th, layers, mol.hf_index)
        x0 = np.random.default_rng(seed).normal(0, 0.05, k)

    def f(th):
        psi = build(th).run()
        e = float(np.real(np.vdot(psi, H @ psi)))
        hist.append(e)
        return e + (pen * float(np.real(np.vdot(psi, D2 @ psi))) if pen else 0.0)

    r = minimize(f, x0, method="BFGS", options={"maxiter": maxiter, "gtol": 1e-8})
    th = r.x
    c = build(th)
    psi = c.run()
    e = float(np.real(np.vdot(psi, H @ psi)))
    out = VQEResult(ansatz, e, th, len(hist), hist, e - mol.e_fci, c,
                    electrons=float(np.real(np.vdot(psi, mol.number @ psi))))
    if noise is not None:
        F = sim.fidelity_estimate(c.counts(), noise)
        out.fidelity = F
        out.noisy_energy = F * e + (1 - F) * float(np.trace(H).real) / H.shape[0]
    return out


def shot_error(mol: Molecule, psi: np.ndarray, shots: int) -> float:
    """Standard error of the energy when every Pauli term is measured with `shots` shots."""
    var = 0.0
    for lab, c in mol.paulis.items():
        if set(lab) == {"I"}:
            continue
        P = sim.pauli_matrix(lab)
        ev = float(np.real(np.vdot(psi, P @ psi)))
        var += c * c * max(0.0, 1 - ev * ev) / shots
    return math.sqrt(var)


def curve(name: str, R_list, ansatz: str = "uccsd", noise: dict | None = None, layers: int = 2) -> list[dict]:
    rows = []
    for R in R_list:
        m = molecule(name, float(R))
        v = vqe(m, ansatz, layers=layers, noise=noise)
        rows.append({"molecule": name, "R_angstrom": float(R), "e_hf": m.e_hf, "e_fci": m.e_fci, "e_vqe": v.energy,
                     "error_mha": 1e3 * v.error, "within_chemical_accuracy": bool(abs(v.error) <= CHEMICAL_ACCURACY),
                     "e_vqe_noisy": v.noisy_energy, "fidelity": v.fidelity, "evaluations": v.evaluations,
                     "electrons": v.electrons,
                     "pauli_terms": len(m.paulis), "cx": v.circuit.counts()["cx_equivalent"]})
    return rows
