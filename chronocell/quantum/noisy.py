"""
A noisy quantum computer, simulated gate by gate, and zero-noise extrapolation (round 2, October 2026).

Density-matrix simulation (up to 10 qubits): every gate U acts as rho -> U rho U^dagger and is followed by a
depolarizing channel on the qubits it touched, rho -> (1 - p) rho + p Tr_Q(rho) (x) I/2^k, with p = p1 after one-qubit
gates and p2 after two-qubit gates (the local noise of today's superconducting chips, as in Qiskit Aer's
depolarizing_error); readout errors flip each measured bit with probability `readout` (diagonal observables only).
This is local noise, unlike the single global-fidelity estimate used elsewhere in the lab.

Zero-noise extrapolation (Temme, Bravyi & Gambetta, PRL 2017; Giurgica-Tiron et al., IEEE QCE 2020): the circuit is
run at amplified noise by unitary folding, C -> C (C^dagger C)^k (noise scale 2k + 1, the same ideal result), and the
measured values at scales 1, 3, 5 are extrapolated back to zero noise (Richardson: the quadratic through the three
points; linear: a least-squares line). Nothing about the ideal answer is used.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import sim

MAX_QUBITS = 10
INVERSE = {"h": "h", "x": "x", "y": "y", "z": "z", "s": "sdg", "sdg": "s", "cx": "cx", "cz": "cz", "swap": "swap",
           "cswap": "cswap"}


def _gate_matrix(g: str, ps: tuple) -> np.ndarray:
    return sim.FIXED[g] if g in sim.FIXED else sim.ROTATION[g](*ps)


def apply_unitary_dm(rho: np.ndarray, U: np.ndarray, qubits, n: int) -> np.ndarray:
    k = len(qubits)
    t = rho.reshape((2,) * (2 * n))
    ar = [n - 1 - q for q in reversed(qubits)]
    Ut = U.reshape((2,) * (2 * k))
    out = np.moveaxis(np.tensordot(Ut, t, axes=(list(range(k, 2 * k)), ar)), list(range(k)), ar)
    ac = [n + a for a in ar]
    out = np.moveaxis(np.tensordot(Ut.conj(), out, axes=(list(range(k, 2 * k)), ac)), list(range(k)), ac)
    return out.reshape(2 ** n, 2 ** n)


def depolarize(rho: np.ndarray, qubits, n: int, p: float) -> np.ndarray:
    """(1 - p) rho + p Tr_Q(rho) (x) I/2^k on the qubits Q."""
    if p <= 0:
        return rho
    t = rho.reshape((2,) * (2 * n))
    m = t
    for q in qubits:
        a = n - 1 - q
        red = np.trace(m, axis1=a, axis2=a + n)                         # removes both axes
        red = np.expand_dims(np.expand_dims(red, a), a + n)              # back to 2n axes, size 1
        eye = (np.eye(2) / 2).reshape([2 if i in (a, a + n) else 1 for i in range(2 * n)])
        m = red * eye
    return (1 - p) * rho + p * m.reshape(2 ** n, 2 ** n)


def run(circ: sim.Circuit, p1: float, p2: float, rho: np.ndarray | None = None) -> np.ndarray:
    n = circ.n
    if n > MAX_QUBITS:
        raise ValueError(f"density-matrix simulation is limited to {MAX_QUBITS} qubits")
    if rho is None:
        rho = np.zeros((2 ** n, 2 ** n), complex)
        rho[0, 0] = 1.0
    for g, qs, ps in circ.gates():
        rho = apply_unitary_dm(rho, _gate_matrix(g, ps), qs, n)
        rho = depolarize(rho, qs, n, p1 if len(qs) == 1 else p2)
    return rho


def readout_probabilities(rho: np.ndarray, n: int, readout: float) -> np.ndarray:
    """Measured bitstring probabilities with each bit flipped with probability `readout`."""
    p = np.real(np.diag(rho)).clip(0)
    if readout <= 0:
        return p
    t = p.reshape((2,) * n)
    M = np.array([[1 - readout, readout], [readout, 1 - readout]])
    for a in range(n):
        t = np.moveaxis(np.tensordot(M, t, axes=([1], [a])), 0, a)
    return t.reshape(-1)


def inverse(circ: sim.Circuit) -> sim.Circuit:
    out = sim.Circuit(circ.n, circ.name + "^dagger")
    for g, qs, ps in reversed(list(circ.gates())):
        if g in sim.ROTATION:
            out.add(g, qs, tuple(-x for x in ps))
        else:
            out.add(INVERSE[g], qs)
    return out


def fold(circ: sim.Circuit, scale: int) -> sim.Circuit:
    """Global unitary folding: C (C^dagger C)^k with scale = 2k + 1 (odd)."""
    if scale < 1 or scale % 2 == 0:
        raise ValueError("scale must be an odd integer >= 1")
    out = sim.Circuit(circ.n, f"{circ.name} x{scale}")
    inv = inverse(circ)
    for c in [circ] + [inv, circ] * ((scale - 1) // 2):
        for g, qs, ps in c.gates():
            out.add(g, qs, ps)
    return out


def extrapolate(scales, values, method: str = "richardson", asymptote: float | None = None) -> float:
    """richardson: the polynomial through the points at 0; linear: least-squares line; exp: E(s) = E_inf + b a^s with
    E_inf the fully mixed state's value (known from the observable alone), log-linear fit."""
    s, v = np.asarray(scales, float), np.asarray(values, float)
    if method == "linear":
        return float(np.polyfit(s, v, 1)[1])
    if method == "exp" and asymptote is not None:
        d = v - asymptote
        if np.all(d > 0) or np.all(d < 0):
            sign = np.sign(d[0])
            slope, icpt = np.polyfit(s, np.log(np.abs(d)), 1)
            return float(asymptote + sign * math.exp(icpt))
        return float(np.polyval(np.polyfit(s, v, len(s) - 1), 0.0))
    return float(np.polyval(np.polyfit(s, v, len(s) - 1), 0.0))          # Richardson: exact polynomial through the points


@dataclass
class ZNEResult:
    ideal: float
    noisy: float
    mitigated: float
    scales: list
    values: list
    method: str

    @property
    def error_noisy(self) -> float:
        return self.noisy - self.ideal

    @property
    def error_mitigated(self) -> float:
        return self.mitigated - self.ideal


def expectation(circ: sim.Circuit, observable, noise: dict, scale: int = 1, keep: np.ndarray | None = None) -> float:
    """Expectation at a noise scale: `observable` is a diagonal (array of 2^n values, read out with readout errors) or
    a Hermitian matrix (Tr rho H; readout errors not modelled for non-diagonal measurements). keep: a 0/1 diagonal
    projector for symmetry verification (rho -> P rho P / Tr P rho P: results with the wrong electron number are
    discarded, as by post-selection)."""
    c = fold(circ, scale) if scale > 1 else circ
    rho = run(c, noise.get("p1", 0.0), noise.get("p2", 0.0))
    if keep is not None:
        P = np.asarray(keep, float)
        rho = (P[:, None] * rho) * P[None, :]
        rho = rho / max(float(np.real(np.trace(rho))), 1e-300)
    obs = np.asarray(observable)
    if obs.ndim == 1:
        return float(readout_probabilities(rho, circ.n, noise.get("readout", 0.0)) @ obs)
    return float(np.real(np.trace(rho @ obs)))


def mixed_value(observable, keep: np.ndarray | None = None) -> float:
    """The observable on the fully mixed state (within the kept sector): the zero-signal limit of strong noise."""
    obs = np.asarray(observable)
    d = obs if obs.ndim == 1 else np.real(np.diag(obs))
    if keep is None:
        return float(d.mean())
    k = np.asarray(keep, bool)
    return float(d[k].mean())


def zne(circ: sim.Circuit, observable, noise: dict, scales=(1, 3, 5), method: str = "richardson",
        ideal: float | None = None, keep: np.ndarray | None = None) -> ZNEResult:
    vals = [expectation(circ, observable, noise, s, keep) for s in scales]
    if ideal is None:
        psi = circ.run()
        obs = np.asarray(observable)
        ideal = float(np.abs(psi) ** 2 @ obs) if obs.ndim == 1 else float(np.real(psi.conj() @ obs @ psi))
    return ZNEResult(ideal, vals[0], extrapolate(scales, vals, method, mixed_value(observable, keep)), list(scales), vals,
                     method)


# ======================================================================================
# Gate-level VQE circuits for the active-space molecules (molecules.py)
# ======================================================================================
def excitation_paulis(A, n: int, tol: float = 1e-10) -> dict:
    """Pauli decomposition of i * A for a real antisymmetric excitation generator A (so exp(theta A) =
    exp(-i theta (i A)) with i A Hermitian); the strings of one excitation commute, so the product of their
    exponentials is exact."""
    M = 1j * np.asarray(A.toarray() if hasattr(A, "toarray") else A)
    return {k: float(np.real(v)) for k, v in sim.pauli_decompose(M, tol).items()}


def vqe_circuit(qm, gens: list, thetas: np.ndarray, name: str = "UCCSD-VQE") -> sim.Circuit:
    """Hartree-Fock preparation (X on the occupied spin orbitals) and the disentangled UCCSD product
    prod_k exp(theta_k A_k), compiled to Pauli exponentials (basis change, CX ladder, RZ)."""
    n = qm.n_qubits
    c = sim.Circuit(n, name)
    for q in range(n):
        if (qm.hf_index >> q) & 1:
            c.x(q)
    for (_, A), th in zip(gens, thetas):
        for label, coef in excitation_paulis(A, n).items():
            sim.pauli_exponential(c, label, th * coef)
    return c
