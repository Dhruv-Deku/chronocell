"""
Quantum kernels and state overlaps.

ZZ feature map (Havlicek et al., Nature 2019; Qiskit's ZZFeatureMap with full entanglement): each feature x_i is a
qubit; per repetition, H on all qubits, then phase 2 phi_i on qubit i and phase 2 phi_ij on the parity of qubits i
and j, with phi_i = s x_i and phi_ij = (pi - s x_i)(pi - s x_j); s is the bandwidth (Shaydulin & Wild, PRA 2022:
without a tuned bandwidth quantum kernels concentrate and learn nothing). Kernel K(x, y) = |<phi(x)|phi(y)>|^2,
the probability of reading all zeros after U(x) then U(y)^dagger (compute-uncompute), estimated with shots on
hardware; here exact, or with simulated shot noise.

QSVM: a support-vector machine on that kernel (scikit-learn, precomputed kernel).

Amplitude encoding and the swap test: a non-negative vector v of length 2^q becomes the q-qubit state v / |v|;
the swap test (an ancilla, H, controlled swaps, H) reads 0 with probability (1 + |<a|b>|^2) / 2. For non-negative
real vectors |<a|b>|^2 is the squared cosine similarity: the circuit measures a classical quantity, shown to teach
how a quantum computer compares two states.
"""

from __future__ import annotations

import math

import numpy as np

from . import sim


def _phases(x: np.ndarray, s: float) -> np.ndarray:
    """Diagonal phase of one ZZ feature-map layer over all 2^q basis states."""
    x = np.asarray(x, float) * s
    q = len(x)
    bits = sim.bits_of(np.arange(2 ** q), q).astype(float)
    ph = 2 * bits @ x
    for i in range(q):
        for j in range(i + 1, q):
            par = (bits[:, i] + bits[:, j]) % 2
            ph += 2 * (math.pi - x[i]) * (math.pi - x[j]) * par
    return ph


def _hadamard_all(psi: np.ndarray, q: int) -> np.ndarray:
    for k in range(q):
        psi = sim.apply_unitary(psi, sim.H, (k,), q)
    return psi


def feature_state(x: np.ndarray, bandwidth: float = 1.0, reps: int = 2) -> np.ndarray:
    q = len(x)
    psi = sim.zero_state(q)
    ph = np.exp(1j * _phases(x, bandwidth))
    for _ in range(reps):
        psi = _hadamard_all(psi, q) * ph
    return psi


def feature_map_circuit(x: np.ndarray, bandwidth: float = 1.0, reps: int = 2) -> sim.Circuit:
    """The same map as gates (RZ in place of the phase gate: equal up to a global phase)."""
    x = np.asarray(x, float) * bandwidth
    q = len(x)
    c = sim.Circuit(q, f"ZZ feature map ({reps} reps)")
    for _ in range(reps):
        for i in range(q):
            c.h(i)
        for i in range(q):
            c.rz(2 * x[i], i)
        for i in range(q):
            for j in range(i + 1, q):
                c.cx(i, j)
                c.rz(2 * (math.pi - x[i]) * (math.pi - x[j]), j)
                c.cx(i, j)
    return c


def states(X: np.ndarray, bandwidth: float = 1.0, reps: int = 2) -> np.ndarray:
    return np.stack([feature_state(x, bandwidth, reps) for x in np.asarray(X, float)])


def kernel(XA: np.ndarray, XB: np.ndarray | None = None, bandwidth: float = 1.0, reps: int = 2,
           shots: int | None = None, seed: int = 0) -> np.ndarray:
    SA = states(XA, bandwidth, reps)
    SB = SA if XB is None else states(XB, bandwidth, reps)
    K = np.abs(SA.conj() @ SB.T) ** 2
    if shots:
        rng = np.random.default_rng(seed)
        K = rng.binomial(int(shots), np.clip(K, 0, 1)) / float(shots)
        if XB is None:
            K = np.triu(K, 1) + np.triu(K, 1).T + np.eye(len(K))
    return K


class QSVM:
    """Support-vector classifier on the ZZ quantum kernel."""

    def __init__(self, C: float = 1.0, bandwidth: float = 1.0, reps: int = 2, shots: int | None = None, seed: int = 0):
        self.C, self.bandwidth, self.reps, self.shots, self.seed = C, bandwidth, reps, shots, seed

    def fit(self, X: np.ndarray, y: np.ndarray) -> "QSVM":
        from sklearn.svm import SVC
        self.X_ = np.asarray(X, float)
        self.svc_ = SVC(C=self.C, kernel="precomputed", class_weight="balanced")
        self.svc_.fit(kernel(self.X_, None, self.bandwidth, self.reps, self.shots, self.seed), np.asarray(y))
        return self

    def decision_function(self, X: np.ndarray) -> np.ndarray:
        K = kernel(np.asarray(X, float), self.X_, self.bandwidth, self.reps, self.shots, self.seed + 1)
        return self.svc_.decision_function(K)

    def predict(self, X: np.ndarray) -> np.ndarray:
        K = kernel(np.asarray(X, float), self.X_, self.bandwidth, self.reps, self.shots, self.seed + 1)
        return self.svc_.predict(K)


def scale_features(train: np.ndarray, *others: np.ndarray) -> tuple:
    """Min-max scale to [0, pi] using the training set's range (clipped for the others)."""
    tr = np.asarray(train, float)
    lo, hi = np.nanmin(tr, axis=0), np.nanmax(tr, axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    f = (lambda a: np.clip((np.asarray(a, float) - lo) / span, 0, 1) * math.pi)
    return (f(tr),) + tuple(f(o) for o in others)


# ======================================================================================
# Amplitude encoding and the swap test
# ======================================================================================
def amplitude_vector(v: np.ndarray, q: int | None = None) -> np.ndarray:
    v = np.nan_to_num(np.abs(np.asarray(v, float)))
    q = int(math.ceil(math.log2(max(2, len(v))))) if q is None else q
    out = np.zeros(2 ** q)
    out[:min(len(v), 2 ** q)] = v[:2 ** q]
    nrm = np.linalg.norm(out)
    if nrm == 0:
        out[0] = 1.0
        nrm = 1.0
    return out / nrm


def _gray(i: int) -> int:
    return i ^ (i >> 1)


def amplitude_prep(c: sim.Circuit, amps: np.ndarray, qubits: list[int]) -> sim.Circuit:
    """Prepare a non-negative real state on `qubits` (qubits[0] = lowest bit) from |0...0> with uniformly controlled
    RY rotations (Mottonen et al., QIC 2005), decomposed into RY and CX with a Gray code."""
    a = np.asarray(amps, float)
    q = len(qubits)
    for k in range(q - 1, -1, -1):                   # target bit k; controls = bits above k
        m = q - 1 - k
        theta = np.zeros(2 ** m)
        for cval in range(2 ** m):
            block = a.reshape(2 ** m, 2, 2 ** k)[cval]
            n0, n1 = np.linalg.norm(block[0]), np.linalg.norm(block[1])
            theta[cval] = 2 * math.atan2(n1, n0)
        tgt = qubits[k]
        if m == 0:
            c.ry(theta[0], tgt)
            continue
        ctrls = [qubits[k + 1 + t] for t in range(m)]        # ctrls[t] holds bit t of the control value
        M = np.array([[(-1) ** bin(j & _gray(i)).count("1") for j in range(2 ** m)] for i in range(2 ** m)], float)
        alpha = M @ theta / 2 ** m
        for i in range(2 ** m):
            c.ry(alpha[i], tgt)
            nxt = (i + 1) % (2 ** m)
            diff = _gray(i) ^ _gray(nxt)
            bit = int(math.log2(diff)) if diff else m - 1
            c.cx(ctrls[bit], tgt)
    return c


def swap_test_circuit(a: np.ndarray, b: np.ndarray) -> sim.Circuit:
    """Ancilla qubit 0; register A on qubits 1..q, register B on q+1..2q; state preparation included."""
    q = int(round(math.log2(len(a))))
    c = sim.Circuit(2 * q + 1, f"swap test ({q}+{q} qubits + ancilla)")
    amplitude_prep(c, a, list(range(1, q + 1)))
    amplitude_prep(c, b, list(range(q + 1, 2 * q + 1)))
    c.barrier("prepared")
    c.h(0)
    for k in range(q):
        c.cswap(0, 1 + k, q + 1 + k)
    c.h(0)
    return c


def swap_test(a: np.ndarray, b: np.ndarray, shots: int = 2000, seed: int = 0, fidelity: float = 1.0) -> dict:
    """Run the swap test on the simulator: exact P(0), overlap estimated from shots, and the exact overlap."""
    a, b = amplitude_vector(a), amplitude_vector(b, int(round(math.log2(len(amplitude_vector(a))))))
    q = int(round(math.log2(len(a))))
    psi = np.kron(np.kron(b, a), np.array([1.0, 0.0])).astype(complex)       # ancilla is the lowest qubit
    circ = sim.Circuit(2 * q + 1, "swap test")
    circ.h(0)
    for k in range(q):
        circ.cswap(0, 1 + k, q + 1 + k)
    circ.h(0)
    out = circ.run(psi)
    p0 = float(np.sum(np.abs(out[0::2]) ** 2))
    rng = np.random.default_rng(seed)
    p0n = fidelity * p0 + (1 - fidelity) * 0.5
    zeros = rng.binomial(int(shots), p0n)
    est = max(0.0, 2 * zeros / shots - 1)
    se = 2 * math.sqrt(p0n * (1 - p0n) / shots)
    return {"overlap_exact": float(np.dot(a, b) ** 2), "p0": p0, "overlap_estimate": float(est), "standard_error": se,
            "shots": int(shots), "qubits": 2 * q + 1}
