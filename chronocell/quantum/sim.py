"""
Statevector simulator for the quantum lab.

Convention (as Qiskit): qubit 0 is the least significant bit of a basis-state index, so the bitstring of index k
read right to left gives qubits 0, 1, 2, ...  Gates follow Qiskit's definitions:
RX(t) = exp(-i t X / 2), RY(t) = exp(-i t Y / 2), RZ(t) = exp(-i t Z / 2), RZZ(t) = exp(-i t Z(x)Z / 2).

Two paths:
* `Circuit` holds a gate list, runs it gate by gate (numpy), counts gates and depth, and writes OpenQASM 2.0 so the
  same circuit can be loaded by Qiskit or uploaded to a quantum computer by the user.
* QAOA on a diagonal cost (any QUBO or higher-order cost given as its 2^n energies) runs on a fast path: the cost
  layer is an element-wise phase and the mixer is RX on every qubit; on the GPU when PyTorch sees one.
  `qaoa_circuit` writes the identical algorithm as gates for a QUBO, and the tests check both paths agree.

Noise (labelled approximate): a global depolarising channel whose survival probability is the product of per-gate
fidelities, F = (1 - p1)^n1 (1 - p2)^n2, plus independent read-out bit flips. Typical rates of current
superconducting devices are the defaults. It is not an emulation of any particular machine, assumes all-to-all
connectivity (no SWAP overhead) and ignores coherent errors.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np

MAX_QUBITS = 24            # 2^24 amplitudes = 256 MB in complex128
NOISE_DEFAULT = {"p1": 1e-3, "p2": 1e-2, "readout": 2e-2}

try:                        # optional GPU path
    import torch
except Exception:           # pragma: no cover - torch is in requirements but stay importable without it
    torch = None


# ======================================================================================
# Gate matrices
# ======================================================================================
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
S = np.array([[1, 0], [0, 1j]], dtype=complex)
SDG = S.conj().T
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}


def rx(t: float) -> np.ndarray:
    c, s = math.cos(t / 2), math.sin(t / 2)
    return np.array([[c, -1j * s], [-1j * s, c]], dtype=complex)


def ry(t: float) -> np.ndarray:
    c, s = math.cos(t / 2), math.sin(t / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def rz(t: float) -> np.ndarray:
    return np.array([[np.exp(-0.5j * t), 0], [0, np.exp(0.5j * t)]], dtype=complex)


def rzz(t: float) -> np.ndarray:
    a, b = np.exp(-0.5j * t), np.exp(0.5j * t)
    return np.diag([a, b, b, a]).astype(complex)


CX = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]], dtype=complex)   # qubits (control, target)
CZ = np.diag([1, 1, 1, -1]).astype(complex)
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)
CSWAP = np.eye(8, dtype=complex)                                                          # qubits (control, a, b)
CSWAP[[3, 5]] = CSWAP[[5, 3]]

FIXED = {"h": H, "x": X, "y": Y, "z": Z, "s": S, "sdg": SDG, "cx": CX, "cz": CZ, "swap": SWAP, "cswap": CSWAP}
ROTATION = {"rx": rx, "ry": ry, "rz": rz, "rzz": rzz}


def apply_unitary(psi: np.ndarray, U: np.ndarray, qubits: tuple[int, ...] | list[int], n: int) -> np.ndarray:
    """Apply a k-qubit unitary. U is little-endian in `qubits` (qubits[0] is the lowest bit of U's index)."""
    k = len(qubits)
    t = psi.reshape((2,) * n)                       # tensor axis a holds qubit n-1-a
    axes = [n - 1 - q for q in reversed(qubits)]    # U's tensor axes run from qubits[k-1] down to qubits[0]
    out = np.tensordot(U.reshape((2,) * (2 * k)), t, axes=(list(range(k, 2 * k)), axes))
    return np.moveaxis(out, list(range(k)), axes).reshape(-1)


def zero_state(n: int) -> np.ndarray:
    psi = np.zeros(2 ** n, dtype=complex)
    psi[0] = 1.0
    return psi


def bits_of(index: int | np.ndarray, n: int) -> np.ndarray:
    """0/1 array of qubits 0..n-1 for one index (shape (n,)) or many (shape (m, n))."""
    idx = np.asarray(index, dtype=np.int64)
    return ((idx[..., None] >> np.arange(n)) & 1).astype(np.int8)


def bitstring(index: int, n: int) -> str:
    """Qiskit-style label: qubit n-1 on the left, qubit 0 on the right."""
    return format(int(index), f"0{n}b")


# ======================================================================================
# Circuits
# ======================================================================================
@dataclass
class Circuit:
    n: int
    name: str = "circuit"
    ops: list = field(default_factory=list)          # (gate, qubits tuple, params tuple)

    def add(self, gate: str, qubits, params=()) -> "Circuit":
        qubits = tuple(int(q) for q in (qubits if isinstance(qubits, (tuple, list)) else (qubits,)))
        if any(q < 0 or q >= self.n for q in qubits) or len(set(qubits)) != len(qubits):
            raise ValueError(f"bad qubits {qubits} for {gate} on {self.n} qubits")
        self.ops.append((gate, qubits, tuple(float(p) for p in params)))
        return self

    # convenience
    def h(self, q): return self.add("h", q)
    def x(self, q): return self.add("x", q)
    def s(self, q): return self.add("s", q)
    def sdg(self, q): return self.add("sdg", q)
    def rx(self, t, q): return self.add("rx", q, (t,))
    def ry(self, t, q): return self.add("ry", q, (t,))
    def rz(self, t, q): return self.add("rz", q, (t,))
    def rzz(self, t, a, b): return self.add("rzz", (a, b), (t,))
    def cx(self, c, t): return self.add("cx", (c, t))
    def cz(self, a, b): return self.add("cz", (a, b))
    def swap(self, a, b): return self.add("swap", (a, b))
    def cswap(self, c, a, b): return self.add("cswap", (c, a, b))

    def barrier(self, label: str = "") -> "Circuit":
        self.ops.append(("barrier", tuple(range(self.n)), (), label))
        return self

    def gates(self):
        for op in self.ops:
            if op[0] != "barrier":
                yield op[0], op[1], op[2]

    def run(self, psi: np.ndarray | None = None) -> np.ndarray:
        if self.n > MAX_QUBITS:
            raise ValueError(f"{self.n} qubits is more than this simulator handles ({MAX_QUBITS}).")
        psi = zero_state(self.n) if psi is None else np.asarray(psi, dtype=complex).copy()
        for g, qs, ps in self.gates():
            U = FIXED[g] if g in FIXED else ROTATION[g](*ps)
            psi = apply_unitary(psi, U, qs, self.n)
        return psi

    def counts(self) -> dict:
        """Gate counts and depth (gates on disjoint qubits share a layer); RZZ counts as 2 CX + 1 RZ when compiled."""
        per = {}
        layer = np.zeros(self.n, dtype=int)
        n1 = n2 = n3 = 0
        for g, qs, _ in self.gates():
            per[g] = per.get(g, 0) + 1
            d = int(layer[list(qs)].max()) + 1
            layer[list(qs)] = d
            if len(qs) == 1:
                n1 += 1
            elif len(qs) == 2:
                n2 += 1
            else:
                n3 += 1
        cx_equiv = per.get("cx", 0) + per.get("cz", 0) + 2 * per.get("rzz", 0) + 3 * per.get("swap", 0) + 8 * per.get("cswap", 0)
        return {"qubits": self.n, "gates": sum(per.values()), "by_gate": per, "one_qubit": n1, "two_qubit": n2,
                "three_qubit": n3, "depth": int(layer.max()) if self.ops else 0, "cx_equivalent": cx_equiv}

    def to_qasm(self, measure: bool = True) -> str:
        """OpenQASM 2.0 (qelib1.inc gates only), loadable by Qiskit and by IBM Quantum's composer."""
        lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"// {self.name} · written by ChronoCell-5D quantum lab",
                 f"qreg q[{self.n}];"]
        if measure:
            lines.append(f"creg c[{self.n}];")
        for op in self.ops:
            g, qs, ps = op[0], op[1], op[2]
            if g == "barrier":
                lines.append("barrier q;")
                continue
            par = f"({','.join(repr(float(p)) for p in ps)})" if ps else ""
            lines.append(f"{g}{par} " + ",".join(f"q[{q}]" for q in qs) + ";")
        if measure:
            lines.append("measure q -> c;")
        return "\n".join(lines) + "\n"


def fidelity_estimate(counts: dict, noise: dict | None = None) -> float:
    """Survival probability of the global depolarising model, from compiled gate counts."""
    nz = {**NOISE_DEFAULT, **(noise or {})}
    n1 = counts["one_qubit"] + counts.get("by_gate", {}).get("rzz", 0)     # each RZZ compiles to 2 CX + 1 RZ
    n2 = counts["cx_equivalent"]
    return float((1 - nz["p1"]) ** n1 * (1 - nz["p2"]) ** n2)


# ======================================================================================
# Sampling and noise
# ======================================================================================
def probabilities(psi) -> np.ndarray:
    if torch is not None and isinstance(psi, torch.Tensor):
        p = (psi.real ** 2 + psi.imag ** 2).double().cpu().numpy()
    else:
        p = np.abs(np.asarray(psi)) ** 2
    s = p.sum()
    return p / s if s > 0 else p


def sample(probs: np.ndarray, shots: int, rng: np.random.Generator, fidelity: float = 1.0, readout: float = 0.0,
           n: int | None = None) -> np.ndarray:
    """Shot outcomes (basis indices). With fidelity < 1 a shot comes from the uniform distribution with probability
    1 - F (global depolarising); with readout > 0 each measured bit flips independently."""
    probs = np.asarray(probs, float)
    dim = len(probs)
    n = int(round(math.log2(dim))) if n is None else n
    cdf = np.cumsum(probs)
    cdf /= cdf[-1]
    out = np.searchsorted(cdf, rng.random(shots), side="right").astype(np.int64)
    out = np.minimum(out, dim - 1)
    if fidelity < 1.0:
        bad = rng.random(shots) > fidelity
        out[bad] = rng.integers(0, dim, int(bad.sum()))
    if readout > 0:
        flips = (rng.random((shots, n)) < readout).astype(np.int64)
        out ^= (flips << np.arange(n)).sum(axis=1)
    return out


def noisy_probabilities(probs: np.ndarray, fidelity: float) -> np.ndarray:
    return fidelity * np.asarray(probs, float) + (1 - fidelity) / len(probs)


# ======================================================================================
# QAOA (fast path on a diagonal cost)
# ======================================================================================
def _device(prefer_gpu: bool = True):
    if torch is None:
        return None
    if prefer_gpu and torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


class _Engine:
    """Applies QAOA layers to a state; torch (GPU or CPU) when available, numpy otherwise."""

    def __init__(self, diag: np.ndarray, prefer_gpu: bool = True):
        self.n = int(round(math.log2(len(diag))))
        if 2 ** self.n != len(diag):
            raise ValueError("the cost vector length must be a power of two")
        if self.n > MAX_QUBITS:
            raise ValueError(f"{self.n} qubits is more than this simulator handles ({MAX_QUBITS}).")
        self.dev = _device(prefer_gpu) if self.n >= 10 else (torch.device("cpu") if torch is not None else None)
        if self.dev is not None:
            self.diag = torch.as_tensor(np.asarray(diag, float), dtype=torch.float64, device=self.dev)
        else:
            self.diag = np.asarray(diag, float)

    @property
    def device_name(self) -> str:
        if self.dev is None:
            return "numpy"
        return "GPU (" + torch.cuda.get_device_name(0) + ")" if self.dev.type == "cuda" else "CPU (torch)"

    def plus(self):
        dim = 2 ** self.n
        if self.dev is not None:
            return torch.full((dim,), 1 / math.sqrt(dim), dtype=torch.complex128, device=self.dev)
        return np.full(dim, 1 / math.sqrt(dim), dtype=complex)

    def cost(self, psi, gamma: float):
        if self.dev is not None:
            return psi * torch.polar(torch.ones_like(self.diag), -gamma * self.diag)
        return psi * np.exp(-1j * gamma * self.diag)

    def mixer(self, psi, beta: float):
        c, s = math.cos(beta), math.sin(beta)
        n = self.n
        for q in range(n):
            v = psi.reshape(2 ** (n - q - 1), 2, 2 ** q)
            a, b = v[:, 0, :], v[:, 1, :]
            if self.dev is not None:
                psi = torch.stack((c * a - 1j * s * b, -1j * s * a + c * b), dim=1).reshape(-1)
            else:
                psi = np.stack((c * a - 1j * s * b, -1j * s * a + c * b), axis=1).reshape(-1)
        return psi

    def state(self, gammas, betas):
        psi = self.plus()
        for g, b in zip(gammas, betas):
            psi = self.mixer(self.cost(psi, g), b)
        return psi

    def probs(self, psi):
        if self.dev is not None:
            return (psi.real ** 2 + psi.imag ** 2)
        return np.abs(psi) ** 2

    def expectation(self, probs) -> float:
        if self.dev is not None:
            return float(torch.dot(probs, self.diag))
        return float(probs @ self.diag)

    def cvar(self, probs, alpha: float) -> float:
        """Mean energy of the lowest-energy alpha fraction of the distribution (Barkoutsos et al., Quantum 2020)."""
        if not hasattr(self, "_order"):
            self._order = torch.argsort(self.diag, stable=True) if self.dev is not None else np.argsort(self.diag, kind="stable")
            self._sorted = self.diag[self._order]
        ps = probs[self._order]
        if self.dev is not None:
            cum = torch.cumsum(ps, 0)
            k = min(int(torch.searchsorted(cum, torch.tensor([alpha], dtype=cum.dtype, device=cum.device))[0]), len(ps) - 1)
            head = float(torch.dot(ps[:k], self._sorted[:k])) if k > 0 else 0.0
            before = float(cum[k - 1]) if k > 0 else 0.0
            return (head + (alpha - before) * float(self._sorted[k])) / alpha
        cum = np.cumsum(ps)
        k = min(int(np.searchsorted(cum, alpha)), len(ps) - 1)
        head = float(ps[:k] @ self._sorted[:k]) if k > 0 else 0.0
        before = float(cum[k - 1]) if k > 0 else 0.0
        return (head + (alpha - before) * float(self._sorted[k])) / alpha


@dataclass
class QAOAResult:
    p: int
    gammas: np.ndarray
    betas: np.ndarray
    objective: str
    expectation: float                  # <E> of the final state (original energy units)
    p_optimal: float                    # probability of measuring a minimum-energy bitstring
    best: int                           # best sampled index
    best_energy: float
    hit: bool                           # a sampled bitstring has the minimum energy
    samples: np.ndarray
    probs: np.ndarray
    history: list
    evaluations: int
    seconds: float
    device: str
    shots: int
    noise_fidelity: float = 1.0
    readout: float = 0.0
    energy_min: float = float("nan")
    energy_max: float = float("nan")

    @property
    def approximation_ratio(self) -> float:
        """(E_max - E_best) / (E_max - E_min): 1 = optimal."""
        span = self.energy_max - self.energy_min
        return float((self.energy_max - self.best_energy) / span) if span > 0 else 1.0


def _ramp(p: int, dt: float = 0.75) -> tuple[np.ndarray, np.ndarray]:
    """Linear-ramp start (annealing-like): gamma rises, beta falls."""
    k = (np.arange(p) + 0.5) / p
    return dt * k, dt * (1 - k)


def _interp(params: np.ndarray, p: int) -> np.ndarray:
    """INTERP growth from depth p to p+1 (Zhou et al., PRX 2020)."""
    out = np.zeros(p + 1)
    for i in range(p + 1):
        a = params[i - 1] if i - 1 >= 0 else 0.0
        b = params[i] if i < p else 0.0
        out[i] = (i / p) * a + ((p - i) / p) * b
    return out


def qaoa(diag: np.ndarray, p: int = 3, shots: int = 2048, objective: str = "cvar", alpha: float = 0.1,
         maxiter: int = 120, seed: int = 0, prefer_gpu: bool = True, noise: dict | None = None,
         fidelity: float | None = None, init: tuple[np.ndarray, np.ndarray] | None = None,
         grow: bool = True, callback=None) -> QAOAResult:
    """Optimise a depth-p QAOA on the diagonal cost `diag` (2^n energies, index = bitstring of qubits) and sample it.

    The cost is rescaled to zero mean and unit standard deviation for the angles (energies are reported unscaled).
    objective: "expectation" (standard) or "cvar" (lowest alpha fraction). Depth grows 1..p with INTERP when `grow`.
    noise: None (ideal) or rates {"p1", "p2", "readout"}; `fidelity` overrides the circuit-derived survival.
    """
    from scipy.optimize import minimize
    t0 = time.perf_counter()
    diag = np.asarray(diag, float)
    mu, sd = float(diag.mean()), float(diag.std()) or 1.0
    z = (diag - mu) / sd
    eng = _Engine(z, prefer_gpu)
    history: list[float] = []
    evals = 0

    def f(x):
        nonlocal evals
        k = len(x) // 2
        pr = eng.probs(eng.state(x[:k], x[k:]))
        v = eng.cvar(pr, alpha) if objective == "cvar" else eng.expectation(pr)
        evals += 1
        history.append(v * sd + mu)
        if callback is not None:
            callback(evals, history[-1])
        return v

    if init is not None:
        g0, b0 = np.asarray(init[0], float), np.asarray(init[1], float)
        depths = [len(g0)] if len(g0) == p else list(range(len(g0), p + 1))
    elif grow:
        g0, b0 = _ramp(1)
        depths = list(range(1, p + 1))
    else:
        g0, b0 = _ramp(p)
        depths = [p]
    x = np.concatenate([g0, b0])
    for d in depths:
        k = len(x) // 2
        if k < d:                                   # grow one layer
            x = np.concatenate([_interp(x[:k], k), _interp(x[k:], k)])
        r = minimize(f, x, method="COBYLA", options={"maxiter": int(maxiter), "rhobeg": 0.3})
        x = r.x
    k = len(x) // 2
    gam, bet = x[:k], x[k:]
    psi = eng.state(gam, bet)
    probs = probabilities(psi)
    emin = float(diag.min())
    opt = np.isclose(diag, emin, rtol=0, atol=1e-9 * max(1.0, abs(emin)))
    rng = np.random.default_rng(seed)
    nz = {**NOISE_DEFAULT, **noise} if noise else None
    F = 1.0
    if fidelity is not None:
        F = float(fidelity)
    elif nz is not None:
        F = fidelity_estimate(qaoa_counts(diag_terms_count(diag), eng.n, k), nz)
    ro = float(nz["readout"]) if nz is not None else 0.0
    smp = sample(probs, shots, rng, F, ro, eng.n)
    e_s = diag[smp]
    b = int(smp[int(np.argmin(e_s))])
    return QAOAResult(p=k, gammas=gam / sd, betas=bet, objective=objective, expectation=float(probs @ diag),
                      p_optimal=float(probs[opt].sum()), best=b, best_energy=float(diag[b]),
                      hit=bool(opt[smp].any()), samples=smp, probs=probs, history=history, evaluations=evals,
                      seconds=time.perf_counter() - t0, device=eng.device_name, shots=int(shots), noise_fidelity=F,
                      readout=ro, energy_min=emin, energy_max=float(diag.max()))


def diag_terms_count(diag: np.ndarray, tol: float = 1e-10) -> dict:
    """Number of Pauli-Z product terms of each order in a diagonal cost (fast Walsh-Hadamard transform)."""
    w = walsh_coefficients(diag)
    n = int(round(math.log2(len(diag))))
    order = np.array([bin(i).count("1") for i in range(len(diag))])
    big = np.abs(w) > tol * max(1.0, float(np.abs(w).max()))
    return {int(o): int(np.sum(big & (order == o))) for o in range(1, n + 1) if np.any(big & (order == o))}


def walsh_coefficients(diag: np.ndarray) -> np.ndarray:
    """Coefficients w_S of diag = sum_S w_S prod_{i in S} Z_i (S encoded as a bitmask index)."""
    return _fwht(np.asarray(diag, float)) / len(diag)


def _fwht(v: np.ndarray) -> np.ndarray:
    a = v.copy()
    h = 1
    n = len(a)
    while h < n:
        a = a.reshape(n // (2 * h), 2, h)
        x, y = a[:, 0, :].copy(), a[:, 1, :].copy()
        a[:, 0, :], a[:, 1, :] = x + y, x - y
        a = a.reshape(n)
        h *= 2
    return a


def qaoa_counts(terms: dict, n: int, p: int) -> dict:
    """Compiled gate counts of a depth-p QAOA whose cost has `terms` Pauli-Z products by order (k: count): an
    order-k term compiles to 2(k-1) CX + 1 RZ; the mixer is n RX per layer; n H to start. Depth is not estimated
    here (it depends on scheduling); a QUBO's gate-level circuit (qaoa_circuit) reports its own."""
    one = n + p * (n + sum(terms.values()))
    cx = p * sum(2 * (k - 1) * c for k, c in terms.items())
    return {"qubits": n, "one_qubit": one, "two_qubit": cx, "three_qubit": 0, "cx_equivalent": cx,
            "by_gate": {"h": n, "rz": p * sum(terms.values()), "rx": p * n, "cx": cx}, "gates": one + cx,
            "depth": None}


def qaoa_circuit(h: np.ndarray, J: dict, gammas, betas, name: str = "QAOA") -> Circuit:
    """Gate-level QAOA for an Ising cost sum_i h_i Z_i + sum_{i<j} J_ij Z_i Z_j (constant dropped):
    exp(-i g H_C) = prod RZ(2 g h_i) prod RZZ(2 g J_ij); mixer RX(2 b) on every qubit."""
    n = len(h)
    c = Circuit(n, name)
    for q in range(n):
        c.h(q)
    for layer, (g, b) in enumerate(zip(gammas, betas)):
        c.barrier(f"layer {layer + 1}")
        for (i, j), v in sorted(J.items()):
            if abs(v) > 1e-12:
                c.rzz(2 * g * v, i, j)
        for i in range(n):
            if abs(h[i]) > 1e-12:
                c.rz(2 * g * h[i], i)
        for i in range(n):
            c.rx(2 * b, i)
    return c


def state_fidelity(a: np.ndarray, b: np.ndarray) -> float:
    """|<a|b>|^2 for normalised states (insensitive to global phase)."""
    a = np.asarray(a, complex)
    b = np.asarray(b, complex)
    return float(abs(np.vdot(a / np.linalg.norm(a), b / np.linalg.norm(b))) ** 2)


def memory_bytes(n: int) -> int:
    """Memory of one complex128 statevector of n qubits."""
    return 16 * 2 ** n


# ======================================================================================
# Pauli operators (small systems)
# ======================================================================================
def pauli_matrix(label: str) -> np.ndarray:
    """Matrix of a Pauli string; label[0] acts on qubit n-1 and label[-1] on qubit 0 (Qiskit order)."""
    m = np.array([[1.0 + 0j]])
    for ch in label:
        m = np.kron(m, PAULI[ch])
    return m


def pauli_decompose(M: np.ndarray, tol: float = 1e-10) -> dict[str, float]:
    """Hermitian 2^n x 2^n matrix -> {Pauli label: real coefficient}."""
    import itertools
    dim = M.shape[0]
    n = int(round(math.log2(dim)))
    out = {}
    for lab in itertools.product("IXYZ", repeat=n):
        lab = "".join(lab)
        c = np.trace(pauli_matrix(lab) @ M) / dim
        if abs(c) > tol:
            out[lab] = float(np.real(c))
    return out


def pauli_exponential(circ: Circuit, label: str, theta: float) -> Circuit:
    """Append exp(-i theta P) for a Pauli string P (Qiskit order): basis change, CX ladder, RZ(2 theta), undo."""
    n = circ.n
    act = [(n - 1 - k, ch) for k, ch in enumerate(label) if ch != "I"]
    if not act:
        return circ
    act.sort()
    for q, ch in act:
        if ch == "X":
            circ.h(q)
        elif ch == "Y":
            circ.sdg(q)
            circ.h(q)
    qs = [q for q, _ in act]
    for a, b in zip(qs[:-1], qs[1:]):
        circ.cx(a, b)
    circ.rz(2 * theta, qs[-1])
    for a, b in reversed(list(zip(qs[:-1], qs[1:]))):
        circ.cx(a, b)
    for q, ch in act:
        if ch == "X":
            circ.h(q)
        elif ch == "Y":
            circ.h(q)
            circ.s(q)
    return circ
