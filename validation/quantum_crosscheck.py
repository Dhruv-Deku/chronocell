"""
Independent check of ChronoCell's statevector simulator (chronocell/quantum/sim.py) against Qiskit.

Every kind of circuit the quantum lab builds (QAOA on a real-data domain QUBO, UCCSD and hardware-efficient VQE,
the ZZ feature map, amplitude preparation + swap test, and a random circuit using every gate) is written as
OpenQASM 2.0, loaded by Qiskit in the isolated environment .chronocell_cache/quantum-venv (qiskit 2.2.1; see
requirements-quantum.txt) and simulated there; the two final states are compared (fidelity |<a|b>|^2, insensitive
to the global phase that RZ and the phase gate differ by). Not a gate: a software check, written to
results_quantum_crosscheck.json.

    python validation/quantum_crosscheck.py
"""

from __future__ import annotations

import datetime as dt
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

from chronocell.quantum import chem as CH, kernels as KQ, problems as PR, qubo as QB, sim  # noqa: E402

VENV_PY = ROOT.parent / ".chronocell_cache" / "quantum-venv" / "Scripts" / "python.exe"
RUNNER = r'''
import json, sys, numpy as np, qiskit
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
out = {"qiskit": qiskit.__version__, "states": {}}
for name, path in json.loads(sys.argv[1]).items():
    qc = QuantumCircuit.from_qasm_file(path)
    sv = Statevector(qc).data
    np.save(path.replace(".qasm", "_qiskit.npy"), sv)
    out["states"][name] = {"qubits": qc.num_qubits, "ops": int(sum(qc.count_ops().values()))}
print(json.dumps(out))
'''


def circuits() -> dict[str, sim.Circuit]:
    rng = np.random.default_rng(0)
    out = {}
    # QAOA on a real-data domain QUBO (GM12878 practice map, a 12-bin window -> 11 qubits)
    import quantum_gateq as QG
    C = PR.coarsen(QG.region_counts("gm12878", "chr1", 100_000_000), 4)
    q, _ = PR.tad_qubo(C, (6, 18), 1.5, 3, 6, boundary_cost=1.5)
    h, J, _ = q.ising()
    r = sim.qaoa(q.energies(), p=2, shots=256, maxiter=40, prefer_gpu=False)
    out["QAOA domain walls (11 qubits, p=2)"] = sim.qaoa_circuit(h, J, r.gammas, r.betas)
    m = CH.molecule("H2", 0.7414)
    out["UCCSD-VQE H2"] = CH.uccsd_circuit(m, CH.vqe(m).thetas)
    out["Hardware-efficient VQE HeH+"] = CH.hea_circuit(4, rng.normal(0, 1, 12), 2, 3)
    out["ZZ feature map (5 qubits, 2 reps)"] = KQ.feature_map_circuit(rng.random(5) * 3, 0.7, 2)
    a, b = KQ.amplitude_vector(rng.random(4)), KQ.amplitude_vector(rng.random(4))
    out["Swap test with state preparation (5 qubits)"] = KQ.swap_test_circuit(a, b)
    c = sim.Circuit(5, "random circuit, every gate")
    gates1 = ["h", "x", "s", "sdg"]
    for _ in range(60):
        k = rng.integers(0, 8)
        qs = rng.choice(5, 3, replace=False)
        if k < 4:
            c.add(gates1[k], int(qs[0]))
        elif k == 4:
            c.add(["rx", "ry", "rz"][rng.integers(0, 3)], int(qs[0]), (float(rng.normal()),))
        elif k == 5:
            c.add(["cx", "cz", "swap"][rng.integers(0, 3)], (int(qs[0]), int(qs[1])))
        elif k == 6:
            c.rzz(float(rng.normal()), int(qs[0]), int(qs[1]))
        else:
            c.cswap(int(qs[0]), int(qs[1]), int(qs[2]))
    out[c.name] = c
    return out


def main() -> None:
    if not VENV_PY.exists():
        raise SystemExit(f"{VENV_PY} not found: create it with requirements-quantum.txt")
    circs = circuits()
    tmp = Path(tempfile.mkdtemp(prefix="chronocell_qx_"))
    paths, mine = {}, {}
    for i, (name, c) in enumerate(circs.items()):
        p = tmp / f"c{i}.qasm"
        p.write_text(c.to_qasm(measure=False), encoding="utf-8")
        paths[name] = str(p)
        mine[name] = c.run()
    res = subprocess.run([str(VENV_PY), "-c", RUNNER, json.dumps(paths)], capture_output=True, text=True, check=True)
    q = json.loads(res.stdout.strip().splitlines()[-1])
    rows = []
    for name, c in circs.items():
        theirs = np.load(paths[name].replace(".qasm", "_qiskit.npy"))
        cnt = c.counts()
        rows.append({"circuit": name, "qubits": c.n, "gates": cnt["gates"], "qiskit_ops": q["states"][name]["ops"],
                     "fidelity": sim.state_fidelity(mine[name], theirs),
                     "max_probability_difference": float(np.max(np.abs(np.abs(mine[name]) ** 2 - np.abs(theirs) ** 2)))})
        print(f"{name}: fidelity {rows[-1]['fidelity']:.15f}", flush=True)
    out = {"made": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "qiskit": q["qiskit"],
           "circuits": rows, "verdict": "agree" if all(r["fidelity"] > 1 - 1e-9 for r in rows) else "DISAGREE"}
    (ROOT / "results_quantum_crosscheck.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(out["verdict"], flush=True)


if __name__ == "__main__":
    main()
