"""
Gate Q9 (October 2026): can error mitigation recover chemical accuracy from a noisy quantum computer?

The UCCSD-VQE circuit of a molecule's 2-electron, 2-orbital active space (4 qubits; molecules.py, compiled to Pauli
exponentials: about 200 gates, 64 CNOT-equivalents; noisy.vqe_circuit) runs on a density-matrix simulator with local
depolarizing noise after every gate (chronocell/quantum/noisy.py). Its noisy energy is compared with the noise-free
energy of the same circuit, before and after mitigation: zero-noise extrapolation (unitary folding at noise scales
1, 3, 5; Richardson, linear or exponential extrapolation) with or without symmetry verification (discarding results
with the wrong electron number). Nothing about the ideal energy is used by the mitigation.

Practice: 15 cases (library molecules at equilibrium and the six Q6 geometries); noise as today's best
superconducting devices (two-qubit error 3e-3, one-qubit 3e-4) and the lab's pessimistic setting (1e-2, 1e-3).
Test (run once; rule in frozen.QUANTUM_MITIGATION): 16 molecule / bond-length cases never run.

    python validation/mitigation_gate.py --practice
    python validation/mitigation_gate.py --test
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

from chronocell.quantum import molecules as M, noisy as N      # noqa: E402

PRACTICE = [(m, 1.0) for m in ("LiH", "HF", "H2O", "NH3", "CH4", "N2", "CO", "HCN", "CH2O")] + \
    [("H2O", 1.5), ("NH3", 1.5), ("N2", 1.5), ("HF", 2.0), ("CH2O", 1.3), ("HCN", 1.3)]
NOISE = {"device": {"p1": 3e-4, "p2": 3e-3}, "pessimistic": {"p1": 1e-3, "p2": 1e-2}}
METHODS = [(m, sv, sc) for m in ("richardson", "linear", "exp") for sv in (False, True) for sc in ((1, 3, 5), (1, 3))]


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def case(name: str, scale: float):
    _, fn, ch = M.LIBRARY[name]
    ints = M.integrals(fn(scale), ch)
    qm = M.qubit_hamiltonian(M.active_space(ints, M.rhf(ints), 2, 2))
    r = M.vqe_uccsd(qm)
    c = N.vqe_circuit(qm, r.generators, r.thetas, f"{name} {scale}x")
    keep = np.zeros(2 ** qm.n_qubits)
    keep[M.sector(qm)] = 1.0
    return c, qm.H.toarray(), keep, r


def evaluate(cases: list, noise: dict, methods: list) -> list[dict]:
    rows = []
    for name, scale in cases:
        c, H, keep, r = case(name, scale)
        row = {"molecule": name, "bond_scale": scale, "ideal": r.energy, "fci": r.fci, "gates": c.counts()["gates"],
               "cx": c.counts()["cx_equivalent"], "methods": {}}
        for meth, sv, sc in methods:
            z = N.zne(c, H, noise, scales=sc, method=meth, ideal=r.energy, keep=keep if sv else None)
            row["methods"][f"{meth}{'+sv' if sv else ''}{''.join(map(str, sc))}"] = {
                "noisy_mEh": 1e3 * z.error_noisy, "mitigated_mEh": 1e3 * z.error_mitigated, "values": z.values}
        rows.append(row)
    return rows


def summarise(rows: list[dict]) -> dict:
    out = {}
    for k in rows[0]["methods"]:
        e = np.array([abs(r["methods"][k]["mitigated_mEh"]) for r in rows])
        n = np.array([abs(r["methods"][k]["noisy_mEh"]) for r in rows])
        out[k] = {"within_chemical_accuracy": int((e <= 1.6).sum()), "cases": len(rows), "median_abs_mEh": float(np.median(e)),
                  "max_abs_mEh": float(e.max()), "median_noisy_mEh": float(np.median(n)),
                  "median_reduction": float(np.median(n / np.maximum(e, 1e-6)))}
    return out


def practice() -> dict:
    out = {"made": _now(), "cases": PRACTICE, "noise": NOISE, "levels": {}}
    for lvl, noise in NOISE.items():
        rows = evaluate(PRACTICE, noise, METHODS)
        out["levels"][lvl] = {"rows": rows, "summary": summarise(rows)}
        for k, v in out["levels"][lvl]["summary"].items():
            print(lvl, k, v["within_chemical_accuracy"], "/", v["cases"], "median", round(v["median_abs_mEh"], 3), flush=True)
    s = out["levels"]["device"]["summary"]
    best = max(s, key=lambda k: (s[k]["within_chemical_accuracy"], -s[k]["median_abs_mEh"]))
    out["choice"] = best
    (ROOT / "results_mitigation_practice.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("choice", best, s[best], flush=True)
    return out


def test() -> dict:
    import frozen as F
    R = F.QUANTUM_MITIGATION
    path = ROOT / "results_mitigation.json"
    if path.exists():
        sys.exit("results_mitigation.json exists: Gate Q9 runs once.")
    meth = (R["method"], R["symmetry_verification"], tuple(R["scales"]))
    rows = evaluate([tuple(c) for c in R["cases"]], R["noise"], [meth])
    key = next(iter(rows[0]["methods"]))
    e = np.array([abs(r["methods"][key]["mitigated_mEh"]) for r in rows])
    n = np.array([abs(r["methods"][key]["noisy_mEh"]) for r in rows])
    within = int((e <= R["tolerance_mEh"]).sum())
    out = {"made": _now(), "rule": R, "rows": rows, "method_key": key, "within": within, "cases": len(rows),
           "median_abs_mEh": float(np.median(e)), "median_noisy_mEh": float(np.median(n)),
           "median_reduction": float(np.median(n / np.maximum(e, 1e-6))),
           "pass": bool(within >= R["min_fraction"] * len(rows))}
    path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("Gate Q9:", "pass" if out["pass"] else "fail", within, "of", len(rows), flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    practice() if a.practice else test()


if __name__ == "__main__":
    main()
