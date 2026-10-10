"""
Re-check every accuracy test from its saved result file, without re-running any test.

    python motion/explainer/recheck.py            # prints the audit (used, recorded, in the explainer video)
    python motion/explainer/recheck.py --pace     # the same, a line at a time (for recording)

Each test was pre-registered (rule frozen in validation/frozen.py before its data were read) and run once, so it is
never run again. What can be checked again is the arithmetic: for the tests whose result files keep the per-case
numbers, this script recomputes the score from those numbers (energy errors from energies, F1 from loop counts,
docking success from pose RMSDs, ...), re-applies the frozen rule, and checks the verdict against the Scoreboard.
"""

from __future__ import annotations

import json
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
VAL = ROOT / "validation"
PACE = "--pace" in sys.argv
OK, BAD = "ok", "MISMATCH"
mismatches = 0


def out(s: str = "", dt: float = 0.035) -> None:
    print(s, flush=True)
    if PACE:
        time.sleep(dt)


def load(name: str) -> dict:
    return json.loads((VAL / name).read_text(encoding="utf-8"))


def check(label: str, recomputed, stored) -> str:
    global mismatches
    same = (abs(recomputed - stored) < 1e-6) if isinstance(recomputed, float) else recomputed == stored
    if not same:
        mismatches += 1
    return f"{label} {OK if same else BAD}"


def head(title: str, rule: str) -> None:
    out()
    out(f"== {title}", 0.12)
    out(f"   rule (frozen before the data were read): {rule}", 0.08)


def f1(m: dict) -> float:
    p, r = m["matched"] / m["calls"], m["matched"] / m["reference"]
    return 2 * p * r / (p + r)


def main() -> None:
    from ui import scoreboard as SB
    out("ChronoCell-5D · re-checking the accuracy tests from their saved result files", 0.2)
    out("rules: validation/frozen.py · results: validation/results_*.json · no test is re-run (each ran once)", 0.2)

    # ------------------------------------------------------------------ the Scoreboard, one line per test
    df = SB.rows()
    out()
    out(f"Scoreboard: {len(df)} tests on held-out real data", 0.15)
    tag = {"pass": "PASS ", "fail": "FAIL ", "other": "OTHER"}
    for i, r in enumerate(df.itertuples(), 1):
        name = r.test if len(r.test) <= 58 else r.test[:57] + "…"
        meas = r.measured if len(r.measured) <= 34 else r.measured[:33] + "…"
        out(f"  {i:>2}  {tag[r.status]}  {name:<58}  {meas}", 0.06)

    # ------------------------------------------------------------------ recomputations
    g6b = load("results_gate6b.json")
    head("Gate 6b · DNA loops vs ENCODE reference loops", "ChronoCell's F1 must beat both published tools on each new cell type")
    ok6 = True
    for cell, v in g6b["sets"].items():
        fc, fs, fm = f1(v["chronocell"]), f1(v["chromosight"]), f1(v["mustache"])
        beats = fc > max(fs, fm)
        ok6 &= beats
        c = v["chronocell"]
        out(f"   {cell:<5} {c['matched']} of {c['reference']} reference loops found ({c['calls']} calls) -> F1 {fc:.3f}"
            f" · chromosight {fs:.3f} · Mustache {fm:.3f}   {'beats both' if beats else 'does not'}", 0.1)
        out(f"         {check('stored F1', round(fc, 6), round(v['chronocell']['f1'], 6))}", 0.05)
    out(f"   verdict: {'PASS' if ok6 else 'FAIL'}   {check('Scoreboard', 'pass' if ok6 else 'fail', _status(df, 'Gate 6b'))}", 0.15)

    r2 = load("results_round2.json")
    q = r2["q6d"]
    tol = q["rule"]["tolerance_mEh"]
    head("Gate Q6d · stretched molecules on a simulated quantum computer (ADAPT-VQE)",
         f"every |error| <= {tol} mHa (chemical accuracy) against the exact lowest singlet")
    errs = []
    for x in q["rows"]:
        e = abs(x["error_vs_singlet_mEh"])
        errs.append(e)
        out(f"   {x['molecule']:<5} bond x{x['bond_scale']:<4} |error| = {e:8.4f} mHa   {'within' if e <= tol else 'OUTSIDE'}", 0.06)
    within = sum(e <= tol for e in errs)
    out(f"   {within} of {len(errs)} within, worst {max(errs):.3f} mHa -> {'PASS' if within == len(errs) else 'FAIL'}   "
        f"{check('stored', (within, round(max(errs), 6)), (q['within'], round(abs(q['max_abs_error_mEh']), 6)))}   "
        f"{check('Scoreboard', 'pass' if within == len(errs) else 'fail', _status(df, 'Gate Q6d'))}", 0.15)

    mit = load("results_mitigation.json")
    R = mit["rule"]
    head("Gate Q9 · error mitigation on a simulated noisy chip",
         f"after mitigation, at least {R['min_fraction']:.0%} of {len(R['cases'])} molecules within {R['tolerance_mEh']} mHa")
    noisy, fixed = [], []
    for x in mit["rows"]:
        m = x["methods"][mit["method_key"]]
        n, f = abs(m["noisy_mEh"]), abs(m["mitigated_mEh"])
        noisy.append(n)
        fixed.append(f)
        out(f"   {x['molecule']:<5} bond x{x['bond_scale']:<4} noisy {n:7.3f} mHa -> mitigated {f:6.3f} mHa   "
            f"{'within' if f <= R['tolerance_mEh'] else 'OUTSIDE'}", 0.06)
    w = sum(f <= R["tolerance_mEh"] for f in fixed)
    red = statistics.median(n / max(f, 1e-6) for n, f in zip(noisy, fixed))
    ok9 = w >= R["min_fraction"] * len(fixed)
    out(f"   {w} of {len(fixed)} within · median {statistics.median(fixed):.2f} mHa · median {red:.0f}x smaller -> "
        f"{'PASS' if ok9 else 'FAIL'}   {check('stored', (w, round(red, 6)), (mit['within'], round(mit['median_reduction'], 6)))}   "
        f"{check('Scoreboard', 'pass' if ok9 else 'fail', _status(df, 'Gate Q9'))}", 0.15)

    q7 = r2["q7b"]
    head("Gate Q7b · quantum docking on 85 new protein-drug complexes",
         "a pose counts if it lies within 2 Å of the crystal; the quantum route must dock at least as often as random search")
    use = [r for r in q7["per_complex"] if r.get("usable")]
    sq = sum(r["qaoa"]["rmsd"] <= 2.0 for r in use)
    sr = sum(r["random"]["rmsd"] <= 2.0 for r in use)
    for r in use[:6]:
        out(f"   {r['cid']:<10} quantum pose {r['qaoa']['rmsd']:6.2f} Å   random-search pose {r['random']['rmsd']:6.2f} Å", 0.05)
    out(f"   … {len(use) - 6} more complexes", 0.05)
    ok7 = q7["qaoa_hit_rate"] >= q7["rule"]["min_hit_rate"] and sq >= sr and sq / len(use) >= q7["rule"]["min_success"]
    out(f"   quantum route {sq}/{len(use)} = {100 * sq / len(use):.1f} %   random search {sr}/{len(use)} = {100 * sr / len(use):.1f} %"
        f" -> {'PASS' if ok7 else 'FAIL'}   {check('stored', round(sq / len(use), 6), round(q7['success_qaoa'], 6))}   "
        f"{check('Scoreboard', 'pass' if ok7 else 'fail', _status(df, 'Gate Q7b'))}", 0.15)

    g = load("results_gateq.json")
    head("Gate Q3 · small molecules, exact energies",
         "|VQE energy − exact (FCI) energy| <= 1.6 mHa for every geometry, with the UCCSD circuit")
    rows3 = [x for x in g["q3"]["rows"] if x.get("ansatz", "UCCSD") == "UCCSD"]
    e3 = [abs(x["e_vqe"] - x["e_fci"]) * 1000 for x in rows3]
    for x, e in zip(rows3, e3):
        out(f"   {x['molecule']:<5} R = {x['R_angstrom']:<6} Å   E_vqe {x['e_vqe']:.6f}   E_exact {x['e_fci']:.6f}   error {e:.1e} mHa", 0.05)
    ok3 = max(e3) <= 1.6
    out(f"   worst {max(e3):.1e} mHa -> {'PASS' if ok3 else 'FAIL'}   "
        f"{check('stored', round(max(e3), 15), round(g['q3']['max_abs_error_mha'], 15))}   "
        f"{check('Scoreboard', 'pass' if ok3 else 'fail', _status(df, 'Gate Q3'))}", 0.15)
    hw = [x for x in g["q3"]["rows"] if x.get("ansatz", "UCCSD") != "UCCSD"]
    if hw:
        miss = sum(abs(x["e_vqe"] - x["e_fci"]) * 1000 > 1.6 for x in hw)
        out(f"   for comparison (not part of the test): a simpler hardware-efficient circuit misses on {miss} of {len(hw)} geometries", 0.1)

    adm = load("results_admet.json")
    head("Gate Q8 · 21-property drug safety profile (ADMET)",
         f"the quantum model within {adm['rule']['margin']} of a classical model on at least {adm['rule']['min_endpoints']} of 21 endpoints")
    met = sum(bool(v["pass"]) for v in adm["endpoints"].values())
    ok8 = met >= adm["rule"]["min_endpoints"]
    out(f"   {met} of {len(adm['endpoints'])} endpoints met -> {'PASS' if ok8 else 'FAIL'}   "
        f"{check('stored', met, adm['passed'])}   {check('Scoreboard', 'pass' if ok8 else 'fail', _status(df, 'Gate Q8'))}", 0.15)

    if (VAL / "results_admet_q8b.json").exists():
        q8b = load("results_admet_q8b.json")
        R = q8b["rule"]
        head(f"Gate Q8b · drug properties with the encoding fixed, {len(q8b['endpoints'])} new endpoints",
             f"quantum >= classical - {R['margin']} and its 95 % interval above chance, on at least {R['min_endpoints']}")
        met = 0
        for n, e in q8b["endpoints"].items():
            m, lo = e["metric"], e["ci95"]["quantum"][0]
            ok = m["quantum"] >= m["rbf8"] - R["margin"] and lo > (0.5 if e["task"] == "cls" else 0.0)
            met += ok
            check("endpoint", ok, bool(e["pass"]))
        ok8b = met >= R["min_endpoints"]
        out(f"   each endpoint's verdict recomputed from its scores and interval: {met} of {len(q8b['endpoints'])} met -> "
            f"{'PASS' if ok8b else 'FAIL'}   {check('stored', met, q8b['passed'])}   "
            f"{check('Scoreboard', 'pass' if ok8b else 'fail', _status(df, 'Gate Q8b'))}", 0.15)

    q1 = g["q1"]
    head("Gate Q1 · quantum optimiser finds the best answer", f"hit rate >= {g['rule']['q1_min_hit_rate']:.0%} of windows")
    ok1 = q1["pooled_hit_rate"] >= g["rule"]["q1_min_hit_rate"]
    out(f"   pooled hit rate {q1['pooled_hit_rate']:.1%} -> {'PASS' if ok1 else 'FAIL'}   "
        f"{check('Scoreboard', 'pass' if ok1 else 'fail', _status(df, 'Gate Q1'))}", 0.15)

    # ------------------------------------------------------------------ summary
    p, f, o = (int((df["status"] == s).sum()) for s in ("pass", "fail", "other"))
    out()
    out(f"Summary: {len(df)} tests · {p} passed · {f} failed (kept on the record) · {o} mixed, blocked or baseline", 0.2)
    out(f"Recomputed scores and verdicts: {'all match the saved results and the Scoreboard' if not mismatches else f'{mismatches} MISMATCHES'}", 0.2)


def _status(df, prefix: str) -> str:
    hit = df[df["test"].str.startswith(prefix + ":") | df["test"].str.startswith(prefix + " ")]
    return hit["status"].iloc[0] if len(hit) else "missing"


if __name__ == "__main__":
    main()
