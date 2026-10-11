"""
Gate 2e (Phase A3): a per-pair reliability score, tested on data no reliability test has touched.

Gate 2c's score (the model's misfit to the pair's own input) ranked held-out pair error weakly with
imaging-derived input and not at all with Hi-C input. Phase A3 candidates, all from the input alone:

  boot_sd    minus the spread (SD) of log median distance over K refits of the model on resampled input
             (imaging: half A's copies drawn with replacement; Hi-C: every count redrawn from a Poisson law)
  evidence   contacts observed around the pair: log(1 + sum of counts in its 3 x 3 neighbourhood)
  combined   the average of the two, as within-stratum ranks
  misfit     Gate 2c's score, for reference

Error and scoring as Gate 2c (validation/reliability.py): the scale-free pair error, Spearman of score
vs minus error within 10 separation strata, averaged. Genome-scale sets are scored per chromosome
unit and pooled (pair-weighted); 95 % intervals resample the units (genome-scale) or the loci.

Test data: the genome-scale sets su_genome and the new su_genome_amanitin, never used for reliability
(Gate 2's six test sets were used twice for this question and are not reused).

    python validation/reliability_v2.py --practice
    python validation/reliability_v2.py --test          # Gate 2e, run once (rule in frozen.RELIABILITY_V2)
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

import phase_a as A             # noqa: E402
import protocol as PR           # noqa: E402
import reliability as R2c       # noqa: E402
from chronocell import population as P                                     # noqa: E402

K_BOOT = 8
CANDIDATES = ("boot_sd", "evidence", "combined", "misfit")


def boot_log_sd(spec, k: int = K_BOOT) -> np.ndarray:
    """SD over k refits on resampled input of each pair's log model median (cached)."""
    path = A.CACHE / f"boot{k}_{spec.name}.npz"
    if path.exists():
        return np.load(path)["sd"]
    from benchmark.run import _hic_input
    u = A.load(spec)
    n = len(u["median"])
    rng = np.random.default_rng(7 + spec.split)
    logs = []
    if spec.input == "imaging":
        _, _, xyz = A._traces(spec)
        a, _ = PR.split(len(xyz), spec.split)
        r_img = 150.0 if u["meta"]["step_bp"] == 30_000 else float(u["meta"]["r_c_nm"])
        for _ in range(k):
            ha = PR.half_stats(xyz[rng.choice(a, len(a), replace=True)], r_img)
            res = A._fit(ha.freq, ha.seen, r_img, n, spec.split)
            logs.append(np.log(np.maximum(res.median_distance_nm, 1e-6)))
    else:
        c = np.asarray(u["counts"], dtype=np.float64)
        iu = np.triu_indices(n, 1)
        for _ in range(k):
            cb = np.zeros((n, n))
            cb[iu] = rng.poisson(np.nan_to_num(c[iu]))
            cb = cb + cb.T
            freq, n_eff, _ = _hic_input(cb)
            res = A._fit(freq, n_eff, float(u["meta"]["r_c_nm"]), n, spec.split)
            logs.append(np.log(np.maximum(res.median_distance_nm, 1e-6)))
    sd = np.std(np.stack(logs), axis=0).astype(np.float32)
    np.savez_compressed(path, sd=sd)
    return sd


def candidate_scores(spec) -> dict[str, np.ndarray]:
    u = A.load(spec)
    n = len(u["median"])
    c = np.nan_to_num(np.asarray(u["counts"], dtype=np.float64))
    pad = np.pad(c, 1)
    nb = sum(pad[1 + di:1 + di + n, 1 + dj:1 + dj + n] for di in (-1, 0, 1) for dj in (-1, 0, 1))
    out = {"boot_sd": -boot_log_sd(spec).astype(np.float64), "evidence": np.log1p(nb)}

    class _Res:                         # the misfit score needs a result object with the fitted spreads
        config = {"r_c_nm": float(u["meta"]["r_c_nm"])}
        median_distance_nm = u["median"].astype(np.float64)
        model = None
    out["misfit"] = R2c.scores(_Res, np.asarray(u["freq"], np.float64), c, np.asarray(u["seen"], np.float64))["misfit"]
    for v in out.values():
        np.fill_diagonal(v, np.nan)
    return out


def _strat_ranks(v: np.ndarray, strata: np.ndarray) -> np.ndarray:
    out = np.full(len(v), np.nan)
    for k in np.unique(strata):
        m = strata == k
        x = np.where(np.isfinite(v[m]), v[m], -1e9)
        out[m] = R2c.PR._rank(x) / max(m.sum() - 1, 1)
    return out


def unit_eval(spec) -> dict:
    """Per candidate: stratified pair-level Spearman (score vs -error) and the pair count, for one unit."""
    u = A.load(spec)
    n = len(u["median"])
    iu = np.triu_indices(n, 1)
    st = R2c._strata(u["sep"], iu)
    err = R2c.pair_error(u["median"].astype(np.float64), u["truth"].astype(np.float64), st, iu)
    sc = {k: v[iu] for k, v in candidate_scores(spec).items()}
    sc["combined"] = (_strat_ranks(sc["boot_sd"], st) + _strat_ranks(sc["evidence"], st)) / 2
    return {k: {"rho": R2c.stratified_spearman(sc[k], err, st), "pairs": int(np.isfinite(err).sum())} for k in CANDIDATES}


def pooled(evals: list[dict], cand: str, boot: int = 1000, seed: int = 0) -> dict:
    """Pair-weighted mean over units, with a 95 % interval from resampling the units."""
    rho = np.array([e[cand]["rho"] for e in evals])
    w = np.array([e[cand]["pairs"] for e in evals], dtype=float)
    ok = np.isfinite(rho)
    rho, w = rho[ok], w[ok]
    est = float(np.average(rho, weights=w))
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(boot):
        pick = rng.integers(0, len(rho), len(rho))
        bs.append(np.average(rho[pick], weights=w[pick]))
    return {"rho": est, "ci95": [float(x) for x in np.percentile(bs, [2.5, 97.5])], "units": int(len(rho))}


def practice(inputs=("imaging", "hic")) -> None:
    """One chunk per input type; each chunk's table is saved on its own, and the practice file is written from
    every chunk saved so far (refits are cached per unit, so an interrupted chunk resumes)."""
    import phase_a_units as U
    from hic_size_calibration import group_of
    for inp in inputs:
        table = {}
        specs = [s for s in U.PRACTICE if s.input == inp and s.thin == 1.0]
        ev = {s: unit_eval(s) for s in specs}
        for g in sorted({group_of(s) for s in specs}):
            evs = [e for s, e in ev.items() if group_of(s) == g]
            table[f"{g} · {inp}"] = {c: pooled(evs, c) for c in CANDIDATES}
            print(f"{g:26s} {inp:8s} " + "  ".join(f"{c} {table[f'{g} · {inp}'][c]['rho']:+.3f}" for c in CANDIDATES),
                  flush=True)
        (ROOT / f"results_reliability_v2_practice_{inp}.json").write_text(json.dumps(
            {"table": table, "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")},
            indent=1, default=float))
    table = {}
    for inp in ("imaging", "hic"):
        part = ROOT / f"results_reliability_v2_practice_{inp}.json"
        if part.exists():
            table.update(json.loads(part.read_text())["table"])
    (ROOT / "results_reliability_v2_practice.json").write_text(json.dumps(
        {"mode": "practice (all candidates; in-sample choice)", "table": table,
         "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}, indent=1, default=float))


TEST_PARTS = [(k, inp) for k in ("su_genome", "su_genome_amanitin") for inp in ("imaging", "hic")]


def test(parts: list[tuple[str, str]] | None = None) -> None:
    """Gate 2e, run once, in chunks: each (set, input) part is saved when done (git-ignored cache); the result
    file is written when every part is there. The score is fixed per input type in frozen.RELIABILITY_V2."""
    from frozen import RELIABILITY_V2 as R
    if (ROOT / "results_reliability_v2.json").exists():
        sys.exit("results_reliability_v2.json exists: Gate 2e runs once.")
    for key, inp in parts or TEST_PARTS:
        part = A.CACHE / f"gate2e_part_{key}_{inp}.json"
        if part.exists():
            continue
        score = R["score"][inp]
        evs = [unit_eval(s) for s in A.genome_units(key, inp)]
        r = pooled(evs, score)
        r["within_rule"] = bool(r["rho"] >= R["min_rho"] and r["ci95"][0] > 0)
        r["score"] = score
        part.write_text(json.dumps({**r, "all_candidates": {c: pooled(evs, c) for c in CANDIDATES}}, default=float))
        print(f"{key:20s} {inp:8s} {score} {r['rho']:+.3f} [{r['ci95'][0]:+.3f}, {r['ci95'][1]:+.3f}]  "
              f"{'yes' if r['within_rule'] else 'no'}", flush=True)
    files = {f"{k} · {i}": A.CACHE / f"gate2e_part_{k}_{i}.json" for k, i in TEST_PARTS}
    if not all(f.exists() for f in files.values()):
        print("parts still to run:", [k for k, f in files.items() if not f.exists()])
        return
    sets = {k: json.loads(f.read_text()) for k, f in files.items()}
    verdict = {inp: ("pass" if all(v["within_rule"] for k, v in sets.items() if k.endswith(inp)) else "fail")
               for inp in ("imaging", "hic")}
    (ROOT / "results_reliability_v2.json").write_text(json.dumps(
        {"mode": "test (run once, in chunks; score and rule frozen)", "rule": R, "sets": sets, "verdict": verdict,
         "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}, indent=1, default=float))
    print("verdict:", verdict)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    ap.add_argument("--input", choices=("imaging", "hic"), help="one input type only (a chunk)")
    ap.add_argument("--set", choices=("su_genome", "su_genome_amanitin"), help="test: one set only (a chunk)")
    a = ap.parse_args()
    if a.practice:
        practice((a.input,) if a.input else ("imaging", "hic"))
    else:
        test([(k, i) for k, i in TEST_PARTS if (a.set in (None, k)) and (a.input in (None, i))])


if __name__ == "__main__":
    main()
