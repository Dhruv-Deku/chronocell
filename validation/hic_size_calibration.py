"""
Gate 1c (Phase A1): calibrated sizes from Hi-C.

Models built from sequencing Hi-C rank distances well but are 0.35-0.57x the imaged size (Gate 1b). A
calibration learned on practice pairs of Hi-C + imaging multiplies every model distance by exp(h), with

    h = f(log10 separation) + optional terms in log10 locus spacing, log10 depth (n_eff) and protocol

The factor depends only on the pair's separation and on properties of the whole input, so it scales all
pairs at one separation alike: the distance pattern beyond the separation trend (the "% of the
reproducible pattern" score) is unchanged by construction; only absolute sizes move.

Forms compared on practice data by leave-one-dataset-out cross-validation (each held-out group: its
units at every depth and protocol); the chosen form is refitted on all practice units and frozen in
chronocell/data/hic_size_calibration.json before the test.

    python validation/hic_size_calibration.py --practice    # compare forms, freeze the chosen one
    python validation/hic_size_calibration.py --test        # Gate 1c, run once (rule in frozen.HIC_SIZE)
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

OUT = ROOT.parent / "chronocell" / "data" / "hic_size_calibration.json"
PROTOCOLS = ("in situ", "intact", "dilution")
FORMS = {   # name -> feature columns (besides the intercept)
    "global": [],
    "sep": ["L"],
    "sep2": ["L", "L2"],
    "sep2+step": ["L", "L2", "step"],
    "sep2+step+depth": ["L", "L2", "step", "depth"],
    "sep2+step+depth+protocol": ["L", "L2", "step", "depth", "intact", "dilution"],
}


def unit_rows(u: dict) -> dict:
    """Pairs i < j of one unit with finite, positive model and truth medians, and their features."""
    n = len(u["truth"])
    iu = np.triu_indices(n, 1)
    m, t, s = u["median"][iu].astype(np.float64), u["truth"][iu].astype(np.float64), u["sep"][iu].astype(np.float64)
    ok = np.isfinite(m) & np.isfinite(t) & (m > 0) & (t > 0) & (s > 0)
    meta = u["meta"]
    return {"model": m[ok], "truth": t[ok], "sep": s[ok], "step": float(meta["step_bp"]), "n_eff": float(meta["n_eff"]),
            "protocol": meta["protocol"]}


def design(rows: dict, cols: list[str]) -> np.ndarray:
    L = np.log10(rows["sep"])
    feats = {"L": L, "L2": L ** 2, "step": np.full_like(L, np.log10(rows["step"])),
             "depth": np.full_like(L, np.log10(max(rows["n_eff"], 1.0))),
             "intact": np.full_like(L, float(rows["protocol"] == "intact")),
             "dilution": np.full_like(L, float(rows["protocol"] == "dilution"))}
    return np.column_stack([np.ones_like(L)] + [feats[c] for c in cols])


def fit(units: list[dict], cols: list[str]) -> np.ndarray:
    """Least squares of log(truth / model) on the features, every unit weighted equally."""
    X = np.concatenate([design(r, cols) for r in units])
    y = np.concatenate([np.log(r["truth"] / r["model"]) for r in units])
    w = np.concatenate([np.full(len(r["model"]), 1.0 / len(r["model"])) for r in units])
    sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    return beta


def apply(rows: dict, cols: list[str], beta: np.ndarray) -> np.ndarray:
    return rows["model"] * np.exp(design(rows, cols) @ beta)


def scores(pred: np.ndarray, truth: np.ndarray, sep: np.ndarray) -> dict:
    s = PR._flat_scores(pred, truth, sep)
    return {"lin_ccc_nm": s["lin_ccc_nm"], "trend_removed_rho": s["spearman_distance_corrected"],
            "raw_spearman": s["spearman"], "size_ratio": float(np.median(pred / truth)), "pairs": int(len(pred))}


def group_of(spec) -> str:
    base = spec.dataset.split(":")[0]
    return {"su_chr2_parm_rep": "su_chr2"}.get(base, base)


def practice() -> None:
    import phase_a_units as U
    specs = [s for s in U.PRACTICE if s.input == "hic"]
    units = [(s, unit_rows(A.load(s))) for s in specs]
    groups = sorted({group_of(s) for s, _ in units})
    table = {}
    for form, cols in FORMS.items():
        per = {}
        for g in groups:
            train = [r for s, r in units if group_of(s) != g]
            held = [(s, r) for s, r in units if group_of(s) == g]
            beta = fit(train, cols)
            pred = np.concatenate([apply(r, cols, beta) for _, r in held])
            truth = np.concatenate([r["truth"] for _, r in held])
            sep = np.concatenate([r["sep"] for _, r in held])
            per[g] = scores(pred, truth, sep)
        table[form] = {"per_group": per, "mean_ccc": float(np.mean([v["lin_ccc_nm"] for v in per.values()])),
                       "min_ccc": float(np.min([v["lin_ccc_nm"] for v in per.values()]))}
        print(f"{form:28s} LODO CCC mean {table[form]['mean_ccc']:.3f} min {table[form]['min_ccc']:.3f}  " +
              "  ".join(f"{g}: {v['lin_ccc_nm']:.2f}/{v['size_ratio']:.2f}" for g, v in per.items()), flush=True)
    uncal = {g: scores(np.concatenate([r["model"] for s, r in units if group_of(s) == g]),
                       np.concatenate([r["truth"] for s, r in units if group_of(s) == g]),
                       np.concatenate([r["sep"] for s, r in units if group_of(s) == g])) for g in groups}
    # choose: highest mean LODO CCC; a more complex form must beat the simpler best by 0.01 to be chosen
    order = list(FORMS)
    best = order[0]
    for f in order[1:]:
        if table[f]["mean_ccc"] > table[best]["mean_ccc"] + 0.01:
            best = f
    beta = fit([r for _, r in units], FORMS[best])
    out = {"form": best, "columns": FORMS[best], "beta": beta.tolist(),
           "features": {"L": "log10 separation (bp)", "L2": "(log10 separation)^2", "step": "log10 locus spacing (bp)",
                        "depth": "log10 n_eff (median adjacent count / p_adjacent)", "intact": "1 for intact Hi-C",
                        "dilution": "1 for dilution Hi-C"},
           "target": "log(imaged median / model median); calibrated = model * exp(X beta)",
           "trained_on": sorted({s.name for s, _ in units}), "groups": groups,
           "selection": "highest mean leave-one-dataset-out CCC; a more complex form needs +0.01 over the simpler best",
           "fitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "note": "Practice data only (validation/phase_a_units.py). Its held-out test is Gate 1c (validation/RESULTS.md)."}
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    (ROOT / "results_hic_size_practice.json").write_text(json.dumps(
        {"mode": "practice (leave-one-dataset-out; in-sample choice)", "forms": table, "uncalibrated": uncal,
         "chosen": best, "beta": beta.tolist(), "units": [s.name for s, _ in units]}, indent=1, default=float))
    print(f"chosen: {best}  beta {np.round(beta, 4).tolist()}")


def test() -> None:
    import phase_a_units as U
    from frozen import HIC_SIZE as R
    cal = json.loads(OUT.read_text(encoding="utf-8"))
    cols, beta = cal["columns"], np.asarray(cal["beta"])
    sets = {"main": {}, "secondary": {}}
    plan = [("main", s.dataset, [s]) for s in U.TEST_HIC_MAIN]
    plan += [("main", "su_genome", A.genome_units("su_genome", "hic"))]
    plan += [("secondary", f"{s.dataset} · {s.source}", [s]) for s in U.TEST_HIC_SECONDARY]
    plan += [("secondary", "su_genome_amanitin", A.genome_units("su_genome_amanitin", "hic"))]
    for tier, name, specs in plan:
        rows = [unit_rows(A.load(s)) for s in specs]
        m = np.concatenate([r["model"] for r in rows])
        cal_pred = np.concatenate([apply(r, cols, beta) for r in rows])
        t = np.concatenate([r["truth"] for r in rows])
        sep = np.concatenate([r["sep"] for r in rows])
        before, after = scores(m, t, sep), scores(cal_pred, t, sep)
        ok = (after["lin_ccc_nm"] >= R["min_ccc"] and R["size_ratio"][0] <= after["size_ratio"] <= R["size_ratio"][1]
              and after["trend_removed_rho"] >= before["trend_removed_rho"] - R["pattern_tolerance"])
        sets[tier][name] = {"uncalibrated": before, "calibrated": after, "within_rule": bool(ok), "units": len(specs)}
        print(f"{tier:9s} {name:42s} CCC {before['lin_ccc_nm']:.3f} -> {after['lin_ccc_nm']:.3f}  size "
              f"{before['size_ratio']:.2f} -> {after['size_ratio']:.2f}  pattern rho {before['trend_removed_rho']:.4f} -> "
              f"{after['trend_removed_rho']:.4f}  {'yes' if ok else 'no'}", flush=True)
    verdict = "pass" if all(v["within_rule"] for v in sets["main"].values()) else "fail"
    (ROOT / "results_hic_size.json").write_text(json.dumps(
        {"mode": "test (run once; calibration frozen)", "rule": R, "calibration": cal, "sets": sets, "verdict": verdict,
         "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}, indent=1, default=float))
    print(f"verdict: {verdict}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    practice() if a.practice else test()


if __name__ == "__main__":
    main()
