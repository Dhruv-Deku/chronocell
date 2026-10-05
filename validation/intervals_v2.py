"""
Gate 2d (Phase A2): honest uncertainty ranges for imaging-derived AND sequencing Hi-C input.

The Gaussian model's stated intervals are too narrow (Gate 2), and with Hi-C input far too narrow, even
after the PIT-histogram recalibration (Gate 2b). Phase A2 states intervals on the log-distance scale:

    range for pair (i, j) at level p = median_ij * exp([q_band((1-p)/2), q_band((1+p)/2)])

with median_ij the model's median distance (Hi-C input: the A1 size-calibrated median when the A1
calibration file exists) and q_band the empirical quantiles of r = log(single-cell distance / median_ij)
over practice units of the same input type, pooled per separation band with every unit weighted equally
(split-conformal intervals per input type and band). Heavier tails than the Maxwell law are therefore
taken from the data, not assumed.

Variants compared on practice data by leave-one-dataset-out (each held-out group's coverage with
quantiles from the other groups): per band vs one pooled band; Hi-C medians calibrated (A1) vs not; and
the existing Maxwell + PIT recalibration as the reference. The chosen variant is frozen in
chronocell/data/intervals_v2.json before the test.

    python validation/intervals_v2.py --practice
    python validation/intervals_v2.py --test          # Gate 2d, run once (rule in frozen.INTERVALS_V2)
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

OUT = ROOT.parent / "chronocell" / "data" / "intervals_v2.json"
BANDS_BP = [0, 100_000, 300_000, 1_000_000, 3_000_000, 10_000_000, 300_000_000]
GRID = np.linspace(-4.0, 4.0, 3201)           # log-ratio histogram edges (0.0025 wide)
LEVELS = (0.5, 0.8, 0.9)
MIN_PAIRS = 200                               # a band of a unit is scored only with at least this many pairs


def band_of(sep: np.ndarray) -> np.ndarray:
    return np.clip(np.searchsorted(BANDS_BP, sep, side="right") - 1, 0, len(BANDS_BP) - 2)


def size_calibrated_median(u: dict) -> np.ndarray | None:
    """The A1 calibrated median (Hi-C units), from the frozen chronocell/data/hic_size_calibration.json."""
    path = ROOT.parent / "chronocell" / "data" / "hic_size_calibration.json"
    if u["meta"]["protocol"] == "imaging" or not path.exists():
        return None
    import hic_size_calibration as HS
    cal = json.loads(path.read_text(encoding="utf-8"))
    n = len(u["median"])
    iu = np.triu_indices(n, 1)
    rows = {"model": u["median"][iu].astype(np.float64), "sep": np.maximum(u["sep"][iu].astype(np.float64), 1.0),
            "step": float(u["meta"]["step_bp"]), "n_eff": float(u["meta"]["n_eff"]), "protocol": u["meta"]["protocol"]}
    out = np.zeros((n, n))
    out[iu] = HS.apply(rows, cal["columns"], np.asarray(cal["beta"]))
    return out + out.T


def log_ratio_hist(spec, median: np.ndarray, rows: int = 16) -> dict[int, np.ndarray]:
    """Per separation band: histogram (on GRID) of log(single-copy distance / median) over half B's copies."""
    u = A.load(spec)
    x = np.asarray(A.single_copies(spec), dtype=np.float32)
    n = x.shape[1]
    band = band_of(u["sep"])
    upper = np.triu(np.ones((n, n), bool), 1)
    med = np.where(upper & (median > 0), median, np.nan).astype(np.float32)
    hists = {b: np.zeros(len(GRID) - 1) for b in range(len(BANDS_BP) - 1)}
    for a in range(0, n, rows):
        e = min(n, a + rows)
        d = np.linalg.norm(x[:, a:e, None, :] - x[:, None, :, :], axis=-1)        # (copies, block, n)
        with np.errstate(divide="ignore"):
            r = np.log(d / med[a:e][None])            # a zero distance (identical coordinates) is -inf: below every interval
        for b in hists:
            m = (band[a:e] == b)[None] & ~np.isnan(r)
            if m.any():
                hists[b] += np.histogram(np.clip(r[m], GRID[0], GRID[-1]), bins=GRID)[0]
    return hists


def quantiles(hist: np.ndarray, qs) -> np.ndarray:
    cdf = np.concatenate([[0.0], np.cumsum(hist)]) / max(hist.sum(), 1.0)
    return np.interp(qs, cdf, GRID)


def coverage(hist: np.ndarray, lo: float, hi: float) -> float:
    cdf = np.concatenate([[0.0], np.cumsum(hist)]) / max(hist.sum(), 1.0)
    return float(np.interp(hi, GRID, cdf) - np.interp(lo, GRID, cdf))


def pooled(hists: list[dict], bands_on: bool) -> dict[int, np.ndarray]:
    """Equal weight per unit (each unit's histogram normalised), per band or all bands pooled."""
    out = {}
    for b in hists[0]:
        parts = [h[b] / h[b].sum() for h in hists if h[b].sum() >= MIN_PAIRS]
        out[b] = np.mean(parts, axis=0) if parts else None
    if not bands_on:
        allp = [h for h in out.values() if h is not None]
        one = np.mean(allp, axis=0)
        out = {b: one for b in out}
    return out


def unit_hists(spec, variant: dict) -> dict[int, np.ndarray]:
    u = A.load(spec)
    med = u["median"].astype(np.float64)
    if variant.get("calibrated") and u["meta"]["protocol"] != "imaging":
        c = size_calibrated_median(u)
        med = c if c is not None else med
    return log_ratio_hist(spec, med)


def evaluate(held: list[dict], q: dict[int, np.ndarray]) -> dict:
    """Coverage per level and band of the pooled held-out histograms with quantiles q."""
    out = {}
    for b in held[0]:
        hs = [h[b] for h in held if h[b].sum() >= MIN_PAIRS]
        if not hs or q.get(b) is None:
            continue
        tot = np.sum(hs, axis=0)
        qq = quantiles(q[b], [p for lv in LEVELS for p in ((1 - lv) / 2, (1 + lv) / 2)])
        out[b] = {str(lv): coverage(tot, qq[2 * k], qq[2 * k + 1]) for k, lv in enumerate(LEVELS)}
        out[b]["width90"] = float(np.exp(qq[5] - qq[4]))
        out[b]["pairs_copies"] = float(tot.sum())
    return out


VARIANTS = {"bands": {"bands": True, "calibrated": False}, "pooled": {"bands": False, "calibrated": False},
            "bands+A1": {"bands": True, "calibrated": True}, "pooled+A1": {"bands": False, "calibrated": True}}


def practice() -> None:
    import phase_a_units as U
    from hic_size_calibration import group_of
    specs = {inp: [s for s in U.PRACTICE if s.input == inp and s.thin == 1.0] for inp in ("imaging", "hic")}
    table = {}
    for name, var in VARIANTS.items():
        for inp, ss in specs.items():
            if var["calibrated"] and inp == "imaging":
                continue
            hs = {s: unit_hists(s, var) for s in ss}
            groups = sorted({group_of(s) for s in ss})
            per = {}
            for g in groups:
                q = pooled([h for s, h in hs.items() if group_of(s) != g], var["bands"])
                per[g] = evaluate([h for s, h in hs.items() if group_of(s) == g], q)
            worst = max(abs(v["0.9"] - 0.9) for g in per.values() for v in g.values())
            table[f"{name}|{inp}"] = {"per_group": per, "worst_abs_err_90": worst}
            print(f"{name:10s} {inp:8s} worst |cov90 - 0.90| = {worst:.3f}  " + "  ".join(
                f"{g}: " + "/".join(f"{100 * v['0.9']:.0f}" for v in per[g].values()) for g in groups), flush=True)
    (ROOT / "results_intervals_v2_practice.json").write_text(json.dumps({"mode": "practice (leave-one-dataset-out)",
                                                                       "variants": table}, indent=1, default=float))
    freeze(table, specs)


TIE = 0.005     # a variant is "as good" as the best if its worst-case error is within this (A1's factor is global, so
                # with and without it the intervals differ only by histogram binning)


def choose(table: dict, inp: str) -> str:
    """Smallest worst-case |coverage at 90 % - 0.90|; within TIE of it, the simplest variant (order of VARIANTS:
    without the A1 calibration before with it, per-band after pooled only if it is better by more than TIE)."""
    cands = [k for k in table if k.endswith("|" + inp)]
    best = min(table[k]["worst_abs_err_90"] for k in cands)
    simple = sorted(cands, key=lambda k: (VARIANTS[k.split("|")[0]]["calibrated"], VARIANTS[k.split("|")[0]]["bands"]))
    return next(k for k in simple if table[k]["worst_abs_err_90"] <= best + TIE)


def freeze(table: dict, specs: dict) -> None:
    chosen = {}
    for inp in ("imaging", "hic"):
        best = choose(table, inp)
        var = VARIANTS[best.split("|")[0]]
        hs = [unit_hists(s, var) for s in specs[inp]]
        q = pooled(hs, var["bands"])
        chosen[inp] = {"variant": best.split("|")[0], **var,
                       "quantiles": {str(b): (quantiles(h, np.linspace(0, 1, 201)).tolist() if h is not None else None)
                                     for b, h in q.items()},
                       "trained_on": [s.name for s in specs[inp]]}
        print(f"chosen for {inp}: {best}")
    out = {"method": "split-conformal log-distance intervals per input type and separation band (Phase A2)",
           "bands_bp": BANDS_BP, "quantile_levels": np.linspace(0, 1, 201).tolist(), "inputs": chosen,
           "selection": f"smallest worst-case |coverage at 90 % - 0.90| over held-out practice groups and bands; within "
                        f"{TIE} of it the simplest variant (no A1 calibration, pooled bands)",
           "fitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "note": "Practice data only. Its held-out test is Gate 2d (validation/RESULTS.md)."}
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")


def refreeze() -> None:
    """Re-run only the choice and the freeze from the saved practice table (no new comparison)."""
    import phase_a_units as U
    table = json.loads((ROOT / "results_intervals_v2_practice.json").read_text())["variants"]
    specs = {inp: [s for s in U.PRACTICE if s.input == inp and s.thin == 1.0] for inp in ("imaging", "hic")}
    freeze(table, specs)


def test() -> None:
    import phase_a_units as U
    from frozen import INTERVALS_V2 as R
    cal = json.loads(OUT.read_text(encoding="utf-8"))
    sets = {}
    plan = [(s.dataset, "imaging", [s]) for s in U.TEST_IMAGING]
    plan += [(s.dataset, "hic", [s]) for s in U.TEST_HIC_MAIN]
    plan += [("su_genome", "imaging", A.genome_units("su_genome", "imaging")),
             ("su_genome", "hic", A.genome_units("su_genome", "hic"))]
    for name, inp, specs in plan:
        c = cal["inputs"][inp]
        qgrid = np.asarray(cal["quantile_levels"])
        q = {int(b): (np.asarray(v) if v is not None else None) for b, v in c["quantiles"].items()}
        held = [unit_hists(s, c) for s in specs]
        res = {}
        for b in held[0]:
            hs = [h[b] for h in held if h[b].sum() >= MIN_PAIRS]
            if not hs or q.get(b) is None:
                continue
            tot = np.sum(hs, axis=0)
            row = {}
            for lv in LEVELS:
                lo, hi = np.interp([(1 - lv) / 2, (1 + lv) / 2], qgrid, q[b])
                row[str(lv)] = coverage(tot, lo, hi)
            lo, hi = np.interp([0.05, 0.95], qgrid, q[b])
            row["width90"] = float(np.exp(hi - lo))
            row["pairs_copies"] = float(tot.sum())
            res[str(b)] = row
        ok = all(R["cov90"][0] <= v["0.9"] <= R["cov90"][1] and R["cov50"][0] <= v["0.5"] <= R["cov50"][1]
                 for v in res.values())
        sets[f"{name} · {inp}"] = {"bands": res, "within_rule": bool(ok)}
        print(f"{name:26s} {inp:8s} " + "  ".join(f"b{b}: {100 * v['0.9']:.0f}/{100 * v['0.5']:.0f}%" for b, v in res.items())
              + f"  {'yes' if ok else 'no'}", flush=True)
    verdict = {inp: ("pass" if all(v["within_rule"] for k, v in sets.items() if k.endswith(inp)) else "fail")
               for inp in ("imaging", "hic")}
    (ROOT / "results_intervals_v2.json").write_text(json.dumps(
        {"mode": "test (run once; quantiles frozen)", "rule": R, "intervals": {k: cal["inputs"][k]["variant"] for k in cal["inputs"]},
         "bands_bp": BANDS_BP, "sets": sets, "verdict": verdict,
         "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}, indent=1, default=float))
    print("verdict:", verdict)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    g.add_argument("--refreeze", action="store_true", help="choose and freeze again from the saved practice table")
    a = ap.parse_args()
    practice() if a.practice else refreeze() if a.refreeze else test()


if __name__ == "__main__":
    main()
