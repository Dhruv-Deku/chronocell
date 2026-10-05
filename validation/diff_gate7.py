"""
Gate 7 (Phase B2): does the differential analysis (chronocell/differential.py) control the false-discovery rate
at the stated level on real replicate maps, and how much of a known change does it find?

Data: four independent untreated intact Hi-C experiments of the same cells (HCT116 RAD21-AID, Aiden lab, ENCODE,
GRCh38 MAPQ >= 30 maps read by region): condition A = ENCSR401WPL (ENCFF011IRM) + ENCSR002OIN (ENCFF058WJC),
condition B = ENCSR579TBL (ENCFF556CFS) + ENCSR697MNL (ENCFF720FKA). No biological change separates A and B.
Spike-ins: in each region, SPIKE pixels drawn at random (fixed seeds) among pixels 3-100 bins apart with a mean of
at least 20 reads get their counts in both B replicates multiplied by a fold (added Poisson reads: plant_changes).
Every significant pixel that was not planted is a false discovery (a conservative count: any real difference
between the experiments also counts as false).

    python validation/diff_gate7.py --practice     # practice region chr21:28-30 Mb: settings
    python validation/diff_gate7.py --test         # run once (rule in frozen.DIFF_GATE7)
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                                  # noqa: E402
from chronocell import differential as DF             # noqa: E402

ENC = "https://www.encodeproject.org/files/{0}/@@download/{0}.hic"
REPLICATES = {"A": [("encode_hct116_rad21aid_ctrl_401WPL", "ENCFF011IRM"), ("encode_hct116_rad21aid_ctrl_002OIN", "ENCFF058WJC")],
              "B": [("encode_hct116_rad21aid_ctrl_579TBL", "ENCFF556CFS"), ("encode_hct116_rad21aid_ctrl_697MNL", "ENCFF720FKA")]}
for _side in REPLICATES.values():
    for _k, _acc in _side:
        D.HIC_SOURCES.setdefault(_k, (ENC.format(_acc), "GRCh38", "HCT116"))
RES = 10_000
WINDOW = 2_000_000
FDR = 0.05
SPIKE = 100
FOLDS = (2.0, 4.0)
SEEDS = (0, 1, 2)
MAX_SEP = 200
PRACTICE = [("chr21", 28_000_000)]
SETTINGS = {"min5_dist": {"min_count": 5.0, "distance_normalise": True},
            "min10_dist": {"min_count": 10.0, "distance_normalise": True},
            "min5_nodist": {"min_count": 5.0, "distance_normalise": False}}


def maps(chrom: str, start: int) -> tuple[list[np.ndarray], list[np.ndarray]]:
    get = lambda k: D.hic_region(k, chrom, start, start + WINDOW, RES).astype(float)        # noqa: E731
    return [get(k) for k, _ in REPLICATES["A"]], [get(k) for k, _ in REPLICATES["B"]]


def candidates(A: list[np.ndarray], B: list[np.ndarray]) -> np.ndarray:
    n = len(A[0])
    mean = np.mean(A + B, axis=0)
    i, j = np.triu_indices(n, 3)
    ok = (j - i <= 100) & (mean[i, j] >= 20)
    return np.stack([i[ok], j[ok]], axis=1)


def run_region(chrom: str, start: int, setting: dict) -> list[dict]:
    A, B = maps(chrom, start)
    out = []
    null = DF.differential_pixels(A, B, RES, FDR, max_sep=MAX_SEP, **setting)
    out.append({"region": f"{chrom}:{start // 1_000_000}-{(start + WINDOW) // 1_000_000} Mb", "fold": 1.0, "seed": None,
                "discoveries": int(len(null.significant)), "false": int(len(null.significant)), "planted": 0,
                "found": 0, "tested": int(len(null.pixels))})
    cand = candidates(A, B)
    for fold in FOLDS:
        for seed in SEEDS:
            rng = np.random.default_rng(1000 * seed + int(fold))
            pick = cand[rng.choice(len(cand), size=min(SPIKE, len(cand)), replace=False)]
            planted = {(int(a), int(b)) for a, b in pick}
            Bp = [DF.plant_changes(m, list(planted), fold, rng) for m in B]
            res = DF.differential_pixels(A, Bp, RES, FDR, max_sep=MAX_SEP, **setting)
            sig = {(int(r.i), int(r.j)) for r in res.significant.itertuples()}
            tp = len(sig & planted)
            out.append({"region": out[0]["region"], "fold": fold, "seed": seed, "discoveries": len(sig), "false": len(sig) - tp,
                        "planted": len(planted), "found": tp, "tested": int(len(res.pixels))})
    return out


def summarise(rows: list[dict]) -> dict:
    spike = [r for r in rows if r["fold"] > 1]
    fdp = np.array([r["false"] / max(r["discoveries"], 1) for r in spike])
    regions = sorted({r["region"] for r in spike})
    rng = np.random.default_rng(7)
    boots = []
    for _ in range(1000):
        pick = rng.choice(regions, len(regions))
        boots.append(np.mean([r["false"] / max(r["discoveries"], 1) for g in pick for r in spike if r["region"] == g]))
    rec = {str(f): float(np.mean([r["found"] / max(r["planted"], 1) for r in spike if r["fold"] == f])) for f in FOLDS}
    null = [r for r in rows if r["fold"] == 1.0]
    return {"mean_fdp": float(fdp.mean()), "mean_fdp_ci95": [float(x) for x in np.percentile(boots, [2.5, 97.5])],
            "recall": rec, "null_discoveries": {r["region"]: r["discoveries"] for r in null}, "runs": len(spike)}


def practice() -> None:
    table = {}
    for name, setting in SETTINGS.items():
        rows = [r for c, s in PRACTICE for r in run_region(c, s, setting)]
        table[name] = {"setting": setting, "summary": summarise(rows), "rows": rows}
        sm = table[name]["summary"]
        print(f"{name:12s} mean FDP {sm['mean_fdp']:.3f}  recall x2 {sm['recall']['2.0']:.2f}  x4 {sm['recall']['4.0']:.2f}  "
              f"null discoveries {list(sm['null_discoveries'].values())}", flush=True)
    (ROOT / "results_gate7_practice.json").write_text(json.dumps({"mode": "practice (chr21:28-30 Mb)", "table": table,
                                                                 "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")},
                                                                indent=1, default=float), encoding="utf-8")


def test() -> None:
    from frozen import DIFF_GATE7 as R
    out = ROOT / "results_gate7.json"
    if out.exists():
        sys.exit("results_gate7.json exists: Gate 7 runs once.")
    rows = []
    for chrom, start in R["regions"]:
        part = D.DATA / "cache" / f"gate7_{chrom}_{start}.json"
        if part.exists():
            rr = json.loads(part.read_text())
        else:
            rr = run_region(chrom, start, R["setting"])
            part.parent.mkdir(parents=True, exist_ok=True)
            part.write_text(json.dumps(rr))
        rows += rr
        print(f"{rr[0]['region']}: null {rr[0]['discoveries']} discoveries; spike-in FDP "
              + " ".join(f"{r['false'] / max(r['discoveries'], 1):.3f}" for r in rr[1:]), flush=True)
    sm = summarise(rows)
    verdict = "pass" if sm["mean_fdp"] <= R["fdr"] else "fail"
    out.write_text(json.dumps({"mode": "test (run once; settings frozen)", "rule": R, "summary": sm, "rows": rows,
                               "verdict": verdict, "not_run": R["not_run"],
                               "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")},
                              indent=1, default=float), encoding="utf-8")
    print(f"Gate 7: mean FDP {sm['mean_fdp']:.3f} (nominal {R['fdr']}) -> {verdict}; recall {sm['recall']}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        practice() if a.practice else test()


if __name__ == "__main__":
    main()
