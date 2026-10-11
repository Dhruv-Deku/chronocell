"""
Gate 6b (October 2026): Gate 6 again, with the settings chosen on more data and tested on new cell lines.

Gate 6 chose ChronoCell's loop-calling settings on three GM12878 windows and failed on K562 (F1 0.40 against
Mustache's 0.49: 315 calls for 156 reference loops, precision 0.30). Gate 6b: the same caller, settings re-chosen on
all nine windows already used (GM12878, K562, IMR-90), with two post-filters added to the grid (a minimum donut
enrichment and a minimum cluster size, as HiCCUPS uses), chromosight and Mustache run on the same windows; test on
two cell lines never used for loops: HMEC (ENCSR711AVS: map ENCFF943JRY, HiCCUPS loops ENCFF999UXN) and HAP-1
(ENCSR390HMC: map ENCFF898HRO, loops ENCFF557WIU), in three 10 Mb regions not used by the other round-2 tests.
Reference, matching and comparators exactly as Gate 6 (validation/loops_gate6.py).

    python validation/loops_gate6b.py --practice
    python validation/loops_gate6b.py --test          # run once (rule in frozen.LOOPS_GATE6B)
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                                                  # noqa: E402
import loops_gate6 as G6                                              # noqa: E402
from chronocell import analysis as AN, normalize as NZ               # noqa: E402

NEW = {"hmec": ("ENCFF943JRY", "ENCFF999UXN", "b38bad690aff7dcb26ec77c9d236d8c0"),
       "hap1": ("ENCFF898HRO", "ENCFF557WIU", "9454561b55b69123cadb53078a533354")}
for _k, (_h, _l, _m) in NEW.items():
    G6.SETS.setdefault(_k, (_h, _l, _m))
    D.HIC_SOURCES.setdefault(f"encode_loops_{_k}", (G6.ENC_HIC.format(_h), "GRCh38", _k))
SEEN = {"gm12878": G6.PRACTICE["gm12878"],
        "k562": [("chr4", 100_000_000), ("chr7", 100_000_000), ("chr11", 60_000_000)],
        "imr90": [("chr4", 100_000_000), ("chr7", 100_000_000), ("chr11", 60_000_000)]}
# FDR widened once (0.001-0.005) when the first practice run's best value sat at 0.01, the bottom of its grid
GRID = {"fdr": [0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.2], "balance": ["none", "kr"], "min_oe": [0.0, 2.5, 3.0],
        "min_cluster": [1, 2]}
TOOLS = ("chromosight", "mustache")


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def raw_calls(i, j, v, n, fdr: float, balance: str) -> list:
    w = NZ.kr_balance(i, j, v, n).weights if balance == "kr" else None
    return AN.call_loops(i, j, v, n, G6.RES, w, max_sep_bp=G6.MAX_SEP, min_sep_bp=G6.MIN_SEP, fdr=fdr)


def to_frame(calls: list, chrom: str, start: int, min_oe: float, min_cluster: int) -> pd.DataFrame:
    keep = [c for c in calls if c.oe_donut >= min_oe and c.cluster_size >= min_cluster]
    return pd.DataFrame([(chrom, start + c.i * G6.RES + G6.RES // 2, start + c.j * G6.RES + G6.RES // 2) for c in keep],
                        columns=["chrom", "a", "b"])


def windows_result(cell: str, windows: list, settings: list[dict], tag: str) -> dict:
    """ChronoCell under every setting, the two tools, and the reference, pooled over the windows of one cell line."""
    ref = G6.reference(cell)
    acc = {"reference": [], **{t: [] for t in TOOLS}, **{k: [] for k in range(len(settings))}}
    notes = []
    for chrom, start in windows:
        i, j, v, n = G6.window_counts(cell, chrom, start)
        acc["reference"].append(G6.in_window(ref, chrom, start))
        raw = {}
        for k, st in enumerate(settings):
            key = (st["fdr"], st["balance"])
            if key not in raw:
                raw[key] = raw_calls(i, j, v, n, *key)
            acc[k].append(to_frame(raw[key], chrom, start, st["min_oe"], st["min_cluster"]))
        G6.WORK.mkdir(parents=True, exist_ok=True)
        cool = G6.write_cool(i, j, v, n, chrom, start, G6.WORK / f"{tag}_{cell}_{chrom}_{start}.mcool")
        for tool in TOOLS:
            got = G6.tool_calls(tool, cool, chrom, start)
            if isinstance(got, str):
                notes.append(f"{chrom}:{start}: {got}")
                acc[tool].append(None)
            else:
                acc[tool].append(got)
        print(cell, chrom, start, "done", flush=True)
    refs = pd.concat(acc["reference"])
    out = {"notes": notes, "settings": [], "tools": {}}
    for k, st in enumerate(settings):
        out["settings"].append({**st, **G6.match(pd.concat(acc[k]), refs)})
    for tool in TOOLS:
        out["tools"][tool] = ({"status": "not run on every window"} if any(x is None for x in acc[tool])
                              else G6.match(pd.concat(acc[tool]), refs))
    return out


def practice() -> dict:
    settings = [{"fdr": f, "balance": b, "min_oe": o, "min_cluster": c} for f in GRID["fdr"] for b in GRID["balance"]
                for o in GRID["min_oe"] for c in GRID["min_cluster"]]
    per = {cell: windows_result(cell, w, settings, "g6b_practice") for cell, w in SEEN.items()}
    rows = []
    for k, st in enumerate(settings):
        f1 = {c: per[c]["settings"][k]["f1"] for c in per}
        tools = {c: max(per[c]["tools"][t]["f1"] for t in TOOLS if "f1" in per[c]["tools"][t]) for c in per}
        margin = {c: f1[c] - tools[c] for c in per}
        rows.append({**st, "f1": f1, "best_tool_f1": tools, "margin": margin, "min_margin": min(margin.values()),
                     "detail": {c: per[c]["settings"][k] for c in per}})
        print(st, {c: round(x, 3) for c, x in f1.items()}, "min margin", round(rows[-1]["min_margin"], 3), flush=True)
    best = max(rows, key=lambda r: (r["min_margin"], -r["fdr"]))
    out = {"made": _now(), "windows": SEEN, "grid": rows, "tools": {c: per[c]["tools"] for c in per},
           "notes": {c: per[c]["notes"] for c in per},
           "choice": {k: best[k] for k in ("fdr", "balance", "min_oe", "min_cluster")}, "choice_scores": best}
    (ROOT / "results_gate6b_practice.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("choice", out["choice"], best["f1"], best["best_tool_f1"], flush=True)
    return out


def test() -> dict:
    from frozen import LOOPS_GATE6B as R
    path = ROOT / "results_gate6b.json"
    if path.exists():
        sys.exit("results_gate6b.json exists: Gate 6b runs once.")
    sets, ok = {}, {}
    for cell, windows in R["test"].items():
        res = windows_result(cell, [tuple(w) for w in windows], [R["choice"]], "g6b_test")
        s = {"chronocell": res["settings"][0], **res["tools"], "notes": res["notes"]}
        others = [s[t]["f1"] for t in TOOLS if "f1" in s[t]]
        ok[cell] = bool(len(others) == len(TOOLS)) and s["chronocell"]["f1"] >= max(others)
        sets[cell] = s
        print(cell, {t: round(s[t]["f1"], 3) for t in ("chronocell",) + TOOLS if "f1" in s[t]}, flush=True)
    verdict = "pass" if all(ok.values()) else "fail"
    out = {"mode": "test (run once; settings frozen)", "rule": R, "sets": sets, "at_least_tools": ok, "verdict": verdict,
           "run_utc": _now()}
    path.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("Gate 6b:", verdict, ok, flush=True)
    return out


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
