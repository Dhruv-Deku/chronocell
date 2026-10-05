"""
Gate 6 (Phase B3) and the loop part of B9: do ChronoCell's loop calls (chronocell.analysis.call_loops) agree with
the reference HiCCUPS calls at least as well as the comparable tools that run here?

Reference: HiCCUPS loops that ENCODE called on the very map we read ("hic-loop-calling-step"; GRCh38; MD5 from the
portal). Practice: GM12878 (ENCSR968KAY: map ENCFF256UOW, loops ENCFF041XLP). Test, held-out cell lines:
K562 (ENCSR545YBD: ENCFF616PUW, ENCFF693XIL) and IMR-90 (ENCSR852KQC: ENCFF188SSH, ENCFF527JOL).
Windows of 10 Mb at 10 kb, read by region; loops 30 kb - 2 Mb apart, both anchors inside the window.
Comparators, run on the same windows (written as .mcool with ICE weights) from an isolated environment
(.chronocell_cache/tools-venv): chromosight 1.6.3 (default loop kernel), Mustache 1.3.3 (defaults; its .hic
reader hic-straw does not build here, so it reads the .mcool). HiCCUPS itself (Juicer tools, Java) is not
installed here, and is the reference.
Matching: a call matches a reference loop when both anchors (midpoints) lie within 25 kb, one-to-one, greedy by
distance. Precision, recall, F1 per test cell line (windows pooled).

    python validation/loops_gate6.py --practice     # chooses FDR and balancing for ChronoCell on GM12878
    python validation/loops_gate6.py --test         # run once (rule in frozen.LOOPS_GATE6)
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                                                  # noqa: E402
from chronocell import analysis as AN, exports as EX, normalize as NZ  # noqa: E402

ENC_HIC = "https://www.encodeproject.org/files/{0}/@@download/{0}.hic"
SETS = {"gm12878": ("ENCFF256UOW", "ENCFF041XLP", "963d23ee4c1096eb07f690aa7adead8c"),
        "k562": ("ENCFF616PUW", "ENCFF693XIL", "ae663464bdbe60998e422254ea0dac2c"),
        "imr90": ("ENCFF188SSH", "ENCFF527JOL", "6e0278ffbe387176c0c008ce3e8bae30")}
for _k, (_h, _, _) in SETS.items():
    D.HIC_SOURCES.setdefault(f"encode_loops_{_k}", (ENC_HIC.format(_h), "GRCh38", _k))
RES = 10_000
WINDOW = 10_000_000
MIN_SEP, MAX_SEP = 30_000, 2_000_000
TOL = 25_000
PRACTICE = {"gm12878": [("chr1", 100_000_000), ("chr2", 100_000_000), ("chr3", 100_000_000)]}
CHOICES = {f"fdr{f}_{b}": {"fdr": f, "balance": b} for f in (0.05, 0.1, 0.2) for b in ("none", "kr")}
VENV = ROOT.parent / ".chronocell_cache" / "tools-venv" / "Scripts"
WORK = D.DATA / "cache" / "gate6"


def reference(cell: str) -> pd.DataFrame:
    _, acc, md5 = SETS[cell]
    p = D.fetch_url(f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bedpe.gz", D.DATA / "encode" / f"{acc}.bedpe.gz", md5)
    rows = []
    with gzip.open(p, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or not f[1].isdigit():
                continue
            c = f[0] if f[0].startswith("chr") else "chr" + f[0]
            rows.append((c, (int(f[1]) + int(f[2])) // 2, (int(f[4]) + int(f[5])) // 2))
    return pd.DataFrame(rows, columns=["chrom", "a", "b"])


def window_counts(cell: str, chrom: str, start: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    m = D.hic_region(f"encode_loops_{cell}", chrom, start, start + WINDOW, RES)
    n = len(m)
    i, j = np.triu_indices(n, 1)
    v = m[i, j]
    ok = v > 0
    return i[ok], j[ok], v[ok].astype(float), n


def in_window(df: pd.DataFrame, chrom: str, start: int) -> pd.DataFrame:
    sel = (df["chrom"] == chrom) & (df["a"] >= start) & (df["b"] < start + WINDOW) & \
          ((df["b"] - df["a"]) >= MIN_SEP) & ((df["b"] - df["a"]) <= MAX_SEP)
    return df[sel]


def match(calls: pd.DataFrame, ref: pd.DataFrame) -> dict:
    pairs = []
    for x, c in enumerate(calls.itertuples(index=False)):
        for y, r in enumerate(ref.itertuples(index=False)):
            if c.chrom == r.chrom and abs(c.a - r.a) <= TOL and abs(c.b - r.b) <= TOL:
                pairs.append((abs(c.a - r.a) + abs(c.b - r.b), x, y))
    used_c, used_r = set(), set()
    for _, x, y in sorted(pairs):
        if x not in used_c and y not in used_r:
            used_c.add(x)
            used_r.add(y)
    tp = len(used_c)
    p = tp / max(len(calls), 1)
    r = tp / max(len(ref), 1)
    return {"calls": int(len(calls)), "reference": int(len(ref)), "matched": tp, "precision": p, "recall": r,
            "f1": 2 * p * r / (p + r) if p + r > 0 else 0.0}


def chronocell_calls(i, j, v, n, chrom, start, fdr: float, balance: str) -> pd.DataFrame:
    w = None
    if balance == "kr":
        w = NZ.kr_balance(i, j, v, n).weights
    calls = AN.call_loops(i, j, v, n, RES, w, max_sep_bp=MAX_SEP, min_sep_bp=MIN_SEP, fdr=fdr)
    return pd.DataFrame([(chrom, start + c.i * RES + RES // 2, start + c.j * RES + RES // 2) for c in calls],
                        columns=["chrom", "a", "b"])


def write_cool(i, j, v, n, chrom: str, start: int, path: Path) -> Path:
    """The window as a one-chromosome .mcool ('chrom' renamed to its window, positions from 0) with ICE weights."""
    b = NZ.ice_balance(i, j, v, n)
    w = np.where(np.isfinite(b.bias), 1.0 / b.bias, np.nan)
    return EX.write_mcool(path, [(chrom, n * RES)], {chrom: (i, j, v)}, RES, factors=(1,), weights={chrom: w})


def tool_calls(tool: str, cool: Path, chrom: str, start: int) -> pd.DataFrame | str:
    out = cool.with_suffix(f".{tool}")
    try:
        if tool == "chromosight":
            subprocess.run([str(VENV / "chromosight.exe"), "detect", "--pattern=loops", f"--min-dist={MIN_SEP}",
                            f"--max-dist={MAX_SEP}", "--no-plotting", f"{cool}::/resolutions/{RES}", str(out)],
                           check=True, capture_output=True, text=True, timeout=3600)
            df = pd.read_csv(str(out) + ".tsv", sep="\t")
            a = (df["start1"] + df["end1"]) // 2
            b = (df["start2"] + df["end2"]) // 2
        else:
            subprocess.run([str(VENV / "mustache.exe"), "-f", str(cool), "-r", str(RES), "-o", str(out) + ".tsv",
                            "-d", str(MAX_SEP), "-ch", chrom, "-norm", "weight"],
                           check=True, capture_output=True, text=True, timeout=3600)
            df = pd.read_csv(str(out) + ".tsv", sep="\t")
            a = (df["BIN1_START"] + df["BIN1_END"]) // 2
            b = (df["BIN2_START"] + df["BIN2_END"]) // 2
    except (subprocess.SubprocessError, OSError, KeyError, pd.errors.EmptyDataError) as exc:
        err = getattr(exc, "stderr", "") or str(exc)
        return f"{tool} failed: {str(err).strip().splitlines()[-1][:200] if str(err).strip() else type(exc).__name__}"
    calls = pd.DataFrame({"chrom": chrom, "a": start + a.to_numpy(), "b": start + b.to_numpy()})
    sep = calls["b"] - calls["a"]
    return calls[(sep >= MIN_SEP) & (sep <= MAX_SEP)]


def practice() -> None:
    ref = reference("gm12878")
    table = {k: {"calls": [], "ref": []} for k in CHOICES}
    for chrom, start in PRACTICE["gm12878"]:
        i, j, v, n = window_counts("gm12878", chrom, start)
        r = in_window(ref, chrom, start)
        for name, c in CHOICES.items():
            table[name]["calls"].append(chronocell_calls(i, j, v, n, chrom, start, c["fdr"], c["balance"]))
            table[name]["ref"].append(r)
    res = {}
    for name in CHOICES:
        res[name] = match(pd.concat(table[name]["calls"]), pd.concat(table[name]["ref"]))
        print(f"{name:14s} P {res[name]['precision']:.3f}  R {res[name]['recall']:.3f}  F1 {res[name]['f1']:.3f}  "
              f"calls {res[name]['calls']} / reference {res[name]['reference']}", flush=True)
    best = max(res, key=lambda k: res[k]["f1"])
    (ROOT / "results_gate6_practice.json").write_text(json.dumps({"mode": "practice (GM12878, 3 windows)", "table": res,
                                                                 "chosen": best, "choice": CHOICES[best]}, indent=1),
                                                     encoding="utf-8")
    print("chosen:", best)


def test() -> None:
    from frozen import LOOPS_GATE6 as R
    out = ROOT / "results_gate6.json"
    if out.exists():
        sys.exit("results_gate6.json exists: Gate 6 runs once.")
    WORK.mkdir(parents=True, exist_ok=True)
    sets = {}
    for cell, windows in R["test"].items():
        ref = reference(cell)
        acc = {"chronocell": [], "chromosight": [], "mustache": [], "reference": []}
        notes = []
        for chrom, start in windows:
            i, j, v, n = window_counts(cell, chrom, start)
            r = in_window(ref, chrom, start)
            acc["reference"].append(r)
            acc["chronocell"].append(chronocell_calls(i, j, v, n, chrom, start, R["choice"]["fdr"], R["choice"]["balance"]))
            cool = write_cool(i, j, v, n, chrom, start, WORK / f"{cell}_{chrom}_{start}.mcool")
            for tool in ("chromosight", "mustache"):
                got = tool_calls(tool, cool, chrom, start)
                if isinstance(got, str):
                    notes.append(f"{chrom}:{start}: {got}")
                    acc[tool].append(None)
                else:
                    acc[tool].append(got)
        refs = pd.concat(acc["reference"])
        sets[cell] = {"notes": notes}
        for tool in ("chronocell", "chromosight", "mustache"):
            if any(x is None for x in acc[tool]):
                sets[cell][tool] = {"status": "not run on every window"}
                continue
            sets[cell][tool] = match(pd.concat(acc[tool]), refs)
        line = "  ".join(f"{t} F1 {sets[cell][t]['f1']:.3f}" for t in ("chronocell", "chromosight", "mustache") if "f1" in sets[cell][t])
        print(f"{cell}: {line}  {notes[:2]}", flush=True)
    ok = {}
    for cell, s in sets.items():
        others = [s[t]["f1"] for t in ("chromosight", "mustache") if "f1" in s[t]]
        ok[cell] = bool(others) and s["chronocell"]["f1"] >= max(others)
    verdict = "pass" if all(ok.values()) and len(ok) == len(R["test"]) else "fail"
    out.write_text(json.dumps({"mode": "test (run once; settings frozen)", "rule": R, "sets": sets, "at_least_tools": ok,
                               "verdict": verdict, "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")},
                              indent=1, default=float), encoding="utf-8")
    print("Gate 6:", verdict, ok)


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
