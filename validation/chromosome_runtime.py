"""
Pillar 6: measured runtime and peak memory of the app's population model on every main chromosome
of every bundled assembly, at the app's default resolution (genome.default_resolution, at most
6,000 beads).

    python validation/chromosome_runtime.py                       # -> validation/chromosome_runtime.json
    python validation/chromosome_runtime.py --assemblies hg38     # one assembly
    python validation/chromosome_runtime.py --resume              # continue an interrupted run (rows are saved as they finish)

Every chromosome runs in a fresh Python process, so peak memory is that run's own peak (Windows: peak
working set; Linux/macOS: max resident set size). Input: the SYNTHETIC reference for that chromosome
(chronocell.synthetic: planted fractal globule, band-informed tracks, Poisson contacts) at its default
resolution, fitted with the app's settings (population.config_for: the frozen whole-chromosome
defaults; v3.3 at <= 400 beads). The input is synthetic: these numbers measure cost, not accuracy.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))


from scale_benchmark import _system_load      # noqa: E402


def run_one(assembly: str, chrom: str) -> dict:
    import numpy as np
    import psutil
    from chronocell import ensemble as E, genome, physics, population as P, synthetic
    from scale_benchmark import _peak_bytes
    res_bp = genome.default_resolution(chrom, assembly)
    ch = genome.chrom(chrom, res_bp, assembly)
    t_build = time.time()
    ref = synthetic.build(ch, seed=7, b0=physics.bond_length_for(ch.resolution))
    build_s = time.time() - t_build
    n = ch.n_bins
    base = psutil.Process().memory_info().rss
    t0 = time.time()
    if n <= P.V33_MAX_BEADS:
        res = E.fit_from_counts(ref.ci, ref.cj, ref.cm, n, ref.valid, b0_nm=physics.bond_length_for(res_bp))
    else:
        res = P.fit_population_from_counts(ref.ci, ref.cj, ref.cm, n, ref.valid, b0_nm=physics.bond_length_for(res_bp),
                                           cfg=P.config_for(n))
    total = time.time() - t0
    return {"assembly": assembly, "chrom": chrom, "resolution_bp": int(res_bp), "bins": int(n),
            "beads": int(np.asarray(ref.valid, bool).sum()), "contacts": int(len(ref.ci)),
            "model": res.config.get("model", "ensemble_v3_3"), "rank": res.config.get("rank", n),
            "iterations": int(res.config.get("iterations", 0) or 0),
            "seconds_fit": round(float(res.config["fit_seconds"]), 2),
            "seconds_sampling": round(float(res.config["sampling_seconds"]), 2), "seconds_total": round(total, 2),
            "seconds_synthetic_input": round(build_s, 2), "peak_memory_mb": round(_peak_bytes() / 2 ** 20, 1),
            "rss_before_fit_mb": round(base / 2 ** 20, 1), "best_misfit": float(res.history["best_loss"][0]),
            "device": res.config.get("device_used", "cpu"), "input": "synthetic reference (planted fractal globule)"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--one", nargs=2, metavar=("ASSEMBLY", "CHROM"))
    ap.add_argument("--assemblies", nargs="*", default=None)
    ap.add_argument("--out", default=str(ROOT / "chromosome_runtime.json"))
    ap.add_argument("--resume", action="store_true", help="keep the finished rows of --out and run only the rest")
    a = ap.parse_args()
    if a.one:
        sys.path.insert(0, str(ROOT))
        print(json.dumps(run_one(*a.one)))
        return
    import psutil
    import torch
    from chronocell import genome
    names = a.assemblies or sorted(genome.assemblies())
    rows, wall, discarded = [], 0.0, []
    if a.resume and Path(a.out).exists():
        old = json.loads(Path(a.out).read_text())
        rows = [r for r in old["rows"] if "error" not in r]
        wall = float(old["meta"].get("wall_seconds", 0.0))
        discarded = old["meta"].get("discarded_rows", [])      # kept on record; those chromosomes are re-measured
    done = {(r["assembly"], r["chrom"]) for r in rows}
    vm = psutil.virtual_memory()
    meta = {"machine": {"platform": platform.platform(), "processor": platform.processor(), "cpus": os.cpu_count(),
                        "ram_gb": round(vm.total / 2 ** 30, 1), "python": platform.python_version(),
                        "torch": torch.__version__, "cuda": torch.cuda.is_available(),
                        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None},
            "assemblies": names, "wall_seconds": wall, "discarded_rows": discarded,
            "planned": [[asm, c] for asm in names for c in genome.main_chromosomes(asm)],
            "note": ("Synthetic input (one planted reference per chromosome at the app's default resolution, at most "
                     "6,000 beads), fitted with the app's settings; cost only, not accuracy. Each row ran in its own "
                     "process, one at a time; the last column is the whole-machine CPU load measured just before the row "
                     "(other work slows timings).")}
    meta["status"] = "partial: interrupted before the last chromosome; continue with --resume"   # until the end
    for asm in names:
        for chrom in genome.main_chromosomes(asm):
            if (asm, chrom) in done:
                continue
            cmd = [sys.executable, str(Path(__file__).resolve()), "--one", asm, chrom]
            load = _system_load()
            t0 = time.time()
            out = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT.parent)
            if out.returncode != 0:
                rows.append({"assembly": asm, "chrom": chrom, "error": out.stderr[-1500:],
                             "wall_seconds": round(time.time() - t0, 1)})
            else:
                rows.append(json.loads(out.stdout.strip().splitlines()[-1]))
            rows[-1]["system_cpu_percent_before"] = load
            meta["wall_seconds"] = round(meta["wall_seconds"] + time.time() - t0, 1)
            print(json.dumps(rows[-1]), flush=True)
            Path(a.out).write_text(json.dumps({"meta": meta, "rows": rows}, indent=1))     # checkpoint every row
    meta["status"] = "complete" if {(r["assembly"], r["chrom"]) for r in rows if "error" not in r} >= {
        tuple(x) for x in meta["planned"]} else "partial: some chromosomes failed; see the error rows"
    order = {(asm, c): k for k, (asm, c) in enumerate((asm, c) for asm in names for c in genome.main_chromosomes(asm))}
    rows.sort(key=lambda r: order.get((r["assembly"], r["chrom"]), 10 ** 6))
    Path(a.out).write_text(json.dumps({"meta": meta, "rows": rows}, indent=1))
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
