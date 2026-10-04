"""
Pillar 1: measured runtime and peak memory of the population model against bead count.

    python validation/scale_benchmark.py                 # -> validation/scale_benchmark.json

Every (method, N) runs in a fresh Python process, so peak memory is that run's own peak (Windows:
peak working set; Linux/macOS: max resident set size). Input: windows of the SYNTHETIC reference
chr22 contact map at 10 kb (chronocell.synthetic, planted fractal globule) up to the whole
chromosome (5,082 beads), converted to probabilities with the app's adjacent-contact assumption.
The input is synthetic: these numbers measure cost, not accuracy.

Methods
- v3.3   chronocell.ensemble.fit_from_counts (dense autograd; the app caps it at 400 beads)
- v4     chronocell.population.fit_population_from_counts with the frozen whole-chromosome settings
         (validation/TUNING.md, section 5)
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

SIZES = (200, 400, 800, 1600, 3200, 5082)
V33_MAX = 1600                       # beyond this the dense v3.3 fit is not attempted (see the JSON note)


def _peak_bytes() -> int:
    import psutil
    info = psutil.Process().memory_info()
    if hasattr(info, "peak_wset"):                         # Windows
        return int(info.peak_wset)
    import resource
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(r) if sys.platform == "darwin" else int(r) * 1024


def run_one(method: str, n: int, iterations: int | None) -> dict:
    import numpy as np
    import psutil
    from chronocell import ensemble as E, genome, physics, population as P, synthetic
    ch = genome.chrom("chr22", 10_000)
    ref = synthetic.build(ch, seed=7, b0=physics.bond_length_for(ch.resolution))
    valid = ref.valid
    # the n-bead window with the most assembled beads (the whole chromosome when n = 5,082)
    lo = 0 if n >= ch.n_bins else int(np.argmax(np.convolve(valid, np.ones(n), "valid")))
    hi = lo + n
    m = (ref.ci >= lo) & (ref.ci < hi) & (ref.cj >= lo) & (ref.cj < hi)
    ci, cj, cm = ref.ci[m] - lo, ref.cj[m] - lo, ref.cm[m]
    base = psutil.Process().memory_info().rss
    t0 = time.time()
    if method == "v3.3":
        cfg = E.EnsembleConfig(iterations=iterations or E.EnsembleConfig().iterations, device="cpu")
        res = E.fit_from_counts(ci, cj, cm, n, valid[lo:hi], b0_nm=50.0, cfg=cfg)
        fit_s, samp_s = res.config["fit_seconds"], res.config["sampling_seconds"]
    else:
        from validation.frozen import WHOLE_CHROMOSOME
        cfg = P.PopulationConfig(**(WHOLE_CHROMOSOME | ({"iterations": iterations} if iterations else {})))
        res = P.fit_population_from_counts(ci, cj, cm, n, valid[lo:hi], b0_nm=50.0, cfg=cfg)
        fit_s, samp_s = res.config["fit_seconds"], res.config["sampling_seconds"]
    total = time.time() - t0
    return {"method": method, "n_beads": n, "window_bins": [lo, hi], "contacts": int(len(ci)),
            "iterations": int(cfg.iterations), "seconds_total": round(total, 2), "seconds_fit": round(fit_s, 2),
            "seconds_sampling": round(samp_s, 2), "peak_memory_mb": round(_peak_bytes() / 2 ** 20, 1),
            "rss_before_fit_mb": round(base / 2 ** 20, 1), "best_misfit": float(res.history["best_loss"][0]),
            "contact_fit_spearman": float(res.contact_fit), "trajectory_frames": int(res.trajectories_nm.shape[0]),
            "replicas": int(res.trajectories_nm.shape[1]), "device": res.config.get("device_used", "cpu"),
            "rank": res.config.get("rank", n), "input": "synthetic reference chr22 (planted fractal globule)"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--one", nargs=2, metavar=("METHOD", "N"))
    ap.add_argument("--iterations", type=int, default=None)
    ap.add_argument("--sizes", type=int, nargs="*", default=list(SIZES))
    ap.add_argument("--out", default=str(ROOT / "scale_benchmark.json"))
    a = ap.parse_args()
    if a.one:
        print(json.dumps(run_one(a.one[0], int(a.one[1]), a.iterations)))
        return
    import torch
    rows = []
    for n in a.sizes:
        for method in ("v3.3", "v4"):
            if method == "v3.3" and n > V33_MAX:
                rows.append({"method": method, "n_beads": n, "skipped": f"dense v3.3 fit not attempted above {V33_MAX} beads"})
                continue
            cmd = [sys.executable, str(Path(__file__).resolve()), "--one", method, str(n)]
            if a.iterations:
                cmd += ["--iterations", str(a.iterations)]
            t0 = time.time()
            out = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT.parent)
            if out.returncode != 0:
                rows.append({"method": method, "n_beads": n, "error": out.stderr[-1500:],
                             "wall_seconds": round(time.time() - t0, 1)})
            else:
                rows.append(json.loads(out.stdout.strip().splitlines()[-1]))
            print(json.dumps(rows[-1]), flush=True)
    meta = {"machine": {"platform": platform.platform(), "processor": platform.processor(), "cpus": os.cpu_count(),
                        "torch": torch.__version__, "cuda": torch.cuda.is_available()},
            "note": "Synthetic input: cost only, not accuracy. Each row ran in its own process."}
    Path(a.out).write_text(json.dumps({"meta": meta, "rows": rows}, indent=1))
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
