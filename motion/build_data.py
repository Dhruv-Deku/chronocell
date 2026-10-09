"""
Everything the film shows as a fact, read from the app and the result files (nothing typed in).

    python motion/build_data.py          # writes motion/data/data.js (window.DATA = {...})

- fold: the app's reference model of chr22 (the same synthetic fractal globule the 3D workspace shows, seed 7),
  downsampled for drawing, with genomic position and the H3K27ac track per bead; the 3D workspace's measurements
  (R_g, scaling exponent, span, contour length); its contact list summed into a 128 x 128 map.
- sv: the 4D workspace's 22q11.2 deletion preset, simulated with the app's defaults (24 frames, 15 sweeps).
- tests: every row of the Scoreboard (ui/scoreboard.rows(), which reads validation/results_*.json).
- highlights: the numbers quoted on screen, each from its result file.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
VAL = ROOT / "validation"


def _load(name: str) -> dict:
    return json.loads((VAL / name).read_text(encoding="utf-8"))


def _r(a, nd=3):
    return np.round(np.asarray(a, float), nd).tolist()


MAP_BINS = 128


def contact_map(ds) -> dict:
    """The dataset's own contact list (simulated Micro-C for the reference model) summed into MAP_BINS x MAP_BINS."""
    n = ds.n
    bi = (np.asarray(ds.ci) * MAP_BINS // n).astype(int)
    bj = (np.asarray(ds.cj) * MAP_BINS // n).astype(int)
    m = np.zeros((MAP_BINS, MAP_BINS))
    np.add.at(m, (bi, bj), np.asarray(ds.cm, float))
    m = np.triu(m) + np.triu(m, 1).T
    lg = np.log1p(m)
    lev = np.round(255 * (lg / (np.percentile(lg[lg > 0], 99.5) or 1.0)).clip(0, 1) ** 0.65).astype(int)   # log, then a display gamma
    return {"bins": MAP_BINS, "bin_kb": round(ds.chrom.size / MAP_BINS / 1000), "total": int(np.asarray(ds.cm).sum()),
            "level": lev.ravel().tolist(),
            "count": np.round(m).astype(int).ravel().tolist()}


def fold_and_sv() -> tuple[dict, dict]:
    from chronocell import agent, genome, physics, scenarios as SC
    from ui.common import load_dataset
    ds = load_dataset("chr22", 7, None, None, "Auto", None, False)
    ch = ds.chrom
    x = ds.frames[0]
    rg = float(physics.radius_of_gyration(x))
    c = x.mean(0)
    scale = float(np.percentile(np.linalg.norm(x - c, axis=1), 98))
    step = 2                                                    # every other bead: 2,541 points, still the same fold
    idx = np.arange(0, ds.n, step)
    epi = np.nan_to_num(ds.epi[idx], nan=0.0)
    epi = (epi - epi.min()) / (np.ptp(epi) or 1.0)
    b0 = physics.bond_length_for(ch.resolution)
    met = agent.compute_metrics(x, ds.epi, ds.valid, b0, ch)          # the numbers the 3D workspace shows
    fold = {"chrom": ch.name, "beads": int(ds.n), "assembled": int(ds.valid.sum()), "resolution_kb": ch.resolution // 1000,
            "size_mb": round(ch.size / 1e6, 2), "contacts": int(len(ds.ci)), "rg_nm": round(rg, 1),
            "nu": round(float(met.nu), 3), "span_nm": round(float(met.span_nm)),
            "contour_um": round(float(np.linalg.norm(np.diff(x, axis=0), axis=1).sum()) / 1000, 2),
            "xyz": _r((x[idx] - c) / scale), "pos": _r(idx / (ds.n - 1)), "valid": ds.valid[idx].astype(int).tolist(),
            "epi": _r(epi, 2), "map": contact_map(ds), "chromosomes": list(genome.MAIN_CHROMOSOMES)}

    p = SC.PRESETS["22q11del"]
    params = dict(SC.preset_region(p, ch))
    for key in ("a", "b", "breakpoint"):
        if key in params:
            params[key] = int(params[key]) - ds.bin0
    traj = SC.simulate(x, ch, p.operation, params, b0, p.title, p.disease, 24, 15)
    f = traj.frames
    n = f.shape[1]
    jdx = np.arange(0, n, 3)
    disp = np.linalg.norm(f[-1] - f[0], axis=1)
    d = disp[jdx] / (np.percentile(disp, 99) or 1.0)
    a, b = params["a"], params["b"]
    sv = {"title": p.title, "disease": p.disease, "frames": int(f.shape[0]), "beads": int(n),
          "deleted_mb": round((b - a) * ch.resolution / 1e6, 2),
          "xyz": [_r((fr[jdx] - c) / scale) for fr in f], "disp": _r(np.clip(d, 0, 1), 2),
          "pos": _r(traj.bead_bin[jdx] / (ch.n_bins - 1)),
          "max_disp_nm": round(float(disp.max()), 1)}
    return fold, sv


def tests() -> dict:
    from ui import scoreboard as S
    df = S.rows()
    def short(name: str) -> str:            # "Gate Q6d: ..." -> "Q6d"; "v3.3 windows, ..." -> "v3.3"
        head = name.split(":")[0].split(",")[0].strip()
        return head.replace("Gate ", "").split(" ")[0]
    rows = [{"area": r["area"], "test": r["test"], "id": short(r["test"]), "measured": r["measured"], "status": r["status"]}
            for _, r in df.iterrows()]
    cnt = df["status"].value_counts().to_dict()
    return {"rows": rows, "n": len(rows), "pass": int(cnt.get("pass", 0)), "fail": int(cnt.get("fail", 0)),
            "other": int(cnt.get("other", 0)),
            "areas": {a: {s: int(((df["area"] == a) & (df["status"] == s)).sum()) for s in ("pass", "fail", "other")}
                      for a in S.AREAS}}


def highlights() -> dict:
    g6b = _load("results_gate6b.json")["sets"]
    loops = {cell: {t: round(v[t]["f1"], 2) for t in ("chronocell", "chromosight", "mustache")} for cell, v in g6b.items()}
    r2 = _load("results_round2.json")
    q6d = r2["q6d"]
    mol = [{"m": x["molecule"], "s": x["bond_scale"], "err": round(abs(x["error_vs_singlet_mEh"]), 3)} for x in q6d["rows"]]
    q6c = [{"m": x["molecule"], "s": x["bond_scale"], "err": round(abs(x["adapt_error_vs_singlet_mEh"]), 3)}
           for x in r2["q6c"]["rows"]]
    mit = _load("results_mitigation.json")
    q9 = [{"m": x["molecule"], "s": x["bond_scale"], "noisy": round(abs(x["methods"][mit["method_key"]]["noisy_mEh"]), 3),
           "fixed": round(abs(x["methods"][mit["method_key"]]["mitigated_mEh"]), 3)} for x in mit["rows"]]
    adm = _load("results_admet.json")
    q7b = r2["q7b"]
    gq = _load("results_gateq.json")
    q1_rate = gq["q1"]["pooled_hit_rate"]
    return {"loops": loops,
            "q6d": {"within": q6d["within"], "cases": q6d["cases"], "worst": round(q6d["max_abs_error_mEh"], 2), "rows": mol},
            "q6c": {"within": r2["q6c"]["within"], "cases": r2["q6c"]["cases"], "rows": q6c},
            "q9": {"within": mit["within"], "cases": mit["cases"], "median": round(mit["median_abs_mEh"], 2),
                   "noisy": round(mit["median_noisy_mEh"], 1), "reduction": round(mit["median_reduction"]), "rows": q9},
            "q8": {"passed": adm["passed"], "of": adm["of"]},
            "q7b": {"qaoa": round(100 * q7b["success_qaoa"]), "random": round(100 * q7b["success_random"]),
                    "complexes": q7b["complexes"]},
            "q1": round(100 * q1_rate),
            "chem_accuracy": 1.6}


def main() -> None:
    fold, sv = fold_and_sv()
    data = {"fold": fold, "sv": sv, "tests": tests(), "hi": highlights()}
    out = HERE / "data" / "data.js"
    out.parent.mkdir(exist_ok=True)
    out.write_text("window.DATA = " + json.dumps(data, separators=(",", ":")) + ";\n", encoding="utf-8")
    t = data["tests"]
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB): fold {len(fold['xyz'])} points, sv {sv['frames']} frames x "
          f"{len(sv['disp'])} points, tests {t['n']} ({t['pass']} pass, {t['fail']} fail, {t['other']} other)")
    print(json.dumps(data["hi"], indent=None)[:1500])


if __name__ == "__main__":
    main()
