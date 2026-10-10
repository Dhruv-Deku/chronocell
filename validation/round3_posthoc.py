"""
Round 3 (October 2026): can the failed structure gates be fixed? Post-hoc studies on data every gate has now seen
(practice units and the old test units of Gates 1c, 2b, 2d, 2e), leave-one-dataset-out, to decide whether a new,
pre-registered retest on new data is worth running. Nothing here changes a verdict: a failed gate stays a fail.

  sizes      Gate 1c's calibration forms, trained on every Hi-C unit but the held-out dataset's. Rule of 1c, per
             held-out dataset: CCC >= 0.8 and median size ratio 0.8-1.25.
  intervals  Hi-C-input distance ranges (Gate 2d's split-conformal method, per separation band) around the richest
             size calibration (sep2+step+depth+protocol), every other dataset weighted equally. Rule of 2d, per
             held-out dataset and band: the 90 % range holds 85-95 %, the 50 % range 40-60 %.

    python validation/round3_posthoc.py            # both; histograms cached in validation/data/cache/phaseA
"""

from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import hic_size_calibration as HS        # noqa: E402
import intervals_v2 as I                 # noqa: E402
import phase_a as A                      # noqa: E402
import phase_a_units as U                # noqa: E402

RICH = "sep2+step+depth+protocol"
OUT = ROOT / "results_round3_posthoc.json"


def group(s) -> str:
    base = s.dataset.split(":")[0]
    return {"su_chr2_parm_rep": "su_chr2", "su_chr21_rep": "su_chr21", "su_genome_amanitin": "su_genome"}.get(base, base)


def hic_specs() -> list:
    specs = [s for s in U.PRACTICE if s.input == "hic"] + U.TEST_HIC_MAIN + U.TEST_HIC_SECONDARY + U.test_genome("hic")
    return [s for s in specs if (A.CACHE / f"{s.name}.npz").exists()]


def sizes(specs: list) -> dict:
    rows = {s: HS.unit_rows(A.load(s)) for s in specs}
    groups = sorted({group(s) for s in specs})
    out = {}
    for form, cols in HS.FORMS.items():
        per = {}
        for g in groups:
            beta = HS.fit([r for s, r in rows.items() if group(s) != g], cols)
            held = [r for s, r in rows.items() if group(s) == g and s.thin == 1.0]
            sc = HS.scores(np.concatenate([HS.apply(r, cols, beta) for r in held]),
                           np.concatenate([r["truth"] for r in held]), np.concatenate([r["sep"] for r in held]))
            per[g] = {"ccc": sc["lin_ccc_nm"], "size_ratio": sc["size_ratio"],
                      "within_1c_rule": bool(sc["lin_ccc_nm"] >= 0.8 and 0.8 <= sc["size_ratio"] <= 1.25)}
        out[form] = per
        print(f"sizes {form:26s} " + "  ".join(f"{g[:12]}:{v['ccc']:.2f}/{v['size_ratio']:.2f}" for g, v in per.items()), flush=True)
    return out


def _calibrated_median(s, cols, beta) -> np.ndarray:
    u = A.load(s)
    n = len(u["median"])
    iu = np.triu_indices(n, 1)
    r = {"model": u["median"][iu].astype(np.float64), "sep": np.maximum(u["sep"][iu].astype(np.float64), 1.0),
         "step": float(u["meta"]["step_bp"]), "n_eff": float(u["meta"]["n_eff"]), "protocol": u["meta"]["protocol"]}
    m = np.zeros((n, n))
    m[iu] = HS.apply(r, cols, beta)
    return m + m.T


def intervals(specs: list) -> dict:
    cols = HS.FORMS[RICH]
    rows = {s: HS.unit_rows(A.load(s)) for s in specs}
    groups = sorted({group(s) for s in specs})
    hists = {}
    for g in groups:
        beta = HS.fit([r for s, r in rows.items() if group(s) != g], cols)
        for s in specs:
            if group(s) != g or s.thin != 1.0:
                continue
            p = A.CACHE / f"r3hist_{s.name}.npz"
            if p.exists():
                z = np.load(p)
                hists[s] = {int(k): z[k] for k in z.files}
                continue
            hists[s] = I.log_ratio_hist(s, _calibrated_median(s, cols, beta))
            np.savez(p, **{str(k): v for k, v in hists[s].items()})
            print("histogram", s.name, flush=True)
    out = {}
    for g in groups:
        q = {}
        for b in range(len(I.BANDS_BP) - 1):          # every other dataset weighted equally (its units pooled first)
            parts = []
            for g2 in groups:
                hs = [h[b] for s, h in hists.items() if group(s) == g2 and g2 != g and h[b].sum() >= I.MIN_PAIRS]
                if hs:
                    parts.append(np.mean([x / x.sum() for x in hs], axis=0))
            q[b] = np.mean(parts, axis=0) if parts else None
        ev = I.evaluate([h for s, h in hists.items() if group(s) == g], q)
        out[g] = {"bands": {I.BANDS_BP[b]: v for b, v in ev.items()},
                  "within_2d_rule": bool(all(0.85 <= v["0.9"] <= 0.95 and 0.40 <= v["0.5"] <= 0.60 for v in ev.values()))}
        print(f"intervals {g:22s} " + "  ".join(f"b{b}:{100 * v['0.9']:.0f}/{100 * v['0.5']:.0f}" for b, v in ev.items()),
              flush=True)
    return out


def main() -> None:
    specs = hic_specs()
    out = {"made": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "note": "Post-hoc on seen data (leave-one-dataset-out); decides whether a retest is worth pre-registering.",
           "units": len(specs), "sizes": sizes(specs), "intervals_form": RICH, "intervals_hic": intervals(specs)}
    OUT.write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
