"""
The two accuracy scores, kept separate everywhere (UI, PDF dossier, JSON report).

* Contact-map fit: how well a model reproduces the contact data it was BUILT FROM. It is measured
  live on the user's own window. It shows the optimisation converged. It is not evidence that
  the 3D model is right: any good optimiser scores high on its own input.
* Microscopy accuracy (benchmark): how well the METHOD predicts distances measured by imaging in
  cells it never saw (validation/, Bintu et al. 2018, held-out test datasets). It is a property of
  the method, read from validation/results.json. It is not measured on the user's data, which has
  no imaging ground truth.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

RESULTS = Path(__file__).resolve().parent.parent / "validation" / "results.json"

CONTACT_FIT_DEFINITION = ("Spearman rank correlation between the model's contacts (or closeness) and the input "
                          "contact map it was built from. Shows convergence; not evidence of accuracy.")
MICROSCOPY_DEFINITION = ("Share of the reproducible 3D folding pattern recovered, relative to how well the experiment "
                         "agrees with itself; trend-removed Spearman rho vs held-out chromatin-tracing distances "
                         "(Bintu et al., Science 2018). Overall = sum(model) / sum(ceiling) over 3 test datasets x 3 splits.")


def _rank(v: np.ndarray) -> np.ndarray:
    return np.argsort(np.argsort(v)).astype(np.float64)


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, np.float64), np.asarray(b, np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    ra, rb = _rank(a[ok]), _rank(b[ok])
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def contact_fit_structure(coords: np.ndarray, ci: np.ndarray, cj: np.ndarray, cm: np.ndarray) -> float:
    """Contact-map fit of a single structure: Spearman(-distance, contact count) over observed pixels."""
    x = np.asarray(coords, dtype=np.float64)
    ci, cj = np.asarray(ci, np.int64), np.asarray(cj, np.int64)
    if ci.size < 3:
        return float("nan")
    d = np.linalg.norm(x[ci] - x[cj], axis=1)
    return spearman(-d, np.asarray(cm, np.float64))


def load_benchmark(path: Path | str = RESULTS) -> dict | None:
    """Summary of the held-out microscopy benchmark, or None if validation has not been run."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        s = data["summary"]
    except (OSError, ValueError, KeyError, TypeError):
        return None
    out = {"source": "validation/results.json", "definition": MICROSCOPY_DEFINITION, "models": {}}
    for model in ("single_structure_v3_2", "ensemble_v3_3"):
        if model not in s.get("overall", {}):
            continue
        out["models"][model] = {
            "overall_percent_of_ceiling": round(float(s["overall"][model]["percent_of_ceiling"]), 1),
            "per_dataset_percent_of_ceiling": {name: round(float(d[model]["percent_of_ceiling"]), 1)
                                               for name, d in s["datasets"].items() if model in d},
        }
    v4 = _gate1_benchmark(Path(path).parent)
    if v4:
        out["models"]["population_v4"] = v4
    pred = _predictor_benchmark(Path(path).parent)
    if pred:
        out["models"]["predicted_sequence_ctcf"] = pred
    return out if out["models"] else None


def _predictor_benchmark(folder: Path) -> dict | None:
    """The no-contact-data predictor's held-out benchmark (Gate 5): % of ceiling per test dataset, all pairs."""
    try:
        r = json.loads((folder / "results_predictor.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    per, num, den = {}, 0.0, 0.0
    for k in r["settings"]["test"]:
        v = r["summary"].get(k, {}).get("sequence + CTCF", {}).get("all_pairs", {}).get("predictor")
        if v:
            per[k] = round(float(v["percent_of_ceiling"]), 1)
            num += float(v["percent_of_ceiling"]) * float(v["ceiling"])
            den += float(v["ceiling"])
    if not per:
        return None
    return {"overall_percent_of_ceiling": round(num / den, 1) if den else None,
            "per_dataset_percent_of_ceiling": per, "verdict": r["verdict"],
            "definition": ("Prediction from sequence + CTCF peaks with no contact data (Gate 5): trend-removed Spearman rho "
                           "vs held-out tracing medians as % of the half-A vs half-B ceiling, all locus pairs; overall = "
                           "sum over test datasets weighted by their ceilings.")}


def _gate1_benchmark(folder: Path) -> dict | None:
    """The v4 whole-chromosome model's held-out benchmark (Gate 1, Su et al. 2020 chr21 tracing)."""
    per, hic = {}, []
    for fname, label in (("results_gate1.json", "IMR-90 chr21"), ("results_gate1_rep.json", "IMR-90 chr21 replicate")):
        try:
            s = json.loads((folder / fname).read_text(encoding="utf-8"))["datasets"][0]["summary"]
        except (OSError, ValueError, KeyError, IndexError):
            continue
        h = s.get("hic_input", {}).get("all_pairs", {}).get("whole_v4[literature_b0]")
        if h:
            per[f"{label} · sequencing Hi-C input (Rao 2014)"] = round(float(h["percent_of_ceiling"]), 1)
            hic.append(float(h["percent_of_ceiling"]))
        per[f"{label} · imaging-derived contacts"] = round(float(s["imaging_input"]["all_pairs"]["whole_v4"]["percent_of_ceiling"]), 1)
    if not per:
        return None
    headline = np.mean(hic) if hic else np.mean(list(per.values()))
    return {"overall_percent_of_ceiling": round(float(headline), 1),
            "per_dataset_percent_of_ceiling": per,
            "definition": ("Trend-removed Spearman rho vs held-out tracing medians (Su et al. 2020, chr21, 651 loci), as % "
                           "of the half-A vs half-B ceiling, all locus pairs, 3 splits. Headline: sequencing Hi-C input, the "
                           "app's input type; with Hi-C input the absolute sizes (nm) are not calibrated.")}


VALIDATION = RESULTS.parent


def interval_evidence() -> dict | None:
    """How well the population's stated 90 % intervals held on held-out single-cell distances, by input
    type, read from the result files: imaging-derived contacts (Gate 2, validation/results_calibration.json)
    and sequencing Hi-C (Gate 2b, validation/results_calibration_hic.json: raw, with its own practice-fitted
    recalibration, with the imaging recalibration, and the pre-registered verdict; before Gate 2b was run,
    the benchmark's raw and imaging-recalibrated coverage). Ranges are min-max over datasets, in %. None
    if no calibration test has been run."""
    try:
        c = json.loads((VALIDATION / "results_calibration.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    raw = [100 * r["coverage"][2] for r in c["rows"]]
    rec = [100 * r["coverage_recalibrated"][2] for r in c["rows"] if r.get("coverage_recalibrated")]
    out = {"imaging_raw_90": (min(raw), max(raw)), "imaging_recalibrated_90": (min(rec), max(rec)) if rec else None,
           "hic_raw_90": None, "hic_recalibrated_90": None, "hic_source": None, "hic_verdict": None,
           "hic_imaging_recalibration_90": None}
    try:                                   # Gate 2b: the pre-registered test of intervals with Hi-C input
        h = json.loads((VALIDATION / "results_calibration_hic.json").read_text(encoding="utf-8"))
        hr = [100 * r["coverage"][2] for r in h["rows"]]
        hc = [100 * r["coverage_recalibrated"][2] for r in h["rows"]]
        hi = [100 * r["coverage_imaging_recalibration"][2] for r in h["rows"] if r.get("coverage_imaging_recalibration")]
        out.update(hic_raw_90=(min(hr), max(hr)), hic_recalibrated_90=(min(hc), max(hc)),
                   hic_imaging_recalibration_90=(min(hi), max(hi)) if hi else None,
                   hic_source="held-out test, Gate 2b", hic_verdict=h["verdict"])
        return out
    except (OSError, ValueError, KeyError):
        pass
    for fname, label in (("benchmark/results.json", "held-out test"), ("benchmark/results_practice.json", "practice")):
        try:
            b = json.loads((VALIDATION / fname).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        hr, hc = [], []
        for ds in b["summary"].values():
            for model, v in ds.get("hic", {}).get("all_pairs", {}).items():
                if model in ("v3_3_windowed", "v4_whole") and "coverage" in v:
                    hr.append(100 * v["coverage"]["90"])
                    if "coverage_recalibrated" in v:
                        hc.append(100 * v["coverage_recalibrated"]["90"])
        if hr:
            out.update(hic_raw_90=(min(hr), max(hr)), hic_imaging_recalibration_90=(min(hc), max(hc)) if hc else None,
                       hic_source=label)
            break
    return out


def interval_input_kind(res) -> str:
    """Which interval test applies to a population result: "hic" (built from sequencing counts),
    "predicted" (fitted to a sequence + CTCF prediction), "imaging" (imaging-derived frequencies, as in
    validation) or "other" (another distance map)."""
    src = str(getattr(res, "config", {}).get("input", ""))
    if src.startswith("predicted"):
        return "predicted"
    if src == "sequencing counts":
        return "hic"
    return "other" if src else "imaging"


def interval_note(ev: dict | None, kind: str) -> str | None:
    """The measured coverage of the stated intervals for this input kind (HTML-safe text, numbers from the
    result files), or None if nothing was measured."""
    if not ev:
        return None
    rng = lambda r: f"{r[0]:.0f}–{r[1]:.0f} %"                              # noqa: E731
    if kind == "hic":
        if ev.get("hic_verdict") == "pass":
            return (f"Held-out test with sequencing Hi-C input (Gate 2b, passed): the recalibrated 90 % interval held "
                    f"{rng(ev['hic_recalibrated_90'])} of real single-cell distances.")
        if ev.get("hic_verdict") == "fail":
            far = ev["hic_recalibrated_90"][1] < 75
            return (f"<b>Coverage {'far ' if far else ''}below nominal with sequencing Hi-C input.</b> On held-out "
                    f"data (Gate 2b) a stated 90 % interval held only {rng(ev['hic_raw_90'])} of real single-cell "
                    f"distances, and {rng(ev['hic_recalibrated_90'])} after a recalibration fitted on practice Hi-C "
                    "(pre-registered pass: 83–97 %; not met). Sizes from Hi-C are not calibrated in nm (Gate 1b). Read "
                    "the range as the model's cell-to-cell spread, not as a 90 % range for real cells.")
        if ev.get("hic_raw_90"):
            return (f"With sequencing Hi-C input, stated 90 % intervals held {rng(ev['hic_raw_90'])} of real single-cell "
                    f"distances ({ev['hic_source']}): far below nominal. Read the range as the model's cell-to-cell "
                    "spread, not as a 90 % range for real cells.")
        return None
    if kind == "predicted":
        return ("These intervals belong to a population fitted to a prediction from sequence + CTCF; their coverage was "
                "not tested (Gate 5 tested the median distances only).")
    if kind == "imaging":
        note = (f"Held-out check of these intervals (Gate 2): with imaging-derived contacts, stated 90 % intervals held "
                f"{rng(ev['imaging_raw_90'])} of real single-cell distances")
        if ev.get("imaging_recalibrated_90"):
            note += f", {rng(ev['imaging_recalibrated_90'])} recalibrated"
        return note + ". Read the interval as a lower bound on cell-to-cell spread."
    return None


def sv_evidence() -> dict | None:
    """The structural-variant test (validation/sv_validation.py, pre-registered) as measured, or None if it
    has not been run. Nothing here is typed in: every number is read from validation/results_sv.json."""
    try:
        r = json.loads((VALIDATION / "results_sv.json").read_text(encoding="utf-8"))
        return {"verdict": r["verdict"], "validated": r["verdict"].startswith("validated"),
                "event": r["settings"]["deleted_source"], "window": r["window"], "pairs": r["spanning_pairs"],
                "spearman": r["spearman"], "ci95": r["ci95_block_bootstrap"], "criteria": r["criteria"]}
    except (OSError, ValueError, KeyError):
        return None


def method_evidence() -> dict:
    """Headline numbers of every held-out test that has been run, read from validation/*.json (measured
    values only; a test that has not been run is reported as such, never filled in)."""
    out: dict = {}
    bench = load_benchmark()
    out["v3_3_microscopy_benchmark"] = bench["models"] if bench else "not run"
    for key, fname in (("gate1_whole_chromosome", "results_gate1.json"),
                       ("gate1_whole_chromosome_replicate", "results_gate1_rep.json")):
        try:
            d = json.loads((VALIDATION / fname).read_text(encoding="utf-8"))["datasets"][0]
            s = d["summary"]
            out[key] = {"dataset": d["dataset"], "loci": d["loci"],
                        **{f"{inp}.{m}.{model}": round(v["percent_of_ceiling"], 1)
                           for inp, rows in s.items() for m, models in rows.items() for model, v in models.items()
                           if m in ("within_tiles", "all_pairs")}}
        except (OSError, ValueError, KeyError, IndexError):
            out[key] = "not run"
    try:
        p = json.loads((VALIDATION / "results_perturbation.json").read_text(encoding="utf-8"))
        out["gate4_cohesin_depletion"] = {"dataset": p["dataset"],
                                          **{f"{k}.change_spearman": round(v["change_spearman"], 3)
                                             for k, v in p["summary"].items() if isinstance(v, dict)
                                             and np.isfinite(v.get("change_spearman", np.nan))}}
    except (OSError, ValueError, KeyError):
        out["gate4_cohesin_depletion"] = "not run"
    sv = sv_evidence()
    out["gate4_structural_variants"] = ({"verdict": sv["verdict"], "spanning_pairs": sv["pairs"],
                                         **{f"{k}.spearman": round(v, 3) for k, v in sv["spearman"].items()}}
                                        if sv else "not run")
    try:
        c = json.loads((VALIDATION / "results_calibration.json").read_text(encoding="utf-8"))
        out["gate2_calibration"] = {r["dataset"]: {"levels": r["levels"], "coverage": [round(v, 3) for v in r["coverage"]],
                                                   "coverage_recalibrated": [round(v, 3) for v in r.get("coverage_recalibrated", [])]}
                                    for r in c["rows"]}
    except (OSError, ValueError, KeyError):
        out["gate2_calibration"] = "not run"
    return out


def two_scores(model: str, contact_fit: float | None) -> dict:
    """The two labelled scores for one displayed model ("single_structure_v3_2" / "ensemble_v3_3" / "input")."""
    bench = load_benchmark()
    micro = None
    if bench and model in bench["models"]:
        m = bench["models"][model]
        micro = {"overall_percent_of_ceiling": m["overall_percent_of_ceiling"],
                 "per_dataset_percent_of_ceiling": m["per_dataset_percent_of_ceiling"],
                 "scope": "method-level benchmark on held-out imaging data; not measured on this structure",
                 "definition": MICROSCOPY_DEFINITION, "source": bench["source"]}
    fit = None if contact_fit is None or not np.isfinite(contact_fit) else round(float(contact_fit), 4)
    return {"contact_map_fit": {"value": fit, "scope": "measured on this window's input contacts",
                                "definition": CONTACT_FIT_DEFINITION},
            "microscopy_accuracy": micro}
