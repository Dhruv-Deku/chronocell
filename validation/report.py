"""
Every numeric table in validation/RESULTS.md, generated from the result files (nothing typed in).

    python validation/report.py            # rewrites the generated blocks of validation/RESULTS.md

RESULTS.md holds hand-written text and generated blocks between markers:
    <!-- BEGIN generated:NAME -->  ...  <!-- END generated:NAME -->
A block whose result file does not exist yet says so ("not run").
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS_MD = ROOT / "RESULTS.md"


def _load(name: str):
    try:
        return json.loads((ROOT / name).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _pct(v, nd=1):
    return "—" if v is None else f"{v:.{nd}f} %"


def _f(v, nd=3):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return "—"
    return f"{x:.{nd}f}" if x == x and abs(x) != float("inf") else "—"


def _ci(ci, nd=1, pct=True):
    if not ci:
        return ""
    return f" [{ci[0]:.{nd}f}, {ci[1]:.{nd}f}]"


# --------------------------------------------------------------------------------------------- Gate 1
def gate1() -> str:
    rows = []
    for fname in ("results_gate1.json", "results_gate1_rep.json"):
        r = _load(fname)
        if not r:
            rows.append(f"_{fname}: not run._")
            continue
        for d in r["datasets"]:
            s = d["summary"]
            rows += [f"**{d['dataset']}** ({d['loci']} loci, {d.get('copies', '—')} copies; mean over "
                     f"{len(d.get('runs', [])) or '—'} splits; {r['mode']})", "",
                     "| Input | Pairs | Model | % of ceiling | Raw ρ | CCC (nm) | Size ratio |", "|---|---|---|---|---|---|---|"]
            for inp, sets in s.items():
                for mname, models in sets.items():
                    for model, v in models.items():
                        rows.append(f"| {inp.replace('_input', '')} | {mname.replace('_', ' ')} | {model} | "
                                    f"{_pct(v['percent_of_ceiling'])} | {_f(v.get('raw_spearman'))} | "
                                    f"{_f(v.get('lin_ccc_nm'))} | {_f(v.get('median_scale'), 2)} |")
            diffs = []
            for run in d.get("runs", []):
                w = run.get("imaging_input", {}).get("within_tiles", {})
                if "whole_v4" in w and "windowed_v3_3" in w:
                    diffs.append(w["whole_v4"]["percent_of_ceiling"] - w["windowed_v3_3"]["percent_of_ceiling"])
            if diffs:
                rows += ["", "Whole − windowed, imaging input, within tiles, per split: "
                         + ", ".join(f"{x:+.2f}" for x in diffs) + " points."]
            rows.append("")
    return "\n".join(rows)


# --------------------------------------------------------------------------------------------- Gate 2
def gate2() -> str:
    c = _load("results_calibration.json")
    if not c:
        return "_results_calibration.json: not run._"
    out = [f"{c['mode']}; recalibration fitted on {', '.join(c['calibration_used']['fitted_on'])}.", "",
           "| Test dataset | Loci | Model | Stated 50 / 80 / 90 %: raw | Recalibrated | Reliability vs error (ρ) |",
           "|---|---|---|---|---|---|"]
    for r in c["rows"]:
        raw = " / ".join(f"{100 * v:.0f}" for v in r["coverage"])
        rec = " / ".join(f"{100 * v:.0f}" for v in r.get("coverage_recalibrated", []))
        out.append(f"| {r['dataset']} | {r['loci']} | {r['model']} | {raw} % | {rec} % | "
                   f"{r['reliability_vs_error_spearman']:+.2f} |")
    bins = sorted({tuple(b["sep_bp"]) for r in c["rows"] for b in r["by_separation"]})
    out += ["", "Coverage of the stated 90 % interval by genomic separation (raw → recalibrated):", "",
            "| Test dataset | " + " | ".join(f"{lo / 1e6:g}–{hi / 1e6:g} Mb" for lo, hi in bins) + " |",
            "|---|" + "---|" * len(bins)]
    for r in c["rows"]:
        have = {tuple(b["sep_bp"]): b for b in r["by_separation"]}
        cells = []
        for key in bins:
            b = have.get(key)
            if b is None:
                cells.append("—")
            elif "coverage_recalibrated" in b:
                cells.append(f"{100 * b['coverage'][2]:.0f} → {100 * b['coverage_recalibrated'][2]:.0f} %")
            else:
                cells.append(f"{100 * b['coverage'][2]:.0f} %")
        out.append(f"| {r['dataset']} | " + " | ".join(cells) + " |")
    return "\n".join(out)


def gate2b() -> str:
    h = _load("results_calibration_hic.json")
    if not h:
        return "_results_calibration_hic.json: not run (python validation/calibration.py --test --input hic)._"
    lo90, hi90 = h["rule"]["pass_90"]
    lo50, hi50 = h["rule"]["pass_50"]
    out = [f"{h['mode']}; Hi-C recalibration fitted on {', '.join(h['calibration_used']['fitted_on'])}. Pass (pre-registered): "
           f"recalibrated 90 % interval holds {100 * lo90:.0f}–{100 * hi90:.0f} % and 50 % interval {100 * lo50:.0f}–"
           f"{100 * hi50:.0f} % on every test dataset.", "",
           "| Test dataset | Loci | Model | Size ratio (model / measured) | Stated 50 / 80 / 90 %: raw | With the imaging "
           "recalibration | With the Hi-C recalibration | Within the rule |", "|---|---|---|---|---|---|---|---|"]
    f = lambda v: " / ".join(f"{100 * x:.0f}" for x in v) + " %"  # noqa: E731
    for r in h["rows"]:
        out.append(f"| {r['dataset']} | {r['loci']} | {r['model']} | {r['median_scale_model_over_real']:.2f} | "
                   f"{f(r['coverage'])} | {f(r.get('coverage_imaging_recalibration', []))} | {f(r['coverage_recalibrated'])} | "
                   f"{'yes' if h['per_dataset_pass'][r['dataset']] else 'no'} |")
    r0 = h["rows"][0]
    out += ["", f"Width of the stated 90 % interval (upper / lower bound): raw {r0['width_90_raw']:.2f}, Hi-C recalibrated "
                f"{r0['width_90_recalibrated']:.2f}. {sum(h['per_dataset_pass'].values())} of {len(h['per_dataset_pass'])} "
                f"test datasets within the rule. Verdict: **{h['verdict']}**."]
    p = _load("results_calibration_hic_practice.json")
    if p:
        out += ["", "Practice (in-sample, after the fit): " + "; ".join(
            f"{r['dataset']} {f(r['coverage_recalibrated'])}" for r in p["rows"]) + "."]
    return "\n".join(out)


def gate2c() -> str:
    r = _load("results_reliability.json")
    p = _load("results_reliability_practice.json")
    out = []
    if p:
        names = [k for k in ("input_se", "misfit", "combined") if k in p["rows"][0]]
        out += ["Practice (all candidates; in-sample choice): pair-level stratified Spearman / bead-level Spearman.", "",
                "| Practice dataset | Input | " + " | ".join(names) + " |", "|---|---|" + "---|" * len(names)]
        for row in p["rows"]:
            out.append(f"| {row['dataset']} | {row['input']} | " + " | ".join(
                f"{row[k]['pair_stratified_spearman']:+.3f} / {row[k]['bead_spearman']:+.3f}" for k in names) + " |")
        out.append("")
    if not r:
        return "\n".join(out + ["_results_reliability.json: not run (python validation/reliability.py --test)._"])
    out += [f"Test (run once): score **{r['rule']['score']}**; pass on every test dataset: Spearman ≥ "
            f"{r['rule']['min_spearman']:.2f} and its 95 % interval above 0.", "",
            "| Test dataset | Input | Pair level: stratified ρ [95 %] | Bead level: ρ [95 %] |", "|---|---|---|---|"]
    for row in r["rows"]:
        v = row[r["rule"]["score"]]
        out.append(f"| {row['dataset']} | {row['input']} | {v['pair_stratified_spearman']:+.3f} "
                   f"[{v['pair_ci95'][0]:+.3f}, {v['pair_ci95'][1]:+.3f}] | {v['bead_spearman']:+.3f} "
                   f"[{v['bead_ci95'][0]:+.3f}, {v['bead_ci95'][1]:+.3f}] |")
    out += ["", "Verdicts: " + "; ".join(f"{k.replace('_', ' ')} **{v['verdict']}** ({v['passing']} of {v['datasets']})"
                                       for k, v in r["verdict"].items()) + "."]
    return "\n".join(out)


# --------------------------------------------------------------------------------------------- Phase A
def gate1c() -> str:
    r = _load("results_hic_size.json")
    p = _load("results_hic_size_practice.json")
    out = []
    if p:
        out += ["Practice (leave-one-dataset-out; CCC of the held-out group / its size ratio):", "",
                "| Form | " + " | ".join(next(iter(p["forms"].values()))["per_group"]) + " | Mean CCC |",
                "|---|" + "---|" * (len(next(iter(p["forms"].values()))["per_group"]) + 1)]
        for f, v in p["forms"].items():
            out.append(f"| {f}{' (chosen)' if f == p['chosen'] else ''} | " + " | ".join(
                f"{g['lin_ccc_nm']:.2f} / {g['size_ratio']:.2f}" for g in v["per_group"].values()) + f" | {v['mean_ccc']:.3f} |")
        out.append("")
    if not r:
        return "\n".join(out + ["_results_hic_size.json: not run (python validation/hic_size_calibration.py --test)._"])
    R = r["rule"]
    out += [f"Test (run once): pass on every main set: CCC ≥ {R['min_ccc']}, size ratio {R['size_ratio'][0]}–"
            f"{R['size_ratio'][1]}, pattern unchanged.", "",
            "| Set | Tier | CCC (nm): uncalibrated → calibrated | Size ratio | Trend-removed ρ | Within the rule |",
            "|---|---|---|---|---|---|"]
    for tier in ("main", "secondary"):
        for name, v in r["sets"][tier].items():
            b, a = v["uncalibrated"], v["calibrated"]
            out.append(f"| {name} | {tier} | {b['lin_ccc_nm']:.3f} → {a['lin_ccc_nm']:.3f} | {b['size_ratio']:.2f} → "
                       f"{a['size_ratio']:.2f} | {a['trend_removed_rho']:.3f} | {'yes' if v['within_rule'] else 'no' if tier == 'main' else '—'} |")
    out += ["", f"Verdict: **{r['verdict']}**."]
    return "\n".join(out)


def gate2d() -> str:
    r = _load("results_intervals_v2.json")
    if not r:
        return "_results_intervals_v2.json: not run (python validation/intervals_v2.py --test)._"
    bands = r["bands_bp"]
    lab = {str(k): f"{bands[k] / 1e6:g}–{bands[k + 1] / 1e6:g} Mb" for k in range(len(bands) - 1)}
    out = [f"Test (run once): variants {', '.join(f'{k}: {v}' for k, v in r['intervals'].items())}; pass on every set and "
           f"band: 90 % ranges hold {100 * r['rule']['cov90'][0]:.0f}–{100 * r['rule']['cov90'][1]:.0f} %, 50 % ranges "
           f"{100 * r['rule']['cov50'][0]:.0f}–{100 * r['rule']['cov50'][1]:.0f} %.", "",
           "| Set · input | Coverage at 90 % / 50 % by separation band (width of the 90 % range, upper / lower) | Within the rule |",
           "|---|---|---|"]
    for name, v in r["sets"].items():
        cells = "; ".join(f"{lab[b]}: {100 * x['0.9']:.0f} / {100 * x['0.5']:.0f} % (×{x['width90']:.1f})"
                          for b, x in v["bands"].items())
        out.append(f"| {name} | {cells} | {'yes' if v['within_rule'] else 'no'} |")
    out += ["", "Verdicts: " + "; ".join(f"{k} input **{v}**" for k, v in r["verdict"].items()) + "."]
    return "\n".join(out)


def gate2e() -> str:
    p = _load("results_reliability_v2_practice.json")
    r = _load("results_reliability_v2.json")
    out = []
    if p:
        cands = list(next(iter(p["table"].values())))
        out += ["Practice (pair-weighted stratified ρ per group):", "", "| Practice group · input | " + " | ".join(cands) + " |",
                "|---|" + "---|" * len(cands)]
        for g, v in p["table"].items():
            out.append(f"| {g} | " + " | ".join(f"{v[c]['rho']:+.3f}" for c in cands) + " |")
        out.append("")
    if not r:
        return "\n".join(out + ["_results_reliability_v2.json: not run (python validation/reliability_v2.py --test)._"])
    sc = r["rule"]["score"]
    sc_txt = ", ".join(f"{k} input **{v}**" for k, v in sc.items()) if isinstance(sc, dict) else f"**{sc}**"
    out += [f"Test (run once): score {sc_txt}; pass per input type: ρ ≥ {r['rule']['min_rho']:.2f} with "
            "the 95 % interval above 0 on every new test set.", "",
            "| New test set · input | Score | ρ [95 %] | Within the rule |", "|---|---|---|---|"]
    for k, v in r["sets"].items():
        out.append(f"| {k} | {v.get('score', sc)} | {v['rho']:+.3f} [{v['ci95'][0]:+.3f}, {v['ci95'][1]:+.3f}] | "
                   f"{'yes' if v['within_rule'] else 'no'} |")
    out += ["", "Verdicts: " + "; ".join(f"{k} input **{v}**" for k, v in r["verdict"].items()) + "."]
    return "\n".join(out)


def gate3b() -> str:
    p = _load("results_learned_correction_practice.json")
    r = _load("results_gate3b.json")
    out = []
    if p:
        out += ["Practice (leave-one-dataset-out; mean held-out trend-removed ρ / CCC):", "",
                "| Candidate | Pattern ρ | CCC |", "|---|---|---|"]
        for c, v in p["table"].items():
            out.append(f"| {c}{' (chosen)' if c == p['chosen'] else ''} | {v['mean_pattern_rho']:+.4f} | {v['mean_ccc']:.3f} |")
        out.append("")
    if not r:
        sys.path.insert(0, str(ROOT))
        try:
            from frozen import LEARNED_CORRECTION as LC
        except ImportError:
            LC = {}
        if LC.get("status") == "not run":
            return "\n".join(out + [f"Test: **not run** — {LC['reason']}. Pre-registered rule (frozen.LEARNED_CORRECTION): "
                                     f"best on ≥ {100 * LC['hic_best_fraction']:.0f} % of Hi-C units and within "
                                     f"{LC['imaging_max_loss_points']:g} point of the current model on every imaging unit."])
        return "\n".join(out + ["_results_gate3b.json: not run (python validation/learned_correction.py --test)._"])
    out += [r["summary_line"], "", "| Unit · input | Learned correction | Best other method | Current model | Within the rule |",
            "|---|---|---|---|---|"]
    for u in r["units"]:
        out.append(f"| {u['unit']} · {u['input']} | {u['learned']:.1f} % | {u['best_other_name']} {u['best_other']:.1f} % | "
                   f"{u['current']:.1f} % | {'yes' if u['ok'] else 'no'} |")
    out += ["", f"Verdict: **{r['verdict']}**."]
    return "\n".join(out)


def gate5b() -> str:
    p = _load("results_predictor_v2_practice.json")
    r = _load("results_predictor_v2.json")
    out = []
    if p:
        best = {}
        for k, v in p["table"].items():
            fs, lam = k.split("|")
            if fs not in best or v["mean"] > best[fs][1]:
                best[fs] = (lam, v["mean"])
        out += ["Practice (leave-one-dataset-out over K562 / HCT116 regions; best ridge per feature set):", "",
                "| Feature set | Ridge | Held-out trend-removed ρ (mean) |", "|---|---|---|"]
        for fs, (lam, m) in best.items():
            out.append(f"| {fs}{' (chosen)' if fs == p['chosen']['feature_set'] else ''} | {lam} | {m:+.3f} |")
        out.append("")
    if not r:
        return "\n".join(out + ["_results_predictor_v2.json: not run (python validation/predictor_v2.py --test)._"])
    R = r["rule"]
    out += [f"Test (run once): pass on ≥ {R['min_sets']} of 5 sets: ≥ {R['min_percent']:.0f} % of the pattern, interval "
            "above 0, and above Gate 5's model.", "", "| Test set | Predictor v2: % of ceiling [95 %] | Gate 5 model | Pass |",
            "|---|---|---|---|"]
    for k, v in r["per_set"].items():
        ci = v["ci95"]
        out.append(f"| {k} | {v['percent_of_ceiling']:.1f} [{ci[0]:.1f}, {ci[1]:.1f}] | {v['gate5_percent']:.1f} | "
                   f"{'yes' if v['passes'] else 'no'} |")
    ctl = r["summary"].get("bintu_hct116_34_37_auxin", {}).get("sequence + CTCF + marks", {}).get("all_pairs", {}).get("predictor_v2")
    if ctl:
        out.append(f"| bintu_hct116_34_37_auxin (control) | {ctl['percent_of_ceiling']:.1f} | — | — |")
    out += ["", f"{r['datasets_passing']} of {len(r['per_set'])} test sets pass. Verdict: **{r['verdict']}**. Not run: "
                + "; ".join(r["not_run"]) + " (reasons in validation/predictor_v2.py)."]
    return "\n".join(out)


# --------------------------------------------------------------------------------------------- Gate 3
def gate3() -> str:
    b = _load("benchmark/results.json")
    if not b:
        return "_validation/benchmark/results.json: not run (python -m validation.benchmark.run)._"
    summ = b["summary"]
    methods = b["meta"]["methods"]
    out = [f"Run {b['meta']['utc']} ({b['meta']['seconds'] / 3600:.1f} h). All pairs, mean over splits; brackets: "
           "95 % interval (resampling half B, split 0). Full tables: `validation/benchmark/RESULTS_TABLE.md`.", ""]
    for inp in ("imaging", "hic"):
        datasets = [ds for ds in summ if inp in summ[ds] and "all_pairs" in summ[ds][inp]]
        if not datasets:
            continue
        cols = [m for m in methods if any(m in summ[ds][inp]["all_pairs"] for ds in datasets)]
        out += [f"**Input: {'imaging-derived contacts' if inp == 'imaging' else 'sequencing Hi-C (Rao et al. 2014)'}** — "
                "% of ceiling (raw ρ)", "",
                "| Dataset | " + " | ".join(cols) + " |", "|---|" + "---|" * len(cols)]
        for ds in datasets:
            cells = []
            for m in cols:
                v = summ[ds][inp]["all_pairs"].get(m)
                if not v:
                    cells.append("—")
                    continue
                ci = v.get("ci95", {}).get("percent_of_ceiling")
                cells.append(f"{v['percent_of_ceiling']:.1f}{_ci(ci)} ({v['raw_spearman']:.3f})")
            out.append(f"| {ds} | " + " | ".join(cells) + " |")
        out.append("")
    for unit, why in (b.get("not_finished") or {}).items():
        out.append(f"- Not finished: **{unit}**: {why}")
    if b.get("not_finished"):
        out.append("")
    out += [f"'Where we lose' entries: {len(b['where_we_lose'])} (listed in RESULTS_TABLE.md).", ""]
    return "\n".join(out)


# --------------------------------------------------------------------------------------------- Gate 4
def gate4() -> str:
    out = []
    p = _load("results_perturbation.json")
    if p:
        out += [f"Cohesin depletion — {p['dataset']}; {p['mode']}.", "",
                "| Prediction | Change agreement (Spearman) | Change RMSE (log) | CCC vs auxin (nm) | Raw ρ vs auxin | "
                "% of auxin ceiling |", "|---|---|---|---|---|---|"]
        for k, v in p["summary"].items():
            if isinstance(v, dict):
                out.append(f"| {k.replace('_', ' ')} | {_f(v['change_spearman'])} | {_f(v['change_rmse_log'])} | "
                           f"{_f(v['auxin_lin_ccc_nm'])} | {_f(v['auxin_raw_spearman'])} | "
                           f"{_pct(v['auxin_trend_removed_pct_of_ceiling'])} |")
        out.append("")
    else:
        out.append("_results_perturbation.json: not run._")
    s = _load("results_sv.json")
    if s:
        sp, ci = s["spearman"], s["ci95_block_bootstrap"]
        out += [f"Structural variants — {s['settings']['deleted_source']}; window {s['window']}; "
                f"{s['spanning_pairs']:,} pairs across the deletions (new separation 50 kb – 2 Mb).", "",
                "| Prediction of K562 counts | Spearman [95 % block-bootstrap interval] |", "|---|---|"]
        for k in ("model", "distance_shift", "no_change", "model_minus_distance_shift", "model_trend_removed"):
            out.append(f"| {k.replace('_', ' ')} | {sp[k]:+.3f} [{ci[k][0]:+.3f}, {ci[k][1]:+.3f}] |")
        chk = s["check_deleted_bins"]
        out += ["", f"Check: K562 reads on the published deleted bins = {chk['k562_deleted_over_kept_coverage']:.4f} "
                    f"of the kept bins (GM12878: {chk['gm12878_deleted_over_kept_coverage']:.2f}). "
                    f"Verdict: **{s['verdict']}**.", ""]
        ph = _load("results_sv_posthoc.json")
        if ph:
            out += ["Post hoc (after the verdict; cannot change it): K562 contacts across each junction as a fraction "
                    "of a contiguous chain at the same separation (1 = joined on every copy):", "",
                    "| Junction (kept end – kept start) | s′ = 1 | 2 | 4 | 8 | 16 bins |", "|---|---|---|---|---|---|"]
            for j in ph["junctions"]:
                r = j["contiguous_chain_ratio"]
                out.append(f"| {j['left_end_bp'] / 1e6:.2f} – {j['right_start_bp'] / 1e6:.2f} Mb | " +
                           " | ".join(f"{r[k]:.3f}" for k in ("1", "2", "4", "8", "16")) + " |")
            out.append("")
    else:
        out.append("_results_sv.json: not run._")
    return "\n".join(out)


# --------------------------------------------------------------------------------------------- Gate 5
def gate5() -> str:
    r = _load("results_predictor.json")
    if not r:
        return "_results_predictor.json: not run._"
    summ = r["summary"]
    out = ["| Dataset | Role | Predictor: % of ceiling [95 %] | raw ρ | CCC | Baseline (trend, no data): raw ρ | Pass |",
           "|---|---|---|---|---|---|---|"]
    keys = list(r["settings"]["test"]) + list(r["settings"]["control"])
    for k in keys:
        s = summ.get(k, {}).get("sequence + CTCF", {}).get("all_pairs", {})
        p, g = s.get("predictor"), s.get("genomic_trend_no_data")
        if not p or not g:
            out.append(f"| {k} | — | not scored | | | | |")
            continue
        ci = p.get("ci95", {}).get("percent_of_ceiling")
        gp = r["gate5"].get(k)
        ok = "—" if gp is None else ("yes" if gp["i"] and gp["ii"] else "no")
        role = "control" if k in r["settings"]["control"] else "test"
        out.append(f"| {k} | {role} | {p['percent_of_ceiling']:.1f}{_ci(ci)} | {p['raw_spearman']:.3f} | "
                   f"{p['lin_ccc_nm']:.3f} | {g['raw_spearman']:.3f} | {ok} |")
    out += ["", f"{r['datasets_passing']} of {len(r['settings']['test'])} test datasets pass both criteria "
                f"(needed: {r['settings']['pass_min_datasets']}). Verdict: **{r['verdict']}**."]
    return "\n".join(out)


# --------------------------------------------------------------------------------------------- Pillar 1 / 6
def scale() -> str:
    r = _load("scale_benchmark.json")
    if not r:
        return "_validation/scale_benchmark.json: not run (python validation/scale_benchmark.py)._"
    out = [_machine(r["meta"]["machine"]) + " " + r["meta"]["note"], "",
           "| Beads | Model | Fit (s) | Sampling (s) | Total (s) | Peak memory (MB) | Device | Rank | CPU load before |",
           "|---|---|---|---|---|---|---|---|---|"]
    for row in r["rows"]:
        if "skipped" in row:
            out.append(f"| {row['n_beads']:,} | {row['method']} | — | — | — | — | — | {row['skipped']} | |")
        elif "error" in row:
            out.append(f"| {row['n_beads']:,} | {row['method']} | failed | | | | | {row['error'][-80:]} | |")
        else:
            out.append(f"| {row['n_beads']:,} | {row['method']} | {row['seconds_fit']:.1f} | {row['seconds_sampling']:.1f} | "
                       f"{row['seconds_total']:.1f} | {row['peak_memory_mb']:,.0f} | {row['device']} | {row['rank']} | "
                       f"{_cpu_load(row)} |")
    return "\n".join(out)


def per_chromosome() -> str:
    r = _load("chromosome_runtime.json")
    if not r:
        return "_validation/chromosome_runtime.json: not run (python validation/chromosome_runtime.py)._"
    out = [_machine(r["meta"]["machine"]) + " " + r["meta"]["note"], "",
           "| Assembly | Chromosome | Resolution | Bins | Assembled beads | Model | Fit (s) | Total (s) | Peak memory (MB) | CPU load before |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for row in r["rows"]:
        if "error" in row:
            out.append(f"| {row['assembly']} | {row['chrom']} | | {row.get('bins', '—')} | failed: {row['error'][-60:]} | | | | | |")
        else:
            out.append(f"| {row['assembly']} | {row['chrom']} | {row['resolution_bp'] // 1000} kb | {row['bins']:,} | "
                       f"{row['beads']:,} | {row['model']} | {row['seconds_fit']:.1f} | {row['seconds_total']:.1f} | "
                       f"{row['peak_memory_mb']:,.0f} | {_cpu_load(row)} |")
    ok = [row for row in r["rows"] if "error" not in row]
    if ok:
        out += ["", f"{len(ok)} chromosomes; total {sum(x['seconds_total'] for x in ok) / 60:.1f} min of fitting; slowest "
                    f"{max(x['seconds_total'] for x in ok):.0f} s; largest peak memory {max(x['peak_memory_mb'] for x in ok):,.0f} MB."]
    planned = [tuple(x) for x in r["meta"].get("planned", [])]
    have = {(x["assembly"], x["chrom"]) for x in ok}
    missing = [f"{a} {c}" for a, c in planned if (a, c) not in have]
    if missing or not r["meta"].get("status", "").startswith("complete"):
        out += ["", f"**Status: {r['meta'].get('status', 'partial')}.** Not measured yet ({len(missing)} of "
                    f"{len(planned)}): {', '.join(missing) or '—'}."]
    for d in r["meta"].get("discarded_rows", []):
        out.append(f"- Discarded: {d['assembly']} {d['chrom']} ({d['seconds_total']:.0f} s): {d['reason']}.")
    return "\n".join(out)


def _cpu_load(row: dict) -> str:
    v = row.get("system_cpu_percent_before")
    return "—" if v is None else f"{v:.0f} %"


def _machine(m: dict) -> str:
    gpu = f", GPU {m['gpu']}" if m.get("gpu") else ", no GPU (CPU only)"
    ram = f", {m['ram_gb']:g} GB RAM" if m.get("ram_gb") else ""
    return (f"Machine: {m['platform']}, {m.get('processor') or 'CPU'}, {m['cpus']} logical CPUs{ram}{gpu}; "
            f"torch {m['torch']}" + (f", Python {m['python']}" if m.get("python") else "") + ".")


def readme_accuracy() -> str:
    """The README's accuracy table: one line per test, every number from a result file."""
    rows = ["| Test (held-out, real data) | Measured | Verdict |", "|---|---|---|"]
    v33 = _load("results.json")
    if v33:
        o = v33["summary"]["overall"]
        rows.append(f"| v3.3 windows, Bintu et al. 2018 tracing (3 test sets, Σ model / Σ ceiling) | population "
                    f"{o['ensemble_v3_3']['percent_of_ceiling']:.1f} % of the reproducible structure; v3.2 single "
                    f"structure {o['single_structure_v3_2']['percent_of_ceiling']:.1f} % | baseline for v4 |")
    within, cross, hic, ccc = [], [], [], []
    for fname in ("results_gate1.json", "results_gate1_rep.json"):
        r = _load(fname)
        for d in (r or {}).get("datasets", []):
            s = d["summary"]
            w = s["imaging_input"]["within_tiles"]
            within.append((d["dataset"], w["whole_v4"]["percent_of_ceiling"], w["windowed_v3_3"]["percent_of_ceiling"]))
            cross.append(s["imaging_input"]["cross_tiles"]["whole_v4"]["percent_of_ceiling"])
            h = s["hic_input"]["all_pairs"]
            for k, v in h.items():
                if k.startswith("whole_v4"):
                    hic.append(v["percent_of_ceiling"])
                    ccc.append(v["lin_ccc_nm"])
    if within:
        gap = max(v - w for _, w, v in within)
        verdict1 = ("matches or beats" if gap <= 0 else f"matches within {gap:.1f} points, slightly below" if gap <= 1.0
                    else f"below the windowed model by up to {gap:.1f} points")
        rows.append("| Gate 1: whole chromosome (v4) vs windows (v3.3), Su et al. 2020 chr21 + replicate, same pairs | "
                    + "; ".join(f"{w:.1f} vs {v:.1f} %" for _, w, v in within)
                    + f"; cross-window pairs (v4 only) {min(cross):.1f}–{max(cross):.1f} % | {verdict1} |")
    if hic:
        rows.append(f"| Gate 1b: sequencing Hi-C (Rao et al. 2014) → imaged distances, all pairs | ranks {min(hic):.1f}–"
                    f"{max(hic):.1f} % of the ceiling; absolute size CCC {min(ccc):.2f}–{max(ccc):.2f} | "
                    f"{'ranks transfer' if min(hic) >= 50 else 'ranks transfer poorly'}, "
                    f"{'nanometres do not' if max(ccc) < 0.5 else 'nanometres partly'} |")
    c = _load("results_calibration.json")
    if c:
        raw = [100 * r["coverage"][2] for r in c["rows"]]
        rec = [100 * r["coverage_recalibrated"][2] for r in c["rows"]]
        rel = [r["reliability_vs_error_spearman"] for r in c["rows"]]
        rows.append(f"| Gate 2: do stated 90 % intervals hold 90 % of real single-cell distances? (6 test sets) | raw "
                    f"{min(raw):.0f}–{max(raw):.0f} %, recalibrated on practice data {min(rec):.0f}–{max(rec):.0f} %; "
                    f"per-bead reliability vs error ρ {min(rel):+.2f} to {max(rel):+.2f} | raw intervals too narrow; "
                    f"recalibrated {'within 7 points of nominal' if min(rec) >= 83 else 'still too narrow'}; "
                    f"{'no usable per-bead reliability' if min(rel) < 0.2 else 'per-bead reliability usable'} |")
    h = _load("results_calibration_hic.json")
    if h:
        raw = [100 * r["coverage"][2] for r in h["rows"]]
        rec = [100 * r["coverage_recalibrated"][2] for r in h["rows"]]
        rows.append(f"| Gate 2b: the same, with sequencing Hi-C input ({len(h['rows'])} test sets) | raw {min(raw):.0f}–"
                    f"{max(raw):.0f} %, recalibrated on practice Hi-C {min(rec):.0f}–{max(rec):.0f} % | "
                    f"{'pass' if h['verdict'] == 'pass' else 'fail: intervals far too narrow for Hi-C input; the app says so'} |")
    rr = _load("results_reliability.json")
    if rr:
        v = rr["verdict"]
        rows.append(f"| Gate 2c: does a per-pair score from the input say which distances are wrong? | test sets reaching "
                    f"the pre-registered bar (ρ ≥ {rr['rule']['min_spearman']:.2f}): "
                    + "; ".join(f"{k.replace('_', ' ').replace('hic', 'Hi-C')} {x['passing']} of {x['datasets']}"
                                for k, x in v.items())
                    + f" | {'usable' if any(x['verdict'] == 'pass' for x in v.values()) else 'no usable reliability score'} |")
    b = _load("benchmark/results.json")
    if b:
        best, pastis, n_ds = [], [], 0
        for ds, inputs in b["summary"].items():
            ap = inputs.get("imaging", {}).get("all_pairs", {})
            ours = [v["percent_of_ceiling"] for k, v in ap.items() if k in ("v3_3_windowed", "v4_whole")]
            pas = [v["percent_of_ceiling"] for k, v in ap.items() if k.startswith("pastis")]
            if ours and pas:
                n_ds += 1
                best.append(max(ours))
                pastis.append(max(pas))
        if n_ds:
            rows.append(f"| Gate 3: benchmark vs PASTIS 0.4.0 and baselines, imaging-derived input ({n_ds} test units) | "
                        f"ChronoCell {min(best):.1f}–{max(best):.1f} %, best PASTIS {min(pastis):.1f}–{max(pastis):.1f} % "
                        f"of the ceiling; {len(b['where_we_lose'])} 'where we lose' entries | see RESULTS.md |")
    p = _load("results_perturbation.json")
    if p:
        f, t = p["summary"]["full_model"], p["summary"]["trend_only"]
        rows.append(f"| Gate 4: cohesin loss (RAD21 degron, Bintu et al. 2018), held-out region | change agreement "
                    f"{f['change_spearman']:.3f} vs {t['change_spearman']:.3f} for a trend-only shift | "
                    f"{'pass, one region' if f['change_spearman'] > t['change_spearman'] and f['auxin_lin_ccc_nm'] > p['summary']['no_change']['auxin_lin_ccc_nm'] else 'fail'} |")
    s = _load("results_sv.json")
    if s:
        sp = s["spearman"]
        rows.append(f"| Gate 4b: structural variants (K562 chr9 deletions vs GM12878, Rao 2014 Hi-C) | model "
                    f"{sp['model']:.3f}, distance shift {sp['distance_shift']:.3f}, no change {sp['no_change']:.3f} | "
                    f"{s['verdict']} |")
    g5 = _load("results_predictor.json")
    if g5:
        rows.append(f"| Gate 5: distances from sequence + CTCF alone (no contact data) | {g5['datasets_passing']} of "
                    f"{len(g5['settings']['test'])} test sets pass the pre-registered rule | {g5['verdict']} |")
    rows += _phase_ab_rows()
    return "\n".join(rows)


def _phase_ab_rows() -> list[str]:
    """README rows of the Phase A / B gates, each from its result file (a gate not run says so)."""
    rows = []
    r = _load("results_hic_size.json")
    if r:
        m = r["sets"]["main"]
        cc = [v["calibrated"]["lin_ccc_nm"] for v in m.values()]
        rows.append(f"| Gate 1c: calibrated sizes from Hi-C ({len(m)} main sets) | CCC {min(cc):.2f}–{max(cc):.2f} after "
                    f"calibration (needed ≥ 0.8) | {r['verdict']} |")
    r = _load("results_intervals_v2.json")
    if r:
        n = {i: sum(1 for k, v in r["sets"].items() if k.endswith(i) and v["within_rule"]) for i in ("imaging", "hic")}
        t = {i: sum(1 for k in r["sets"] if k.endswith(i)) for i in ("imaging", "hic")}
        rows.append(f"| Gate 2d: conformal distance ranges hold 85–95 % / 40–60 % in every band | imaging input "
                    f"{n['imaging']} of {t['imaging']} sets, Hi-C input {n['hic']} of {t['hic']} | "
                    + ", ".join(f"{k} {v}" for k, v in r["verdict"].items()) + " |")
    r = _load("results_reliability_v2.json")
    if r:
        rows.append("| Gate 2e: per-pair reliability on untouched genome-scale sets (ρ ≥ 0.30) | "
                    + "; ".join(f"{k} {v['rho']:+.2f}" for k, v in r["sets"].items()) + " | "
                    + ", ".join(f"{k} {v}" for k, v in r["verdict"].items()) + " |")
    rows.append("| Gate 3b: a learned correction on the population model | practice chose no correction | not run |")
    r = _load("results_cohesin_hic.json")
    if r:
        rows.append(f"| Gate 4c: cohesin loss vs RAD21-degron Hi-C on 6 held-out regions | {r['regions_beating_trend']} of "
                    f"{len(r['rule']['regions'])} regions beat trend only | {r['verdict']} |")
    rows.append("| Gate 4d: SV effects on new events with Hi-C before and after | two usable events found, three needed | "
                "blocked; variant engine stays a mechanism simulator |")
    r = _load("results_predictor_v2.json")
    if r:
        rows.append(f"| Gate 5b: prediction with cohesin peaks | {sum(v['passes'] for v in r['per_set'].values())} of "
                    f"{len(r['per_set'])} test sets | {r['verdict']} |")
    r = _load("results_predictor_mouse.json")
    if r:
        rows.append("| Gate 5m: the human predictor on mouse ES-cell tracing (4DN) | "
                    + "; ".join(f"{v['percent_of_ceiling']:.1f} % of the ceiling" for v in r["gate5m"].values())
                    + f" | {r['verdict']} (modest) |")
    r = _load("results_gate6.json")
    if r:
        rows.append("| Gate 6: loop calls vs ENCODE HiCCUPS calls, held-out cell lines | "
                    + "; ".join(f"{c}: F1 {v['chronocell']['f1']:.2f} (chromosight {v['chromosight'].get('f1', float('nan')):.2f}, "
                                f"Mustache {v['mustache'].get('f1', float('nan')):.2f})" for c, v in r["sets"].items())
                    + f" | {r['verdict']} |")
    r = _load("results_gate6b.json")
    if r:
        rows.append("| Gate 6b: loop calls, settings re-chosen on three cell lines, two new cell lines | "
                    + "; ".join(f"{c}: F1 {v['chronocell']['f1']:.2f} (chromosight {v['chromosight'].get('f1', float('nan')):.2f}, "
                                f"Mustache {v['mustache'].get('f1', float('nan')):.2f})" for c, v in r["sets"].items())
                    + f" | {r['verdict']} |")
    r = _load("results_gate7.json")
    if r:
        rows.append(f"| Gate 7: false discoveries of the differential analysis (real replicates + planted changes) | "
                    f"mean FDP {r['summary']['mean_fdp']:.3f} at nominal {r['rule']['fdr']}; recall ×2 "
                    f"{r['summary']['recall']['2.0']:.2f}, ×4 {r['summary']['recall']['4.0']:.2f} | {r['verdict']} |")
    return rows


def gate4c() -> str:
    p = _load("results_cohesin_hic_practice.json")
    r = _load("results_cohesin_hic.json")
    out = []
    if p:
        out += ["Practice (the Gate 4 practice region; nothing fitted):", "",
                "| Pair | Region | Model ρ | Trend-only ρ | Difference [95 %] |", "|---|---|---|---|---|"]
        for x in p["rows"]:
            out.append(f"| {x['pair']} | {x['region']} | {x['model']:+.3f} | {x['trend_only']:+.3f} | {x['difference']:+.3f} "
                       f"[{x['difference_ci95'][0]:+.3f}, {x['difference_ci95'][1]:+.3f}] |")
        out.append("")
    if not r:
        return "\n".join(out + ["_results_cohesin_hic.json: not run (python validation/cohesin_hic.py --test)._"])
    R = r["rule"]
    out += [f"Test (run once): pass if ≥ {R['min_regions']} of {len(R['regions'])} held-out regions beat trend only on the "
            f"main pair ({R['main_pair']}), difference > 0 with its 95 % interval above 0.", "",
            "| Pair | Region | Pairs | Model ρ | Trend-only ρ | Difference [95 %] | Beats trend only |",
            "|---|---|---|---|---|---|---|"]
    for x in r["rows"]:
        out.append(f"| {x['pair']} ({x['tier']}) | {x['region']} | {x['pairs']:,} | {x['model']:+.3f} | {x['trend_only']:+.3f} | "
                   f"{x['difference']:+.3f} [{x['difference_ci95'][0]:+.3f}, {x['difference_ci95'][1]:+.3f}] | "
                   f"{'yes' if x['beats_trend'] else 'no'} |")
    out += ["", f"{r['regions_beating_trend']} of {len(R['regions'])} main-pair regions beat trend only. "
                f"Verdict: **{r['verdict']}**."]
    return "\n".join(out)


def gate5m() -> str:
    r = _load("results_predictor_mouse.json")
    if not r:
        return "_results_predictor_mouse.json: not run (python validation/predictor_mouse.py --test)._"
    s = r["summary"]
    out = ["Test (run once; the human model unchanged): pass if both test loci have % of ceiling > 0 with its 95 % "
           "interval above 0 and raw Spearman above the training-trend baseline.", "",
           "| Set | Role | Traces | Loci | Predictor: % of ceiling (3 splits) | Trend baseline | Raw ρ: predictor / trend |",
           "|---|---|---|---|---|---|---|"]
    for k, v in r["sets"].items():
        a = s.get(k, {}).get("sequence + CTCF", {}).get("all_pairs", {})
        p, g = a.get("predictor", {}), a.get("genomic_trend_no_data", {})
        ci = r["gate5m"].get(k, {}).get("ci95")
        ci_txt = f" [{ci[0]:.1f}, {ci[1]:.1f}]" if ci else ""
        out.append(f"| {k} ({v['description']}) | {v['role']} | {v['copies']:,} | {v['loci']} | "
                   f"{_f(p.get('percent_of_ceiling'), 1)}{ci_txt} | {_f(g.get('percent_of_ceiling'), 1)} | "
                   f"{_f(p.get('raw_spearman'))} / {_f(g.get('raw_spearman'))} |")
    out += ["", f"Verdict: **{r['verdict']}** (95 % interval from split 0)."]
    return "\n".join(out)


def gate6() -> str:
    p = _load("results_gate6_practice.json")
    r = _load("results_gate6.json")
    out = []
    if p:
        out += ["Practice (GM12878, three 10 Mb windows; ChronoCell settings):", "", "| Setting | Precision | Recall | F1 |",
                "|---|---|---|---|"]
        for k, v in p["table"].items():
            out.append(f"| {k}{' (chosen)' if k == p['chosen'] else ''} | {v['precision']:.3f} | {v['recall']:.3f} | "
                       f"{v['f1']:.3f} |")
        out.append("")
    if not r:
        return "\n".join(out + ["_results_gate6.json: not run (python validation/loops_gate6.py --test)._"])
    out += ["Test (run once): pass if ChronoCell's F1 is at least the best of chromosight and Mustache on both held-out "
            "cell lines (reference: ENCODE HiCCUPS calls on the same maps).", "",
            "| Cell line | Method | Calls | Reference loops | Matched | Precision | Recall | F1 |",
            "|---|---|---|---|---|---|---|---|"]
    for cell, v in r["sets"].items():
        for tool in ("chronocell", "chromosight", "mustache"):
            m = v.get(tool, {})
            if "f1" not in m:
                out.append(f"| {cell} | {tool} | — | — | — | — | — | {m.get('status', 'not run')} |")
                continue
            out.append(f"| {cell} | {tool} | {m['calls']} | {m['reference']} | {m['matched']} | {m['precision']:.3f} | "
                       f"{m['recall']:.3f} | {m['f1']:.3f} |")
    out += ["", f"Verdict: **{r['verdict']}**."]
    return "\n".join(out)


def gate8_marks() -> str:
    r = _load("results_gate8_marks_practice.json")
    if not r:
        return "_results_gate8_marks_practice.json: not run (python validation/drug_gate8.py --marks)._"
    best = sorted(r["rows"], key=lambda x: x["perm_p"])[:6]
    out = [f"Round-2 practice check (Gate 8's data, all seen): other IMR-90 marks as the drugs' targets, the simplest "
           f"model (a pair changes with the two loci's mark coverage). {r['comparisons']} drug × allele × mark comparisons; "
           f"{r['below_0.05']} had a permutation p ≤ 0.05 ({r['expected_below_0.05_by_chance']:.1f} expected by chance; "
           "the comparisons overlap: the marks are correlated and every drug shares the untreated traces). Smallest p values:",
           "", "| Drug | Allele | Mark | ρ | Permutation p |", "|---|---|---|---|---|"]
    for x in best:
        out.append(f"| {x['drug']} | {x['allele']} | {x['mark']} | {x['rho']:+.2f} | {x['perm_p']:.3f} |")
    out += ["", "The directions differ by drug and would have to be fitted on these data; no fresh tracing after these drugs "
                "exists to test such a model, so Gate 8 is not retested."]
    return "\n".join(out)


def gate6b() -> str:
    p = _load("results_gate6b_practice.json")
    r = _load("results_gate6b.json")
    out = []
    if p:
        b = p["choice_scores"]
        out += [f"Practice (nine windows: GM12878, K562, IMR-90; {len(p['grid'])} settings). Chosen: {p['choice']}.", "",
                "| Cell line | ChronoCell calls | Reference | Matched | Precision | Recall | F1 | chromosight F1 | Mustache F1 |",
                "|---|---|---|---|---|---|---|---|---|"]
        for c, m in b["detail"].items():
            t = p["tools"][c]
            out.append(f"| {c} | {m['calls']} | {m['reference']} | {m['matched']} | {m['precision']:.3f} | {m['recall']:.3f} | "
                       f"{m['f1']:.3f} | {_f(t['chromosight'].get('f1'))} | {_f(t['mustache'].get('f1'))} |")
        out.append("")
    if not r:
        return "\n".join(out + ["_results_gate6b.json: not run (python validation/loops_gate6b.py --test)._"])
    out += ["Test (run once; HMEC and HAP-1, three new windows each):", "",
            "| Cell line | Method | Calls | Reference loops | Matched | Precision | Recall | F1 |", "|---|---|---|---|---|---|---|---|"]
    for cell, v in r["sets"].items():
        for tool in ("chronocell", "chromosight", "mustache"):
            m = v.get(tool, {})
            if "f1" not in m:
                out.append(f"| {cell} | {tool} | — | — | — | — | — | {m.get('status', 'not run')} |")
                continue
            out.append(f"| {cell} | {tool} | {m['calls']} | {m['reference']} | {m['matched']} | {m['precision']:.3f} | "
                       f"{m['recall']:.3f} | {m['f1']:.3f} |")
    out += ["", f"Verdict: **{r['verdict']}**."]
    return "\n".join(out)


def gate7() -> str:
    p = _load("results_gate7_practice.json")
    r = _load("results_gate7.json")
    out = []
    if p:
        out += ["Practice (chr21:28-30 Mb):", "", "| Setting | Mean FDP | Recall ×2 | Recall ×4 | No-change discoveries |",
                "|---|---|---|---|---|"]
        for k, v in p["table"].items():
            sm = v["summary"]
            out.append(f"| {k} | {sm['mean_fdp']:.3f} | {sm['recall']['2.0']:.2f} | {sm['recall']['4.0']:.2f} | "
                       f"{sum(sm['null_discoveries'].values())} |")
        out.append("")
    if not r:
        return "\n".join(out + ["_results_gate7.json: not run (python validation/diff_gate7.py --test)._"])
    sm, R = r["summary"], r["rule"]
    out += [f"Test (run once): pass if the mean false-discovery proportion over {sm['runs']} spike-in runs is ≤ {R['fdr']}.",
            "", "| Region | No-change discoveries | Spike-in FDP (×2 / ×4) | Recall ×2 / ×4 |", "|---|---|---|---|"]

    def fdp(xs):
        return sum(x["false"] / max(x["discoveries"], 1) for x in xs) / max(len(xs), 1)

    def rec(xs):
        return sum(x["found"] / max(x["planted"], 1) for x in xs) / max(len(xs), 1)
    for reg in dict.fromkeys(x["region"] for x in r["rows"]):
        rows = [x for x in r["rows"] if x["region"] == reg]
        null = next(x for x in rows if x["fold"] == 1.0)
        f2 = [x for x in rows if x["fold"] == 2.0]
        f4 = [x for x in rows if x["fold"] == 4.0]
        out.append(f"| {reg} | {null['discoveries']} | {fdp(f2):.3f} / {fdp(f4):.3f} | {rec(f2):.2f} / {rec(f4):.2f} |")
    ci = sm["mean_fdp_ci95"]
    out += ["", f"Mean FDP {sm['mean_fdp']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}] at nominal {R['fdr']}; recall ×2 "
                f"{sm['recall']['2.0']:.2f}, ×4 {sm['recall']['4.0']:.2f}. Not run: "
            + "; ".join(f"{k} ({v})" for k, v in r["not_run"].items()) + f". Verdict: **{r['verdict']}**."]
    return "\n".join(out)


def b9tools() -> str:
    r = _load("results_tools_b9.json")
    if not r:
        return "_results_tools_b9.json: not written (python validation/tools_b9.py)._"
    out = ["| Task | Tool | Status | Where the comparison is |", "|---|---|---|---|"]
    for x in r["tools"]:
        out.append(f"| {x['task']} | {x['tool']} | {x['status']} | {x.get('result', '—')} |")
    return "\n".join(out)


def summary_ab() -> str:
    """Summary rows of the Phase A / B gates for RESULTS.md (the same rows as the README table)."""
    return "\n".join(["| Test (held-out, real data) | Measured | Verdict |", "|---|---|---|"] + _phase_ab_rows())


# --------------------------------------------------------------------------------------------- Gate Q (quantum lab)
def _q_rows() -> list[str]:
    """Summary rows of Gate Q (README and RESULTS), from results_gateq.json."""
    r = _load("results_gateq.json")
    if not r or "overall" not in r:
        return ["| Gate Q: quantum lab (simulated quantum) | not run (python validation/quantum_gateq.py --test all) | — |"]
    v = lambda k: "pass" if r["overall"][k] else "fail"
    q1, q2, q3, q4 = r["q1"], r["q2"], r["q3"], r["q4"]
    f1 = "; ".join(f"{c}: QAOA {rows['QAOA (simulated quantum)']['f1']:.2f} vs insulation "
                   f"{rows['insulation (classical)']['f1']:.2f}, TopDom-like {rows['topdom (classical)']['f1']:.2f}"
                   for c, rows in q2["rows"].items())
    return [f"| Gate Q1: QAOA (simulator) finds the optimum of the domain QUBO, {q1['windows']} held-out windows | "
            f"{100 * q1['pooled_hit_rate']:.0f} % of windows (needed ≥ {100 * r['rule']['q1_min_hit_rate']:.0f} %) | {v('q1')} |",
            f"| Gate Q2: quantum domain calls vs classical callers (F1 vs ENCODE Arrowhead) | {f1} | {v('q2')} |",
            f"| Gate Q3: VQE within chemical accuracy (H2, HeH+) | worst error {q3['max_abs_error_mha']:.2g} mHa "
            f"(needed ≤ 1.6) | {v('q3')} |",
            f"| Gate Q4: quantum-kernel gene classifier vs RBF-SVM, GM12878 → IMR-90 | AUC {q4['auc']['qsvm']:.3f} vs "
            f"{q4['auc']['rbf']:.3f} | {v('q4')} |"]


def summary_q() -> str:
    return "\n".join(["| Test (held-out, real data; simulated quantum) | Measured | Verdict |", "|---|---|---|"] + _q_rows())


def gateq_practice() -> str:
    t = _load("results_gateq_practice_tad.json")
    a = _load("results_gateq_practice_qaoa.json")
    k = _load("results_gateq_practice_qsvm.json")
    c = _load("results_gateq_practice_chem.json")
    out = []
    if t:
        ch, best = t["choice"], t["choice_scores"]
        out += [f"Domain QUBO, practice GM12878 ({best['windows']} windows of 20 bins; {len(t['grid'])} settings, exact "
                f"optimum of each): chosen {ch['res'] // 1000} kb bins, minimum domain {ch['min_size']} bins, "
                f"gamma {ch['gamma']}, boundary cost {ch['boundary_cost']}, {ch['weight']} weights: F1 {best['f1']:.3f} "
                f"(precision {best['precision']:.3f}, recall {best['recall']:.3f}). Classical callers on the same windows, "
                "best of their grids:", "", "| Caller | Setting | Precision | Recall | F1 |", "|---|---|---|---|---|",
                f"| Domain QUBO (exact optimum) | chosen | {best['precision']:.3f} | {best['recall']:.3f} | {best['f1']:.3f} |"]
        for m, v in t["classical"].items():
            b = v["best"]
            out.append(f"| {m} | {b['param']} | {b['precision']:.3f} | {b['recall']:.3f} | {b['f1']:.3f} |")
        out.append("")
    if a:
        out += ["QAOA settings, practice (same windows; 4,096 shots; noise column: approximate hardware-noise model):", "",
                "| Depth p | Objective | Optimum found | Mean P(optimum) | Uniform P(optimum) | With noise | Simulated "
                "annealing | SQA | Mean QAOA s | Mean CX |", "|---|---|---|---|---|---|---|---|---|---|"]
        for r in a["grid"]:
            out.append(f"| {r['p']} | {r['objective']} | {100 * r['hit_rate']:.0f} % | {r['mean_p_optimal']:.4f} | "
                       f"{r['mean_uniform_p_optimal']:.2g} | {100 * r['noisy_hit_rate']:.0f} % | {100 * r['sa_hit_rate']:.0f} % | "
                       f"{100 * r['sqa_hit_rate']:.0f} % | {r['mean_qaoa_seconds']:.1f} | {r['mean_cx']:.0f} |")
        out += ["", f"Chosen: p = {a['choice']['p']}, {a['choice']['objective']}.", ""]
    if k:
        b = k["best_cv_auc"]
        out += [f"Gene classifier, practice: {k['genes']} GM12878 genes ({100 * k['expressed_fraction']:.0f} % expressed), "
                f"5-fold cross-validated AUC: quantum kernel {b['qsvm']:.3f} ({k['choice']['qsvm']}), RBF-SVM "
                f"{b['rbf']:.3f} ({k['choice']['rbf']}), logistic regression {b['logistic']:.3f}.", ""]
    if c:
        out += ["Chemistry anchors (Szabo & Ostlund): " + "; ".join(f"{kk}: {vv:.5f}" for kk, vv in c["anchors"].items()) + "."]
    return "\n".join(out) if out else "_Gate Q practice: not run._"


def gateq() -> str:
    r = _load("results_gateq.json")
    if not r or "overall" not in r:
        return "_results_gateq.json: not run (python validation/quantum_gateq.py --test all)._"
    R = r["rule"]
    out = ["Q1 (run once, held-out K562 and IMR-90 windows; 19 qubits each):", "",
           "| Cell line | Windows | QAOA found the optimum | Mean P(optimum) | Random guessing found it | With noise | "
           "Simulated annealing | SQA | Mean QAOA s | Mean exact DP s |", "|---|---|---|---|---|---|---|---|---|---|"]
    for cell in ("k562", "imr90"):
        x = r["q1"][cell]
        out.append(f"| {cell} | {x['windows']} | {100 * x['hit_rate']:.0f} % | {x['mean_p_optimal']:.4f} | "
                   f"{100 * x['random_hit_rate']:.0f} % | {100 * x['noisy_hit_rate']:.0f} % | {100 * x['sa_hit_rate']:.0f} % | "
                   f"{100 * x['sqa_hit_rate']:.0f} % | {x['mean_qaoa_seconds']:.1f} | {x['mean_dp_seconds']:.4f} |")
    out += ["", f"Pooled: {100 * r['q1']['pooled_hit_rate']:.1f} % of {r['q1']['windows']} windows (needed ≥ "
                f"{100 * R['q1_min_hit_rate']:.0f} %). Q1: **{'pass' if r['q1']['pass'] else 'fail'}**.", "",
            f"Q2 (reference: ENCODE Arrowhead domains on the same maps; a call matches within one bin; pass if QAOA's F1 ≥ "
            f"the better classical caller's − {R['q2_margin']} on both cell lines):", "",
            "| Cell line | Method | Precision | Recall | F1 |", "|---|---|---|---|---|"]
    for cell, rows in r["q2"]["rows"].items():
        for m, v in rows.items():
            out.append(f"| {cell} | {m} | {_f(v['precision'])} | {_f(v['recall'])} | {_f(v['f1'])} |")
    out += ["", f"Q2: **{'pass' if r['q2']['pass'] else 'fail'}** ("
            + ", ".join(f"{c} {'pass' if ok else 'fail'}" for c, ok in r["q2"]["pass_by_cell"].items()) + ").", "",
            "Q3 (UCCSD-VQE vs exact diagonalisation; noise column: hardware-efficient ansatz under the noise model):", "",
            "| Molecule | Bond length (Å) | Ansatz | Hartree–Fock | FCI | VQE | Error (mHa) | With noise |",
            "|---|---|---|---|---|---|---|---|"]
    for x in r["q3"]["rows"]:
        out.append(f"| {x['molecule']} | {x['R_angstrom']:.3f} | {x['ansatz']} | {x['e_hf']:.5f} | {x['e_fci']:.5f} | "
                   f"{x['e_vqe']:.5f} | {x['error_mha']:+.2g} | {_f(x['e_vqe_noisy'], 4)} |")
    q4 = r["q4"]
    ci = q4["ci95"]
    out += ["", f"Q3: **{'pass' if r['q3']['pass'] else 'fail'}** (worst UCCSD error {r['q3']['max_abs_error_mha']:.2g} mHa).", "",
            f"Q4 (trained on {q4['train_genes']} GM12878 genes, tested on {q4['test_genes']} IMR-90 genes, "
            f"{100 * q4['test_expressed_fraction']:.0f} % expressed; 95 % intervals from 1,000 resamples of genes):", "",
            "| Classifier | Test AUC |", "|---|---|",
            f"| Quantum-kernel SVM (simulated, exact kernel) | {q4['auc']['qsvm']:.3f}{_ci(ci['qsvm'], 3)} |",
            f"| Quantum-kernel SVM, kernel from {R['qsvm_shots']} shots per entry | {q4['qsvm_auc_with_shots']:.3f} |",
            f"| RBF-kernel SVM (classical) | {q4['auc']['rbf']:.3f}{_ci(ci['rbf'], 3)} |",
            f"| Logistic regression (classical) | {q4['auc']['logistic']:.3f}{_ci(ci['logistic'], 3)} |",
            "", f"Quantum minus RBF: {q4['auc']['qsvm'] - q4['auc']['rbf']:+.3f}{_ci(ci['qsvm_minus_rbf'], 3)}. Q4: "
                f"**{'pass' if q4['pass'] else 'fail'}** (needed: within {R['q4_margin']} of the RBF-SVM and the lower "
                "95 % bound above 0.5)."]
    return "\n".join(out)


def quantum_crosscheck() -> str:
    r = _load("results_quantum_crosscheck.json")
    if not r:
        return "_results_quantum_crosscheck.json: not written (python validation/quantum_crosscheck.py)._"
    out = [f"Qiskit {r['qiskit']} statevectors of the exported OpenQASM circuits vs ChronoCell's simulator:", "",
           "| Circuit | Qubits | Gates | State fidelity |", "|---|---|---|---|"]
    for x in r["circuits"]:
        out.append(f"| {x['circuit']} | {x['qubits']} | {x['gates']} | {x['fidelity']:.12f} |")
    out += ["", f"Verdict: {r['verdict']}."]
    return "\n".join(out)


# --------------------------------------------------------------------------------------------- Gate 8 and Q5-Q7 (drug tabs)
def gate8() -> str:
    pr = _load("results_gate8_practice.json")
    r = _load("results_gate8.json")
    out = []
    if pr:
        out += ["Practice (no treated test data): untreated traces split in two halves (no drug, so any agreement is noise) "
                "and alpha-amanitin (transcription-inhibitor class):", "",
                "| Comparison | Allele | Drug-lab class | Spearman rho | 95 % interval | Shuffled-target mean | Permutation p |",
                "|---|---|---|---|---|---|---|"]
        for allele, v in pr["split_halves"].items():
            for k, x in v.items():
                out.append(f"| untreated half vs half | {allele} | {k} | {x['rho']:+.3f} | {_ci(x['ci95'], 3, False)} | "
                           f"{x['perm_mean']:+.3f} | {x['perm_p']:.3f} |")
        for allele, x in pr["practice_drug"].items():
            out.append(f"| alpha-amanitin vs untreated | {allele} | {x['drug_class']} | {x['rho']:+.3f} | {_ci(x['ci95'], 3, False)} | "
                       f"{x['perm_mean']:+.3f} | {x['perm_p']:.3f} |")
        out.append("")
    if not r:
        return "\n".join(out + ["_results_gate8.json: not run (python validation/drug_gate8.py --test)._"])
    R = r["rule"]
    out += [f"Test (run once; pass per drug on {R['allele']}: rho >= {R['min_rho']}, 95 % interval above 0, permutation p <= "
            f"{R['max_perm_p']}; gate: at least {R['min_drugs']} of {len(r['drugs'])} drugs):", "",
            "| Drug (class) | Allele | Spearman rho | 95 % interval | Shuffled-target mean | Permutation p | Measured global log change | Traces (untreated / treated) | Pass |",
            "|---|---|---|---|---|---|---|---|---|"]
    for cond, v in r["drugs"].items():
        for allele, x in v.items():
            out.append(f"| {cond} | {allele} | {x['rho']:+.3f} | {_ci(x['ci95'], 3, False)} | {x['perm_mean']:+.3f} | "
                       f"{x['perm_p']:.3f} | {x['measured_global_log_change']:+.3f} | {x['n_untreated']} / {x['n_treated']} | "
                       f"{('yes' if x.get('passes') else 'no') if allele == R['allele'] else '—'} |")
    out += ["", f"{r['passing_drugs']} of {len(r['drugs'])} drugs pass. Verdict: **{r['verdict']}**."]
    return "\n".join(out)


def qdrug_practice() -> str:
    q5, q6, q7 = (_load(f"results_qdrug_practice_{k}.json") for k in ("q5", "q6", "q7"))
    out = []
    if q5:
        a = q5["featuriser_vs_rdkit"]
        out += [f"Q5 practice: TDC hERG, {q5['compounds']} compounds read ({q5['unreadable']} unreadable), "
                f"{100 * q5['blockers']:.0f} % blockers. SMILES descriptors against RDKit (share exact or within tolerance; "
                "Pearson r): " + "; ".join(f"{k} {100 * v['exact_or_within_tol']:.0f} % (r {v['pearson']:.3f})" for k, v in a.items()) + ".",
                "", "| Features | Quantum-kernel SVM (CV AUC) | RBF-SVM | Logistic regression |", "|---|---|---|---|"]
        for fs, v in q5["feature_sets"].items():
            out.append(f"| {fs} | {v['qsvm_best']['cv_auc']:.3f} | {v['rbf_best']['cv_auc']:.3f} | {v['logistic']:.3f} |")
        out += ["", f"Chosen: {q5['choice']}.", ""]
    if q6:
        so = q6["szabo_ostlund"]
        out += ["Q6 practice: STO-3G Hartree-Fock against Szabo & Ostlund (Table 3.13, printed to 1 mEh): "
                + "; ".join(f"{k} {v['diff_mEh']:+.2f} mEh" for k, v in so.items()) + ". Against OpenFermion's H2 data (HF / FCI, mEh): "
                + "; ".join(f"{k.split('_')[-1].replace('.hdf5', '')} A {v['diff_hf_mEh']:+.1e} / {v['diff_fci_mEh']:+.1e}"
                            for k, v in q6["openfermion_h2"].items()) + ".", ""]
    if q7:
        out += [f"Q7 practice ({q7['qaoa'][next(iter(q7['qaoa']))]['usable']} usable complexes; settings {q7['settings']}):", "",
                "| QAOA | Max clique found | Mean P(best clique) | Uniform | Docked (QAOA route) | Docked (classical cliques) | Docked (random search) | Mean QAOA s |",
                "|---|---|---|---|---|---|---|---|"]
        for k, v in q7["qaoa"].items():
            out.append(f"| {k} | {100 * v['qaoa_hit_rate']:.0f} % | {v['mean_qaoa_p_optimal']:.4f} | {v['mean_uniform_p_optimal']:.1e} | "
                       f"{100 * v['success_qaoa']:.0f} % | {100 * v['success_classical_clique']:.0f} % | "
                       f"{100 * v['success_random_search']:.0f} % | {v['mean_qaoa_seconds']:.1f} |")
        out += ["", f"Chosen: {q7['choice']}."]
    return "\n".join(out) if out else "_Q5-Q7 practice: not run._"


def qdrug() -> str:
    r = _load("results_qdrug.json")
    if not r or "overall" not in r:
        return "_results_qdrug.json: not run (python validation/quantum_drug_gates.py --test all)._"
    R = r["rule"]
    q5, q6, q7 = r["q5"], r["q6"], r["q7"]
    ci = q5["ci95"]
    out = [f"Q5 (trained on {q5['train']} TDC hERG compounds, tested on {q5['test']} hERG_Karim compounds not in the training "
           f"set ({q5['test_removed_overlap']} removed), {100 * q5['test_blockers']:.0f} % blockers):", "",
           "| Classifier | Test AUC (95 % interval) |", "|---|---|",
           f"| Quantum-kernel SVM (simulated) | {q5['auc']['qsvm']:.3f}{_ci(ci['qsvm'], 3)} |",
           f"| RBF-kernel SVM (classical) | {q5['auc']['rbf']:.3f}{_ci(ci['rbf'], 3)} |",
           f"| Logistic regression (classical) | {q5['auc']['logistic']:.3f}{_ci(ci['logistic'], 3)} |", "",
           f"Quantum minus RBF {q5['auc']['qsvm'] - q5['auc']['rbf']:+.3f}{_ci(ci['qsvm_minus_rbf'], 3)}. Q5: "
           f"**{'pass' if q5['pass'] else 'fail'}** (needed: within {R['q5']['margin']} of the RBF-SVM, lower bound above 0.5).", "",
           f"Q6 (independent reference: OpenFermion's stored data; tolerance {R['q6']['reference_tol_mEh']} mEh):", "",
           "| Reference | HF (ours) | HF (reference) | Difference (mEh) | FCI (ours) | FCI (reference) | Difference (mEh) |",
           "|---|---|---|---|---|---|---|"]
    for x in q6["reference"]:
        out.append(f"| {x['file']} | {x['hf']:.6f} | {x['ref_hf']:.6f} | {x['diff_hf_mEh']:+.1e} | {x['fci']:.6f} | "
                   f"{x['ref_fci']:.6f} | {x['diff_fci_mEh']:+.1e} |")
    out += ["", "Active-space UCCSD-VQE on stretched molecules (chemical accuracy 1.6 mHa):", "",
            "| Molecule | Bond × equilibrium | Active space | Qubits | Parameters | HF | Active-space FCI | VQE | Error (mHa) |",
            "|---|---|---|---|---|---|---|---|---|"]
    for x in q6["vqe"]:
        out.append(f"| {x['molecule']} | {x['bond_scale']} | {x['active'][0]}e, {x['active'][1]}o | {x['qubits']} | {x['parameters']} | "
                   f"{x['e_hf']:.5f} | {x['e_cas_fci']:.5f} | {x['e_vqe']:.5f} | {x['error_mEh']:+.3f} |")
    out += ["", f"Q6: **{'pass' if q6['pass'] else 'fail'}** (reference {'within' if q6['reference_pass'] else 'outside'} "
                f"tolerance; VQE {'within' if q6['vqe_pass'] else 'outside'} chemical accuracy).", "",
            f"Q7 ({q7['usable']} usable test complexes of {q7['complexes']}; 20-qubit interaction graphs):", "",
            "| Measure | Value |", "|---|---|",
            f"| QAOA found the maximum-weight clique | {100 * q7['qaoa_hit_rate']:.0f} % (needed ≥ {100 * R['q7']['min_hit_rate']:.0f} %) |",
            f"| Mean probability of the best clique (uniform guess) | {q7['mean_qaoa_p_optimal']:.4f} ({q7['mean_uniform_p_optimal']:.1e}) |",
            f"| Docked within 2 A: QAOA route | {100 * q7['success_qaoa']:.1f} % |",
            f"| Docked within 2 A: classical cliques, same graph | {100 * q7['success_classical_clique']:.1f} % |",
            f"| Docked within 2 A: random search, same score, 1,000 poses | {100 * q7['success_random_search']:.1f} % |",
            f"| Mean QAOA time per complex | {q7['mean_qaoa_seconds']:.1f} s |", "",
            f"Q7: **{'pass' if q7['pass'] else 'fail'}** (solver {'pass' if q7['solver_pass'] else 'fail'}, docking "
            f"{'pass' if q7['docking_pass'] else 'fail'}: QAOA route at least as good as random search)."]
    return "\n".join(out)


def _qd_rows() -> list[str]:
    rows = []
    r = _load("results_gate8.json")
    if r:
        rows.append(f"| Gate 8: Drug lab vs chromatin tracing after real drug treatment (IMR-90 chrX, 4 drugs) | "
                    f"{r['passing_drugs']} of {len(r['drugs'])} drugs met the rule | {r['verdict']} |")
    else:
        rows.append("| Gate 8: Drug lab vs chromatin tracing after real drug treatment | not run | — |")
    q = _load("results_qdrug.json")
    if q and "overall" in q:
        v = lambda k: "pass" if q["overall"][k] else "fail"
        rows += [f"| Gate Q5: quantum-kernel hERG screen vs RBF-SVM (TDC hERG → hERG_Karim) | AUC {q['q5']['auc']['qsvm']:.3f} vs "
                 f"{q['q5']['auc']['rbf']:.3f} | {v('q5')} |",
                 f"| Gate Q6: molecule energies vs OpenFermion; stretched-molecule VQE | worst reference difference "
                 f"{max(max(abs(x['diff_hf_mEh']), abs(x['diff_fci_mEh'])) for x in q['q6']['reference']):.1e} mHa; worst VQE error "
                 f"{max(abs(x['error_mEh']) for x in q['q6']['vqe']):.2f} mHa | {v('q6')} |",
                 f"| Gate Q7: QAOA max-clique docking (PoseBusters) | docked {100 * q['q7']['success_qaoa']:.0f} % vs random search "
                 f"{100 * q['q7']['success_random_search']:.0f} %; clique found {100 * q['q7']['qaoa_hit_rate']:.0f} % | {v('q7')} |"]
    else:
        rows.append("| Gates Q5-Q7: quantum drug tabs | not run | — |")
    return rows


def summary_qd() -> str:
    return "\n".join(["| Test (held-out, real data) | Measured | Verdict |", "|---|---|---|"] + _qd_rows())


# --------------------------------------------------------------------------------------------- round 2 (quantum)
def round2_practice() -> str:
    out = []
    q2 = _load("results_round2_practice_q2b.json")
    if q2:
        c, s = q2["choice"], q2["choice_scores"]
        out += [f"Q2b practice: {len(q2['grid'])} QUBO settings, exact optimum on the nine seen windows (GM12878, K562, "
                f"IMR-90). Chosen: {c}. F1 of the QUBO / better classical caller (tuned on the same windows, same "
                "resolution): " + "; ".join(f"{k} {s['f1'][k]:.3f} / {s['classical_best'][k]:.3f}" for k in s["f1"])
                + f". Classical settings: {q2['classical_choice']}.", ""]
    q4 = _load("results_round2_practice_q4b.json")
    if q4:
        out += [f"Q4b practice: {sum(q4['genes'].values())} genes ({', '.join(f'{k} {v}' for k, v in q4['genes'].items())}); "
                "ENCODE RNA-seq labels agree with Q4's GTEx labels for "
                + ", ".join(f"{k} {100 * v['agree']:.0f} %" for k, v in q4["agreement_with_gtex_labels"].items())
                + " of genes. Leave-one-cell-line-out AUC:", "",
                "| Features | Quantum-kernel SVM | RBF-SVM | Logistic regression |", "|---|---|---|---|"]
        for k, v in q4["feature_sets"].items():
            out.append(f"| {k} ({len(v['features'])}) | {v['qsvm_best']['mean']:.3f} | {v['rbf_best']['mean']:.3f} | "
                       f"{v['logistic']['mean']:.3f} |")
        out += ["", f"Chosen: {q4['choice']}.", ""]
    q6 = _load("results_round2_practice_q6b.json")
    if q6:
        out += ["Q6b practice (error vs active-space FCI, mHa; chemical accuracy 1.6):", "",
                "| Molecule | Bond × eq. | Active space | Qubits | Fixed UCCSD | ADAPT, occupied→virtual pool | ADAPT, generalized pool (operators) |",
                "|---|---|---|---|---|---|---|"]
        sd, gsd = q6["grid"]["sd"]["rows"], q6["grid"]["gsd"]["rows"]
        for a, b in zip(sd, gsd):
            out.append(f"| {a['molecule']} | {a['bond_scale']} | {a['active'][0]}e, {a['active'][1]}o | {a['qubits']} | "
                       f"{a['uccsd_error_mEh']:.3f} | {a['adapt_error_mEh']:.3f} | {b['adapt_error_mEh']:.3f} ({b['adapt_operators']}) |")
        out += ["", f"Chosen: {q6['choice']} (the operator cap was raised to 150 in the frozen rule).", ""]
    q7 = _load("results_round2_practice_q7b.json")
    if q7:
        out += [f"Q7b practice ({q7['complexes']} PoseBusters complexes of Q7, all practice now; QAOA {q7['qaoa']}):", "",
                "| Setting | Usable | Docked: QAOA route | Classical cliques | Random search (same score, same refinement) | QAOA unrefined |",
                "|---|---|---|---|---|---|"]
        for k, v in q7["grid"].items():
            lab = (f"clique poses seed a local search ({v['settings']['seeded']} seeds), {v['settings']['refine_top']} refined"
                   if v["settings"].get("seeded") else f"{v['settings']['refine_top']} poses refined")
            out.append(f"| {lab} | {v['usable']} | {100 * v['success_qaoa']:.1f} % | "
                       f"{100 * v['success_classical']:.1f} % | {100 * v['success_random']:.1f} % | {100 * v['success_qaoa_unrefined']:.1f} % |")
        vh = q7.get("virtual_hydrogens")
        if vh:
            same = q7.get("same_subset_given_hydrogens", {})
            out += ["", f"The chosen setting with the structures' hydrogens removed and virtual polar hydrogens added (as the "
                        f"test structures need), on the {vh['usable']} usable Q7-practice complexes: QAOA route "
                        f"{100 * vh['success_qaoa']:.1f} %, classical cliques {100 * vh['success_classical']:.1f} %, random "
                        f"search {100 * vh['success_random']:.1f} % (with the structures' own hydrogens: "
                        f"{100 * same.get('success_qaoa', float('nan')):.1f} %, {100 * same.get('success_classical', float('nan')):.1f} %, "
                        f"{100 * same.get('success_random', float('nan')):.1f} %)."]
        out += ["", f"Chosen: {q7.get('choice')}."]
    return "\n".join(out) if out else "_Round 2 practice: not run._"


def round2() -> str:
    r = _load("results_round2.json")
    if not r:
        return "_results_round2.json: not run (python validation/quantum_round2.py --test q2b|q4b|q6b|q7b)._"
    out = []
    if "q6b" in r:
        q = r["q6b"]
        out += [f"Q6b ({q['cases']} cases never run before; settings {q['rule']['settings']}):", "",
                "| Molecule | Bond × eq. | Active space | Qubits | HF | Active-space FCI | ADAPT-VQE | ADAPT error (mHa) | Operators | Fixed UCCSD error (mHa) |",
                "|---|---|---|---|---|---|---|---|---|---|"]
        for x in q["rows"]:
            out.append(f"| {x['molecule']} | {x['bond_scale']} | {x['active'][0]}e, {x['active'][1]}o | {x['qubits']} | {x['e_hf']:.5f} | "
                       f"{x['e_cas_fci']:.5f} | {x['adapt_energy']:.5f} | {x['adapt_error_mEh']:+.3f} | {x['adapt_operators']} | "
                       f"{x['uccsd_error_mEh']:+.3f} |")
        out += ["", f"Q6b: **{'pass' if q['pass'] else 'fail'}**: ADAPT-VQE within chemical accuracy in {q['within']} of "
                    f"{q['cases']} (worst {q['max_abs_error_mEh']:.3f} mHa); fixed-order UCCSD in {q['uccsd_within']} of {q['cases']}.", ""]
        ph = _load("results_round2_q6b_posthoc.json")
        if ph:
            trip = [x for x in ph["rows"] if x["sector_ground_S2"] > 1e-3]
            within = sum(abs(x["adapt_error_vs_singlet_mEh"]) <= 1.6 for x in ph["rows"])
            uw = sum(abs(x["uccsd_error_vs_singlet_mEh"]) <= 1.6 for x in ph["rows"])
            out += ["Found after the test (not part of the verdict): " + "; ".join(
                        f"{x['molecule']} at {x['bond_scale']}x: the sector's lowest state is a triplet (S(S+1) "
                        f"{x['sector_ground_S2']:.1f}), {x['singlet_minus_sector_mEh']:.3f} mHa below the lowest singlet; "
                        f"ADAPT-VQE is {x['adapt_error_vs_singlet_mEh']:+.4f} mHa from that singlet" for x in trip)
                    + f". Against the lowest singlet, ADAPT-VQE is within chemical accuracy in {within} of {len(ph['rows'])} "
                      f"cases, fixed-order UCCSD in {uw}. The app now reports the spin of the exact state and the lowest "
                      "singlet next to it.", ""]
    if "q6c" in r:
        q = r["q6c"]
        out += [f"Q6c (reference corrected to the lowest singlet; {q['cases']} cases never run; settings {q['rule']['settings']}):", "",
                "| Molecule | Bond × eq. | Active space | Qubits | Lowest singlet | Sector's lowest state (S(S+1)) | ADAPT-VQE | ADAPT error vs singlet (mHa) | Operators | Fixed UCCSD error vs singlet (mHa) |",
                "|---|---|---|---|---|---|---|---|---|---|"]
        for x in q["rows"]:
            out.append(f"| {x['molecule']} | {x['bond_scale']} | {x['active'][0]}e, {x['active'][1]}o | {x['qubits']} | "
                       f"{x['lowest_singlet']:.5f} | {x['e_cas_fci']:.5f} ({x['sector_ground_S2']:.1f}) | {x['adapt_energy']:.5f} | "
                       f"{x['adapt_error_vs_singlet_mEh']:+.3f} | {x['adapt_operators']} | {x['uccsd_error_vs_singlet_mEh']:+.3f} |")
        out += ["", f"Q6c: **{'pass' if q['pass'] else 'fail'}**: ADAPT-VQE within chemical accuracy of the lowest singlet in "
                    f"{q['within']} of {q['cases']} (worst {q['max_abs_error_mEh']:.3f} mHa); fixed-order UCCSD in "
                    f"{q['uccsd_within']} of {q['cases']}; the sector's lowest state was a triplet in {q['triplet_ground']}.", ""]
    if "q4b" in r:
        q = r["q4b"]
        ci = q["ci95"]
        out += [f"Q4b (trained on {q['train_genes']} genes of GM12878, K562 and IMR-90; tested on {q['test_genes']} HMEC genes, "
                f"{100 * q['test_expressed_fraction']:.0f} % active):", "", "| Classifier | Test AUC (95 % interval) |", "|---|---|",
                f"| Quantum-kernel SVM (simulated) | {q['auc']['qsvm']:.3f}{_ci(ci['qsvm'], 3)} |",
                f"| RBF-kernel SVM (classical) | {q['auc']['rbf']:.3f}{_ci(ci['rbf'], 3)} |",
                f"| Logistic regression (classical) | {q['auc']['logistic']:.3f}{_ci(ci['logistic'], 3)} |",
                f"| Quantum-kernel SVM, kernel from 1,000 shots per entry | {q['qsvm_auc_with_shots']:.3f} |", "",
                f"Quantum minus RBF {q['auc']['qsvm'] - q['auc']['rbf']:+.3f}{_ci(ci['qsvm_minus_rbf'], 3)}.", "",
                f"Q4b: **{'pass' if q['pass'] else 'fail'}** (needed: within {q['rule']['margin']} of the RBF-SVM, lower bound "
                "above 0.5).", ""]
    if "q2b" in r:
        q = r["q2b"]
        out += [f"Q2b (new cell lines; settings {q['rule']['tad']}):", "",
                "| Cell line | Method | Precision | Recall | F1 |", "|---|---|---|---|---|"]
        for cell, v in q["per_cell"].items():
            for m, s in v["rows"].items():
                out.append(f"| {cell} | {m} | {_f(s['precision'], 2)} | {_f(s['recall'], 2)} | {_f(s['f1'], 3)} |")
        out += ["", "QAOA found the QUBO's optimum in " + "; ".join(f"{c} {100 * v['hit_rate']:.0f} % of {v['windows']} windows"
                                                              for c, v in q["q1_like"].items()) + ".",
                f"Q2b: **{'pass' if q['pass'] else 'fail'}** (" + "; ".join(f"{c} {'pass' if v['pass'] else 'fail'}"
                                                                         for c, v in q["per_cell"].items()) + ").", ""]
    if "q7b" in r:
        q = r["q7b"]
        out += [f"Q7b ({q['usable']} usable Astex Diverse complexes of {q['complexes']}):", "", "| Measure | Value |", "|---|---|",
                f"| QAOA found the maximum-weight clique | {100 * q['qaoa_hit_rate']:.0f} % |",
                f"| Docked within 2 A: QAOA route, Vina-like score, refined | {100 * q['success_qaoa']:.1f} % |",
                f"| Docked within 2 A: classical cliques, same | {100 * q['success_classical']:.1f} % |",
                f"| Docked within 2 A: random search, same score and refinement | {100 * q['success_random']:.1f} % |",
                f"| Docked within 2 A: QAOA route without refinement | {100 * q['success_qaoa_unrefined']:.1f} % |", "",
                f"Q7b: **{'pass' if q['pass'] else 'fail'}** (clique found in ≥ {100 * q['rule']['min_hit_rate']:.0f} %: "
                f"{'yes' if q['solver_pass'] else 'no'}; QAOA route at least as good as random search: "
                f"{'yes' if q['vs_random_pass'] else 'no'}; QAOA route docks ≥ {100 * q['rule']['min_success']:.0f} %: "
                f"{'yes' if q['improvement_pass'] else 'no'}).", ""]
    return "\n".join(out)


def _r2_rows() -> list[str]:
    r = _load("results_round2.json") or {}
    rows = []
    if "q2b" in r:
        q = r["q2b"]
        f1 = "; ".join(f"{c}: QAOA {v['rows']['QAOA (simulated quantum)']['f1']:.2f} vs insulation "
                       f"{v['rows']['insulation (classical)']['f1']:.2f}, TopDom-like {v['rows']['topdom (classical)']['f1']:.2f}"
                       for c, v in q["per_cell"].items())
        rows.append(f"| Gate Q2b: quantum domain calls vs classical callers, new cell lines | {f1} | {'pass' if q['pass'] else 'fail'} |")
    if "q4b" in r:
        q = r["q4b"]
        rows.append(f"| Gate Q4b: quantum-kernel gene classifier vs RBF-SVM, three cell lines → HMEC | AUC {q['auc']['qsvm']:.3f} vs "
                    f"{q['auc']['rbf']:.3f} | {'pass' if q['pass'] else 'fail'} |")
    if "q6b" in r:
        q = r["q6b"]
        rows.append(f"| Gate Q6b: ADAPT-VQE on {q['cases']} new stretched molecules | worst error {q['max_abs_error_mEh']:.2f} mHa "
                    f"(needed ≤ 1.6); {q['within']} of {q['cases']} within | {'pass' if q['pass'] else 'fail'} |")
    if "q6c" in r:
        q = r["q6c"]
        rows.append(f"| Gate Q6c: ADAPT-VQE on {q['cases']} new stretched molecules, against the lowest singlet | worst error "
                    f"{q['max_abs_error_mEh']:.2f} mHa (needed ≤ 1.6); {q['within']} of {q['cases']} within | "
                    f"{'pass' if q['pass'] else 'fail'} |")
    if "q7b" in r:
        q = r["q7b"]
        rows.append(f"| Gate Q7b: QAOA docking with Vina-like score and refinement (Astex Diverse) | docked "
                    f"{100 * q['success_qaoa']:.0f} % vs random search {100 * q['success_random']:.0f} % | {'pass' if q['pass'] else 'fail'} |")
    return rows or ["| Round 2 (Q2b, Q4b, Q6b, Q7b) | not run | — |"]


def summary_r2() -> str:
    return "\n".join(["| Test (held-out, real data; simulated quantum) | Measured | Verdict |", "|---|---|---|"] + _r2_rows())


BLOCKS = {"gate1": gate1, "gate2": gate2, "gate2b": gate2b, "gate2c": gate2c, "gate3": gate3, "gate4": gate4,
          "gate5": gate5, "gate1c": gate1c, "gate2d": gate2d, "gate2e": gate2e, "gate3b": gate3b, "gate5b": gate5b,
          "gate4c": gate4c, "gate5m": gate5m, "gate6": gate6, "gate7": gate7, "b9tools": b9tools, "summary_ab": summary_ab, "scale": scale,
          "per_chromosome": per_chromosome, "readme_accuracy": readme_accuracy, "gateq_practice": gateq_practice,
          "gateq": gateq, "quantum_crosscheck": quantum_crosscheck, "summary_q": summary_q, "gate8": gate8,
          "qdrug_practice": qdrug_practice, "qdrug": qdrug, "summary_qd": summary_qd,
          "round2_practice": round2_practice, "round2": round2, "summary_r2": summary_r2, "gate6b": gate6b,
          "gate8_marks": gate8_marks}
TARGETS = (RESULTS_MD, ROOT.parent / "README.md")


def render(text: str) -> str:
    for name, fn in BLOCKS.items():
        pat = re.compile(rf"(<!-- BEGIN generated:{name} -->\n)(.*?)(<!-- END generated:{name} -->)", re.S)
        if pat.search(text):
            text = pat.sub(lambda m: m.group(1) + fn().rstrip() + "\n" + m.group(3), text)
    return text


def main() -> None:
    for path in TARGETS:
        text = path.read_text(encoding="utf-8")
        new = render(text)
        path.write_text(new, encoding="utf-8", newline="\n")
        print(f"{path.name}: {sum(f'BEGIN generated:{n}' in new for n in BLOCKS)} generated blocks refreshed")
    if "--print" in sys.argv:
        for n, fn in BLOCKS.items():
            print(f"\n=== {n}\n{fn()}")


if __name__ == "__main__":
    main()
