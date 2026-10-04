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
    m = r["meta"]["machine"]
    out = [f"Machine: {m['platform']}, {m['cpus']} logical CPUs, torch {m['torch']}, CUDA {m['cuda']}. "
           f"{r['meta']['note']}", "",
           "| Beads | Model | Fit (s) | Sampling (s) | Total (s) | Peak memory (MB) | Device | Rank |", "|---|---|---|---|---|---|---|---|"]
    for row in r["rows"]:
        if "skipped" in row:
            out.append(f"| {row['n_beads']:,} | {row['method']} | — | — | — | — | — | {row['skipped']} |")
        elif "error" in row:
            out.append(f"| {row['n_beads']:,} | {row['method']} | failed | | | | | {row['error'][-80:]} |")
        else:
            out.append(f"| {row['n_beads']:,} | {row['method']} | {row['seconds_fit']:.1f} | {row['seconds_sampling']:.1f} | "
                       f"{row['seconds_total']:.1f} | {row['peak_memory_mb']:,.0f} | {row['device']} | {row['rank']} |")
    return "\n".join(out)


def per_chromosome() -> str:
    r = _load("chromosome_runtime.json")
    if not r:
        return "_validation/chromosome_runtime.json: not run (python validation/chromosome_runtime.py)._"
    out = [r["meta"]["note"], "", "| Assembly | Chromosome | Bins | Beads fitted | Fit (s) | Total (s) | Peak memory (MB) |",
           "|---|---|---|---|---|---|---|"]
    for row in r["rows"]:
        if "error" in row:
            out.append(f"| {row['assembly']} | {row['chrom']} | {row.get('bins', '—')} | failed: {row['error'][-60:]} | | | |")
        else:
            out.append(f"| {row['assembly']} | {row['chrom']} | {row['bins']:,} | {row['beads']:,} | {row['seconds_fit']:.1f} | "
                       f"{row['seconds_total']:.1f} | {row['peak_memory_mb']:,.0f} |")
    return "\n".join(out)


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
    return "\n".join(rows)


BLOCKS = {"gate1": gate1, "gate2": gate2, "gate2b": gate2b, "gate2c": gate2c, "gate3": gate3, "gate4": gate4,
          "gate5": gate5, "scale": scale,
          "per_chromosome": per_chromosome, "readme_accuracy": readme_accuracy}
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
