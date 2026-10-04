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


BLOCKS = {"gate1": gate1, "gate2": gate2, "gate3": gate3, "gate4": gate4, "gate5": gate5, "scale": scale,
          "per_chromosome": per_chromosome}


def render(text: str) -> str:
    for name, fn in BLOCKS.items():
        pat = re.compile(rf"(<!-- BEGIN generated:{name} -->\n)(.*?)(<!-- END generated:{name} -->)", re.S)
        text = pat.sub(lambda m: m.group(1) + fn().rstrip() + "\n" + m.group(3), text)
    return text


def main() -> None:
    text = RESULTS_MD.read_text(encoding="utf-8")
    new = render(text)
    RESULTS_MD.write_text(new, encoding="utf-8", newline="\n")
    print(f"{RESULTS_MD}: {sum(f'BEGIN generated:{n}' in new for n in BLOCKS)} generated blocks refreshed")
    if "--print" in sys.argv:
        for n, fn in BLOCKS.items():
            print(f"\n=== {n}\n{fn()}")


if __name__ == "__main__":
    main()
