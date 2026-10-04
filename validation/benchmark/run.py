"""
Benchmark harness (Pillar 3), one command:

    python -m validation.benchmark.run                 # test datasets  -> validation/benchmark/results.json + RESULTS_TABLE.md
    python -m validation.benchmark.run --practice      # practice datasets -> results_practice.json + RESULTS_TABLE_practice.md
    python -m validation.benchmark.run --datasets A B --out part_x   # a subset of the plan -> part_x.json
    python -m validation.benchmark.run --merge part_x part_y [--not-finished '{"unit": "why"}']
                                                       # parts -> results.json + RESULTS_TABLE.md

Every method (validation/benchmark/methods.py) gets the same input for a split and is scored the same
way against the same held-out truth (validation/protocol.py):
  - distance-pattern recovery beyond the genomic-separation trend: trend-removed Spearman as % of the
    ceiling (half A vs half B), within the windowed model's tiles and on all pairs;
  - rank agreement: raw Spearman;
  - absolute size: Lin's concordance correlation in nm, and the median model / measured ratio;
  - calibration (population models): share of half-B single-copy distances inside the stated 50 / 80
    / 90 % intervals, raw and with the recalibration frozen in chronocell/data/calibration.json.
95 % intervals come from resampling half B's copies (truth noise; split 0, pair subsample).

Datasets: the registry roles in validation/datasets.py decide practice vs test; the harness refuses to
mix them. Inputs: imaging-derived contacts for every dataset, and sequencing Hi-C (Rao et al. 2014,
same cell line) where it exists (Bintu IMR-90 regions, Su chr21 / chr2 / genome-scale; K562 on practice).
Hi-C settings are the frozen Gate 1 choices (validation/frozen.py).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                       # noqa: E402
import protocol as PR                      # noqa: E402
from chronocell import ensemble as E      # noqa: E402
from validation.benchmark import methods as M   # noqa: E402

TEST = [("bintu_imr90_28_30", 3), ("bintu_imr90_18_20", 3), ("bintu_a549_28_30", 3), ("bintu_hct116_34_37_auxin", 3),
        ("su_chr21", 1), ("su_chr21_rep", 1), ("su_genome", 1)]
PRACTICE = [("bintu_k562_28_30", 3), ("bintu_hct116_28_30", 3), ("bintu_hct116_28_30_auxin", 3),
            ("bintu_hct116_34_37", 3), ("su_chr2", 1), ("su_chr2_parm_rep", 1)]
CACHE = ROOT / "data" / "cache"
BOOT = {"bintu_csv": 100, "su_trace": 30, "su_genome": 30}


def _cal_pits():
    try:
        cal = json.loads((ROOT.parent / "chronocell" / "data" / "calibration.json").read_text())
        return PR.isotonic_recalibration(np.asarray(cal["pit_histogram"]))
    except (OSError, ValueError, KeyError):
        return None


def _stats(key: str, xyz: np.ndarray, split: int, r_c, tag: str = "") -> dict:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"bench_{key}{tag}_split{split}_{'adj' if r_c is None else int(r_c)}.npz"
    if path.exists():
        z = np.load(path)
        return {k: z[k] for k in z.files}
    a, b = PR.split(len(xyz), split)
    ha = PR.half_stats(xyz[a], r_c)
    rc = ha.adjacent_median if r_c is None else float(r_c)
    hb = PR.half_stats(xyz[b], rc)
    out = {"freq_a": ha.freq, "seen_a": ha.seen, "median_a": ha.median, "median_b": hb.median, "r_c": np.array(rc),
           "b_idx": b}
    np.savez(path, **out)
    return out


def _hic_input(counts: np.ndarray):
    from frozen import GATE1
    c = np.asarray(counts, dtype=float).copy()
    np.fill_diagonal(c, 0.0)
    p_adj = float(GATE1["p_adjacent"])
    p = E.counts_to_probability(c, p_adj)
    adj = np.diag(c, 1)
    n_eff = float(np.median(adj[adj > 0])) / p_adj
    if GATE1.get("hic_zeros") == "unobserved":
        p = np.where(c > 0, p, np.nan)
    return p, n_eff, c


def _score_pred(pred: M.Prediction, st: dict, sep: np.ndarray, masks: dict, xyz_b: np.ndarray, boot: int,
                pits, do_boot: bool, seed: int) -> dict:
    out = {"seconds": round(pred.seconds, 2), "contact_map_fit": None if not np.isfinite(pred.contact_fit)
           else round(float(pred.contact_fit), 4), "notes": pred.notes}
    for mname, m in masks.items():
        if not m.any() or np.isfinite(pred.median_nm[m]).mean() < 0.99:
            continue
        ceil = PR.scores(st["median_a"], st["median_b"], sep, m)["spearman_distance_corrected"]
        s = PR.scores(pred.median_nm, st["median_b"], sep, m)
        row = {"percent_of_ceiling": 100 * s["spearman_distance_corrected"] / ceil, "trend_removed_rho": s["spearman_distance_corrected"],
               "ceiling": ceil, "raw_spearman": s["spearman"], "lin_ccc_nm": s["lin_ccc_nm"],
               "median_scale": s["median_scale_model_over_real"], "pairs": s["pairs"]}
        if do_boot and boot:
            row["ci95"] = _bootstrap(pred.median_nm, st["median_a"], xyz_b, sep, m, boot, seed)
        if pred.sigma_nm is not None:
            cov = PR.coverage(xyz_b, pred.sigma_nm, m)
            row["coverage"] = dict(zip(["50", "80", "90"], [round(v, 4) for v in cov["coverage"]]))
            if pits is not None:
                rc = PR.coverage_with_pits(xyz_b, pred.sigma_nm, m, pits)
                row["coverage_recalibrated"] = dict(zip(["50", "80", "90"], [round(v, 4) for v in rc["coverage"]]))
        out[mname] = row
    return out


def _bootstrap(pred: np.ndarray, median_a: np.ndarray, xyz_b: np.ndarray, sep: np.ndarray, mask: np.ndarray,
               reps: int, seed: int, max_pairs: int = 20_000) -> dict:
    rng = np.random.default_rng(seed)
    iu = np.argwhere(mask)
    if len(iu) > max_pairs:
        iu = iu[rng.choice(len(iu), max_pairs, replace=False)]
    i, j = iu[:, 0], iu[:, 1]
    d = np.linalg.norm(np.asarray(xyz_b, np.float32)[:, i] - np.asarray(xyz_b, np.float32)[:, j], axis=-1)
    p, a, s = pred[i, j], median_a[i, j], sep[i, j]
    pct, raw, ccc = [], [], []
    for _ in range(reps):
        pick = rng.integers(0, len(d), len(d))
        t = PR.nanmedian0(d[pick])
        sp = PR._flat_scores(p, t, s)
        sc = PR._flat_scores(a, t, s)
        pct.append(100 * sp["spearman_distance_corrected"] / sc["spearman_distance_corrected"])
        raw.append(sp["spearman"])
        ccc.append(sp["lin_ccc_nm"])
    q = lambda v: [round(float(x), 4) for x in np.nanpercentile(v, [2.5, 97.5])]  # noqa: E731
    return {"percent_of_ceiling": q(pct), "raw_spearman": q(raw), "lin_ccc_nm": q(ccc), "reps": reps,
            "pairs_sampled": int(len(i))}


def _masks(n: int, tiles: list) -> dict:
    within = np.zeros((n, n), bool)
    for a, b in tiles:
        within[a:b, a:b] = True
    upper = np.triu(np.ones((n, n), bool), 1)
    return {"within_tiles": upper & within, "all_pairs": upper} if len(tiles) > 1 else {"all_pairs": upper}


def run_units(key: str, splits: int, methods: list[str], pits, log=print) -> list[dict]:
    """Units = (dataset or chromosome, input, split). Returns one row per unit and method."""
    e = D.REGISTRY[key]
    rows = []
    if e.kind == "su_genome":
        tr = D.load(key)
        H, hch, hst = D.load_su_hic(D.REGISTRY[e.paired_hic])
        for c in [c for c in dict.fromkeys(tr.chrom) if c not in ("chr2", "chrY")]:     # chr2 overlaps practice
            cols = tr.chrom == c
            own = np.isfinite(tr.xyz[:, cols, 0]).mean(axis=1) >= 0.5
            xyz = tr.xyz[own][:, cols].astype(np.float64)
            starts = tr.starts[cols]
            hsel = (hch == c)
            if not np.array_equal(hst[hsel], starts):
                log(f"{c}: Hi-C loci do not match the imaged loci; Hi-C input skipped")
                hic = None
            else:
                hic = H[np.ix_(hsel, hsel)]
            if len(xyz) < 50 or cols.sum() < 8:
                continue
            rows += _run_one(f"{key}:{c}", e, xyz, starts, hic, 1, methods, pits, None, log, tag=f"_{c}")
        return rows
    tr = D.load(key)
    hic = None
    if e.kind == "bintu_csv":
        try:
            hic, _ = D.bintu_hic(key)
        except KeyError:
            hic = None
        r_c = 150.0
    else:
        H, _, hst = D.load_su_hic(D.REGISTRY[e.paired_hic])
        hic = H if np.array_equal(hst, tr.starts) else None
        r_c = None
    return _run_one(key, e, tr.xyz, tr.starts, hic, splits, methods, pits, r_c, log)


def _run_one(name, e, xyz, starts, hic, splits, methods, pits, r_c, log, tag: str = "") -> list[dict]:
    rows = []
    n = xyz.shape[1]
    sep = PR.separation(n, starts)
    tiles = PR.tiles(n, 400)
    masks = _masks(n, tiles)
    for split in range(splits):
        st = _stats(e.key, xyz, split, r_c, tag)
        xyz_b = xyz[st["b_idx"]]
        rc = float(st["r_c"])
        inputs = {"imaging": M.Input("imaging", n, st["freq_a"], st["seen_a"], st["freq_a"] * np.maximum(st["seen_a"], 1),
                                     rc, st["median_a"], sep, tiles, split)}
        if hic is not None:
            p, n_eff, c = _hic_input(hic)
            # same contact-radius anchor convention as the app: adjacent beads sit at the literature b0
            from chronocell import physics
            from frozen import GATE1
            b0 = physics.bond_length_for(e.step_bp)
            r_h = b0 / float(E.gaussian_median_distance(float(GATE1["p_adjacent"]), 1.0))
            inputs["hic"] = M.Input("hic", n, p, n_eff, c, r_h, st["median_a"], sep, tiles, split)
        for kind, inp in inputs.items():
            for mname in methods:
                if mname == "v3_2_single" and n > 2000:
                    continue
                if mname == "genomic_distance_only" and kind == "hic":
                    continue                                   # identical to the imaging row (uses half A medians)
                t0 = time.time()
                try:
                    pred = M.METHODS[mname](inp)
                    sc = _score_pred(pred, st, sep, masks, xyz_b, BOOT.get(e.kind, 30), pits, split == 0, split)
                    status = "ok"
                except Exception as exc:          # a method that fails is reported, never dropped silently
                    sc, status = {"error": f"{type(exc).__name__}: {exc}"}, "error"
                rows.append({"dataset": name, "role": e.role, "input": kind, "split": split, "method": mname,
                             "loci": n, "status": status, **sc})
                ap = sc.get("all_pairs", {})
                log(f"{name:28s} {kind:7s} s{split} {mname:22s} "
                    + (f"all {ap.get('percent_of_ceiling', float('nan')):6.1f}% raw {ap.get('raw_spearman', float('nan')):.3f} "
                       f"ccc {ap.get('lin_ccc_nm', float('nan')):.3f} ({time.time() - t0:.0f}s)" if status == "ok"
                       else sc["error"][:120]))
    return rows


def summarise(rows: list[dict]) -> dict:
    """Mean over splits per (dataset, input, method, pair set)."""
    out: dict = {}
    for r in rows:
        if r["status"] != "ok":
            continue
        for mname in ("within_tiles", "all_pairs"):
            if mname not in r:
                continue
            key = (r["dataset"], r["input"], r["method"], mname)
            out.setdefault(key, []).append(r[mname])
    summ = {}
    for (ds, inp, meth, mname), vals in out.items():
        agg = {k: float(np.nanmean([v[k] for v in vals])) for k in ("percent_of_ceiling", "raw_spearman", "lin_ccc_nm",
                                                                     "median_scale", "ceiling")}
        agg["splits"] = len(vals)
        ci = next((v["ci95"] for v in vals if "ci95" in v), None)
        if ci:
            agg["ci95"] = ci
        for cov in ("coverage", "coverage_recalibrated"):
            cs = [v[cov] for v in vals if cov in v]
            if cs:
                agg[cov] = {lv: float(np.mean([c[lv] for c in cs])) for lv in ("50", "80", "90")}
        summ.setdefault(ds, {}).setdefault(inp, {}).setdefault(mname, {})[meth] = agg
    return summ


def where_we_lose(summ: dict) -> list[str]:
    """Every (dataset, input, pair set, metric) where a non-ChronoCell method or baseline beats the best
    ChronoCell population model, or where a stated interval misses its coverage by > 5 points."""
    out = []
    for ds, inputs in summ.items():
        for inp, msets in inputs.items():
            for mname, meths in msets.items():
                ours = {k: v for k, v in meths.items() if k in M.OURS}
                if not ours:
                    continue
                for metric in ("percent_of_ceiling", "raw_spearman", "lin_ccc_nm"):
                    best_name, best = max(((k, v[metric]) for k, v in ours.items()), key=lambda kv: kv[1])
                    for other, v in meths.items():
                        if other in M.OURS or not np.isfinite(v.get(metric, np.nan)):
                            continue
                        if v[metric] > best:
                            out.append(f"{ds} · {inp} · {mname} · {metric}: {other} {v[metric]:.3f} > {best_name} {best:.3f}")
                for name, v in ours.items():
                    for cov in ("coverage", "coverage_recalibrated"):
                        for lv, val in v.get(cov, {}).items():
                            if abs(val - int(lv) / 100) > 0.05:
                                out.append(f"{ds} · {inp} · {mname} · {name} {cov} at {lv} %: {100 * val:.1f} %")
    return out


def table_md(summ: dict, title: str, losses: list[str], not_run: dict, meta: dict) -> str:
    lines = [f"# {title}", "", f"Generated {meta['utc']} by `python -m validation.benchmark.run{' --practice' if meta['practice'] else ''}` "
             f"in {meta['seconds']:.0f} s. Every number is computed by the harness from the data; nothing is typed in.", "",
             "Columns: % of ceiling = trend-removed Spearman / ceiling (distance-pattern recovery beyond the "
             "genomic-separation trend); raw ρ = Spearman; CCC = Lin's concordance in nm (absolute size); "
             "coverage = share of held-out single-copy distances inside the stated 50 / 80 / 90 % intervals "
             "(raw → recalibrated). Brackets: 95 % interval from resampling half B's copies (split 0).", ""]
    for ds, inputs in summ.items():
        for inp, msets in inputs.items():
            for mname, meths in msets.items():
                lines += [f"## {ds} · input: {inp} · {mname.replace('_', ' ')}", "",
                          "| Method | % of ceiling | raw ρ | CCC (nm) | size ratio | coverage 50/80/90 % |",
                          "|---|---|---|---|---|---|"]
                for meth, v in sorted(meths.items(), key=lambda kv: -kv[1]["percent_of_ceiling"]):
                    ci = v.get("ci95", {})
                    pct = f"{v['percent_of_ceiling']:.1f}" + (f" [{ci['percent_of_ceiling'][0]:.1f}, {ci['percent_of_ceiling'][1]:.1f}]" if ci else "")
                    raw = f"{v['raw_spearman']:.3f}" + (f" [{ci['raw_spearman'][0]:.3f}, {ci['raw_spearman'][1]:.3f}]" if ci else "")
                    ccc = f"{v['lin_ccc_nm']:.3f}" + (f" [{ci['lin_ccc_nm'][0]:.3f}, {ci['lin_ccc_nm'][1]:.3f}]" if ci else "")
                    cov = ""
                    if "coverage" in v:
                        cov = "/".join(f"{100 * v['coverage'][k]:.0f}" for k in ("50", "80", "90"))
                        if "coverage_recalibrated" in v:
                            cov += " → " + "/".join(f"{100 * v['coverage_recalibrated'][k]:.0f}" for k in ("50", "80", "90"))
                    lines.append(f"| {meth} | {pct} | {raw} | {ccc} | {v['median_scale']:.2f} | {cov} |")
                lines.append("")
    lines += ["## Where we lose", ""] + ([f"- {x}" for x in losses] or ["- nowhere on these data"]) + [""]
    lines += ["## Methods not run", ""] + [f"- **{k}**: {v}" for k, v in not_run.items()] + [""]
    return "\n".join(lines)


def merge(parts: list[str], practice: bool, not_finished: dict[str, str]) -> None:
    """Combine runs made on subsets of the plan (--datasets / --methods, each with its own --out) into
    results.json and RESULTS_TABLE.md, recomputing the summary from their rows. not_finished records units
    that were started but stopped before finishing (unit -> reason)."""
    rows, metas = [], []
    for stem in parts:
        d = json.loads((HERE / f"{stem}.json").read_text())
        if d["meta"]["practice"] != practice:
            sys.exit(f"{stem} is a {'practice' if d['meta']['practice'] else 'test'} run; refusing to mix roles")
        rows += d["rows"]
        metas.append({"part": stem, **d["meta"]})
    summ = summarise(rows)
    losses = where_we_lose(summ)
    meta = {"utc": max(m["utc"] for m in metas), "practice": practice, "seconds": sum(m["seconds"] for m in metas),
            "methods": list(dict.fromkeys(x for m in metas for x in m["methods"])),
            "plan": [p for m in metas for p in m["plan"]], "parts": metas, "not_finished": not_finished}
    stem = "results_practice" if practice else "results"
    (HERE / f"{stem}.json").write_text(json.dumps({"meta": meta, "summary": summ, "where_we_lose": losses,
                                                   "not_run": M.NOT_RUN, "not_finished": not_finished, "rows": rows},
                                                  indent=1, default=float))
    md = table_md(summ, "Benchmark: practice datasets (tuning allowed; in-sample)" if practice else
                  "Benchmark: held-out test datasets", losses, dict(M.NOT_RUN, **not_finished), meta)
    (HERE / ("RESULTS_TABLE_practice.md" if practice else "RESULTS_TABLE.md")).write_text(md, encoding="utf-8")
    print(f"-> {HERE / (stem + '.json')} from {len(parts)} parts; {len(losses)} 'where we lose' entries")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--practice", action="store_true")
    ap.add_argument("--datasets", nargs="*")
    ap.add_argument("--methods", nargs="*", default=list(M.METHODS))
    ap.add_argument("--splits", type=int, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--merge", nargs="*", help="combine the named part runs (their --out stems) into results.json")
    ap.add_argument("--not-finished", default="{}", help='with --merge: JSON {"unit": "reason"} for stopped units')
    a = ap.parse_args()
    if a.merge:
        merge(a.merge, a.practice, json.loads(a.not_finished))
        return
    plan = PRACTICE if a.practice else TEST
    if a.datasets:
        plan = [(k, s) for k, s in plan if k in a.datasets]
    role = "practice" if a.practice else "test"
    for k, _ in plan:
        if D.REGISTRY[k].role != role:
            sys.exit(f"{k} is a {D.REGISTRY[k].role} dataset; refusing to mix roles")
    pits = _cal_pits()
    t0 = time.time()
    rows = []
    for k, s in plan:
        rows += run_units(k, a.splits or s, a.methods, pits)
    summ = summarise(rows)
    losses = where_we_lose(summ)
    meta = {"utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "practice": a.practice,
            "seconds": time.time() - t0, "methods": a.methods, "plan": plan}
    stem = a.out or ("results_practice" if a.practice else "results")
    (HERE / f"{stem}.json").write_text(json.dumps({"meta": meta, "summary": summ, "where_we_lose": losses,
                                                   "not_run": M.NOT_RUN, "rows": rows}, indent=1, default=float))
    md = table_md(summ, "Benchmark: practice datasets (tuning allowed; in-sample)" if a.practice else
                  "Benchmark: held-out test datasets", losses, M.NOT_RUN, meta)
    md_name = f"{a.out}.md" if a.out else ("RESULTS_TABLE_practice.md" if a.practice else "RESULTS_TABLE.md")
    (HERE / md_name).write_text(md, encoding="utf-8")
    print(f"-> {HERE / (stem + '.json')} ({meta['seconds']:.0f} s); {len(losses)} 'where we lose' entries")


if __name__ == "__main__":
    main()
