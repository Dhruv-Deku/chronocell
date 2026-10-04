"""
02 4D dynamics · 04 Variant impact: what a structural variant changes in the fitted population
(contacts, genes, enhancer-promoter pairs), with 90 % intervals from refits on resampled counts.

The prediction is chronocell.perturb's covariance-space construction applied to a population model of
a window containing the variant (built in 01 Structure, or here). Its standing on real data is read
from validation/results_sv.json and shown beside it; no number on this panel is typed in.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
import streamlit as st

from chronocell import accuracy as ACC, genes as G, viz
from chronocell import theme as T
from ui.common import Dataset, banner, esc, fmt, html, readout, telemetry_row

FIT_P_ADJACENT = 0.3      # adjacent contact probability assumed for fits made here (the Gate 1 / SV-test setting)
FLANK_BEADS = 200         # beads on each side of the variant when a window is fitted here
MAX_FIT_BEADS = 2000      # largest window fitted here
MAX_IMPACT_BEADS = 1500   # beads used for the impact: the exact marginal of the fitted window around the variant
MIN_P = 0.01              # list a change only if the contact probability is at least this before or after
MIN_LOG2 = 1.0            # ... and it changes at least two-fold
BOOT_REPS = 8             # refits on Poisson-resampled counts for the 90 % intervals
CACHE_ENTRIES = 4


def extent(op: str, params: dict) -> tuple[int, int]:
    """Half-open bead span a variant touches (local bead coordinates)."""
    if op == "translocation":
        bp = int(params["breakpoint"])
        return bp - 1, bp + 1
    if "segments" in params:
        return min(a for a, _ in params["segments"]), max(b for _, b in params["segments"])
    return int(params["a"]), int(params["b"])


def covers(op: str, params: dict, lo: int, hi: int) -> bool:
    """Does window [lo, hi) hold the variant with at least one flanking bead where one is needed?"""
    a, b = extent(op, params)
    if op == "duplication":
        return lo <= a and b <= hi
    if op == "translocation":
        return lo + 2 <= int(params["breakpoint"]) <= hi - 1
    return lo < a and b < hi


def crop(op: str, params: dict, lo: int, hi: int) -> tuple[int, int] | None:
    """At most MAX_IMPACT_BEADS of [lo, hi) centred on the variant; None if the variant alone is larger."""
    if hi - lo <= MAX_IMPACT_BEADS:
        return lo, hi
    a, b = extent(op, params)
    if b - a + 2 > MAX_IMPACT_BEADS:
        return None
    room = (MAX_IMPACT_BEADS - (b - a)) // 2
    c_lo = max(lo, a - room)
    c_hi = min(hi, c_lo + MAX_IMPACT_BEADS)
    c_lo = max(lo, c_hi - MAX_IMPACT_BEADS)
    return c_lo, c_hi


def populations(ds: Dataset, frame: int) -> list[tuple[int, int, str, object]]:
    """Population models of this dataset and frame built this session: (lo, hi, key, result), largest first."""
    out = []
    for key, res in st.session_state.get("ensembles", {}).items():
        try:
            dkey, fr, lo, hi = key.rsplit(":", 3)
            if dkey == ds.key and int(fr) == int(frame):
                out.append((int(lo), int(hi), key, res))
        except ValueError:
            continue
    return sorted(out, key=lambda t: t[0] - t[1])


def _window_counts(ds: Dataset, lo: int, hi: int):
    m = (ds.ci >= lo) & (ds.ci < hi) & (ds.cj >= lo) & (ds.cj < hi)
    return ds.ci[m] - lo, ds.cj[m] - lo, ds.cm[m]


def _fit_counts(ds: Dataset, lo: int, hi: int, b0: float, p_adj: float, ci, cj, cm, progress=None):
    """The app's population model for a window: v3.3 up to POP.V33_MAX_BEADS beads, v4 above."""
    from chronocell import ensemble as ENS, population as POP
    n = hi - lo
    if n > POP.V33_MAX_BEADS:
        return POP.fit_population_from_counts(ci, cj, cm, n, ds.valid[lo:hi], b0_nm=b0, p_adjacent=p_adj,
                                              cfg=POP.config_for(n), progress=progress)
    return ENS.fit_from_counts(ci, cj, cm, n, ds.valid[lo:hi], b0_nm=b0, p_adjacent=p_adj, progress=progress)


def _status() -> None:
    ev = ACC.sv_evidence()
    if ev is None:
        banner("<b>Mechanism simulator, not validated.</b> The structural-variant test "
               "(<code>validation/sv_validation.py</code>) has not been run on this install.")
        return
    s, ci = ev["spearman"], ev["ci95"]
    txt = (f"Pre-registered test on real data ({esc(ev['event'])}; {esc(ev['window'])}; {ev['pairs']:,} pairs across "
           f"the deletions): Spearman with the measured contacts {s['model']:.3f} (95 % CI {ci['model'][0]:.3f} to "
           f"{ci['model'][1]:.3f}), against {s['distance_shift']:.3f} for a genomic-distance shift and "
           f"{s['no_change']:.3f} for 'no change'.")
    if ev["validated"]:
        banner(f"<b>Checked on one real event only.</b> {txt} Details: <code>validation/RESULTS.md</code> (Gate 4).", "info")
    else:
        banner(f"<b>Mechanism simulator, not validated.</b> {txt} It assumes the variant list fully describes how "
               "the kept pieces are joined. Details: <code>validation/RESULTS.md</code> (Gate 4).")


def _remember(key: str, value) -> None:
    cache = st.session_state.setdefault("vi_cache", {})
    cache[key] = value
    while len(cache) > CACHE_ENTRIES:
        cache.pop(next(iter(cache)))


def _fit_here(ds: Dataset, frame: int, b0: float, op: str, params: dict):
    """Offer (and on request run) a population fit of a window around the variant; returns
    (lo, hi, key, result) once fitted, else None."""
    a, b = extent(op, params)
    if b - a + 2 * 20 > MAX_FIT_BEADS:
        html(f'<p class="cc-note">This variant spans {b - a:,} beads; a window around it would exceed '
             f'{MAX_FIT_BEADS:,} beads. Use a coarser resolution, or fit the window in 01 Structure.</p>')
        return None
    flank = min(FLANK_BEADS, (MAX_FIT_BEADS - (b - a)) // 2)
    lo, hi = max(0, a - flank), min(ds.n, b + flank)
    ci, cj, cm = _window_counts(ds, lo, hi)
    ch = ds.chrom
    html(f'<p class="cc-note">No population model of this dataset covers the variant yet. Fit one for '
         f'{ch.name}:{float(ch.bin_start(ds.bin0 + lo)) / 1e6:.2f}–{float(ch.bin_end(ds.bin0 + hi - 1)) / 1e6:.2f} Mb '
         f'({hi - lo:,} beads, adjacent contact probability {FIT_P_ADJACENT} assumed, as in the held-out tests), '
         'or build one in 01 Structure.</p>')
    if len(cm) < 20:
        html('<p class="cc-note">Too few contacts in that window for a population model.</p>')
        return None
    if not st.button("Fit population around the variant", key="vi_fit", icon=":material/scatter_plot:", width="stretch"):
        return None
    stop = st.empty()
    stop.button("Stop", key="vi_fit_stop", icon=":material/stop_circle:",
                help="Stops the fit at its next progress update; nothing is saved.")
    bar = st.progress(0.0, text="Fitting the population…")
    t0 = time.time()

    def on_step(it: int, tot: int, row: dict) -> None:
        if it % 50 == 0 or it == tot:
            bar.progress(min(it / max(tot, 1), 1.0), text=f"Fitting the population · step {it}/{tot} · "
                                                         f"misfit {row['loss']:.4f} · {time.time() - t0:.0f} s")
    try:
        res = _fit_counts(ds, lo, hi, b0, FIT_P_ADJACENT, ci, cj, cm, on_step)
    except (ValueError, FloatingPointError) as exc:
        html(f'<p class="cc-note">The fit failed: {esc(str(exc))}</p>')
        return None
    finally:
        stop.empty()
        bar.empty()
    c_fit = ACC.spearman(res.contact_probability[ci, cj], cm)
    res.config["window_contact_fit"] = c_fit
    key = f"{ds.key}:{frame}:{lo}:{hi}"
    st.session_state.setdefault("ensembles", {})[key] = res
    model = "v4 population" if res.config.get("model") == "population_v4" else "v3.3 population"
    st.session_state.setdefault("telemetry", []).append(
        telemetry_row(model, ds.bin0 + lo, ds.bin0 + hi, hi - lo, "ensemble fit (variant impact)",
                      res.config["fit_seconds"], res.config.get("device_used", "cpu"), res.history["best_loss"][0], c_fit))
    return lo, hi, key, res


def _gene_rows(ds: Dataset, imp, c_lo: int, op: str, params: dict) -> pd.DataFrame:
    ch = ds.chrom
    res = ch.resolution
    n = len(imp.log2_fc)
    lo_bin = ds.bin0 + c_lo
    tab = G.in_region(ch.name, int(ch.bin_start(lo_bin)), int(ch.bin_end(lo_bin + n - 1)), ch.assembly)
    if tab.empty:
        return pd.DataFrame(columns=["Gene", "TSS (Mb)", "Effect", "Strongest contact change (log₂)", "with (Mb)"])
    a, b = (np.array(extent(op, params)) - c_lo)
    segs = [(x - c_lo, y - c_lo) for x, y in params.get("segments", [])] or [(a, b)]
    pmax = np.fmax(np.nan_to_num(imp.p_before, nan=0.0), np.nan_to_num(imp.p_after, nan=0.0))
    rows = []
    for _, g in tab.iterrows():
        gs, ge = int(g["start"]) // res - lo_bin, (int(g["end"]) - 1) // res - lo_bin
        t = int(g["tss"]) // res - lo_bin
        effect = ""
        if op == "deletion":
            hit = [(x, y) for x, y in segs if gs < y and ge >= x]
            if hit:
                effect = "deleted" if any(x <= gs and ge < y for x, y in hit) else "partly deleted"
        elif op == "duplication" and gs < b and ge >= a:
            effect = "duplicated" if a <= gs and ge < b else "partly duplicated"
        elif op == "inversion" and gs < b and ge >= a:
            effect = "inverted" if a <= gs and ge < b else "breakpoint inside the gene"
        elif op == "translocation":
            bp = int(params["breakpoint"]) - c_lo
            if gs < bp <= ge:
                effect = "breakpoint inside the gene"
            elif gs >= bp:
                effect = "beyond the breakpoint (not on this derivative chromosome)"
        strongest, partner, iv = np.nan, np.nan, ""
        if 0 <= t < n:
            row = np.where((np.abs(np.arange(n) - t) >= 2) & (pmax[t] >= MIN_P), imp.log2_fc[t], np.nan)
            if np.isfinite(row).any():
                k = int(np.nanargmax(np.abs(row)))
                strongest, partner = float(row[k]), (lo_bin + k) * res / 1e6
                if imp.interval is not None:
                    iv = f"{imp.interval['lo'][t, k]:+.2f} to {imp.interval['hi'][t, k]:+.2f}"
                if not effect and abs(strongest) >= MIN_LOG2:
                    effect = "contacts changed"
        if effect:
            rows.append({"Gene": g["name"], "TSS (Mb)": round(int(g["tss"]) / 1e6, 3), "Effect": effect,
                         "Strongest contact change (log₂)": None if not np.isfinite(strongest) else round(strongest, 2),
                         "with (Mb)": None if not np.isfinite(partner) else round(float(partner), 3),
                         **({"90 % interval": iv} if imp.interval is not None else {})})
    order = {"deleted": 0, "partly deleted": 1, "duplicated": 0, "partly duplicated": 1, "inverted": 0,
             "breakpoint inside the gene": 1, "beyond the breakpoint (not on this derivative chromosome)": 2,
             "contacts changed": 3}
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.assign(_o=out["Effect"].map(order),
                         _f=-out["Strongest contact change (log₂)"].abs().fillna(0)).sort_values(["_o", "_f"])
        out = out.drop(columns=["_o", "_f"]).reset_index(drop=True)
    return out


def _ep_rows(ds: Dataset, imp, c_lo: int) -> pd.DataFrame | str:
    from chronocell import perturb as PT
    if ds.signal_is_placeholder:
        return ("Enhancer–promoter pairs need a measured activity track (H3K27ac or ATAC, Data → Tracks); the "
                "current signal is the synthetic placeholder.")
    ch = ds.chrom
    res = ch.resolution
    n = len(imp.log2_fc)
    lo_bin = ds.bin0 + c_lo
    sig = np.asarray(ds.epi[c_lo:c_lo + n], dtype=float)
    if not np.isfinite(sig).any():
        return "The loaded activity track has no values in this window."
    tab = G.in_region(ch.name, int(ch.bin_start(lo_bin)), int(ch.bin_end(lo_bin + n - 1)), ch.assembly)
    promoters = [(str(g["name"]), int(g["tss"]) // res - lo_bin) for _, g in tab.iterrows()]
    promoters = [(name, t) for name, t in promoters if 0 <= t < n]
    thr = float(np.nanpercentile(sig, 90))
    p_beads = {t for _, t in promoters}
    enh = [int(k) for k in np.flatnonzero(np.nan_to_num(sig, nan=-np.inf) >= thr) if int(k) not in p_beads]
    pairs = PT.enhancer_promoter_pairs(imp, promoters, enh, min_p=MIN_P, min_abs_log2=MIN_LOG2)
    rows = []
    for r in pairs[:200]:
        row = {"Gene": r["gene"], "Promoter (Mb)": round((lo_bin + r["promoter_bead"]) * res / 1e6, 3),
               "Enhancer (Mb)": round((lo_bin + r["enhancer_bead"]) * res / 1e6, 3), "Change": r["change"],
               "log₂ change": None if not np.isfinite(r["log2_fc"]) else round(r["log2_fc"], 2),
               "P before": round(r["p_before"], 3), "P after": round(r["p_after"], 3)}
        if "log2_fc_90ci" in r:
            row["90 % interval"] = f"{r['log2_fc_90ci'][0]:+.2f} to {r['log2_fc_90ci'][1]:+.2f}"
        rows.append(row)
    return pd.DataFrame(rows, columns=["Gene", "Promoter (Mb)", "Enhancer (Mb)", "Change", "log₂ change", "P before",
                                       "P after"] + (["90 % interval"] if imp.interval is not None else []))


def render(ds: Dataset, frame: int, b0: float, op: str, params: dict, label: str) -> None:
    _status()
    try:
        from chronocell import perturb as PT
    except ImportError:
        html('<p class="cc-note">Needs PyTorch (<code>pip install torch</code>).</p>')
        return
    if not ds.has_contacts:
        html('<p class="cc-note">Needs contacts (Data → Graph): the impact is computed from a population model of '
             'the contact map.</p>')
        return
    file_dels = st.session_state.get("sc_file_deletions", {}).get(ds.key) or []
    if op == "deletion" and len(file_dels) > 1:
        if st.toggle(f"Apply all {len(file_dels)} deletions in the file on {ds.chrom.name} together", value=False,
                     key="vi_all_dels"):
            params = {"segments": [tuple(s) for s in file_dels]}
            label = f"{len(file_dels)} deletions from the variant file"
    if ds.is_reference:
        banner("Input is the <b>synthetic</b> reference contact map (planted model): these changes illustrate the "
               "mechanism only.", "info")
    cover = [(lo, hi, key, res) for lo, hi, key, res in populations(ds, frame) if covers(op, params, lo, hi)]
    if not cover:
        fitted = _fit_here(ds, frame, b0, op, params)
        if fitted is None:
            return
        cover = [fitted]
    w_lo, w_hi, pop_key, res = cover[0]
    cw = crop(op, params, w_lo, w_hi)
    if cw is None:
        html(f'<p class="cc-note">The variant spans more than {MAX_IMPACT_BEADS:,} beads; use a coarser resolution.</p>')
        return
    c_lo, c_hi = cw
    ch = ds.chrom
    is_v4 = res.config.get("model") == "population_v4"
    p_adj = float(res.config.get("p_adjacent_assumed", FIT_P_ADJACENT))
    html(f'<p class="cc-note">Population: {"v4" if is_v4 else "v3.3"} model of {ch.name}:'
         f'{float(ch.bin_start(ds.bin0 + w_lo)) / 1e6:.2f}–{float(ch.bin_end(ds.bin0 + w_hi - 1)) / 1e6:.2f} Mb '
         f'({w_hi - w_lo:,} beads; adjacent contact probability assumed {p_adj:g}'
         + ("" if abs(p_adj - FIT_P_ADJACENT) < 1e-9 else f", the held-out tests used {FIT_P_ADJACENT:g}")
         + f'). Impact computed on {c_hi - c_lo:,} beads around {esc(label)}.</p>')
    key = f"{pop_key}|{op}|{sorted((k, str(v)) for k, v in params.items())}|{c_lo}:{c_hi}"
    cache = st.session_state.setdefault("vi_cache", {})
    imp = cache.get(key)
    if imp is None:
        t0 = time.time()
        with st.spinner("Applying the variant to the population…"):
            try:
                imp = PT.variant_impact(res, op, PT.shift_params(params, w_lo), window=(c_lo - w_lo, c_hi - w_lo))
            except ValueError as exc:
                html(f'<p class="cc-note">{esc(str(exc))}</p>')
                return
        imp.notes.append(f"computed in {time.time() - t0:.1f} s")
        _remember(key, imp)

    # ---- uncertainty: refits on resampled counts --------------------------------------------------
    if imp.interval is None:
        if st.button(f"Estimate uncertainty ({BOOT_REPS} refits on resampled counts)", key="vi_boot",
                     icon=":material/query_stats:", width="stretch",
                     help="Refits the population on Poisson-resampled counts and re-applies the variant each time; "
                          "gives 90 % intervals of every change. Takes about as long as the fit times "
                          f"{BOOT_REPS}."):
            stop = st.empty()
            stop.button("Stop", key="vi_boot_stop", icon=":material/stop_circle:",
                        help="Stops at the next refit; nothing is saved.")
            bar = st.progress(0.0, text=f"Refit 0 / {BOOT_REPS}…")
            t0 = time.time()
            ci, cj, cm = _window_counts(ds, w_lo, w_hi)

            def fit_fn(a, b, c):
                return _fit_counts(ds, w_lo, w_hi, b0, p_adj, a, b, c)

            def on_rep(done: int, total: int) -> None:
                bar.progress(done / total, text=f"Refit {done} / {total} · {time.time() - t0:.0f} s")
            try:
                iv = PT.bootstrap_impact(fit_fn, (ci, cj, cm), op, PT.shift_params(params, w_lo), reps=BOOT_REPS,
                                         seed=0, level=0.9, progress=on_rep, window=(c_lo - w_lo, c_hi - w_lo))
            finally:
                stop.empty()
                bar.empty()
            imp.interval = iv
            st.session_state.setdefault("telemetry", []).append(
                telemetry_row("v4 population" if is_v4 else "v3.3 population", ds.bin0 + w_lo, ds.bin0 + w_hi,
                              w_hi - w_lo, f"variant-impact bootstrap ({BOOT_REPS} refits)", time.time() - t0,
                              res.config.get("device_used", "cpu"), None, None))

    # ---- read-outs ---------------------------------------------------------------------------------
    fc = imp.log2_fc
    n = len(fc)
    lo_bin = ds.bin0 + c_lo
    res_bp = ch.resolution
    pmax = np.fmax(np.nan_to_num(imp.p_before, nan=0.0), np.nan_to_num(imp.p_after, nan=0.0))
    iu = np.triu_indices(n, 2)
    v = fc[iu]
    ok = np.isfinite(v) & (pmax[iu] >= MIN_P)
    absent = int((~np.isfinite(imp.p_after).any(axis=1)).sum())
    rows = [("Pairs gaining ≥ 2×", f"{int((v[ok] >= MIN_LOG2).sum()):,}", f"of {int(ok.sum()):,} with P ≥ {MIN_P:g}"),
            ("Pairs losing ≥ 2×", f"{int((v[ok] <= -MIN_LOG2).sum()):,}", ""),
            ("Beads absent after the variant", f"{absent:,}", "deleted, or moved off this chromosome")]
    if imp.interval is not None:
        lo_i, hi_i = imp.interval["lo"][iu], imp.interval["hi"][iu]
        sure = ok & (np.abs(v) >= MIN_LOG2) & ((lo_i > 0) | (hi_i < 0))
        rows.append(("Of those, 90 % interval excludes no change", f"{int(sure.sum()):,}",
                     f"{imp.interval['reps']} refits"))
    else:
        rows.append(("Uncertainty", "not estimated", "use the button above"))
    readout(rows)
    if op == "translocation":
        html('<p class="cc-note">Contacts with the partner segment are not in this chromosome\'s coordinates; the '
             'partner is modelled as a homogeneous chain (no partner data), so only losses appear here.</p>')

    t_map, t_pairs, t_genes, t_ep = st.tabs(["Map", "Changed contacts", "Genes", "Enhancer–promoter"])
    with t_map:
        st.plotly_chart(viz.fold_change_map(fc, lo_bin, res_bp), theme=None, width="stretch", config=T.PLOT_CONFIG,
                        key="vi_map")
        html('<p class="cc-note">log₂ (contact probability after / before), averaged into pixels; blank where a bead is '
             'absent after the variant.</p>')
    with t_pairs:
        order = np.argsort(-np.abs(np.where(ok, v, 0.0)))[:25]
        prow = []
        for t in order:
            if not ok[t]:
                continue
            i, j = int(iu[0][t]), int(iu[1][t])
            r = {"Locus A (Mb)": round((lo_bin + i) * res_bp / 1e6, 3), "Locus B (Mb)": round((lo_bin + j) * res_bp / 1e6, 3),
                 "log₂ change": round(float(fc[i, j]), 2), "P before": round(float(imp.p_before[i, j]), 3),
                 "P after": round(float(imp.p_after[i, j]), 3)}
            if imp.interval is not None:
                r["90 % interval"] = f"{imp.interval['lo'][i, j]:+.2f} to {imp.interval['hi'][i, j]:+.2f}"
            prow.append(r)
        pairs_df = pd.DataFrame(prow)
        st.dataframe(pairs_df, hide_index=True, width="stretch", height=min(38 + 35 * max(len(prow), 1), 320),
                     key="vi_pairs")
    with t_genes:
        genes_df = _gene_rows(ds, imp, c_lo, op, params)
        html(f'<p class="cc-note">{len(genes_df):,} genes affected ({G.source_note(ch.assembly)}): inside the variant, '
             f'or with a contact changing at least two-fold (P ≥ {MIN_P:g}).</p>')
        st.dataframe(genes_df, hide_index=True, width="stretch", height=min(38 + 35 * max(len(genes_df), 1), 320),
                     key="vi_genes")
    with t_ep:
        ep = _ep_rows(ds, imp, c_lo)
        if isinstance(ep, str):
            html(f'<p class="cc-note">{esc(ep)}</p>')
            ep_df = pd.DataFrame()
        else:
            ep_df = ep
            html(f'<p class="cc-note">{len(ep_df):,} pairs. Candidate enhancers: beads in the top 10 % of the loaded '
                 f'{esc(ds.signal_label)} signal in this window (an activity proxy, not an enhancer annotation); '
                 'promoters: RefSeq TSS beads.</p>')
            st.dataframe(ep_df, hide_index=True, width="stretch", height=min(38 + 35 * max(len(ep_df), 1), 320),
                         key="vi_ep")
    stem = f"{ch.name}_{op}_impact"
    c1, c2, c3 = st.columns(3)
    c1.download_button("Contacts (CSV)", pairs_df.to_csv(index=False), f"{stem}_contacts.csv", "text/csv",
                       width="stretch", icon=":material/download:", key="vi_dl_pairs")
    c2.download_button("Genes (CSV)", genes_df.to_csv(index=False), f"{stem}_genes.csv", "text/csv",
                       width="stretch", icon=":material/download:", key="vi_dl_genes")
    c3.download_button("E–P pairs (CSV)", ep_df.to_csv(index=False), f"{stem}_ep_pairs.csv", "text/csv",
                       width="stretch", icon=":material/download:", key="vi_dl_ep", disabled=ep_df.empty)
