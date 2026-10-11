"""
01 3D structure · 03 Model & convergence · population model from sequence + CTCF (Pillar 5).

A prediction for windows with no contact data, made with the frozen Gate 5 model
(chronocell/data/predictor.json): CTCF ChIP-seq peaks of the cell type, oriented by the JASPAR CTCF
motif, and GC per bead give each pair's median distance; a population model is then fitted to that
map (population.fit_population_from_medians) so every page can use it. It is labelled "predicted"
wherever it appears, and its held-out record (Gate 5) is read from the result files. When contact
data exist the data-driven model is the right choice; blending the two was not tested.
"""

from __future__ import annotations

import hashlib
import time

import numpy as np
import streamlit as st

from chronocell import accuracy as ACC, genome, predict as PD, theme as T, viz
from ui.common import CACHE_DIR, Dataset, banner, esc, html, readout, telemetry_row

FASTA_DIR = CACHE_DIR / "fasta"
ENCODE_DIR = CACHE_DIR / "encode"
INPUT_LABEL = "predicted from sequence + CTCF (no contact data)"


def is_predicted(res) -> bool:
    return str(getattr(res, "config", {}).get("input", "")).startswith("predicted")


def _fasta_path(ch: genome.Chrom):
    p = FASTA_DIR / f"{ch.assembly}_{ch.name}.fa.gz"
    return p if p.exists() else None


@st.cache_resource(show_spinner="Reading the chromosome sequence…", max_entries=2)
def _sequence_cached(path: str, mtime: float) -> bytes:
    return PD.read_fasta(path)


def _sequence(ch: genome.Chrom) -> bytes | None:
    """The chromosome's sequence (upper-case bytes), from the local cache; None if not downloaded yet."""
    p = _fasta_path(ch)
    return None if p is None else _sequence_cached(str(p), p.stat().st_mtime)


def _status(meta: dict) -> None:
    ev = ACC.load_benchmark()
    pred = (ev or {}).get("models", {}).get("predicted_sequence_ctcf")
    if pred is None:
        banner("<b>Predicted, not measured.</b> The no-contact-data predictor's held-out test (Gate 5) has not been "
               "run on this install.")
        return
    per = " · ".join(f"{esc(k)} {v:.0f} %" for k, v in pred["per_dataset_percent_of_ceiling"].items())
    with_data = []
    for name, label in (("ensemble_v3_3", "v3.3"), ("population_v4", "v4")):
        m = (ev or {}).get("models", {}).get(name)
        if m and m.get("overall_percent_of_ceiling") is not None:
            with_data.append(f"{label} {m['overall_percent_of_ceiling']:.0f} %")
    banner(f"<b>Predicted from sequence + CTCF, with no contact data.</b> Held-out test (Gate 5, verdict: "
           f"{esc(pred['verdict'])}): it recovers {pred['overall_percent_of_ceiling']:.1f} % of the reproducible distance "
           f"pattern beyond the separation trend ({per})"
           + (f"; population models built from contact data: {', '.join(with_data)} (headline benchmarks)" if with_data else "")
           + ". Details: validation/RESULTS.md. Read it as a prior, not a measurement.", "info")
    caveat = control_caveat(pred)
    if caveat:
        html(f'<p class="cc-note">{caveat}</p>')


def control_caveat(pred: dict | None) -> str | None:
    """What the cohesin-depletion control says about the predictor's signal, from the Gate 5 result file."""
    ctl = (pred or {}).get("control_percent_of_ceiling") or {}
    if not ctl:
        return None
    tests = list(pred["per_dataset_percent_of_ceiling"].values())
    vals = " · ".join(f"{esc(k)} {v:.1f} %" for k, v in ctl.items())
    as_high = max(ctl.values()) >= float(pred.get("overall_percent_of_ceiling") or max(tests))
    return (f"<b>Control, cohesin depleted</b> (RAD21 degraded, so CTCF-anchored loops are gone): {vals}, against "
            f"{min(tests):.1f}–{max(tests):.1f} % on untreated test cells (overall "
            f"{pred.get('overall_percent_of_ceiling', float('nan')):.1f} %). "
            + ("It scores as high without loops, so its signal is compartment / insulation level (GC and CTCF density "
               "between loci), <b>not CTCF loops</b>. Do not read loops from it." if as_high else
               "It scores lower without loops, as expected if part of its signal is CTCF loops."))


PEAK_SOURCES = ("Upload or paste", "ENCODE, by cell type")


def _encode_peaks() -> tuple[bytes | None, str, dict]:
    """CTCF IDR peaks of one ENCODE cell line (data/validation_sources.json), downloaded once into the cache
    and checked against the portal's MD5: (file bytes, source name, record of the accession)."""
    srcs = [s for s in PD.encode_ctcf_sources() if s["assembly"] in ("GRCh38", "hg38")]
    if not srcs:
        html('<p class="cc-note">No ENCODE CTCF sources are listed in chronocell/data/validation_sources.json.</p>')
        return None, "", {}
    by_cell = {s["cell_line"]: s for s in srcs}
    cell = st.selectbox("Cell type (ENCODE CTCF ChIP-seq, GRCh38 IDR thresholded peaks)", list(by_cell), key="pred_encode_cell",
                        help="Use the cell type you want to model: CTCF binding differs between cell types.")
    s = by_cell[cell]
    path = ENCODE_DIR / f"{s['accession']}.bed.gz"
    html(f'<p class="cc-note">ENCODE file <a href="{esc(s["page"])}" target="_blank">{esc(s["accession"])}</a> '
         f'(experiment {esc(s["experiment"])}). {esc(s["citation"])} Downloaded on demand into '
         '<code>.chronocell_cache/encode/</code>, checked against the MD5 the portal publishes, never redistributed.</p>')
    if not path.exists():
        if not st.button(f"Download the {cell} CTCF peaks from ENCODE", key="pred_encode_dl", icon=":material/download:"):
            return None, "", {}
        bar = st.progress(0.0, text="Downloading…")

        def on_bytes(done: int, total: int | None) -> None:
            if total:
                bar.progress(min(done / total, 1.0), text=f"Downloading {done / 1e6:,.2f} of {total / 1e6:,.2f} MB")
        try:
            PD.fetch_encode_peaks(s["accession"], s["md5"], ENCODE_DIR, on_bytes)
        except OSError as exc:
            bar.empty()
            html(f'<p class="cc-note">Download failed: {esc(str(exc))}</p>')
            return None, "", {}
        bar.empty()
    meta = {"encode_accession": s["accession"], "encode_experiment": s["experiment"], "encode_cell_line": cell,
            "encode_md5": s["md5"], "encode_citation": s["citation"]}
    return path.read_bytes(), f"ENCODE {s['accession']} ({cell} CTCF)", meta


def _peaks(ch: genome.Chrom) -> tuple[dict | None, str, str, dict]:
    """CTCF peaks of this chromosome from ENCODE, an upload or pasted lines: (peaks, source name, SHA-256,
    extra record fields)."""
    options = PEAK_SOURCES if ch.assembly == "hg38" else PEAK_SOURCES[:1]      # the ENCODE list is human (GRCh38) only
    src = st.segmented_control("CTCF peaks", options, default=options[0], required=True, key="pred_peaks_src")
    meta: dict = {}
    if src == PEAK_SOURCES[1]:
        data, name, meta = _encode_peaks()
        if data is None:
            return None, "", "", {}
    else:
        up = st.file_uploader("CTCF ChIP-seq peaks of the same cell type (narrowPeak or BED; .gz accepted)",
                              type=["bed", "narrowpeak", "gz", "txt", "tsv"], key="pred_peaks")
        text = st.text_area("…or paste peak lines (chrom, start, end[, …, summit])", key="pred_peaks_text", height=70)
        if up is not None:
            data, name = up.getvalue(), up.name
        elif text.strip():
            data, name = text.encode(), "pasted peaks"
        else:
            return None, "", "", {}
    pk = PD.read_peaks(data, ch.name, lambda c: genome.normalize_chrom(c, ch.assembly))
    return pk, name, hashlib.sha256(data).hexdigest(), meta


def render(ds: Dataset, lo: int, hi: int, fit_key: str, b0: float) -> None:
    """Prediction input for the population model; stores the fitted result in st.session_state.ensembles."""
    ss = st.session_state
    loaded = PD.load_model()
    if loaded is None:
        html('<p class="cc-note">The predictor file (chronocell/data/predictor.json) is missing.</p>')
        return
    model, meta = loaded
    _status(meta)
    ch = ds.chrom
    if not PD.assembly_supported(ch.assembly, meta):
        html(f'<p class="cc-note">The predictor was trained and tested on human {esc(meta.get("assembly", ""))} only; it is '
             f'not offered for {esc(ch.genome.display)} (validation/RESULTS.md, Gates 5 and 5m).</p>')
        return
    if ch.assembly != meta.get("assembly"):
        from chronocell import accuracy as ACC
        ev = ACC.mouse_predictor_evidence()
        loci = "; ".join(f"{v['region']}: {v['percent_of_ceiling']:.1f} % of the reproducible pattern (95 % interval "
                         f"{v['ci95'][0]:.1f} to {v['ci95'][1]:.1f})" for v in ev["loci"].values())
        banner(f"<b>Mouse: the human model, unchanged.</b> Pre-registered mouse test (Gate 5m, ORCA tracing in mouse ES "
               f"cells, {esc(ev['assembly_tested'])}): passed on both loci, {esc(loci)}. Raw gains over the separation "
               "trend are small: a prior, not a substitute for contacts. Give mouse CTCF peaks of your cells (upload "
               "or paste); the ENCODE list is human only.", "info")
    from chronocell import population as POP
    n = hi - lo
    if n > POP.MAX_BEADS:
        html(f'<p class="cc-note">This window has {n:,} beads; the population model handles up to {POP.MAX_BEADS:,}.</p>')
        return
    pk, pk_name, pk_sha, pk_meta = _peaks(ch)
    if pk is not None:
        html(f'<p class="cc-note">{len(pk["starts"]):,} peaks on {ch.name} read from {esc(pk_name)}.</p>')
    seq = _sequence(ch)
    if seq is None:
        html(f'<p class="cc-note">GC content and the CTCF motif orientation need the {ch.name} sequence '
             f'({esc(ch.genome.display)}): downloaded once from UCSC into <code>.chronocell_cache/fasta/</code> and '
             'checked against UCSC\'s MD5.</p>')
        if st.button(f"Download the {ch.name} sequence from UCSC", key="pred_fasta_dl", icon=":material/download:"):
            bar = st.progress(0.0, text="Downloading…")

            def on_bytes(done: int, total: int | None) -> None:
                if total:
                    bar.progress(min(done / total, 1.0), text=f"Downloading {done / 1e6:,.1f} of {total / 1e6:,.1f} MB")
            try:
                PD.fetch_chromosome_fasta(ch.assembly, ch.name, FASTA_DIR, on_bytes)
            except OSError as exc:
                html(f'<p class="cc-note">Download failed: {esc(str(exc))}</p>')
                return
            bar.empty()
            seq = _sequence(ch)
    ready = seq is not None and pk is not None and len(pk["starts"]) > 0
    if not st.button("Predict and build population model", key="pred_go", width="stretch", disabled=not ready,
                     help="Needs the peaks and the sequence. Uses the frozen Gate 5 model; nothing is fitted to your data."):
        return
    stop = st.empty()
    stop.button("Stop", key="pred_stop", icon=":material/stop_circle:",
                help="Stops at the next progress update; nothing is saved.")
    bar = st.progress(0.0, text="Predicting distances from sequence + CTCF…")
    t0 = time.time()
    g = ds.gbin(np.arange(lo, hi))
    d, info = PD.predict_window(ch.bin_start(g), ch.bin_end(g), seq, pk, model, meta["settings"])
    t_pred = time.time() - t0
    t1 = time.time()

    def on_step(it: int, tot: int, row: dict) -> None:
        if it % 50 == 0 or it == tot:
            bar.progress(min(it / max(tot, 1), 1.0), text=f"Fitting the population to the prediction · step {it}/{tot} · "
                                                         f"misfit {row['loss']:.4f} · {time.time() - t1:.0f} s")
    try:
        res = POP.fit_population_from_medians(d, progress=on_step)
    except (ValueError, FloatingPointError) as exc:
        html(f'<p class="cc-note">The fit failed: {esc(str(exc))}</p>')
        return
    finally:
        stop.empty()
        bar.empty()
    res.config.update({"input": INPUT_LABEL, "window_contact_fit": None,
                       "predictor": {"model_file": "chronocell/data/predictor.json", "trained_on": model.trained_on,
                                     "motif": meta["settings"]["motif"], "ridge_lambda": model.ridge_lambda,
                                     "validation": meta.get("validation")},
                       "prediction_inputs": {**info, "peaks_file": pk_name, "peaks_sha256": pk_sha, **pk_meta,
                                             "sequence": f"UCSC {ch.assembly} {ch.name}"},
                       "anchor_b0_nm": float(b0)})
    ss.setdefault("ensembles", {})[fit_key] = res
    ss.setdefault("predicted_ensembles", {})[fit_key] = res        # kept for the side-by-side view
    model_name = "predicted · " + ("v4 population" if res.config.get("model") == "population_v4" else "v3.3 population")
    ss.setdefault("telemetry", []).append(
        telemetry_row("sequence + CTCF predictor", ds.bin0 + lo, ds.bin0 + hi, n, "distance prediction (Gate 5 model)",
                      t_pred, "cpu", None, None))
    ss.telemetry.append(telemetry_row(model_name, ds.bin0 + lo, ds.bin0 + hi, n, "population fit to the prediction",
                                      res.config["fit_seconds"], res.config.get("device_used", "cpu"),
                                      res.history["best_loss"][0], None, res.spread_cv))
    ss.show_fit = True
    ss.structure_pending = (fit_key, "Population model")     # applied before the Structure switch is drawn
    st.rerun()


def _block_mean(mat: np.ndarray, max_px: int = 300) -> tuple[np.ndarray, int]:
    """k x k block means of a square map, so that it has at most max_px pixels a side."""
    n = len(mat)
    k = max(1, int(np.ceil(n / max_px)))
    m = int(np.ceil(n / k))
    pad = np.full((m * k, m * k), np.nan)
    pad[:n, :n] = mat
    with np.errstate(invalid="ignore"):
        return np.nanmean(pad.reshape(m, k, m, k), axis=(1, 3)), k


def render_comparison(ds: Dataset, lo: int, hi: int, fit_key: str) -> None:
    """Side by side, when this window has both: the population built from the measured contacts and the one
    predicted from sequence + CTCF, with their rank agreement. Shown, never blended."""
    ss = st.session_state
    pred = ss.get("predicted_ensembles", {}).get(fit_key)
    meas = ss.get("contact_ensembles", {}).get(fit_key)
    if pred is None or meas is None or len(pred.median_distance_nm) != hi - lo or len(meas.median_distance_nm) != hi - lo:
        return
    html('<p class="cc-eyebrow" style="margin-top:14px">Predicted vs built from measured contacts · this window</p>')
    cmp = PD.compare_maps(meas.median_distance_nm, pred.median_distance_nm)
    ss.setdefault("prediction_comparison", {})[fit_key] = cmp
    g_lo = ds.gbin(lo)
    left, right = st.columns(2)
    for col, res, title, key in ((left, meas, "From measured contacts (population model)", "cmp_measured"),
                                 (right, pred, "Predicted · sequence + CTCF, no contacts used", "cmp_predicted")):
        mat, k = _block_mean(np.asarray(res.median_distance_nm, dtype=np.float64))
        col.markdown(f'<p class="cc-note"><b>{esc(title)}</b> · median distance (nm)</p>', unsafe_allow_html=True)
        col.plotly_chart(viz.matrix_chart(mat, g_lo, k, "distance", height=280, resolution=ds.chrom.resolution), theme=None,
                         width="stretch", config=T.PLOT_CONFIG, key=key)
    readout([
        ("Rank agreement, predicted vs from contacts<small>Spearman ρ over all pairs"
         + (" (200,000 sampled)" if cmp["sampled"] else "") + "; mostly the shared separation trend</small>",
         f"{cmp['spearman']:.3f}", "ρ"),
        ("Agreement beyond the separation trend<small>each map divided by its mean at the same separation</small>",
         f"{cmp['spearman_trend_removed']:.3f}", "ρ"),
        ("Size, predicted / from contacts<small>median ratio over pairs</small>", f"{cmp['median_ratio']:.2f}", "×"),
    ])
    html('<p class="cc-note">Agreement between two models, not accuracy: neither map is a distance measurement on this '
         'window. The two are shown side by side and never blended: combining a prediction with contact data was not '
         'tested (validation/RESULTS.md, Gate 5). Where both exist, the model built from contacts is the validated one.</p>')
    if ds.is_reference:
        html('<p class="cc-note"><b>These contacts are the synthetic reference model</b>, not measurements: here the '
             'agreement only compares a planted toy fold with a prediction from real sequence and says nothing about real '
             'folding.</p>')
