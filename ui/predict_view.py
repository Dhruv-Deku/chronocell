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

from chronocell import accuracy as ACC, genome, predict as PD
from ui.common import CACHE_DIR, Dataset, banner, esc, html, telemetry_row

FASTA_DIR = CACHE_DIR / "fasta"
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


def _peaks(ch: genome.Chrom) -> tuple[dict | None, str, str]:
    """CTCF peaks of this chromosome from the upload or pasted lines: (peaks, source name, SHA-256)."""
    up = st.file_uploader("CTCF ChIP-seq peaks of the same cell type (narrowPeak or BED; .gz accepted)",
                          type=["bed", "narrowpeak", "gz", "txt", "tsv"], key="pred_peaks")
    text = st.text_area("…or paste peak lines (chrom, start, end[, …, summit])", key="pred_peaks_text", height=70)
    if up is not None:
        data, name = up.getvalue(), up.name
    elif text.strip():
        data, name = text.encode(), "pasted peaks"
    else:
        return None, "", ""
    pk = PD.read_peaks(data, ch.name, lambda c: genome.normalize_chrom(c, ch.assembly))
    return pk, name, hashlib.sha256(data).hexdigest()


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
    if ch.assembly != meta.get("assembly"):
        html(f'<p class="cc-note">The predictor was trained and tested on human {esc(meta.get("assembly", ""))} only; it is '
             f'not offered for {esc(ch.genome.display)}.</p>')
        return
    from chronocell import population as POP
    n = hi - lo
    if n > POP.MAX_BEADS:
        html(f'<p class="cc-note">This window has {n:,} beads; the population model handles up to {POP.MAX_BEADS:,}.</p>')
        return
    pk, pk_name, pk_sha = _peaks(ch)
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
                       "prediction_inputs": {**info, "peaks_file": pk_name, "peaks_sha256": pk_sha,
                                             "sequence": f"UCSC {ch.assembly} {ch.name}"},
                       "anchor_b0_nm": float(b0)})
    ss.setdefault("ensembles", {})[fit_key] = res
    model_name = "predicted · " + ("v4 population" if res.config.get("model") == "population_v4" else "v3.3 population")
    ss.setdefault("telemetry", []).append(
        telemetry_row("sequence + CTCF predictor", ds.bin0 + lo, ds.bin0 + hi, n, "distance prediction (Gate 5 model)",
                      t_pred, "cpu", None, None))
    ss.telemetry.append(telemetry_row(model_name, ds.bin0 + lo, ds.bin0 + hi, n, "population fit to the prediction",
                                      res.config["fit_seconds"], res.config.get("device_used", "cpu"),
                                      res.history["best_loss"][0], None, res.spread_cv))
    ss.show_fit = True
    ss[f"structure_{fit_key}"] = "Population model"
    st.rerun()
