"""v4 Pillar 8: the reproducibility record (chronocell_audit_log.json and PDF), input hashes, the method
evidence block (measured values only), the REST API's whole-window model, and the Export controls."""

from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pytest

from chronocell import accuracy as ACC, provenance as PROV


def _record(**over):
    kw = dict(software="4.0", method="maximum-entropy population model", parameters={"b0_nm": 50.0, "alpha": 3.0},
              inputs=[{"role": "structure", "name": "x.npy", "sha256": hashlib.sha256(b"abc").hexdigest()}],
              metrics={"rg_nm": 412.5, "nu": float("nan")}, accuracy=ACC.two_scores("ensemble_v3_3", 0.9),
              telemetry=[{"Stage": "ensemble fit", "Time (s)": 1.2}], genome={"assembly": "GRCh38/hg38", "chrom": "chr22"})
    kw.update(over)
    return PROV.build(**kw)


def test_record_has_statement_inputs_equations_and_sources():
    rec = _record()
    assert "not a clinical, regulatory or certified audit" in rec["statement"]
    assert rec["inputs"][0]["sha256"] == hashlib.sha256(b"abc").hexdigest()
    assert "population_objective" in rec["equations"] and rec["software"]["environment"]["numpy"] == np.__version__
    assert rec["measured_metrics"]["nu"] is None                       # NaN is never written as a number
    text = PROV.to_json(rec)                                           # strict JSON (allow_nan=False)
    assert json.loads(text)["record_sha256"] == rec["record_sha256"] and len(rec["record_sha256"]) == 64
    srcs = rec["validation_datasets"]
    assert srcs and all(s.get("license") and s.get("citation") for s in srcs)
    assert _record(inputs=[{"role": "structure", "name": "x.npy", "sha256": "0" * 64}])["record_sha256"] != rec["record_sha256"]


def test_record_pdf_renders():
    pdf = PROV.to_pdf(_record())
    assert pdf[:4] == b"%PDF" and len(pdf) > 3000


def test_method_evidence_reports_measured_values_or_not_run():
    ev = ACC.method_evidence()
    for key in ("v3_3_microscopy_benchmark", "gate1_whole_chromosome", "gate2_calibration", "gate4_cohesin_depletion"):
        assert key in ev
        assert ev[key] == "not run" or isinstance(ev[key], dict)
    if isinstance(ev["gate1_whole_chromosome"], dict):                 # the stored test result, not a constant
        stored = json.loads((ACC.VALIDATION / "results_gate1.json").read_text())
        v = stored["datasets"][0]["summary"]["imaging_input"]["all_pairs"]["whole_v4"]["percent_of_ceiling"]
        assert ev["gate1_whole_chromosome"]["imaging_input.all_pairs.whole_v4"] == round(v, 1)


def test_dataset_records_input_hashes():
    from ui.common import load_dataset
    from chronocell import synthetic, genome
    ref = synthetic.build(genome.chrom("chr22"), seed=7)
    buf = io.BytesIO()
    np.save(buf, ref.coords[:300].astype(np.float32))
    data = buf.getvalue()
    ds = load_dataset("chr22", 7, None, ("upload", "coords.npy", data), "nm", None, False)
    assert ds.inputs == (("structure", "coords.npy", hashlib.sha256(data).hexdigest()),)


def test_api_whole_window_population_model(tmp_path):
    from chronocell import api
    from tests.test_v33 import _population
    d = _population(n=420, cells=3000, seed=2)
    f = (d < 1.0).mean(0)
    i, j = np.triu_indices(420, 1)
    cnt = np.random.default_rng(0).poisson(f[i, j] * 300).astype(float)
    keep = cnt > 0
    contacts = {"i": i[keep].tolist(), "j": j[keep].tolist(), "count": cnt[keep].tolist()}
    log = api.AuditLog(tmp_path / "runs.jsonl")
    with pytest.raises(api.RequestError):                               # v3.3 model keeps its 400-bead limit
        api.reconstruct({"contacts": contacts, "n_beads": 420, "model": "population"}, log=log)
    out = api.reconstruct({"contacts": contacts, "n_beads": 420, "model": "population_v4"}, log=log)
    assert out["model"] == "population_v4" and np.asarray(out["coords_nm"]).shape == (420, 3)
    assert out["accuracy"]["contact_map_fit"]["value"] > 0.5
    with pytest.raises(api.RequestError):
        api.reconstruct({"contacts": {"i": [0], "j": [1], "count": [1]}, "n_beads": 7000, "model": "population_v4"}, log=log)


def test_export_offers_the_audit_record(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    import ui.common as C
    import ui.states_panel as SP
    monkeypatch.setattr(C, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "DEMO_ROOT", tmp_path / "demo")
    at = AppTest.from_file(str(Path(__file__).resolve().parent.parent / "app.py"), default_timeout=900)
    at.run()
    assert not at.exception
    labels = [d.label for d in at.get("download_button")]
    assert "Audit log (JSON)" in labels
    at.button(key="audit_pdf_build").click().run()
    assert not at.exception, [e.value for e in at.exception]
    assert "Audit record (PDF)" in [d.label for d in at.get("download_button")]
