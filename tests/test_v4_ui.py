"""v4 UI: every new control, end to end under Streamlit's AppTest. Each test checks that the default
view is what it was before v4 and that the new control, once used, shows computed values."""

from __future__ import annotations

from pathlib import Path

import pytest

APP = str(Path(__file__).resolve().parent.parent / "app.py")


@pytest.fixture()
def app(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    import ui.common as C
    import ui.states_panel as SP
    monkeypatch.setattr(C, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "DEMO_ROOT", tmp_path / "demo")
    at = AppTest.from_file(APP, default_timeout=900)
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at


def _text(at) -> str:
    return "".join(m.value for m in at.markdown)


def _ok(at) -> bool:
    return not at.exception and "cc-card-warn" not in _text(at)


def test_compare_keeps_its_page_and_adds_the_state_evaluator_tab(app):
    at = app
    at.segmented_control(key="workspace").set_value("Compare").run()
    assert not at.exception, [e.value for e in at.exception]
    labels = [t.label for t in at.tabs]
    assert labels[:2] == ["Side by side", "Self-Math PDB State Evaluator"]
    # nothing is computed before Evaluate is pressed
    assert "pe_last" not in at.session_state
    at.button(key="pe_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    ev = at.session_state["pe_last"]["ev"]
    assert ev.state in ("Normal", "Diseased", "Senescent", "Indeterminate")
    assert ev.rg_nm > 0 and ev.n > 1000 and ev.hull_volume_um3 > 0
    txt = _text(at)
    assert "Radius of gyration" in txt and "Asphericity" in txt and "not a diagnosis" in txt
    assert any(d.label == "Evaluation (JSON)" for d in at.get("download_button"))
    # a window of the structure
    at.toggle(key="pe_whole").set_value(False).run()
    at.button(key="pe_go").click().run()
    assert not at.exception and at.session_state["pe_last"]["ev"].n < ev.n


def test_population_probe_interval_overlay_and_telemetry(app):
    at = app
    at.session_state["custom_window"] = (2000, 2150)
    at.segmented_control(key="region_choice").set_value("custom").run()
    assert _ok(at), [e.value for e in at.exception]
    next(b for b in at.button if (b.label or "").startswith("Build population model")).click().run()
    assert _ok(at), [e.value for e in at.exception]
    row = next(r for r in at.session_state["telemetry"] if r["Model"] == "v3.3 population")
    for col in ("Contact-map fit", "Microscopy accuracy", "Ensemble consistency (CV)", "Energy (final objective)",
                "Time (s)", "Stage"):
        assert col in row
    assert row["Microscopy accuracy"] is None                      # never filled in for a user's window
    # default view unchanged: no overlay
    assert "population spread" not in _text(at)
    at.toggle(key="probe_on").set_value(True).run()
    txt = _text(at)
    assert "Middle 50 % of cells (model)" in txt and "Mean ± SD across cells (model)" in txt
    tip = next(d for d in _chart(at, "viewport")[1]["data"] if d["type"] == "mesh3d")["hovertemplate"]
    assert "(population model)" in tip and "middle 50 %" in tip            # hover: distance to bead A with interval
    at.select_slider(key="probe_level").set_value("90 %").run()
    assert "Middle 90 % of cells (model)" in _text(at)
    at.toggle(key="unc_on").set_value(True).run()
    assert _ok(at) and "population spread (RMSF" in _text(at)


def test_whole_window_population_uses_v4_above_400_beads(app):
    at = app
    at.session_state["custom_window"] = (2000, 2450)              # 450 beads: beyond the v3.3 window limit
    at.segmented_control(key="region_choice").set_value("custom").run()
    assert _ok(at), [e.value for e in at.exception]
    assert "Population model (v4)" in _text(at)
    btn = next(b for b in at.button if (b.label or "").startswith("Build whole-window population model"))
    btn.click().run()
    assert _ok(at), [e.value for e in at.exception]
    ens = list(at.session_state["ensembles"].values())[0]
    assert ens.config["model"] == "population_v4" and ens.representative_nm.shape == (450, 3)
    assert any(r["Model"] == "v4 population" for r in at.session_state["telemetry"])


VCF_TEXT = ("chr22\t20000000\td1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20100000\n"
            "chr9\t100\tx\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000\n")


def _table(at, column: str):
    """The first dataframe on the page with this column (AppTest does not index dataframes by key)."""
    return next((d.value for d in at.dataframe if column in d.value.columns), None)


def test_variant_file_scenario_and_impact_panel(app, monkeypatch):
    import ui.variant_impact as VI
    monkeypatch.setattr(VI, "BOOT_REPS", 2)                        # two refits keep the test short
    at = app
    at.segmented_control(key="workspace").set_value("4D dynamics").run()
    assert _ok(at), [e.value for e in at.exception]
    sel = at.selectbox(key="sc_choice_chr22")
    assert sel.value == "22q11del"                                 # default scenario unchanged
    assert "From a variant file (VCF / BEDPE)" in sel.options
    assert "vi_cache" not in at.session_state                      # nothing computed before it is asked for
    sel.set_value("file").run()
    at.text_area(key="sc_vtext").input(VCF_TEXT).run()
    assert _ok(at), [e.value for e in at.exception]
    txt = _text(at)
    assert "2 variants read, 1 on the loaded window" in txt and "on chr9, not chr22" in txt
    assert "Mechanism simulator, not validated" in txt             # the measured SV-test verdict, read from JSON
    assert "synthetic" in txt                                      # reference input is labelled
    at.button(key="vi_fit").click().run()
    assert _ok(at), [e.value for e in at.exception]
    assert any(r["Stage"] == "ensemble fit (variant impact)" for r in at.session_state["telemetry"])
    txt = _text(at)
    assert "Pairs gaining ≥ 2×" in txt and "Uncertainty" in txt
    genes = _table(at, "Effect")
    assert genes is not None
    pairs = _table(at, "Locus A (Mb)")
    assert len(pairs) > 0 and (pairs["log₂ change"].abs() >= 1).any()
    assert "90 % interval" not in pairs.columns
    at.button(key="vi_boot").click().run()
    assert _ok(at), [e.value for e in at.exception]
    assert "90 % interval" in _table(at, "Locus A (Mb)").columns
    assert "Of those, 90 % interval excludes no change" in _text(at)
    assert any("bootstrap" in r["Stage"] for r in at.session_state["telemetry"])


def _chart(at, key: str):
    import json
    for c in at.get("plotly_chart"):
        if c.proto.id.endswith(f"-{key}"):
            return c.proto, json.loads(c.proto.spec)
    return None, None


def test_bead_tooltips_click_to_pick_and_linked_gene(app):
    at = app
    proto, spec = _chart(at, "viewport")
    assert proto is not None and list(proto.selection_mode) == []            # default: clicks do nothing, as before
    mesh = next(d for d in spec["data"] if d["type"] == "mesh3d")
    assert "bead %{customdata[5]" in mesh["hovertemplate"]                    # bead number in the tooltip
    assert any(t for t in mesh["text"])                                       # gene names at the beads
    assert "to bead" not in mesh["hovertemplate"]
    at.toggle(key="pick_on").set_value(True).run()
    assert _ok(at), [e.value for e in at.exception]
    assert list(_chart(at, "viewport")[0].selection_mode)                     # 3D view now reports clicked beads
    assert list(_chart(at, "matrix")[0].selection_mode)                       # and so does the map
    at.toggle(key="probe_on").set_value(True).run()
    _, spec = _chart(at, "viewport")
    mesh = next(d for d in spec["data"] if d["type"] == "mesh3d")
    assert "to bead 1 (this structure)" in mesh["hovertemplate"]
    heat = _chart(at, "matrix")[1]["data"]
    assert len(heat) == 2                                                      # the probe pair is marked on the map
    # linked selection: a gene picked on the Genes page is marked in the 3D view
    at.segmented_control(key="workspace").set_value("Genes").run()
    at.selectbox(key="genes_pick_chr22").set_value("BCR").run()
    assert _ok(at), [e.value for e in at.exception]
    assert at.session_state["linked_gene"][0] == "BCR"
    at.segmented_control(key="workspace").set_value("3D structure").run()
    assert _ok(at), [e.value for e in at.exception]
    _, spec = _chart(at, "viewport")
    assert any(d.get("text") == ["  BCR"] for d in spec["data"])
