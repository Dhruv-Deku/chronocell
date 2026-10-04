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
