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
