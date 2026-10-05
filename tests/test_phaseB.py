"""Phase B: product features. Each test checks that the default app is what it was before Phase B and that
the new feature works when switched on."""

from __future__ import annotations

from pathlib import Path

import pytest

APP = str(Path(__file__).resolve().parent.parent / "app.py")


@pytest.fixture()
def app(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    import ui.common as C
    import ui.settings as SET
    import ui.states_panel as SP
    monkeypatch.setattr(C, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "DEMO_ROOT", tmp_path / "demo")
    monkeypatch.setattr(SET, "PATH", tmp_path / "settings.json")
    at = AppTest.from_file(APP, default_timeout=900)
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at


def _text(at) -> str:
    return "".join(m.value for m in at.markdown)


def test_settings_default_and_save(tmp_path, monkeypatch):
    import ui.settings as SET
    monkeypatch.setattr(SET, "PATH", tmp_path / "s.json")
    assert SET.load() == {"research_mode_default": True}               # unchanged app by default
    SET.save(research_mode_default=False, unknown_key=1)
    assert SET.load() == {"research_mode_default": False}
    assert "unknown_key" not in (tmp_path / "s.json").read_text()


def test_is_synthetic():
    from types import SimpleNamespace as NS
    from ui.common import SYNTHETIC_PREFIX, is_synthetic
    ref = NS(structure_label="x", inputs=(), is_reference=True)
    demo = NS(structure_label="Healthy · demo/synthetic_demo_healthy.pdb", inputs=(), is_reference=False)
    real = NS(structure_label="patient.pdb", inputs=(("structure", "patient.pdb", "0" * 64),), is_reference=False)
    assert is_synthetic(ref) and is_synthetic(demo) and not is_synthetic(real)
    assert is_synthetic(None, SYNTHETIC_PREFIX + "anything") and not is_synthetic(None, "anything")


def test_research_mode_switch_and_synthetic_labels(app):
    at = app
    assert at.toggle(key="research_mode").value is True                # default: everything as before
    assert any("Drug lab" in o for o in at.segmented_control(key="workspace").options)
    at.toggle(key="research_mode").set_value(False).run()
    assert not at.exception, [e.value for e in at.exception]
    assert not any("Drug lab" in o for o in at.segmented_control(key="workspace").options)
    at.segmented_control(key="workspace").set_value("Compare").run()
    assert any(str(o).startswith("SYNTHETIC · Reference model") for o in at.selectbox(key="pe_source").options)
    at.button(key="pe_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    txt = _text(at)
    assert "Research mode is off" in txt and "computed on SYNTHETIC data" in txt
    assert "Senescent: at least 2 of its 3 criteria" not in txt          # the rule table is hidden
    at.toggle(key="research_mode").set_value(True).run()
    assert "Senescent: at least 2 of its 3 criteria" in _text(at)
