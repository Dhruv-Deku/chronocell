"""
Local app settings (Phase B8), kept in .chronocell_cache/settings.json on this machine (git-ignored).

research_mode_default  whether the app opens in Research mode (all pages, including the mechanism
                       simulators and the rule-based state labels). True keeps the app as it always was.
"""

from __future__ import annotations

import json

from ui.common import CACHE_DIR

PATH = CACHE_DIR / "settings.json"
DEFAULTS = {"research_mode_default": True}


def load() -> dict:
    try:
        data = json.loads(PATH.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        data = {}
    return {**DEFAULTS, **{k: v for k, v in data.items() if k in DEFAULTS}}


def save(**changes) -> dict:
    data = {**load(), **{k: v for k, v in changes.items() if k in DEFAULTS}}
    PATH.parent.mkdir(parents=True, exist_ok=True)
    PATH.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    return data
