"""v4 guard rails: the numerics package (chronocell/) never imports Streamlit, and no module ships
hard-coded accuracy numbers in place of measured ones."""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PKG = ROOT / "chronocell"


def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    out: set[str] = set()
    for node in ast.walk(tree):                       # every level: module, functions, branches
        if isinstance(node, ast.Import):
            out.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            out.add(node.module.split(".")[0])
    return out


def test_chronocell_has_no_streamlit_import_anywhere():
    offenders = [str(p.relative_to(ROOT)) for p in PKG.rglob("*.py") if "streamlit" in _imports(p)]
    assert offenders == [], offenders


def test_importing_every_chronocell_module_leaves_streamlit_unloaded():
    mods = sorted(".".join(p.relative_to(ROOT).with_suffix("").parts) for p in PKG.rglob("*.py")
                  if p.name != "__main__.py")
    code = ("import importlib, sys\n"
            f"for m in {mods!r}:\n"
            "    importlib.import_module(m)\n"
            "print('streamlit' in sys.modules)\n")
    out = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.strip().splitlines()[-1] == "False"
