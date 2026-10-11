"""
Record a real terminal session for the explainer: run a command, keep every output line with the time it appeared.

    python motion/explainer/record_terminal.py pytest      # the project's test suite (tests/), verbose
    python motion/explainer/record_terminal.py recheck     # motion/explainer/recheck.py: every score recomputed

Writes motion/explainer/data/term_<name>.json = {"cmd": shown command, "lines": [[t_seconds, text], ...], "secs": total}.
The explainer replays these lines at their real relative times (sped up, and says so); nothing is typed in by hand.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
HOME = re.compile(r"[A-Za-z]:\\Users\\[^\\\s]+", re.I)
ROOT = HERE.parent.parent
CMDS = {
    "pytest": ("python -m pytest tests -v", [sys.executable, "-u", "-m", "pytest", "tests", "-v", "-p", "no:cacheprovider",
                                             "--color=no", "-o", "console_output_style=progress"]),
    "recheck": ("python motion/explainer/recheck.py", [sys.executable, "-u", str(HERE / "recheck.py"), "--pace"]),
}


def record(name: str) -> Path:
    shown, cmd = CMDS[name]
    env = dict(os.environ, PYTHONIOENCODING="utf-8", COLUMNS="110")
    t0 = time.time()
    p = subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env, text=True,
                         encoding="utf-8", errors="replace", bufsize=1)
    lines = []
    for line in p.stdout:
        lines.append([round(time.time() - t0, 3), HOME.sub("~", line.rstrip("\n"))])     # no user name in the files
        print(line, end="", flush=True)
    p.wait()
    out = HERE / "data" / f"term_{name}.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({"cmd": shown, "lines": lines, "secs": round(time.time() - t0, 2), "code": p.returncode},
                              ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out} ({len(lines)} lines, {time.time() - t0:.0f} s, exit {p.returncode})")
    return out


if __name__ == "__main__":
    record(sys.argv[1] if len(sys.argv) > 1 else "pytest")
