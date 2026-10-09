"""
Pack the explainer's inputs into motion/explainer/data/explainer_data.js (window.EX), so the page works from file://.

    python motion/explainer/build.py

- clips: each screen recording's length and click times (clips/*.json, from record_app.py)
- terms: the recorded terminal sessions (data/term_*.json, from record_terminal.py). The only change made to the output
  is cosmetic: this computer's home folder is shown as "~" so the video does not show a user name.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
HOME = re.compile(r"[A-Za-z]:\\Users\\[^\\\s]+", re.I)


def main() -> None:
    clips = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted((HERE / "clips").glob("c*.json"))}
    terms = {}
    for p in sorted((HERE / "data").glob("term_*.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        d["lines"] = [[t, HOME.sub("~", s)] for t, s in d["lines"]]
        terms[p.stem[5:]] = d
    out = HERE / "data" / "explainer_data.js"
    out.write_text("window.EX = " + json.dumps({"clips": clips, "terms": terms}, ensure_ascii=False) + ";\n", encoding="utf-8")
    print(f"wrote {out}: {len(clips)} clips ({sum(c['secs'] for c in clips.values()):.0f} s), terminals {list(terms)}")


if __name__ == "__main__":
    main()
