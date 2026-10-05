"""
Phase B9: which third-party tools run on this machine, and where each comparison is. Checks are made when the
script runs (executables on PATH, the isolated tools environment .chronocell_cache/tools-venv) and written to
validation/results_tools_b9.json; RESULTS.md shows the table through report.py. Tools that could not run are
recorded with the reason, as validation/benchmark/methods.py does for Gate 3.

    python validation/tools_b9.py
"""

from __future__ import annotations

import datetime as dt
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT.parent / ".chronocell_cache" / "tools-venv" / "Scripts" / "python.exe"


def _venv_version(pkg: str) -> str | None:
    if not VENV.exists():
        return None
    r = subprocess.run([str(VENV), "-m", "pip", "show", pkg], capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.startswith("Version:"):
            return line.split(":", 1)[1].strip()
    return None


def main() -> None:
    java, r_exe = shutil.which("java"), shutil.which("Rscript") or shutil.which("R")
    no_java = "not run: needs Java (Juicer tools); no Java runtime on this machine" if not java else "Java present"
    no_r = "not run: R / Bioconductor package; R is not installed on this machine" if not r_exe else "R present"
    cs, mu = _venv_version("chromosight"), _venv_version("mustache-hic")
    gate6 = "Gate 6 (validation/results_gate6.json)"
    tools = [
        {"task": "loops", "tool": "HiCCUPS", "status": no_java + "; ENCODE's HiCCUPS calls are Gate 6's reference",
         "result": gate6},
        {"task": "loops", "tool": "chromosight", "status": f"run ({cs}, isolated environment)" if cs else "not installed",
         "result": gate6},
        {"task": "loops", "tool": "Mustache", "status": f"run ({mu}, isolated environment with NumPy 1.26; reads .mcool "
                                                        "because hic-straw does not build here)" if mu else "not installed",
         "result": gate6},
        {"task": "TADs", "tool": "TopDom", "status": no_r, "result": "ChronoCell's TopDom-like caller (analysis suite); no "
                                                                      "pre-registered boundary gate"},
        {"task": "TADs", "tool": "Arrowhead", "status": no_java, "result": "ChronoCell's Arrowhead-like caller; no gate"},
        {"task": "TADs", "tool": "insulation (cooltools)", "status": "not run: cooltools has no Windows build (needs a C "
                                                                      "compiler)", "result": "ChronoCell's insulation score"},
        {"task": "compartments", "tool": "dcHiC", "status": no_r, "result": "—"},
        {"task": "differential", "tool": "diffHic", "status": no_r, "result": "Gate 7 (validation/results_gate7.json)"},
        {"task": "differential", "tool": "multiHiCcompare", "status": no_r, "result": "Gate 7"},
        {"task": "differential", "tool": "CHESS", "status": "not run: pip install chess-hic fails building pysam (a FAN-C "
                                                            "dependency) on Windows", "result": "Gate 7"},
        {"task": "SV detection", "tool": "HiNT", "status": "not run: not installable here (needs R and BWA / samtools)",
         "result": "—"},
        {"task": "SV detection", "tool": "hic_breakfinder", "status": "not run: C++ (Eigen, BamTools) with no Windows "
                                                                       "binary; no compiler here", "result": "—"},
        {"task": "SV detection", "tool": "EagleC", "status": "not run: no held-out SV truth set with Hi-C here (Gate 4d "
                                                             "is blocked), so there is nothing to score it on",
         "result": "—"},
        {"task": "SV impact", "tool": "Akita (basenji)", "status": "not run: TensorFlow model code not packaged for pip; "
                                                                   "no SV truth set (Gate 4d)", "result": "—"},
        {"task": "SV impact", "tool": "Orca", "status": "not run: needs selene-sdk (no Windows build) and large weights; "
                                                        "no SV truth set (Gate 4d)", "result": "—"},
    ]
    out = {"machine": {"java": java, "R": r_exe, "tools_venv": str(VENV) if VENV.exists() else None}, "tools": tools,
           "written_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
    (ROOT / "results_tools_b9.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for t in tools:
        print(f"{t['task']:14s} {t['tool']:24s} {t['status'][:90]}")


if __name__ == "__main__":
    main()
