"""
Per-sample and per-variant reports (Phase B6): a self-contained HTML file (no scripts, no external requests,
so it opens offline and sends nothing anywhere) and, on request, a PDF with the same content (fpdf2).

Each report carries: every measured value of the run (summary.json), the standing of the method behind it read
from the validation result files (a gate that has not run or not passed is labelled so; mechanism simulators are
labelled "not validated"), the result tables, the list of exported files, and the reproducibility record
(run.json: inputs with SHA-256, software versions, parameters).
"""

from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd

VALIDATION = Path(__file__).resolve().parent.parent / "validation"
GATES = {"4c": ("results_cohesin_hic.json", "Cohesin loss predicted for held-out Hi-C regions"),
         "4d": ("results_sv_v2.json", "Structural-variant effects on new events (Hi-C before and after)"),
         "6": ("results_gate6.json", "Loop calls against reference calls on held-out data"),
         "6b": ("results_gate6b.json", "Loop calls, settings re-chosen on three cell lines, two new cell lines"),
         "7": ("results_gate7.json", "False-discovery control of the differential analysis")}
KIND_GATES = {"analyze": ["6", "6b"], "diff": ["7"], "impact": ["4c", "4d"]}
CSS = """body{font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;margin:24px auto;max-width:1100px;padding:0 16px;
color:#1d2433;background:#fff}h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 8px;border-bottom:1px solid
#dfe3ea;padding-bottom:4px}table{border-collapse:collapse;width:100%;font-size:12.5px;margin:6px 0}td,th{border:1px solid
#dfe3ea;padding:3px 6px;text-align:left}th{background:#f4f6f9}.tag{display:inline-block;padding:1px 8px;border-radius:10px;
font-size:12px;margin-right:6px}.ok{background:#e3f4ea;color:#185c37}.warn{background:#fdf0dc;color:#7a4b00}.no{background:
#f1f2f5;color:#4b5363}code,pre{font:12px Consolas,monospace;background:#f6f8fa}pre{padding:8px;overflow-x:auto}
.note{color:#4b5363;font-size:13px}"""


def gate_standing(gate: str) -> tuple[str, str]:
    """(label, css class) for a gate from its result file."""
    fname, what = GATES[gate]
    path = VALIDATION / fname
    if not path.exists():
        return f"Gate {gate} ({what}): not run", "no"
    try:
        v = json.loads(path.read_text(encoding="utf-8")).get("verdict", "unknown")
    except (OSError, ValueError):
        return f"Gate {gate}: result file unreadable", "no"
    v = v if isinstance(v, str) else json.dumps(v)
    return f"Gate {gate} ({what}): {v}", "ok" if v == "pass" else "warn"


def _table(df: pd.DataFrame, rows: int = 50) -> str:
    if df is None or not len(df):
        return '<p class="note">No rows.</p>'
    more = f'<p class="note">First {rows} of {len(df):,} rows; the full table is in the folder.</p>' if len(df) > rows else ""
    return df.head(rows).to_html(index=False, escape=True, border=0, float_format=lambda x: f"{x:.4g}") + more


def _kind(folder: Path) -> str:
    if (folder / "variant_ranking.csv").exists():
        return "impact"
    if (folder / "differential_pixels_all.tsv").exists():
        return "diff"
    return "analyze"


def build_html(folder: Path) -> str:
    folder = Path(folder)
    kind = _kind(folder)
    summ = json.loads((folder / "summary.json").read_text(encoding="utf-8")) if (folder / "summary.json").exists() else {}
    run = json.loads((folder / "run.json").read_text(encoding="utf-8")) if (folder / "run.json").exists() else {}
    title = {"analyze": "Contact-map analysis", "diff": "Differential analysis", "impact": "Variant impact"}[kind]
    parts = [f"<h1>{html.escape(title)}</h1>", f'<p class="note">{html.escape(folder.name)} · ChronoCell-5D · research use '
             'only, not a diagnosis.</p>', "<h2>Standing of the method</h2>"]
    standing = summ.get("standing") or (summ.get("summary") or {}).get("standing")
    if standing:
        parts.append(f"<p>{html.escape(standing)}</p>")
    for g in KIND_GATES[kind]:
        lab, cls = gate_standing(g)
        parts.append(f'<span class="tag {cls}">{html.escape(lab)}</span>')
    if kind == "impact":
        parts.append('<p><span class="tag warn">Mechanism simulator, not validated</span>Every change below is a prediction of '
                     'the model, not a measurement.</p>')
    parts.append("<h2>Measured values</h2>")
    s = summ.get("summary", summ)
    if kind == "impact":
        parts.append(_table(pd.DataFrame(summ.get("variants", [])).drop(columns=["standing"], errors="ignore")))
    else:
        flat = {k: v for k, v in s.items() if not isinstance(v, (dict, list)) and k != "standing"}
        parts.append(_table(pd.DataFrame([flat])))
    notes = summ.get("notes") or []
    if notes:
        parts.append("<ul>" + "".join(f"<li>{html.escape(str(n))}</li>" for n in notes) + "</ul>")
    tables = {"analyze": [("Loops", "loops.bedpe", "bedpe"), ("Domains", "domains.bed", "bed"), ("Boundaries", "boundaries.bed", "bed")],
              "diff": [("Differential pixels (significant)", "differential_pixels.bedpe", "bedpe"),
                       ("Loops", "differential_loops.csv", "csv"), ("Boundaries", "differential_boundaries.csv", "csv"),
                       ("Compartment switches", "differential_compartments.csv", "csv")],
              "impact": [("Variant ranking (heuristic, not validated)", "variant_ranking.csv", "csv"),
                         ("Genes affected", "genes.csv", "csv"), ("Enhancer-promoter pairs", "ep_pairs.csv", "csv"),
                         ("Domain boundaries", "boundaries.csv", "csv")]}[kind]
    for name, fname, fmt in tables:
        p = folder / fname
        parts.append(f"<h2>{html.escape(name)}</h2>")
        if not p.exists() or p.stat().st_size == 0:
            parts.append('<p class="note">None.</p>')
            continue
        try:
            df = pd.read_csv(p) if fmt == "csv" else pd.read_csv(p, sep="\t", header=None)
        except (pd.errors.EmptyDataError, pd.errors.ParserError):
            parts.append('<p class="note">None.</p>')
            continue
        parts.append(_table(df))
    files = sorted(x.name for x in folder.iterdir() if x.is_file() and not x.name.startswith("report"))
    parts.append("<h2>Files</h2><p>" + ", ".join(f"<code>{html.escape(f)}</code>" for f in files) + "</p>")
    parts.append('<p class="note">BED / BEDPE / bedGraph open in IGV; loops_juicebox_2d.txt loads in Juicebox as 2D '
                 'annotations; a .mcool export opens in HiGlass.</p>')
    parts.append("<h2>Reproducibility record</h2><pre>" + html.escape(json.dumps(run, indent=1, default=str)) + "</pre>")
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,"
            f"initial-scale=1'><title>{html.escape(title)}</title><style>{CSS}</style></head><body>" + "\n".join(parts)
            + "</body></html>")


def _latin1(text: str) -> str:
    rep = {"·": "-", "–": "-", "—": "-", "≥": ">=", "≤": "<=", "→": "->", "×": "x", "₂": "2", "≈": "~", "…": "..."}
    for a, b in rep.items():
        text = text.replace(a, b)
    return text.encode("latin-1", "replace").decode("latin-1")


def build_pdf(folder: Path) -> bytes:
    from fpdf import FPDF
    folder = Path(folder)
    kind = _kind(folder)
    summ = json.loads((folder / "summary.json").read_text(encoding="utf-8")) if (folder / "summary.json").exists() else {}
    pdf = FPDF()
    pdf.set_auto_page_break(True, 12)
    pdf.add_page()
    pdf.set_font("Helvetica", "B", 15)
    pdf.cell(0, 9, _latin1({"analyze": "Contact-map analysis", "diff": "Differential analysis",
                             "impact": "Variant impact"}[kind]), new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 9)
    pdf.multi_cell(0, 5, _latin1(f"{folder.name} - ChronoCell-5D - research use only, not a diagnosis."), new_x="LMARGIN", new_y="NEXT")
    for g in KIND_GATES[kind]:
        pdf.multi_cell(0, 5, _latin1(gate_standing(g)[0]), new_x="LMARGIN", new_y="NEXT")
    if kind == "impact":
        pdf.multi_cell(0, 5, "Mechanism simulator, not validated: every change is a prediction, not a measurement.", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(2)
    pdf.set_font("Helvetica", "B", 11)
    pdf.cell(0, 7, "Measured values", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Courier", "", 8)
    s = summ.get("summary", summ)
    body = json.dumps(summ.get("variants", s), indent=1, default=str)
    for line in body.splitlines()[:300]:
        pdf.multi_cell(0, 4, _latin1(line[:120]), new_x="LMARGIN", new_y="NEXT")
    run = folder / "run.json"
    if run.exists():
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 11)
        pdf.cell(0, 7, "Reproducibility record", new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Courier", "", 7)
        for line in run.read_text(encoding="utf-8").splitlines()[:400]:
            pdf.multi_cell(0, 3.5, _latin1(line[:140]), new_x="LMARGIN", new_y="NEXT")
    return bytes(pdf.output())


def report_folder(folder: Path, pdf: bool = False) -> Path:
    folder = Path(folder)
    out = folder / "report.html"
    out.write_text(build_html(folder), encoding="utf-8")
    if pdf:
        (folder / "report.pdf").write_bytes(build_pdf(folder))
    return out
