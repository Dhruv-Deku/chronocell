"""
Reproducibility record of one analysis: chronocell_audit_log.json and its PDF rendering.

The record says what was run (method, equations, parameters, software versions), on which inputs
(the SHA-256 of the exact bytes of every input file), which public datasets the method's accuracy
claims rest on (source, licence, citation), and every value measured in the session (structure
metrics, the two accuracy scores kept apart, execution telemetry). It is a reproducibility record:
it is not a clinical, regulatory or certified audit, and the file says so.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np

STATEMENT = ("Reproducibility record. It lets someone rerun this analysis on the same inputs and check the same "
             "numbers. It is not a clinical, regulatory or certified audit trail, and ChronoCell-5D is research "
             "software, not a diagnostic.")
SOURCES_PATH = Path(__file__).with_name("data") / "validation_sources.json"

EQUATIONS = {
    "contact_to_target_distance (v3.2 single structure)": "d*_ij = b0 (M_ij / M_ref)^(-1/alpha), alpha > 0",
    "single_structure_objective (v3.2)": "L = L_contact + lambda1 L_smooth + lambda2 L_steric",
    "contact_frequency_to_pair_variance (population)": "f_ij = P(|r_ij| < r_c) = F_Maxwell(r_c / sigma_ij); s_ij = sigma_ij^2",
    "population_objective": "L = (1/P) sum_{i<j} w_ij (log s_ij - log s*_ij)^2, w_ij = N f_ij / (1 - f_ij)",
    "population_model_v4": "s_ij = |a_i - a_j|^2 + sum_{k=i}^{j-1} v_k  (rank-r chain plus independent-bond random walk)",
    "pair_distance_distribution": "|r_ij| ~ Maxwell(sigma_ij): median 1.5382 sigma, mean 1.5958 sigma, sd 0.6734 sigma",
    "langevin": "dx = -P x dt + sqrt(2) dW, P = covariance^+ (exact Ornstein-Uhlenbeck propagation, mode by mode)",
    "radius_of_gyration": "R_g = sqrt((1/N) sum_i |x_i - x_cm|^2)",
    "scaling_exponent": "R(s) ~ s^nu: OLS of log RMS distance on log s, s in [4, N/10]",
    "contact_decay": "P(s) ~ s^-gamma",
    "microscopy_accuracy (benchmark only)": "trend-removed Spearman rho vs held-out imaging medians / the same for half A vs half B",
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def software_versions(extra: dict | None = None) -> dict:
    out = {"python": sys.version.split()[0], "platform": platform.platform(), "numpy": np.__version__}
    for mod in ("torch", "streamlit", "plotly", "pandas", "scipy"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:                 # optional packages
            out[mod] = None
    out.update(extra or {})
    return out


def validation_sources() -> list[dict]:
    """Datasets behind the method's benchmark (written from validation/datasets.py by the reproduce script)."""
    try:
        return json.loads(SOURCES_PATH.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def build(*, software: str, method: str, parameters: dict, inputs: list[dict], metrics: dict,
          accuracy: dict | None, telemetry: list[dict], genome: dict, extra: dict | None = None) -> dict:
    """The record as a JSON-ready dict. Every value passed in must be a measured or configured value."""
    rec = {
        "record_type": "ChronoCell-5D reproducibility record",
        "statement": STATEMENT,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "software": {"name": "ChronoCell-5D", "version": software, "environment": software_versions()},
        "genome": genome,
        "method": method,
        "equations": EQUATIONS,
        "parameters": parameters,
        "inputs": inputs,
        "measured_metrics": metrics,
        "accuracy": {"note": "Two separate scores, never combined: the contact-map fit is measured on this input; the "
                             "microscopy accuracy is the method's held-out benchmark and is not measured on this input.",
                     "scores": accuracy},
        "execution_telemetry": telemetry,
        "validation_datasets": validation_sources(),
    }
    if extra:
        rec["additional"] = extra
    rec = _clean(rec)
    rec["record_sha256"] = sha256_bytes(json.dumps(rec, sort_keys=True).encode())
    return rec


def to_json(record: dict) -> str:
    return json.dumps(record, indent=2, sort_keys=False, allow_nan=False)


def to_pdf(record: dict) -> bytes:
    """A readable PDF of the same record (fpdf2), with the statement on the first page and in the footer."""
    from .pdf_report import _Doc
    d = _Doc()
    pdf = d.pdf
    pdf.add_page()
    d.font(18, True)
    pdf.multi_cell(pdf.epw, 8, d.clean("ChronoCell-5D · Reproducibility record"), new_x="LMARGIN", new_y="NEXT")
    d.para(f"Generated {record['generated_utc']} · software {record['software']['name']} {record['software']['version']} · "
           f"record SHA-256 {record.get('record_sha256', '')[:16]}…", size=8.5)
    d.banner(record["statement"])

    def kv_table(title: str, items: dict | list, widths=(60, 120)):
        d.heading(title, size=11)
        rows = [["Field", "Value"]]
        if isinstance(items, dict):
            for k, v in items.items():
                rows.append([str(k), json.dumps(v, default=str) if isinstance(v, (dict, list)) else str(v)])
        d.table(rows, widths, size=7.6)

    kv_table("Genome and method", {"genome": record.get("genome"), "method": record.get("method")})
    kv_table("Software environment", record["software"]["environment"])
    d.heading("Inputs (SHA-256 of the exact bytes)", size=11)
    rows = [["Role", "File", "SHA-256"]] + [[i.get("role", ""), i.get("name", ""), i.get("sha256", "")]
                                            for i in record.get("inputs", [])]
    if len(rows) == 1:
        rows.append(["reference model", "synthetic (generated in the app)", "n/a"])
    d.table(rows, (30, 50, 100), size=7.2)
    kv_table("Parameters", record.get("parameters", {}))
    kv_table("Equations", record.get("equations", {}), (70, 110))
    kv_table("Measured metrics", record.get("measured_metrics", {}))
    acc = record.get("accuracy", {})
    d.heading("Accuracy: two separate scores", size=11)
    d.para(acc.get("note", ""), size=8.3)
    if acc.get("scores"):
        d.para(json.dumps(acc["scores"], indent=1, default=str)[:3500], size=7.2)
    tel = record.get("execution_telemetry") or []
    if tel:
        d.heading("Execution telemetry (measured in this session)", size=11)
        cols = list(tel[0].keys())[:8]
        d.table([cols] + [[("" if r.get(c) is None else str(r.get(c))) for c in cols] for r in tel],
                tuple([180 / len(cols)] * len(cols)), size=6.8)
    srcs = record.get("validation_datasets") or []
    if srcs:
        d.heading("Datasets behind the method's benchmark", size=11)
        d.table([["Dataset", "Role", "Licence", "Citation"]] +
                [[s.get("key", ""), s.get("role", ""), s.get("license", "")[:90], s.get("citation", "")[:120]] for s in srcs],
                (34, 16, 60, 70), size=6.6)
    return d.output()
