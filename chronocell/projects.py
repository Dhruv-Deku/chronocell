"""
Saved projects (Phase B7): a named folder under .chronocell_cache/projects/ holding what is needed to reopen a
session after the app restarts:

    project.json   settings (plain values of the app's widgets), the list of data files with SHA-256, notes
    data/          copies of the data files the session used (kept on this computer only)
    results.pkl    results built in the session (population models, telemetry), written and read only by this
                   app on this computer; a project folder from anywhere else is not loaded (see `load`)
    reports/       reports and exports saved with the project

Nothing is uploaded; patient-like data stay on this machine.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import pickle
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / ".chronocell_cache" / "projects"
MARK = "ChronoCell-5D project v1"


def _safe(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9 _.-]+", "", (name or "").strip())[:60].strip()
    if not s:
        raise ValueError("Give the project a name (letters, digits, space, _ . -).")
    return s


def _plain(v) -> bool:
    return isinstance(v, (str, int, float, bool, type(None))) or (isinstance(v, (list, tuple)) and all(_plain(x) for x in v))


def _picklable(results: dict) -> tuple[dict, list[str]]:
    """The result items that can be serialised (each entry of a dict-valued item checked on its own)."""
    keep, skipped = {}, []
    for key, value in (results or {}).items():
        items = value.items() if isinstance(value, dict) else [(None, value)]
        out = {}
        for k, v in items:
            try:
                pickle.dumps(v, protocol=pickle.HIGHEST_PROTOCOL)
            except Exception as exc:                 # noqa: BLE001 - any object the pickler refuses is skipped
                skipped.append(f"{key}{'' if k is None else f' / {k}'}: {type(exc).__name__}")
                continue
            if k is None:
                keep[key] = v
            else:
                out[k] = v
        if isinstance(value, dict):
            keep[key] = out
    return keep, skipped


def save(name: str, settings: dict, files: dict[str, bytes] | None = None, results: dict | None = None,
         notes: str = "", root: Path | None = None) -> Path:
    """Save (or overwrite) a project; only plain settings values are kept. Results that cannot be serialised are
    skipped and listed in project.json ('skipped_results'). The project is written to a temporary folder first and
    moved into place, so a failed save never leaves a half-written project."""
    base_dir = Path(root if root is not None else ROOT)
    d = base_dir / _safe(name)
    tmp = base_dir / f".{_safe(name)}.saving"
    if tmp.exists():
        shutil.rmtree(tmp)
    (tmp / "data").mkdir(parents=True, exist_ok=True)
    listed = []
    for fname, data in (files or {}).items():
        base = Path(fname).name
        (tmp / "data" / base).write_bytes(data)
        listed.append({"name": base, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)})
    keep, skipped = _picklable(results or {})
    meta = {"format": MARK, "name": _safe(name), "saved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "settings": {k: v for k, v in settings.items() if _plain(v)}, "files": listed, "notes": notes,
            "has_results": bool(keep), "skipped_results": skipped}
    if keep:
        with (tmp / "results.pkl").open("wb") as fh:
            pickle.dump(keep, fh, protocol=pickle.HIGHEST_PROTOCOL)
    (tmp / "project.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    if d.exists():
        shutil.rmtree(d)
    tmp.replace(d)
    return d


def list_projects(root: Path | None = None) -> list[dict]:
    out = []
    for p in sorted(Path(root if root is not None else ROOT).glob("*/project.json")):
        try:
            m = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if m.get("format") == MARK:
            out.append({"name": m["name"], "saved_utc": m["saved_utc"], "files": len(m["files"]), "has_results": m["has_results"]})
    return out


def load(name: str, root: Path | None = None) -> dict:
    """Settings, data files (checked against their saved SHA-256) and results of a project saved by this app
    under `root`. Pickled results are read only from that folder, which only this app writes."""
    d = Path(root if root is not None else ROOT) / _safe(name)
    meta = json.loads((d / "project.json").read_text(encoding="utf-8"))
    if meta.get("format") != MARK:
        raise ValueError("Not a ChronoCell project folder.")
    files = {}
    for f in meta["files"]:
        data = (d / "data" / f["name"]).read_bytes()
        if hashlib.sha256(data).hexdigest() != f["sha256"]:
            raise ValueError(f"{f['name']} changed since the project was saved (SHA-256 differs).")
        files[f["name"]] = data
    results = None
    if meta.get("has_results") and (d / "results.pkl").exists():
        with (d / "results.pkl").open("rb") as fh:
            results = pickle.load(fh)                    # noqa: S301 - written by this app in its own cache folder
    return {"settings": meta["settings"], "files": files, "results": results, "notes": meta.get("notes", ""),
            "saved_utc": meta["saved_utc"]}


def delete(name: str, root: Path | None = None) -> None:
    d = Path(root if root is not None else ROOT) / _safe(name)
    if (d / "project.json").exists():
        shutil.rmtree(d)
