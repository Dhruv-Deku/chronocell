"""
Selection plumbing shared by the 3D view, the contact / distance map and the gene table.

Plotly selection events reach the app through st.plotly_chart(on_select=...) as dictionaries in
session state. These helpers turn them into beads (local indices of the view) and apply each event
only once, so a selection that Streamlit keeps between reruns does not override later choices.
"""

from __future__ import annotations

import math

import numpy as np

from chronocell import genes as G, genome


def _points(state) -> list[dict]:
    try:
        sel = state["selection"] if isinstance(state, dict) else state.selection
        pts = sel["points"] if isinstance(sel, dict) else sel.points
        return list(pts or [])
    except (KeyError, AttributeError, TypeError):
        return []


def _rows(state) -> list[int]:
    try:
        sel = state["selection"] if isinstance(state, dict) else state.selection
        rows = sel["rows"] if isinstance(sel, dict) else sel.rows
        return [int(r) for r in (rows or [])]
    except (KeyError, AttributeError, TypeError, ValueError):
        return []


def fresh(ss, key: str, kind: str = "points") -> list:
    """The points (or table rows) of widget `key`'s selection if it changed since last seen, else []."""
    state = ss.get(key)
    items = _rows(state) if kind == "rows" else _points(state)
    sig = repr(items)
    seen = ss.setdefault("_selection_seen", {})
    if not items or seen.get(key) == sig:
        seen[key] = sig
        return []
    seen[key] = sig
    return items


def bead_from_point(point: dict, g_lo: int, n_view: int) -> int | None:
    """Local bead of a clicked 3D point: every pickable trace carries the genomic bin as customdata[0]."""
    cd = point.get("customdata")
    if isinstance(cd, (list, tuple)) and len(cd):
        cd = cd[0]
    try:
        gb = float(cd)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(gb):
        return None
    k = int(round(gb)) - int(g_lo)
    return k if 0 <= k < n_view else None


def pair_from_map_point(point: dict, g_lo: int, n_view: int, resolution: int) -> tuple[int, int] | None:
    """Local beads (i, j) of a clicked contact / distance map pixel (axes in Mb at pixel centres)."""
    try:
        x, y = float(point["x"]), float(point["y"])
    except (KeyError, TypeError, ValueError):
        return None
    i = int(math.floor(x * 1e6 / resolution)) - int(g_lo)
    j = int(math.floor(y * 1e6 / resolution)) - int(g_lo)
    if not (0 <= i < n_view and 0 <= j < n_view) or i == j:
        return None
    return (min(i, j), max(i, j))


def gene_labels(chrom: genome.Chrom, bin0: int, n: int, max_names: int = 2) -> np.ndarray:
    """Genes (RefSeq, the chromosome's assembly) overlapping each bead: '' or 'A, B +2'."""
    tab = G.on_chromosome(chrom.name, chrom.assembly)
    res = chrom.resolution
    names: list[list[str]] = [[] for _ in range(n)]
    lo_bp, hi_bp = bin0 * res, (bin0 + n) * res
    sel = tab[(tab["start"] < hi_bp) & (tab["end"] > lo_bp)]
    for name, s, e in zip(sel["name"], sel["start"], sel["end"]):
        a = max(0, int(s) // res - bin0)
        b = min(n, (int(e) - 1) // res - bin0 + 1)
        for k in range(a, b):
            names[k].append(str(name))
    out = []
    for v in names:
        out.append("" if not v else ", ".join(v[:max_names]) + (f" +{len(v) - max_names}" if len(v) > max_names else ""))
    return np.array(out, dtype=object)
