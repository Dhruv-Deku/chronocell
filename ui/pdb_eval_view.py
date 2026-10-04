"""
03 Compare -> "Self-Math PDB State Evaluator" sub-tab (chronocell/analytics/pdb_evaluator.py).

Sources: any structure loaded in the workstation (whole, or a Mb window), the members of a population
model built this session, or an uploaded PDB file. Nothing is computed until "Evaluate" is pressed
(Streamlit renders every tab on each run, so the Compare page itself costs the same as before).
Every number shown is computed from the coordinates; every state threshold is shown with its source.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import streamlit as st

from chronocell import theme as T, viz
from chronocell.analytics import pdb_evaluator as PE
from ui.common import Dataset, banner, esc, fmt, html, readout, warning_card

ss = st.session_state
STATE_TAG = {"Normal": "ok", "Senescent": "warn", "Diseased": "warn", "Indeterminate": ""}


@st.cache_resource(show_spinner="Evaluating the structure…", max_entries=16)
def _evaluate(key: str, _x: np.ndarray) -> PE.Evaluation:
    return PE.evaluate(_x)


@st.cache_resource(show_spinner="Evaluating every member of the population…", max_entries=8)
def _evaluate_population(key: str, _members: np.ndarray) -> dict:
    return PE.evaluate_population(_members)


def _sources(options: list[tuple[str, Dataset]], populations: dict) -> list[tuple[str, str]]:
    out = [(f"structure:{i}", lab) for i, (lab, _) in enumerate(options)]
    for k, res in populations.items():
        n = res.representative_nm.shape[0]
        kind = (" · predicted from sequence + CTCF (no contact data)"
                if str(res.config.get("input", "")).startswith("predicted") else "")
        out.append((f"population:{k}", f"Population model{kind} · {n} beads · {res.frames_nm.shape[0]} members (this session)"))
    out.append(("upload", "Upload a PDB file"))
    return out


def render(options: list[tuple[str, Dataset]], populations: dict | None = None) -> None:
    populations = populations or {}
    html('<p class="cc-note">Pure mathematics on 3D coordinates: radius of gyration, gyration-tensor shape '
         '(asphericity, acylindricity), packing density in the convex hull, and the power law of distance against '
         'separation. A rule set then names the geometry <b>Normal</b>, <b>Diseased</b>, <b>Senescent</b> or '
         '<b>Indeterminate</b>. No labels or stored answers are used. The rules are not validated against '
         'labelled disease or senescence structures, so the label describes the geometry: it is not a diagnosis.</p>')
    srcs = _sources(options, populations)
    keys = [k for k, _ in srcs]
    labels = dict(srcs)
    c1, c2 = st.columns([1.4, 1])
    pick = c1.selectbox("Structure to evaluate", keys, format_func=lambda k: labels[k], key="pe_source")
    coords, name, members = None, labels[pick], None
    if pick.startswith("structure:"):
        _, ds = options[int(pick.split(":")[1])]
        ch = ds.chrom
        mb0, mb1 = float(ch.bin_start(ds.bin0)) / 1e6, float(ch.bin_end(ds.bin0 + ds.n - 1)) / 1e6
        whole = c2.toggle("Whole loaded structure", True, key="pe_whole")
        if whole:
            coords = ds.frames[0]
        else:
            a = c2.number_input("From (Mb)", min_value=round(mb0, 2), max_value=round(mb1, 2), value=round(mb0, 2),
                                step=0.5, key="pe_a")
            b = c2.number_input("To (Mb)", min_value=round(mb0, 2), max_value=round(mb1, 2),
                                value=round(min(mb1, mb0 + 5.0), 2), step=0.5, key="pe_b")
            lo = int(np.clip(int(a * 1e6 / ch.resolution) - ds.bin0, 0, ds.n - 8))
            hi = int(np.clip(int(np.ceil(b * 1e6 / ch.resolution)) - ds.bin0, lo + 8, ds.n))
            coords = ds.frames[0][lo:hi]
            name = f"{labels[pick]} · {ch.name}:{float(ch.bin_start(ds.bin0 + lo)) / 1e6:.2f}–{float(ch.bin_end(ds.bin0 + hi - 1)) / 1e6:.2f} Mb"
        unit_note = "nm (workstation coordinates)"
        resolution = ch.resolution
    elif pick.startswith("population:"):
        res = populations[pick.split(":", 1)[1]]
        coords, members = res.representative_nm, res.frames_nm
        unit_note, resolution = "nm (population model)", None
    else:
        up = c2.file_uploader("PDB file", type=["pdb", "ent"], key="pe_upload",
                              help="Any PDB: chromatin bead models (ChronoCell exports are read in nm) or atomic "
                                   "structures (Angstrom, converted to nm). Every ATOM/HETATM record is one point, in "
                                   "file order; several MODELs are evaluated as a population.")
        resolution = None
        unit_note = ""
        if up is not None:
            try:
                models, unit_note = PE.coords_from_pdb_text(up.getvalue().decode("utf-8", "replace"))
                coords = models[0]
                members = np.stack(models) if len(models) > 1 and len({len(m) for m in models}) == 1 else None
                name = up.name
            except ValueError as exc:
                warning_card("That PDB file could not be read", str(exc))
                return
    run = st.button("Evaluate", key="pe_go", type="primary", disabled=coords is None)
    if run and coords is not None:
        key = f"{pick}:{len(coords)}:{hash(np.asarray(coords, np.float64).tobytes())}"
        try:
            ev = _evaluate(key, np.asarray(coords, np.float64))
        except ValueError as exc:
            warning_card("This structure cannot be evaluated", str(exc))
            return
        pop = _evaluate_population(key, members) if members is not None and len(members) > 1 else None
        ss.pe_last = {"key": key, "name": name, "unit": unit_note, "ev": ev, "pop": pop, "resolution": resolution}
    last = ss.get("pe_last")
    if not last:
        html('<p class="cc-note">Choose a structure and press <b>Evaluate</b>.</p>')
        return
    _show(last)


def _show(last: dict) -> None:
    ev: PE.Evaluation = last["ev"]
    html(f'<p class="cc-eyebrow" style="margin-top:10px">{esc(last["name"])}</p>'
         f'<p class="cc-meta"><span class="cc-tag {STATE_TAG.get(ev.state, "")}">{ev.state}</span> '
         f'&nbsp;{ev.n:,} points · coordinates in {esc(last["unit"] or "nm")}</p>')
    left, right = st.columns([1, 1.1], gap="large")
    with left:
        readout([
            ("Radius of gyration R<sub>g</sub><small>√((1/N) Σ|r<sub>i</sub> − r<sub>cm</sub>|²)</small>", fmt(ev.rg_nm), "nm"),
            ("Gyration-tensor eigenvalues λ₁ ≥ λ₂ ≥ λ₃", " · ".join(f"{v:,.0f}" for v in ev.eigenvalues_nm2), "nm²"),
            ("Asphericity Δ<small>(3/2) Σ(λ<sub>k</sub> − λ̄)² / (Tr S)² · 0 sphere, 1 rod</small>", fmt(ev.asphericity, 3), ""),
            ("Acylindricity c = λ₂ − λ₃", fmt(ev.acylindricity_nm2, 0), "nm²"),
            ("Convex-hull volume", fmt(ev.hull_volume_um3, 3), "µm³"),
            ("Packing density N / V<sub>hull</sub>", fmt(ev.packing_density_per_um3, 0), "beads per µm³"),
            ("Median bond b", fmt(ev.bond_nm, 1), "nm"),
            (f"Distance-decay slope γ<sub>d</sub><small>OLS log d<sub>ij</sub> on log |i−j|, window {PE.WINDOW_MIN}–"
             f"{max(2 * PE.WINDOW_MIN, ev.n // 10)} · R² {fmt(ev.gamma_d_r2, 3)}</small>", fmt(ev.gamma_d, 3), ""),
            ("Contact-scale exponent γ<sub>c</sub> = α γ<sub>d</sub><small>α = 3, the contact law d ∝ I<sup>−1/α</sup> "
             "used across ChronoCell</small>", fmt(ev.gamma_c, 3), ""),
            (f"Short-range γ<sub>c</sub><small>pairs 2–{PE.SHORT_MAX} beads apart</small>", fmt(ev.gamma_c_short, 3), ""),
            ("Scaling break<small>|γ<sub>d</sub>(lower half) − γ<sub>d</sub>(upper half)| of the window</small>",
             fmt(ev.scaling_break, 3), ""),
            ("Density spikes · anomaly runs", f"{100 * ev.spike_fraction:.1f} % · {len(ev.anomaly_runs)}", ""),
        ])
    with right:
        crit = pd.DataFrame([{"State": c.state, "Criterion": c.name, "Value": round(c.value, 3) if np.isfinite(c.value) else None,
                              "Threshold": c.threshold, "Met": "yes" if c.met else "no", "Threshold source": c.source}
                             for c in ev.criteria])
        st.dataframe(crit, hide_index=True, width="stretch", height=290)
        html('<p class="cc-note">Senescent: at least 2 of its 3 criteria. Diseased: at least 2 of its 3. Normal: γ<sub>c</sub> '
             'in range and neither rule set met. Otherwise Indeterminate. Thresholds marked <i>specification</i> come from '
             'the ChronoCell v4 brief; <i>assumption</i> marks values chosen here because the brief gave none. The '
             'references for R<sub>g</sub> and density are a compact globule of the same N and bond length at the volume '
             'fraction of ChronoCell\'s reference globule (0.19).</p>')
    c1, c2 = st.columns(2, gap="large")
    with c1:
        html('<p class="cc-eyebrow">Distance against separation</p>')
        st.plotly_chart(viz.decay_fit_chart(ev.decay_s, ev.decay_d, ev.gamma_d, (PE.WINDOW_MIN, max(2 * PE.WINDOW_MIN, ev.n // 10)),
                                            ev.bond_nm, last.get("resolution")),
                        theme=None, width="stretch", config=T.PLOT_CONFIG, key="pe_decay")
    with c2:
        html('<p class="cc-eyebrow">Local density along the chain</p>')
        st.plotly_chart(viz.density_z_chart(ev.density_z, PE.Thresholds().spike_z, ev.anomaly_runs), theme=None,
                        width="stretch", config=T.PLOT_CONFIG, key="pe_density")
    pop = last.get("pop")
    if pop:
        html(f'<p class="cc-eyebrow" style="margin-top:6px">Population · {pop["members"]} members evaluated one by one</p>')
        rows = [{"Metric": k, "Median": round(v["median"], 3), "10th pct": round(v["p10"], 3), "90th pct": round(v["p90"], 3)}
                for k, v in pop.items() if isinstance(v, dict) and "median" in v]
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        html('<p class="cc-note">States across members: ' + ", ".join(f"{k} {v}" for k, v in pop["states"].items())
             + '. This is ensemble spread (how consistent the population is), not accuracy.</p>')
    st.download_button("Evaluation (JSON)", json.dumps(ev.summary(), indent=1, default=float),
                       "chronocell_state_evaluation.json", "application/json", key="pe_json", icon=":material/download:")
