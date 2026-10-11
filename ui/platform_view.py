"""
Sidebar panels for the local platform (Phase B7):
    Projects   save the session (settings, population models and telemetry built this session, the last analysis /
               differential / variant reports) under a name, and reopen it after the app restarts;
    Jobs       the local job queue (chronocell.jobs): submit a `chronocell` command as a GPU or CPU job, start the
               worker, watch progress (the job's log), stop a job.
Everything stays on this computer.
"""

from __future__ import annotations

import json
import shlex
import shutil
from pathlib import Path

import streamlit as st

from ui.common import esc, html

ss = st.session_state
SETTINGS_KEYS = ("assembly_choice", "chrom_choice", "region_choice", "custom_window", "workspace", "research_mode",
                 "disp_style", "disp_colour", "disp_radius", "clip_on", "clip_axis", "clip_pos", "unc_on", "probe_on",
                 "probe_level", "pop_input")
REPORT_KEYS = (("an_last", "analysis"), ("df_last", "differential"), ("ve_last", None))


def _open_project(name: str) -> None:
    """on_click callback: runs before the next script run, so widget values can still be set."""
    from chronocell import projects as PRJ
    try:
        p = PRJ.load(name)
    except (OSError, ValueError) as exc:
        ss.prj_msg = f"Could not open {name}: {exc}"
        return
    for k, v in p["settings"].items():
        if k in SETTINGS_KEYS:
            ss[k] = tuple(v) if k == "custom_window" and isinstance(v, list) else v
    res = p["results"] or {}
    if res.get("ensembles"):
        ss["ensembles"] = res["ensembles"]
    if res.get("telemetry"):
        ss["telemetry"] = res["telemetry"]
    ss.prj_msg = f"Opened {name} (saved {p['saved_utc']})."


def projects_panel() -> None:
    from chronocell import projects as PRJ
    with st.expander("Projects", expanded=False):
        name = st.text_input("Project name", key="prj_name", placeholder="e.g. HCT116 cohesin loss")
        if st.button("Save this session", key="prj_save", icon=":material/save:", disabled=not name, width="stretch"):
            settings = {k: (list(ss[k]) if isinstance(ss.get(k), tuple) else ss[k]) for k in SETTINGS_KEYS if k in ss}
            results = {"ensembles": ss.get("ensembles", {}), "telemetry": ss.get("telemetry", [])}
            try:
                d = PRJ.save(name, settings, None, results, notes="Saved from the app.")
                skipped = json.loads((d / "project.json").read_text(encoding="utf-8")).get("skipped_results", [])
                if skipped:
                    ss.prj_msg = (f"Saved {d.name}; {len(skipped)} result(s) could not be stored and were skipped "
                                  f"({'; '.join(skipped[:3])}). Rebuild them after opening.")
                rep = d / "reports"
                for key, sub in REPORT_KEYS:
                    out = (ss.get(key) or {}).get("out")
                    if out and Path(out).exists():
                        shutil.copytree(out, rep / (sub or "variant_impact"), dirs_exist_ok=True)
                st.toast(f"Saved {d.name}.")
            except (OSError, ValueError, TypeError) as exc:
                ss.prj_msg = f"Not saved: {exc}"
        projects = PRJ.list_projects()
        if projects:
            pick = st.selectbox("Saved projects", [p["name"] for p in projects], key="prj_pick",
                                format_func=lambda n: f"{n} · {next(p['saved_utc'][:16] for p in projects if p['name'] == n)}")
            st.button("Open", key="prj_open", on_click=_open_project, args=(pick,), icon=":material/folder_open:",
                      width="stretch")
        if ss.get("prj_msg"):
            html(f'<p class="cc-note">{esc(ss.prj_msg)}</p>')
        html('<p class="cc-note">Projects live in <code>.chronocell_cache/projects/</code> on this computer. Data files '
             'from the sidebar states stay where they are; uploads are not copied.</p>')


def jobs_panel() -> None:
    from chronocell import jobs as JQ
    q = JQ.Queue()
    with st.expander("Jobs", expanded=False):
        dev = st.segmented_control("Device", ["gpu", "cpu"], default="gpu", required=True, key="job_dev",
                                   help="At most one GPU job runs at a time; CPU jobs run in parallel up to a limit.")
        cmd = st.text_input("chronocell command", key="job_cmd",
                            placeholder="analyze data.mcool --region chr21:28000000-30000000 --res 10000")
        if st.button("Submit", key="job_submit", disabled=not cmd.strip(), icon=":material/play_arrow:", width="stretch"):
            q.submit(shlex.split(cmd, posix=False), dev)
            JQ.start_worker()
            st.toast("Submitted; the worker runs it in the background.")
        jobs = q.jobs()[-8:]
        for j in reversed(jobs):
            html(f'<p class="cc-meta"><b>{esc(j["state"])}</b> · {esc(j["device"])} · {esc(j["label"])}</p>')
        running = [j for j in jobs if j["state"] in ("queued", "running")]
        if running:
            stop = st.selectbox("Stop a job", [j["id"] for j in running], key="job_stop_pick",
                                format_func=lambda i: next(j["label"] for j in running if j["id"] == i))
            if st.button("Stop", key="job_stop", icon=":material/stop_circle:", width="stretch"):
                q.stop(stop)
        if jobs:
            st.code(q.log(jobs[-1]["id"], 12) or "(no output yet)", language="text")
        if st.button("Start the worker", key="job_worker", icon=":material/memory:", width="stretch",
                     help="Runs queued jobs; jobs left running when the app stopped are resumed."):
            st.toast(f"Worker running (PID {JQ.start_worker()}).")
