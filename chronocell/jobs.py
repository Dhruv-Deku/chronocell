"""
Local job queue for heavy work (Phase B7): fits, variant refits, differential analysis, batch runs. Everything
runs on this computer; nothing is sent anywhere.

A job is a `chronocell` command line (chronocell.cli) with a device class:
    gpu  at most one at a time (one population fit or refit series on the GPU)
    cpu  up to `max_cpu` in parallel
Jobs live in .chronocell_cache/jobs/<id>/: job.json (state, command, times), log.txt (the command's output),
out/ (its results, the cache). States: queued -> running -> done | failed | stopped.

    python -m chronocell.jobs submit gpu -- impact data.mcool --region chr9:... --res 10000 --variants sv.vcf
    python -m chronocell.jobs worker          # runs queued jobs; restart it any time
    python -m chronocell.jobs list
    python -m chronocell.jobs stop <id>

Resume after a restart: a job left "running" whose process is gone is put back in the queue (attempt + 1) when
the next worker starts; commands that write per-item results (batch) skip what is already done.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / ".chronocell_cache" / "jobs"
MAX_CPU = max(1, (os.cpu_count() or 4) // 4)


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _alive(pid: int | None) -> bool:
    if not pid:
        return False
    try:
        if os.name == "nt":
            out = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"], capture_output=True, text=True).stdout
            return str(pid) in out
        os.kill(pid, 0)
        return True
    except (OSError, subprocess.SubprocessError):
        return False


class Queue:
    def __init__(self, root: Path = ROOT, max_cpu: int = MAX_CPU):
        self.root = Path(root)
        self.max_cpu = max_cpu

    def _path(self, jid: str) -> Path:
        return self.root / jid / "job.json"

    def get(self, jid: str) -> dict:
        return json.loads(self._path(jid).read_text(encoding="utf-8"))

    def _put(self, job: dict) -> None:
        p = self._path(job["id"])
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(".tmp")
        tmp.write_text(json.dumps(job, indent=1), encoding="utf-8")
        tmp.replace(p)

    def submit(self, args: list[str], device: str = "cpu", label: str = "") -> dict:
        if device not in ("cpu", "gpu"):
            raise ValueError("device is 'cpu' or 'gpu'")
        jid = dt.datetime.now().strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
        out = self.root / jid / "out"
        job = {"id": jid, "label": label or " ".join(args[:2]), "args": list(args), "device": device, "state": "queued",
               "created": _now(), "started": None, "finished": None, "pid": None, "attempt": 0, "returncode": None,
               "out": str(out)}
        self._put(job)
        return job

    def jobs(self) -> list[dict]:
        if not self.root.exists():
            return []
        out = []
        for p in sorted(self.root.glob("*/job.json")):
            try:
                out.append(json.loads(p.read_text(encoding="utf-8")))
            except (OSError, ValueError):
                continue
        return out

    def stop(self, jid: str) -> dict:
        job = self.get(jid)
        (self.root / jid / "STOP").write_text(_now(), encoding="utf-8")
        if job["state"] == "queued":
            job.update(state="stopped", finished=_now())
            self._put(job)
        return job

    def log(self, jid: str, tail: int = 40) -> str:
        p = self.root / jid / "log.txt"
        return "\n".join(p.read_text(encoding="utf-8", errors="replace").splitlines()[-tail:]) if p.exists() else ""

    def recover(self) -> list[str]:
        """Jobs left running by a worker that is gone go back to the queue."""
        back = []
        for job in self.jobs():
            if job["state"] == "running" and not _alive(job.get("pid")):
                job.update(state="queued", pid=None, attempt=job["attempt"] + 1)
                self._put(job)
                back.append(job["id"])
        return back

    def _launch(self, job: dict) -> subprocess.Popen:
        d = self.root / job["id"]
        args = list(job["args"])
        if "--out" not in args and args and args[0] in ("analyze", "diff", "impact", "batch"):
            args += ["--out", job["out"]]
        log = (d / "log.txt").open("a", encoding="utf-8")
        log.write(f"--- attempt {job['attempt'] + 1} {_now()}: chronocell {' '.join(args)}\n")
        log.flush()
        env = dict(os.environ, PYTHONUNBUFFERED="1")
        if job["device"] == "cpu":
            env["CUDA_VISIBLE_DEVICES"] = ""
        proc = subprocess.Popen([sys.executable, "-m", "chronocell.cli"] + args, stdout=log, stderr=subprocess.STDOUT,
                                cwd=str(Path(__file__).resolve().parent.parent), env=env)
        job.update(state="running", started=_now(), pid=proc.pid)
        self._put(job)
        return proc

    def step(self, running: dict) -> dict:
        """One scheduling pass: reap finished jobs, honour Stop, start what the limits allow."""
        for jid, proc in list(running.items()):
            stop = (self.root / jid / "STOP").exists()
            if stop and proc.poll() is None:
                proc.terminate()
                proc.wait(timeout=30)
            rc = proc.poll()
            if rc is not None:
                job = self.get(jid)
                job.update(state="stopped" if stop else ("done" if rc == 0 else "failed"), finished=_now(), returncode=rc,
                           pid=None)
                self._put(job)
                del running[jid]
        busy_gpu = any(self.get(j)["device"] == "gpu" for j in running)
        busy_cpu = sum(self.get(j)["device"] == "cpu" for j in running)
        for job in sorted((j for j in self.jobs() if j["state"] == "queued"), key=lambda j: j["created"]):
            if (self.root / job["id"] / "STOP").exists():
                job.update(state="stopped", finished=_now())
                self._put(job)
                continue
            if job["device"] == "gpu" and not busy_gpu:
                running[job["id"]] = self._launch(job)
                busy_gpu = True
            elif job["device"] == "cpu" and busy_cpu < self.max_cpu:
                running[job["id"]] = self._launch(job)
                busy_cpu += 1
        return running

    def work(self, poll: float = 2.0, until_empty: bool = False) -> None:
        self.recover()
        running: dict = {}
        while True:
            running = self.step(running)
            if until_empty and not running and not any(j["state"] == "queued" for j in self.jobs()):
                return
            time.sleep(poll)


def start_worker(root: Path = ROOT) -> int:
    """Start a background worker process (used by the app's Jobs panel); returns its PID."""
    root.mkdir(parents=True, exist_ok=True)
    pidfile = root / "worker.pid"
    if pidfile.exists():
        try:
            pid = int(pidfile.read_text())
            if _alive(pid):
                return pid
        except ValueError:
            pass
    flags = 0x00000008 if os.name == "nt" else 0          # DETACHED_PROCESS on Windows
    proc = subprocess.Popen([sys.executable, "-m", "chronocell.jobs", "worker"], cwd=str(Path(__file__).resolve().parent.parent),
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, creationflags=flags)
    pidfile.write_text(str(proc.pid))
    return proc.pid


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("submit")
    s.add_argument("device", choices=("cpu", "gpu"))
    s.add_argument("args", nargs=argparse.REMAINDER)
    sub.add_parser("worker")
    sub.add_parser("list")
    p = sub.add_parser("stop")
    p.add_argument("id")
    a = ap.parse_args(argv)
    q = Queue()
    if a.cmd == "submit":
        args = a.args[1:] if a.args and a.args[0] == "--" else a.args
        print(json.dumps(q.submit(args, a.device), indent=1))
    elif a.cmd == "worker":
        q.work()
    elif a.cmd == "list":
        for j in q.jobs():
            print(f"{j['id']}  {j['device']}  {j['state']:8s} {j['label']}")
    else:
        print(json.dumps(q.stop(a.id), indent=1))


if __name__ == "__main__":
    main()
