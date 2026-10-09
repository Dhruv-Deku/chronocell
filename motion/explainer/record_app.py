"""
Live screen recordings of the running app for the explainer (motion/explainer/clips/*.mp4).

    streamlit run app.py --server.port 8601 --server.headless true      # in another terminal
    python motion/explainer/record_app.py [--only ws01,ws07] [--port 8601]

Chrome (Playwright) drives the real app the way a person would: the mouse glides, clicks, drags the 3D views, scrolls,
moves sliders and presses Run. Chrome's own screencast (CDP Page.startScreencast) records what is on screen, frame by
frame with timestamps, and ffmpeg turns each take into a 30 fps clip. A drawn cursor (headless Chrome has none) follows
the real mouse; clicks are logged so the explainer can add a click sound. Nothing is staged: every number on screen is
computed by the app during the take.
"""

from __future__ import annotations

import argparse
import base64
import json
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

from playwright.sync_api import Page, sync_playwright

HERE = Path(__file__).resolve().parent
OUT = HERE / "clips"
sys.path.insert(0, str(HERE.parent))
from render import ffmpeg_exe  # noqa: E402

VW, VH, DSF = 1600, 900, 1.2                    # 1920 x 1080 device pixels
HIDE_HEADER = 'header[data-testid="stHeader"]{display:none!important} [data-testid="stToolbar"]{display:none!important}'
HIDE_SIDEBAR = 'section[data-testid="stSidebar"]{display:none!important}'
CURSOR_JS = """() => {
  if (window.__cur) return;
  const c = document.createElement('div');
  c.innerHTML = '<svg width="30" height="38" viewBox="0 0 23 29"><path d="M2 1.5 L2 23 L7.6 18 L11.4 26.6 L15 25 L11.3 16.6 L18.6 16.4 Z" fill="#fff" stroke="#111" stroke-width="1.6" stroke-linejoin="round"/></svg>';
  Object.assign(c.style, {position:'fixed', left:'0', top:'0', zIndex:2147483647, pointerEvents:'none', transformOrigin:'2px 2px',
    filter:'drop-shadow(0 3px 4px rgba(0,0,0,.35))', transition:'none'});
  document.documentElement.appendChild(c);
  window.__cur = (x, y, s) => { c.style.transform = `translate(${x - 2}px, ${y - 2}px) scale(${s || 1})`; };
  window.__ripple = (x, y) => {
    const r = document.createElement('div');
    Object.assign(r.style, {position:'fixed', left:(x - 22) + 'px', top:(y - 22) + 'px', width:'44px', height:'44px', borderRadius:'50%',
      border:'3px solid #3340D1', zIndex:2147483646, pointerEvents:'none', opacity:'0.9', transform:'scale(0.3)',
      transition:'transform .45s ease-out, opacity .45s ease-out'});
    document.documentElement.appendChild(r);
    requestAnimationFrame(() => { r.style.transform = 'scale(1.6)'; r.style.opacity = '0'; });
    setTimeout(() => r.remove(), 600);
  };
}"""


class Take:
    """One recorded clip: screencast frames with timestamps, plus the clicks."""

    def __init__(self, pg: Page):
        self.pg = pg
        self.cdp = pg.context.new_cdp_session(pg)
        self.frames: list[tuple[float, bytes]] = []
        self.clicks: list[dict] = []
        self.on = False
        self.cdp.on("Page.screencastFrame", self._frame)

    def _frame(self, ev):
        if self.on:
            self.frames.append((time.time(), base64.b64decode(ev["data"])))
        try:
            self.cdp.send("Page.screencastFrameAck", {"sessionId": ev["sessionId"]})
        except Exception:
            pass

    def start(self):
        self.frames, self.clicks = [], []
        self.t0 = time.time()
        self.on = True
        self.cdp.send("Page.startScreencast", {"format": "jpeg", "quality": 92, "maxWidth": 1920, "maxHeight": 1080,
                                               "everyNthFrame": 1})
        self.pg.wait_for_timeout(300)

    def stop(self, name: str) -> Path:
        self.pg.wait_for_timeout(250)
        t1 = time.time()
        self.on = False
        self.cdp.send("Page.stopScreencast")
        tmp = OUT / f"_{name}"
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True)
        lines = []
        fr = self.frames
        for i, (t, data) in enumerate(fr):
            p = tmp / f"f{i:05d}.jpg"
            p.write_bytes(data)
            dur = (fr[i + 1][0] if i + 1 < len(fr) else t1) - t
            lines.append(f"file '{p.name}'\nduration {max(dur, 0.001):.4f}")
        lines.append(f"file 'f{len(fr) - 1:05d}.jpg'")            # concat demuxer: the last duration needs a closing entry
        (tmp / "list.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
        out = OUT / f"{name}.mp4"
        subprocess.run([ffmpeg_exe(), "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", str(tmp / "list.txt"),
                        "-vf", "fps=30,scale=1920:1080,format=yuv420p", "-c:v", "libx264", "-preset", "medium", "-crf", "16",
                        "-g", "6", "-movflags", "+faststart", str(out)], check=True)
        start = fr[0][0] if fr else self.t0
        meta = {"name": name, "secs": round(t1 - start, 3), "frames": len(fr),
                "clicks": [dict(c, t=round(c["t"] - start, 3)) for c in self.clicks]}
        (OUT / f"{name}.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
        shutil.rmtree(tmp, ignore_errors=True)
        print(f"   {out.name}: {meta['secs']:.1f} s, {len(fr)} frames, {len(self.clicks)} clicks", flush=True)
        return out


# --------------------------------------------------------------------------------------------- a person's hand
class Hand:
    def __init__(self, pg: Page, take: Take):
        self.pg, self.take = pg, take
        self.x, self.y = VW * 0.62, VH * 0.75

    def cursor(self):
        self.pg.evaluate(CURSOR_JS)
        self.pg.evaluate("([x, y]) => window.__cur(x, y)", [self.x, self.y])

    def glide(self, x: float, y: float, dur: float = 0.6, held: bool = False):
        self.pg.evaluate(CURSOR_JS)                     # idempotent: survives page reloads
        n = max(6, int(dur * 45))
        x0, y0 = self.x, self.y
        for k in range(1, n + 1):
            p = k / n
            e = p * p * (3 - 2 * p)
            xi, yi = x0 + (x - x0) * e, y0 + (y - y0) * e
            self.pg.mouse.move(xi, yi)
            self.pg.evaluate("([x, y]) => window.__cur && window.__cur(x, y)", [xi, yi])
            self.pg.wait_for_timeout(dur * 1000 / n)
        self.x, self.y = x, y

    def to(self, loc, dur: float = 0.6, dx: float = 0.5, dy: float = 0.5) -> bool:
        if not loc.count():
            return False
        loc.first.scroll_into_view_if_needed()
        b = loc.first.bounding_box()
        if not b:
            return False
        self.glide(b["x"] + b["width"] * dx, b["y"] + b["height"] * dy, dur)
        return True

    def click(self, loc=None, dur: float = 0.6, dx: float = 0.5, dy: float = 0.5, after: float = 0.4) -> bool:
        if loc is not None and not self.to(loc, dur, dx, dy):
            return False
        self.pg.evaluate(CURSOR_JS)
        self.take.clicks.append({"t": time.time(), "x": round(self.x * DSF), "y": round(self.y * DSF)})
        self.pg.evaluate("([x, y]) => { window.__cur(x, y, 0.85); window.__ripple(x, y); }", [self.x, self.y])
        self.pg.mouse.down()
        self.pg.wait_for_timeout(90)
        self.pg.mouse.up()
        self.pg.evaluate("([x, y]) => window.__cur(x, y, 1)", [self.x, self.y])
        self.pg.wait_for_timeout(after * 1000)
        return True

    def drag(self, dx: float, dy: float, dur: float = 1.6):
        self.pg.mouse.down()
        self.glide(self.x + dx, self.y + dy, dur, held=True)
        self.pg.mouse.up()

    def wheel(self, total: float, dur: float = 1.5):
        n = max(4, int(dur * 30))
        for _ in range(n):
            self.pg.mouse.wheel(0, total / n)
            self.pg.wait_for_timeout(dur * 1000 / n)

    def wait(self, s: float):
        self.pg.wait_for_timeout(s * 1000)


# --------------------------------------------------------------------------------------------- app helpers
def idle(pg: Page, extra: float = 0.8, limit: float = 300):
    pg.wait_for_timeout(500)
    t0 = time.time()
    while time.time() - t0 < limit:
        if not (pg.locator('[data-testid="stStatusWidget"]').count() + pg.locator('[data-testid="stSpinner"]').count()):
            break
        pg.wait_for_timeout(300)
    pg.wait_for_timeout(extra * 1000)


def style(pg: Page, sidebar: bool):
    css = HIDE_HEADER + ("" if sidebar else HIDE_SIDEBAR)
    pg.evaluate("css => { document.querySelectorAll('style[data-cc-shot]').forEach(s => s.remove());"
                " const s = document.createElement('style'); s.dataset.ccShot = '1'; s.textContent = css; document.head.appendChild(s); }", css)


def scroll_top(pg: Page):
    pg.evaluate("""() => { const c = [...document.querySelectorAll('*')].filter(e => e.scrollHeight > e.clientHeight + 50
        && getComputedStyle(e).overflowY.match(/auto|scroll/)); c.sort((a, b) => b.clientHeight - a.clientHeight);
        (c[0] || document.scrollingElement).scrollTo(0, 0); }""")
    pg.wait_for_timeout(400)


def ws(pg: Page, label: str):
    return pg.locator(".st-key-seg_workspace button").filter(has_text=label)


def plots(pg: Page):
    return pg.locator(".js-plotly-plot")


def plotly_button(pg: Page, text: str):
    return pg.locator(".js-plotly-plot .updatemenu-button").filter(has_text=text)


def btn(pg: Page, text: str):
    return pg.locator("button").filter(has_text=text)


# --------------------------------------------------------------------------------------------- the takes
def open_ws(pg: Page, h: Hand, label: str, sidebar: bool = False):
    h.click(ws(pg, label), 0.7, after=0.2)
    idle(pg, 1.0)
    style(pg, sidebar)
    h.cursor()


def t_launch(pg: Page, h: Hand, take: Take):
    style(pg, True)
    scroll_top(pg)
    h.cursor()
    take.start()
    h.wait(0.8)
    for label in ("3D structure", "4D dynamics", "Compare", "Drug lab", "Genes", "Guide", "Quantum lab", "Scoreboard"):
        h.to(ws(pg, label), 0.38)
        h.wait(0.12)
    h.to(ws(pg, "3D structure"), 0.6)
    h.wait(0.6)
    take.stop("c01_launch")


def t_ws01(pg: Page, h: Hand, take: Take):
    open_ws(pg, h, "3D structure")
    scroll_top(pg)
    take.start()
    h.wait(0.6)
    p = plots(pg)
    if p.count():
        p.first.scroll_into_view_if_needed()
        h.wait(0.5)
        h.to(p.first, 0.8, 0.45, 0.45)
        h.drag(380, 40, 2.2)
        h.wait(0.3)
        h.click(plotly_button(pg, "Turntable"), 0.7, after=3.2)
        h.click(plotly_button(pg, "Pause"), 0.4, after=0.4)
    h.glide(VW * 0.82, VH * 0.5, 0.6)
    h.wheel(700, 1.6)
    h.wait(1.4)
    take.stop("c02_structure")


def t_ws02(pg: Page, h: Hand, take: Take):
    open_ws(pg, h, "4D dynamics")
    scroll_top(pg)
    take.start()
    h.wait(0.5)
    p = plots(pg)
    if p.count():
        p.first.scroll_into_view_if_needed()
        h.wait(0.5)
        h.click(plotly_button(pg, "Play"), 0.8, after=3.6)
        h.to(p.first, 0.6, 0.45, 0.4)
        h.drag(-300, 30, 1.8)
    h.wait(0.8)
    take.stop("c03_dynamics")


def t_ws03(pg: Page, h: Hand, take: Take):
    style(pg, True)
    demo = pg.get_by_text("Load demo patients (synthetic)")
    if demo.count():
        demo.first.click()
        idle(pg, 2)
    open_ws(pg, h, "Compare")
    scroll_top(pg)
    take.start()
    h.wait(0.6)
    fr = pg.locator("iframe")
    target = fr if fr.count() else plots(pg)
    if target.count():
        target.first.scroll_into_view_if_needed()
        h.wait(0.6)
        h.to(target.first, 0.8, 0.25, 0.5)
        h.drag(260, 20, 2.0)
        h.wait(0.6)
        h.to(target.first, 0.6, 0.75, 0.5)
        h.drag(-240, -20, 1.8)
    h.wait(0.8)
    take.stop("c04_compare")


def t_ws04(pg: Page, h: Hand, take: Take):
    open_ws(pg, h, "Drug lab")
    scroll_top(pg)
    take.start()
    h.wait(0.6)
    p = plots(pg)
    if p.count():
        p.first.scroll_into_view_if_needed()
        h.wait(0.6)
        h.click(plotly_button(pg, "Play"), 0.8, after=4.0)
    h.glide(VW * 0.7, VH * 0.55, 0.5)
    h.wheel(500, 1.4)
    h.wait(1.0)
    take.stop("c05_druglab")


def t_ws05(pg: Page, h: Hand, take: Take):
    open_ws(pg, h, "Genes")
    scroll_top(pg)
    take.start()
    h.wait(0.6)
    h.glide(VW * 0.6, VH * 0.6, 0.6)
    h.wheel(650, 1.8)
    h.wait(0.6)
    p = plots(pg)
    if p.count():
        h.to(p.first, 0.6, 0.45, 0.45)
        h.drag(300, 0, 1.8)
    h.wait(0.6)
    h.wheel(600, 1.6)
    h.wait(1.0)
    take.stop("c06_genes")


def t_ws06(pg: Page, h: Hand, take: Take):
    open_ws(pg, h, "Guide")
    scroll_top(pg)
    take.start()
    h.wait(1.0)
    h.glide(VW * 0.6, VH * 0.6, 0.5)
    for _ in range(3):
        h.wheel(520, 1.5)
        h.wait(0.9)
    take.stop("c07_guide")


def t_ws07(pg: Page, h: Hand, take: Take):
    open_ws(pg, h, "Quantum lab")
    scroll_top(pg)
    h.click(btn(pg, "Molecule (VQE)"), 0.5, after=0.2)
    idle(pg, 1.0)
    lab = pg.get_by_text("Bond length (Å)")
    if lab.count():                                  # the molecule controls near the top of the window
        lab.first.scroll_into_view_if_needed()
        b = lab.first.bounding_box()
        if b:
            pg.mouse.move(VW * 0.5, VH * 0.5)
            pg.mouse.wheel(0, b["y"] - 120)
            pg.wait_for_timeout(600)
    h.cursor()
    take.start()
    h.wait(0.6)
    h.click(btn(pg, "HeH+"), 0.8, after=0.3)
    idle(pg, 0.4)
    h.cursor()
    sl = pg.locator('[role="slider"]:visible')
    if sl.count():                                   # drag the bond length out a little
        h.to(sl.last, 0.7)
        h.drag(70, 0, 0.9)
        idle(pg, 0.4)
        h.cursor()
    h.click(btn(pg, "Run VQE"), 0.8, after=0.2)
    idle(pg, 2.2, limit=180)
    h.cursor()
    h.glide(VW * 0.86, VH * 0.7, 0.8)
    h.wait(1.6)
    h.click(pg.get_by_text("Add hardware noise"), 0.8, after=0.3)
    idle(pg, 0.4)
    h.cursor()
    h.click(btn(pg, "Run VQE"), 0.7, after=0.2)
    idle(pg, 2.2, limit=180)
    h.cursor()
    h.glide(VW * 0.6, VH * 0.7, 0.5)
    h.wheel(260, 0.8)
    h.glide(VW * 0.5, VH * 0.8, 0.8)
    h.wait(2.0)
    take.stop("c08_quantum")


def t_ws08(pg: Page, h: Hand, take: Take):
    open_ws(pg, h, "Scoreboard")
    scroll_top(pg)
    take.start()
    h.wait(1.0)
    h.glide(VW * 0.6, VH * 0.6, 0.5)
    h.wheel(560, 1.6)
    h.wait(0.6)
    for tab in ("Molecules (quantum chemistry)", "DNA loops", "Error mitigation", "Docking"):
        t = pg.get_by_role("tab", name=tab)
        if t.count():
            h.click(t, 0.6, after=1.6)
    t = pg.get_by_role("tab", name="All tests")
    if t.count():
        h.click(t, 0.6, after=0.8)
        h.glide(VW * 0.6, VH * 0.6, 0.4)
        h.wheel(900, 2.2)
    h.wait(1.0)
    take.stop("c09_scoreboard")


TAKES = {"launch": t_launch, "ws01": t_ws01, "ws02": t_ws02, "ws03": t_ws03, "ws04": t_ws04, "ws05": t_ws05,
         "ws06": t_ws06, "ws07": t_ws07, "ws08": t_ws08}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="")
    ap.add_argument("--port", type=int, default=8601)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    todo = [k for k in TAKES if not a.only or k in a.only.split(",")]
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", args=["--enable-gpu", "--ignore-gpu-blocklist", "--hide-scrollbars"])
        pg = b.new_page(viewport={"width": VW, "height": VH}, device_scale_factor=DSF)
        pg.goto(f"http://localhost:{a.port}", timeout=180_000)
        pg.wait_for_selector("text=ChronoCell-5D", timeout=180_000)
        idle(pg, 2.5)
        take = Take(pg)
        hand = Hand(pg, take)
        for k in todo:
            print(k, flush=True)
            TAKES[k](pg, hand, take)
        b.close()


if __name__ == "__main__":
    main()
