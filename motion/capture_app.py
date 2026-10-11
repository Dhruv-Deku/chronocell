"""
Screenshots of the running app for the film and the Guide (motion/assets/app/*.png, 3200 x 1800).

    streamlit run app.py --server.port 8601 --server.headless true      # in another terminal
    python motion/capture_app.py [--only ws01,ws07] [--port 8601]

Each shot is a real page of the app driven by Chrome (Playwright): a workspace opened, buttons pressed, the page
scrolled to a heading or chart. Nothing is drawn on top. Streamlit's header bar is hidden; the sidebar is hidden in
the close-ups and shown in the wide shots.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from playwright.sync_api import Page, sync_playwright

OUT = Path(__file__).resolve().parent / "assets" / "app"
HIDE_HEADER = 'header[data-testid="stHeader"]{display:none!important} [data-testid="stToolbar"]{display:none!important}'
HIDE_SIDEBAR = 'section[data-testid="stSidebar"]{display:none!important}'


def idle(pg: Page, extra: float = 1.5, limit: float = 600) -> None:
    """Wait until Streamlit has finished running the script and every spinner is gone."""
    time.sleep(0.8)
    t0 = time.time()
    while time.time() - t0 < limit:
        busy = pg.locator('[data-testid="stStatusWidget"]').count() + pg.locator('[data-testid="stSpinner"]').count()
        if not busy:
            break
        time.sleep(0.4)
    time.sleep(extra)


def style(pg: Page, sidebar: bool) -> None:
    pg.evaluate("() => document.querySelectorAll('style[data-cc-shot]').forEach(s => s.remove())")
    css = HIDE_HEADER + ("" if sidebar else HIDE_SIDEBAR)
    pg.evaluate("css => { const s = document.createElement('style'); s.dataset.ccShot = '1'; s.textContent = css; "
                "document.head.appendChild(s); }", css)


def scroller(pg: Page) -> None:
    pg.evaluate("""() => { const c = [...document.querySelectorAll('*')].filter(e => e.scrollHeight > e.clientHeight + 50
        && getComputedStyle(e).overflowY.match(/auto|scroll/)); c.sort((a, b) => b.clientHeight - a.clientHeight);
        window.__sc = c[0] || document.scrollingElement; }""")


def scroll_y(pg: Page, y: float) -> None:
    scroller(pg)
    pg.evaluate(f"() => window.__sc.scrollTo(0, {y})")
    time.sleep(0.9)


def scroll_to(pg: Page, locator, offset: int = 90) -> None:
    """Scroll the main column so the element's top sits `offset` px below the window top."""
    scroller(pg)
    locator.first.scroll_into_view_if_needed()
    box = locator.first.bounding_box()
    if box:
        pg.evaluate(f"() => window.__sc.scrollBy(0, {box['y'] - offset})")
    time.sleep(1.0)


def workspace(pg: Page, label: str) -> None:
    pg.locator(".st-key-seg_workspace button").filter(has_text=label).first.click()
    idle(pg, 2.5)


def shot(pg: Page, name: str, el=None, pad: int = 0) -> None:
    path = OUT / f"{name}.png"
    if el is not None:
        if not el.first.is_visible() or (el.first.bounding_box() or {}).get("height", 0) < 40:
            print("   (hidden, skipped)", name, flush=True)          # a chart inside a closed expander or popover
            return
        el.first.scroll_into_view_if_needed()
        time.sleep(0.6)
        if pad:
            b = el.first.bounding_box()
            pg.screenshot(path=str(path), clip={"x": max(b["x"] - pad, 0), "y": max(b["y"] - pad, 0),
                                                 "width": b["width"] + 2 * pad, "height": b["height"] + 2 * pad})
        else:
            el.first.screenshot(path=str(path))
    else:
        pg.screenshot(path=str(path))
    print("  ", path.name, flush=True)


def button(pg: Page, text: str):
    return pg.locator("button").filter(has_text=text)


def plots(pg: Page):
    return pg.locator(".js-plotly-plot")


# --------------------------------------------------------------------------------------------- the shots
def ws01(pg: Page) -> None:
    workspace(pg, "3D structure")
    style(pg, sidebar=True)
    scroll_y(pg, 0)
    shot(pg, "ws01_wide")
    style(pg, sidebar=False)
    scroll_y(pg, 0)
    shot(pg, "ws01_top")
    scroll_to(pg, plots(pg), 40)
    shot(pg, "ws01_fold")
    shot(pg, "fold3d", plots(pg))
    n = plots(pg).count()
    for i in range(1, n):
        shot(pg, f"ws01_plot{i}", plots(pg).nth(i))


def ws02(pg: Page) -> None:
    workspace(pg, "4D dynamics")
    style(pg, sidebar=False)
    scroll_y(pg, 0)
    shot(pg, "ws02_top")
    scroll_to(pg, plots(pg), 40)
    shot(pg, "ws02_fold")
    shot(pg, "sv3d", plots(pg))
    for i in range(1, plots(pg).count()):
        shot(pg, f"ws02_plot{i}", plots(pg).nth(i))


def ws03(pg: Page) -> None:
    style(pg, sidebar=True)
    demo = pg.get_by_text("Load demo patients (synthetic)")
    if demo.count():
        demo.first.click()
        idle(pg, 3)
    workspace(pg, "Compare")
    style(pg, sidebar=False)
    scroll_y(pg, 0)
    shot(pg, "ws03_top")
    if plots(pg).count():
        scroll_to(pg, plots(pg), 40)
        shot(pg, "ws03_fold")
        for i in range(plots(pg).count()):
            shot(pg, f"ws03_plot{i}", plots(pg).nth(i))


def ws04(pg: Page) -> None:
    workspace(pg, "Drug lab")
    style(pg, sidebar=False)
    scroll_y(pg, 0)
    shot(pg, "ws04_top")
    for i in range(plots(pg).count()):
        shot(pg, f"ws04_plot{i}", plots(pg).nth(i))
    scroll_to(pg, plots(pg).nth(0), 60)
    shot(pg, "ws04_mid")


def ws05(pg: Page) -> None:
    workspace(pg, "Genes")
    style(pg, sidebar=False)
    scroll_y(pg, 0)
    shot(pg, "ws05_top")
    if plots(pg).count():
        scroll_to(pg, plots(pg), 40)
        shot(pg, "ws05_fold")
        shot(pg, "genes3d", plots(pg))


def ws06(pg: Page) -> None:
    workspace(pg, "Guide")
    style(pg, sidebar=False)
    scroll_y(pg, 0)
    shot(pg, "ws06_top")
    for name, text in (("ws06_pages", "What each page does"), ("ws06_quantum", "Quantum lab, in plain words"),
                       ("ws06_new", "New: the scoreboard")):
        loc = pg.get_by_text(text)
        if loc.count():
            scroll_to(pg, loc, 60)
            shot(pg, name)


def _qtab(pg: Page, label: str) -> None:
    pg.locator("button").filter(has_text=label).first.click()
    idle(pg, 2)


def ws07(pg: Page) -> None:
    workspace(pg, "Quantum lab")
    style(pg, sidebar=False)
    scroll_y(pg, 0)
    shot(pg, "ws07_top")
    runs = (("TAD boundaries", "q_tad"), ("Molecule (VQE)", "q_vqe"), ("Noise & mitigation", "q_noise"),
            ("Lattice fold", "q_lattice"), ("Quantum walk", "q_walk"), ("Drug combination", "q_combo"))
    for tab, name in runs:
        _qtab(pg, tab)
        go = pg.locator("button").filter(has_text="Run").filter(has_text="simulator")
        if not go.count():
            go = pg.locator('button[kind="primary"], button[data-testid="stBaseButton-primary"]')
        if go.count():
            go.first.click()
            idle(pg, 3, limit=900)
        n = plots(pg).count()
        print(f"   {tab}: {n} charts", flush=True)
        if n:
            scroll_to(pg, plots(pg), 60)
            shot(pg, f"{name}_view")
            for i in range(min(n, 4)):
                shot(pg, f"{name}_plot{i}", plots(pg).nth(i))


def ws08(pg: Page) -> None:
    workspace(pg, "Scoreboard")
    style(pg, sidebar=False)
    scroll_y(pg, 0)
    shot(pg, "ws08_top")
    tabs = ("Before → after", "Molecules (quantum chemistry)", "DNA loops", "Drug safety (ADMET)", "Error mitigation",
            "Docking")
    for k, tab in enumerate(tabs):
        t = pg.get_by_role("tab", name=tab)
        if t.count():
            t.first.click()
            time.sleep(2.0)
            vis = pg.locator('[role="tabpanel"]:not([hidden]) .js-plotly-plot')
            if vis.count():
                shot(pg, f"sb_{k}", vis.first, pad=8)
    t = pg.get_by_role("tab", name="All tests")
    if t.count():
        t.first.click()
        time.sleep(1.5)
        scroll_to(pg, pg.get_by_role("tab", name="All tests"), 60)
        shot(pg, "ws08_list")


SHOTS = {"ws01": ws01, "ws02": ws02, "ws03": ws03, "ws04": ws04, "ws05": ws05, "ws06": ws06, "ws07": ws07, "ws08": ws08}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="", help="comma-separated shot groups, e.g. ws01,ws07")
    ap.add_argument("--port", type=int, default=8601)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    todo = [k for k in SHOTS if not a.only or k in a.only.split(",")]
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", args=["--enable-gpu", "--ignore-gpu-blocklist"])
        pg = b.new_page(viewport={"width": 1600, "height": 900}, device_scale_factor=2)
        pg.goto(f"http://localhost:{a.port}", timeout=180_000)
        pg.wait_for_selector("text=ChronoCell-5D", timeout=180_000)
        idle(pg, 3)
        for k in todo:
            t = time.time()
            print(k, flush=True)
            SHOTS[k](pg)
            print(f"   {time.time() - t:.0f} s", flush=True)
        b.close()


if __name__ == "__main__":
    main()
