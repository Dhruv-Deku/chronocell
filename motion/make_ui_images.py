"""
Pictures for the app's Guide and the README, made from the film and the app's own screenshots.

    python motion/make_ui_images.py

- docs/images/hero.jpg          the chr22 model from the film's 3D layer (no text), cropped to a banner
- docs/images/film_poster.jpg   the film's title frame (poster for the video player)
- docs/images/tour/wsNN.jpg     each workspace, from motion/assets/app (motion/capture_app.py), 1280 x 720
"""

from __future__ import annotations

from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
DOCS = HERE.parent / "docs" / "images"
APP = HERE / "assets" / "app"
TOUR = {"ws01": "ws01_wide", "ws02": "ws02_fold", "ws03": "ws03_top", "ws04": "ws04_top", "ws05": "ws05_fold",
        "ws06": "q_lattice_view", "ws07": "ws08_top", "ws08": "ws06_top"}     # app order: 06 Quantum lab, 07 Scoreboard, 08 Guide


def stills() -> None:
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        b = p.chromium.launch(channel="chrome", args=["--enable-gpu", "--ignore-gpu-blocklist", "--use-angle=d3d11"])
        for name, query, crop in (("hero", "bare&t=12.6", (0, 140, 1920, 880)), ("film_poster", "t=22.6", None)):
            pg = b.new_page(viewport={"width": 1920, "height": 1080})
            pg.goto((HERE / "index.html").as_uri() + "?" + query, timeout=120_000)
            pg.wait_for_function("window.READY !== undefined")
            pg.evaluate("async () => { await window.READY; }")
            pg.wait_for_timeout(800)
            tmp = DOCS / f"_{name}.png"
            pg.screenshot(path=str(tmp))
            im = Image.open(tmp).convert("RGB")
            if crop:
                im = im.crop(crop)
            dest = DOCS / f"{name}.jpg"
            im.save(dest, quality=88, optimize=True, progressive=True)
            tmp.unlink()
            pg.close()
            print("  ", dest, im.size)
        b.close()


def tour() -> None:
    out = DOCS / "tour"
    out.mkdir(parents=True, exist_ok=True)
    for k, src in TOUR.items():
        im = Image.open(APP / f"{src}.png").convert("RGB")
        im = im.resize((1280, round(1280 * im.height / im.width)), Image.LANCZOS)
        im.save(out / f"{k}.jpg", quality=84, optimize=True, progressive=True)
        print("  ", out / f"{k}.jpg", im.size)


if __name__ == "__main__":
    DOCS.mkdir(parents=True, exist_ok=True)
    stills()
    tour()
