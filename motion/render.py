"""
Render the film (motion/index.html) to MP4, frame by frame, in Chrome.

    python motion/render.py --preview              # stills every 2 s + a contact sheet (motion/out/preview/)
    python motion/render.py --preview --at 16.2,30  # chosen stills
    python motion/render.py --share                 # a compact copy (..._share.mp4) for messaging apps
    python motion/render.py                         # motion/out/ChronoCell-5D_film.mp4 (1920x1080, 60 fps, sound,
                                                    #   motion blur: 120 drawn frames per second, averaged in pairs)
    python motion/render.py --page explainer        # motion/out/ChronoCell-5D_explainer.mp4 (30 fps, motion blur)

The page draws any time t on request (window.seek), so every frame is exact: no screen recording, no dropped
frames. Frames are piped straight into ffmpeg (H.264, yuv420p, faststart). The soundtrack (motion/sound.py) is
synthesised from the film's own cue list (window.CUES) and muxed in.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
PAGES = {"film": (HERE / "index.html", "ChronoCell-5D_film.mp4", 60),
         "explainer": (HERE / "explainer" / "index.html", "ChronoCell-5D_explainer.mp4", 30)}


def ffmpeg_exe() -> str:
    import shutil
    exe = shutil.which("ffmpeg")
    if exe:
        return exe
    import imageio_ffmpeg
    return imageio_ffmpeg.get_ffmpeg_exe()


def open_film(p, scale: float = 1.0, page: Path = PAGES["film"][0]):
    b = p.chromium.launch(channel="chrome", args=["--enable-gpu", "--ignore-gpu-blocklist", "--enable-webgl",
                                                  "--use-angle=d3d11", "--disable-gpu-vsync"])
    pg = b.new_page(viewport={"width": 1920, "height": 1080}, device_scale_factor=scale)
    errors = []
    pg.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
    pg.on("pageerror", lambda e: errors.append(str(e)))
    pg.goto(page.as_uri(), timeout=120_000)
    pg.wait_for_function("window.READY !== undefined", timeout=60_000)
    pg.evaluate("async () => { await window.READY; }")
    if errors:
        print("page errors:", *errors[:10], sep="\n  ", flush=True)
    return b, pg, errors


def cue_file(pg, stem: str) -> Path:
    """The page's sound cues (and, if it sets them, its music sections) for sound.py."""
    path = OUT / f"{stem}_cues.json"
    path.write_text(json.dumps({"duration": pg.evaluate("window.DURATION"), "cues": pg.evaluate("window.CUES"),
                                "sections": pg.evaluate("window.SECTIONS || null")}, indent=1), encoding="utf-8")
    return path


def preview(times: list[float], page: str = "film") -> None:
    from PIL import Image, ImageDraw
    from playwright.sync_api import sync_playwright
    d = OUT / "preview"
    d.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as p:
        b, pg, errors = open_film(p, page=PAGES[page][0])
        cue_file(pg, Path(PAGES[page][1]).stem)
        paths = []
        for t in times:
            pg.evaluate("t => window.seek(t)", t)
            path = d / f"t{t:06.2f}.jpg"
            pg.screenshot(path=str(path), type="jpeg", quality=88)
            paths.append((t, path))
        if errors:
            print("errors while seeking:", *errors[-5:], sep="\n  ")
        b.close()
    W, H, cols = 480, 270, 4
    rows = (len(paths) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * W, rows * (H + 18)), "black")
    dr = ImageDraw.Draw(sheet)
    for i, (t, path) in enumerate(paths):
        im = Image.open(path).resize((W, H))
        x, y = (i % cols) * W, (i // cols) * (H + 18)
        sheet.paste(im, (x, y + 18))
        dr.text((x + 4, y + 3), f"t = {t:.2f} s", fill="white")
    sheet.save(d / "sheet.jpg", quality=88)
    print("wrote", d / "sheet.jpg", f"({len(paths)} stills)")


def render(fps: int, crf: int, out: Path, t_from: float, t_to: float | None, sound: bool, blur: bool = True,
           page: str = "film") -> None:
    """blur: draw two frames per output frame and average them (motion blur over half the frame time)."""
    from playwright.sync_api import sync_playwright
    OUT.mkdir(exist_ok=True)
    with sync_playwright() as p:
        b, pg, errors = open_film(p, page=PAGES[page][0])
        dur = float(pg.evaluate("window.DURATION"))
        cues = cue_file(pg, out.stem)
        t_to = min(t_to or dur, dur)
        sub = 2 if blur else 1
        rate = fps * sub
        n0, n1 = int(round(t_from * rate)), int(round(t_to * rate))
        video = out.with_suffix(".video.mp4") if sound else out
        vf = "scale=in_range=pc:out_range=tv:out_color_matrix=bt709,format=yuv420p"
        if blur:
            vf = f"tmix=frames=2:weights='1 1',fps={fps}," + vf
        cmd = [ffmpeg_exe(), "-y", "-loglevel", "error", "-f", "image2pipe", "-framerate", str(rate), "-c:v", "mjpeg",
               "-i", "-", "-vf", vf,
               "-c:v", "libx264", "-preset", "slow", "-crf", str(crf), "-colorspace", "bt709", "-color_primaries", "bt709",
               "-color_trc", "bt709", "-color_range", "tv", "-movflags", "+faststart", str(video)]
        ff = subprocess.Popen(cmd, stdin=subprocess.PIPE)
        t0 = time.time()
        for k, n in enumerate(range(n0, n1)):
            pg.evaluate("t => window.seek(t)", n / rate)
            ff.stdin.write(pg.screenshot(type="jpeg", quality=95))
            if k % (300 * sub) == 0:
                el = time.time() - t0
                print(f"frame {n}/{n1} (t = {n / rate:.1f} s) · {el:.0f} s elapsed · "
                      f"{(n1 - n0 - k) * el / max(k, 1):.0f} s to go", flush=True)
        ff.stdin.close()
        ff.wait()
        if errors:
            print("page errors:", *errors[-5:], sep="\n  ")
        b.close()
    print(f"video: {video} ({video.stat().st_size / 1e6:.1f} MB, {time.time() - t0:.0f} s)")
    if sound:
        sys.path.insert(0, str(HERE))
        import sound as S
        wav = OUT / f"{out.stem}_soundtrack.wav"
        S.synth(cues, wav)
        cmd = [ffmpeg_exe(), "-y", "-loglevel", "error", "-i", str(video), "-ss", f"{t_from:.3f}", "-i", str(wav),
               "-c:v", "copy", "-c:a", "aac", "-b:a", "256k", "-shortest", "-movflags", "+faststart", str(out)]
        subprocess.run(cmd, check=True)
        video.unlink()
        print(f"film: {out} ({out.stat().st_size / 1e6:.1f} MB)")


def share(src: Path) -> None:
    """A smaller copy for messaging apps and uploads (same picture size and frame rate, stronger compression)."""
    dst = src.with_name(src.stem + "_share.mp4")
    subprocess.run([ffmpeg_exe(), "-y", "-loglevel", "error", "-i", str(src), "-c:v", "libx264", "-preset", "slow",
                    "-crf", "24", "-pix_fmt", "yuv420p", "-c:a", "aac", "-b:a", "160k", "-movflags", "+faststart", str(dst)],
                   check=True)
    print(f"share copy: {dst} ({dst.stat().st_size / 1e6:.1f} MB)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--preview", action="store_true")
    ap.add_argument("--at", default="", help="comma-separated times for --preview")
    ap.add_argument("--every", type=float, default=2.0)
    ap.add_argument("--page", choices=list(PAGES), default="film", help="film (the reel) or explainer")
    ap.add_argument("--fps", type=int, default=None, help="default: 60 for the film, 30 for the explainer")
    ap.add_argument("--crf", type=int, default=18)
    ap.add_argument("--from", dest="t_from", type=float, default=0.0)
    ap.add_argument("--to", dest="t_to", type=float, default=None)
    ap.add_argument("--no-sound", action="store_true")
    ap.add_argument("--no-blur", action="store_true", help="one drawn frame per output frame (twice as fast)")
    ap.add_argument("--share", action="store_true", help="only make a compact copy of the rendered film")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    a.out = a.out or str(OUT / PAGES[a.page][1])
    a.fps = a.fps or PAGES[a.page][2]
    if a.share:
        share(Path(a.out))
        return
    if a.preview:
        if a.at:
            times = [float(x) for x in a.at.split(",")]
        else:
            times, t = [], 0.3
            while t < (100 if a.page == "film" else 300):
                times.append(round(t, 2))
                t += a.every
        preview(times, a.page)
    else:
        render(a.fps, a.crf, Path(a.out), a.t_from, a.t_to, not a.no_sound, not a.no_blur, a.page)


if __name__ == "__main__":
    main()
