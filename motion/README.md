# ChronoCell-5D — the film

A 92-second motion-graphics film of the app (1920 × 1080, 60 fps, with a synthesised soundtrack), made from the app
itself. The film is code: an HTML page whose every frame is a pure function of time, rendered frame by frame in
Chrome and encoded with ffmpeg.

```
python motion/render.py                 # -> motion/out/ChronoCell-5D_film.mp4 (about 12 minutes)
python motion/render.py --preview       # stills every 2 s + a contact sheet in motion/out/preview/
open motion/index.html?play             # watch it live in a browser (?t=42.5 shows one frame)
```

## What is on screen, and where it comes from

| Scene | Time | Source |
|---|---|---|
| DNA particles form a double helix; "2 metres of DNA", "10 micrometres" | 0–7.5 s | artwork; textbook figures |
| The helix unwinds and the fibre folds, bead by bead, into chr22 | 7.5–16 s | the app's reference model of chr22 (seed 7), its own coordinates (`build_data.py`); the chips (5,082 beads, 10 kb, 50.82 Mb, 70,466 contacts) are read from it |
| Title | 16–21 s | the app's logo mark |
| The eight workspaces | 21–48 s | real screenshots of the running app (`capture_app.py`); each sentence is the app's own one-line purpose (`app.py` `WORKSPACES`) |
| 4D dynamics | 24.4–27.8 s | the 4D workspace's 22q11.2 deletion preset simulated with the app's defaults (24 frames); colour = how far each bead moved |
| Mosaic | 48–52 s | the same screenshots |
| "Every claim is a test" | 52–57 s | the project's test rules (`validation/frozen.py`) |
| 37 tests: 14 passed, 16 failed, 7 mixed or blocked | 57–63 s | the Scoreboard's rows (`ui/scoreboard.rows()`, read from `validation/results_*.json`), one dot per test |
| DNA loops (Gate 6b), molecules (Gate Q6d), error mitigation (Gate Q9) | 63–75 s | `results_gate6b.json`, `results_round2.json`, `results_mitigation.json` |
| The honest fail: docking (Gate Q7b) | 75–78 s | `results_round2.json` |
| Feature burst | 78–84 s | screenshots |
| Outro | 84–92 s | the chr22 model again; the test count from the Scoreboard |

Nothing on screen is typed in: `build_data.py` writes `data/data.js` from the app and the result files, and the film
reads it. Re-run it after a new test and the numbers in the film follow.

## Files

| File | What it does |
|---|---|
| `build_data.py` | chr22 model, 22q11.2 trajectory, Scoreboard rows and quoted results → `data/data.js` |
| `capture_app.py` | drives the running app in Chrome and saves the screenshots in `assets/app/` |
| `index.html`, `film.css`, `film.js` | the film's 2D layer: titles, app windows, charts, transitions, the cue list for the sound |
| `gl.js` | the 3D layer (three.js r128 + bloom): dust, helix, the chr22 fold, the 22q11.2 deletion |
| `sound.py` | the soundtrack, synthesised from the film's cue list (120 bpm, A minor; no samples) |
| `render.py` | frame-by-frame capture (Playwright + Chrome) piped into ffmpeg; muxes the soundtrack |
| `make_ui_images.py` | the Guide's banner, the film poster and the tour pictures in `docs/images/` |
| `vendor/` | three.js r128 and its bloom passes (MIT) |

## Rebuilding everything

```
pip install playwright imageio-ffmpeg       # uses the installed Google Chrome; ffmpeg comes with imageio-ffmpeg
streamlit run app.py --server.port 8601 --server.headless true   # for the screenshots only
python motion/capture_app.py
python motion/build_data.py
python motion/make_ui_images.py
python motion/render.py
```

The rendered film (`out/`) is not committed: it is about 100 MB and is rebuilt by the last command.
