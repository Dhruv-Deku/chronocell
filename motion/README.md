# ChronoCell-5D — the film

A 100-second motion-graphics reel of the app (1920 × 1080, 60 fps, motion blur, synthesised soundtrack), made from
the app itself. The film is code: an HTML page whose every frame is a pure function of time, rendered frame by frame
in Chrome and encoded with ffmpeg. Its look follows a reference reel the user supplied: full-bleed brand colours, chunky
extruded type, a cursor that clicks the UI, rolling counters, confetti, a bar race, a bento grid, word slams and a
sunburst ending. The colours are the app's own (`chronocell/theme.py`: cobalt, terracotta, ochre, violet, paper, ink).

```
python motion/render.py                 # -> motion/out/ChronoCell-5D_film.mp4 (about 20-30 minutes)
python motion/render.py --share         # -> ..._share.mp4, a smaller copy for messaging apps
python motion/render.py --preview       # stills every 2 s + a contact sheet in motion/out/preview/
open motion/index.html?play             # watch it live in a browser (?t=42.5 shows one frame)
```

## What is on screen, and where it comes from

| # | Scene | Time | Source |
|---|---|---|---|
| 01 | The workstation: the app's name pill, a click on "Open chr22" | 0–4.5 s | chromosome names from `chronocell/genome.py` |
| 02 | The fold: chr22 draws itself bead by bead; 50.82 Mb, 5,082 beads, 10 kb, R<sub>g</sub> 580.9 nm | 4.5–14.5 s | the app's reference model of chr22 (seed 7), its own coordinates and measurements (`build_data.py`) |
| 03 | The contacts: pixels fly in to build the contact map; 70,466 contact pairs | 14.5–20 s | the reference model's own contact list (simulated Micro-C), summed into 128 × 128 squares; log scale |
| — | Title | 20–24.5 s | the app's logo mark |
| 04 | 3D structure: the model turning; R<sub>g</sub>, ν, span, contour length | 24.5–29 s | the 3D workspace's measurements (`chronocell/agent.py compute_metrics`) |
| 05 | 4D dynamics: the 22q11.2 deletion playing, largest displacement 368 nm | 29–34 s | the 4D workspace's 22q11.2 preset, simulated with the app's defaults (24 frames × 15 sweeps) |
| 06–09 | Compare, Drug lab, Genes + Guide, Quantum lab | 34–52.5 s | real screenshots of the running app (`capture_app.py`); each sentence is the app's own one-line purpose |
| 10 | Everything: six feature tiles | 52.5–58 s | the same data and screenshots |
| 11 | The scoreboard: 37 tests burst out as a network, then sort into 14 passed, 16 failed, 7 mixed | 58–67 s | `ui/scoreboard.rows()` (read from `validation/results_*.json`), one node per test, labelled with its gate |
| 12 | DNA loops: a bar race; ChronoCell first on both new cell types | 67–72.5 s | `results_gate6b.json` (F1 vs ENCODE loops) |
| 13 | Molecules: 14 / 14 within chemical accuracy, worst 1.16 mHa | 72.5–76.5 s | `results_round2.json` (Gate Q6d) |
| 14 | Noisy chips: 16 molecules fall from noisy to mitigated error | 76.5–80 s | `results_mitigation.json` (Gate Q9), per molecule |
| 15 | The honest fail: quantum docking 35 % vs random search 58 % | 80–84 s | `results_round2.json` (Gate Q7b) |
| 16 | In one word: Fold. Time. Dose. Genes. Qubits. Proof. | 84–90 s | — |
| — | Outro: 37 tests · 14 passed · 16 failed · 7 mixed | 90–100 s | the Scoreboard |

Nothing on screen is typed in: `build_data.py` writes `data/data.js` from the app and the result files, and the film
reads it (`tests/test_round3.py` checks the numbers against the Scoreboard). Re-run it after a new test and the film
follows. The bars grow at one speed in the loop race, so the order you see emerge is the order of the final scores.

## Files

| File | What it does |
|---|---|
| `build_data.py` | chr22 model, its measurements and contact map, 22q11.2 trajectory, Scoreboard rows, quoted results → `data/data.js` |
| `capture_app.py` | drives the running app in Chrome and saves the screenshots in `assets/app/` |
| `index.html`, `film.css`, `film.js` | the film: scenes, cursor, counters, confetti, transitions, the chapter badges, the cue list for the sound |
| `gl.js` | the 3D layer (three.js r128): the chr22 tube that draws itself, and the 22q11.2 deletion as a moving tube |
| `sound.py` | the soundtrack, synthesised from the film's cue list (120 bpm, A minor; no samples) |
| `render.py` | frame-by-frame capture (Playwright + Chrome) piped into ffmpeg; motion blur; muxes the soundtrack |
| `make_ui_images.py` | the Guide's banner, the film poster and the tour pictures in `docs/images/` |
| `vendor/` | three.js r128 (MIT) |

## Rebuilding everything

```
pip install playwright imageio-ffmpeg       # uses the installed Google Chrome; ffmpeg comes with imageio-ffmpeg
streamlit run app.py --server.port 8601 --server.headless true   # for the screenshots only
python motion/capture_app.py
python motion/build_data.py
python motion/render.py
python motion/render.py --share
python motion/make_ui_images.py
```

The rendered film (`out/`) is not committed: it is about 100 MB and is rebuilt by `render.py`.
