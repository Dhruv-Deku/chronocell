# ChronoCell-5D — the explainer

A 5-minute explainer video (1920 × 1080, 30 fps, motion blur, synthesised soundtrack, captions): what the app is,
the app running live, and how its accuracy is tested. It is built like the reel next door (`../`): an HTML page whose
every frame is a function of time, rendered frame by frame in Chrome. Unlike the reel it is made of **real
recordings**: screen recordings of the running app, and terminal sessions of real test runs.

```
python motion/render.py --page explainer            # -> motion/out/ChronoCell-5D_explainer.mp4 (about 45 minutes)
python motion/render.py --page explainer --preview  # stills every 2 s + a contact sheet
open motion/explainer/index.html?play               # watch it live in a browser
```

## What is in it

| Part | Time | What is on screen | Where it comes from |
|---|---|---|---|
| Title | 0:00 | name, what the video covers | — |
| 01 What it is | 0:09 | 2 m of DNA in a nucleus; the chr22 fold draws itself; the contact map builds itself; 3D / 4D / 5D | the app's reference model and its contact list (`../build_data.py`); the Guide's and README's own wording |
| 02 The app, live | 0:41 | nine screen recordings: the eight workspaces in use, in the app's tab order (the Guide last), including a live VQE calculation checked against the exact energy and the Scoreboard's fails with their "later passed as" notes | `record_app.py`: Chrome drives the running app; every number on screen was computed during the take. A badge shows the playback speed when a take is sped up. |
| 03 How it is tested | 2:47 | the four steps (practice data, frozen rule, run once, keep the result); the test suite running; `recheck.py` recomputing the scores | `record_terminal.py`: real output with its real timing (the suite is shown about 60× faster; the re-check is paced for reading) |
| 04 Where it stands | 4:12 | 15 passed, 16 failed (7 later passed as a retest, 9 still open), 7 mixed or blocked; the highlights and the honest fails | the Scoreboard (`ui/scoreboard.rows()`) and the result files |
| The pass / fail table | 4:26 | all 38 tests in four columns: passed; failed, then passed as a retest (each lights up with its retest); failed and still open; no pass / fail verdict | the Scoreboard rows; a short name and result per test in `../build_data.py` (`TABLE`), whose every number `tests/test_round3.py` checks against the row |
| Outro | 4:46 | the three commands to run the app | README |

**No accuracy test is re-run.** Each of the 38 tests was pre-registered and run once on held-out data; re-running it
would break that. `recheck.py` instead recomputes each score from the per-case numbers saved in
`validation/results_*.json` (energy errors from energies, F1 from loop counts, docking success from pose RMSDs, ...),
re-applies the frozen rule from `validation/frozen.py`, and checks the verdict against the Scoreboard
(`tests/test_round3.py` runs it too). The only edit to recorded terminal output is cosmetic: the home folder is shown
as `~`.

## Files

| File | What it does |
|---|---|
| `record_app.py` | screen-records the running app (Chrome screencast, a drawn cursor, logged clicks) → `clips/*.mp4` + `clips/*.json` |
| `record_terminal.py` | runs a command and keeps every output line with its time → `data/term_*.json` |
| `recheck.py` | recomputes every score it can from the saved results and checks each verdict |
| `build.py` | packs the clip lengths, click times and terminal sessions into `data/explainer_data.js` |
| `index.html`, `explainer.css`, `explainer.js` | the video: scenes, captions, screen frames, terminal replays, the chapter tags |

## Rebuilding

```
streamlit run app.py --server.port 8601 --server.headless true    # in another terminal
python motion/explainer/record_app.py                               # about 5 minutes
python motion/explainer/record_terminal.py pytest                   # the full test suite, about 17 minutes
python motion/explainer/record_terminal.py recheck
python motion/explainer/build.py
python motion/render.py --page explainer
```

The recordings (`clips/*.mp4`, about 190 MB) and the full-size video are not committed. The compact copy
(`motion/out/ChronoCell-5D_explainer_share.mp4`, about 70 MB) is, so the Guide plays the explainer straight after a
download; re-render with `--share` and commit the new copy when the explainer changes.
