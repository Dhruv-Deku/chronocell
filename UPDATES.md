# ChronoCell-5D — development log

Project: ChronoCell-5D, a contact-guided 3D/4D chromatin reconstruction and interpretation
workstation for human chromosomes (GRCh38).

Repository: `github.com/WolframNU13/Team-NU-13`
Context: HackBlitz 2.0, ALLEN Global / Global City International School, Bengaluru — team HBZ-13
("Team NU — Null Pointers"), 27 September 2026.

---

## About this document

**This log was reconstructed on 28 September 2026 from git history and filesystem timestamps. It
was not written contemporaneously.** That matters if it is ever relied on as a record, so the
evidence behind every entry is stated below rather than presented as memory.

**Evidence classes used:**

| Class | Strength | What it proves |
|---|---|---|
| Git commit timestamps (3 commits) | Strong | The listed files existed in that state at that moment |
| File modification times (mtime) | Weaker | When a file was **last written** — not when it was created |

**Known limits of the timestamps — read before citing any of them:**

- **mtime is last-write, not creation.** A file created at 10:30 and edited at 14:10 shows only
  14:10. Several files below were certainly started earlier than the time shown.
- **This repository lives in a OneDrive folder.** Cloud sync can rewrite mtimes. Git commit times
  are the only timestamps here not subject to that.
- **Commit times are when work was committed,** which trails when it was written.
- **Times are local (IST, UTC+05:30).**

**Attribution is missing and must be added by hand.** All three commits are authored
`Shivoham Pandey`, because one machine did the committing. The git record therefore does **not**
reflect who contributed which idea or which module. Every `Contributors:` line below is left blank
deliberately — fill them in together, from memory, while it is still fresh. Do not guess, and do not
let one person fill in the other's entries.

**Contributors to be named:**

- Shivoham Pandey — <pandey.sharma@gmail.com>
- _(co-originator of the concept — name, email, and contribution to be filled in)_
- Dhruv Kumar, Kaanishka S, Harismitha Srinath — HBZ-13 teammates; individual contributions to be
  recorded below where applicable

---

## 26 September 2026 — before the event

Work exists on this date, the evening before the hackathon build window. It is recorded here because
it is what the evidence shows.

**20:17 — 20:18 · Reference material collected**
`AI-Powered Codon Optimization for Vaccines.pdf`, and a second copy.
The project's starting brief was codon optimization, not chromatin. The pivot is visible in the
commit history below.
Contributors: _______________

**20:40 · First prototype**
`legacy/app_v1.py` — the version later retired into `legacy/`.
Contributors: _______________

**21:32 · `chronocell` package created**
`chronocell/__init__.py`.
Contributors: _______________

**21:48 · First tests**
`tests/test_core.py`.
Contributors: _______________

**22:27 · Benchmark harness**
`chronocell/benchmark.py`.
Contributors: _______________

**23:00 · Architecture written**
`ARCHITECTURE.md`. The system design predates the build window.
Contributors: _______________

---

## 27 September 2026 — event day

Event schedule for reference: report 08:15 · build 10:00–13:00 and 14:00–16:30 · judging 17:30–19:00.

### 09:39:21 — Commit `326c3d6` "Initial commit: ChronoCell - AI-Powered Codon Optimization platform"

25 files. **The commit message names codon optimization, but the contents are already the chromatin
system** — `egnn.py`, `genome.py`, `physics.py`, `build_graph.py`, `features.py`, `formats.py`,
`synthetic.py`, `train.py`, `viz.py`, `theme.py`, `benchmark.py`, plus `ARCHITECTURE.md`,
`AUDIT.md`, `README.md`, `app.py`, `pytest.ini`, and the two reference PDFs.

This is the clearest single piece of evidence for when the pivot from codon optimization to chromatin
folding happened: **before this commit**, i.e. on or before the evening of 26 September.

Contributors: _______________
Conception note — whose idea was the pivot, and when was it discussed: _______________

### Morning build block, 10:00–13:00

**10:37 · Reference genome data + synthetic generator** — `chronocell/data/hg38.json`, `synthetic.py`
**10:38 · Feature extraction** — `features.py`
**10:44 · Disease scenarios** — `scenarios.py`
**10:49 · File formats (PDB/wwPDB export path)** — `formats.py`
**10:52 · UI split begins** — `ui/__init__.py`; the monolithic `app.py` starts becoming a package
**10:58 · Graph construction + Colab bridge** — `build_graph.py`, `colab.py`
**10:59 · Colab code packer** — `colab/pack_code.py`
**11:02 · E(3)-equivariant graph network** — `egnn.py`
**11:06 · GPU notebook + v3 tests** — `colab/ChronoCell5D_Colab.ipynb`, `tests/test_v3.py`
**11:08 · EGNN tests** — `tests/test_egnn.py`
**11:17 · Training loop** — `train.py`
**12:26 · Polymer physics analytics** — `physics.py` (R_g, ν, contact-decay exponent γ, packing)
**12:29 · Genome coordinate handling** — `genome.py`
**12:46 · Demo state generator** — `demo_states.py` (synthetic healthy / cancer / senescent)

Contributors, morning block: _______________

### 13:31:02 — Commit `24830e4` "Update codebase with new features and tests"

31 files. Adds the ChronoAgent interpreter (`agent.py`), the biological-state system (`states.py`),
disease scenarios, the Colab GPU pipeline, and the first real UI package — `ui/common.py`,
`ui/four_d.py`, `ui/agent_panel.py`, `ui/states_panel.py` — plus `tests/test_v3.py` and
`tests/test_v31.py` and `APP_GUIDE.md`.

Contributors: _______________

### Afternoon build block, 14:00–16:30

**13:39 · Gene annotation dataset** — `chronocell/data/genes_hg38.json.gz`
**13:43 · TAD / domain detection** — `domains.py`

**13:45 · `chronocell/therapy.py` — the drug-lab model** ⭐
The virtual treatment model: each drug class reduced to *where* it acts (target weights per bead from
the signal track and 3D crowding) and *which way* it pushes (+1 open / −1 compact / 0 re-loop); dose
scales displacement; direction-gated movement toward a healthy baseline; every treated conformation
relaxed back to valid polymer geometry (bond lengths → b₀, excluded volume); restoration scored as
`100 × (1 − RMSD(treated, healthy) / RMSD(untreated, healthy))`.
**This module is the strongest candidate for novel technical subject matter in the project.** Record
its conception carefully.
Contributors: _______________
Conception note — who proposed the two-axis drug reduction, the direction gating, and the restoration
metric, and when: _______________

**13:47 · Session snapshot/restore** — `snapshot.py`
**13:49 · PDF dossier generator** — `pdf_report.py`
**13:51 · Biological states reworked** — `states.py`
**13:52 · Shared UI/dataset assembly** — `ui/common.py`
**13:54 · Streamlit config** — `.streamlit/config.toml`
**13:55 · Theme + visualisation** — `theme.py`, `viz.py`
**13:56 · Compare workspace** — `ui/compare.py` (side-by-side states, linked rotation)
**13:58 · Genes workspace** — `ui/genes_view.py`
**13:59 · ChronoAgent + plain-language guide** — `agent.py`, `ui/guide.py`
**14:10 · Drug lab UI, agent panel, 4D workspace** — `ui/drug_lab.py`, `ui/agent_panel.py`,
`ui/four_d.py`. The drug lab pre-computes all 11 doses so the dose slider animates in-browser.
**14:11 · States panel** — `ui/states_panel.py`
**14:12 · Data ingest** — `ingest.py`
**14:14 — 14:16 · Test suites** — `tests/test_v31.py`, `tests/test_v32.py`
**14:17 · Gene lookup** — `genes.py`
**14:22 · Dependencies + sync view** — `requirements.txt`, `ui/sync_view.py`
**14:24 — 14:25 · Documentation** — `JUDGES_GUIDE.md`, `APP_GUIDE.md`, `AUDIT.md`,
`coordinates/README.md`
**14:30 · Application shell finalised** — `app.py`

Contributors, afternoon block: _______________

### 15:00 — 15:02 · Validation against real microscopy ⭐

**15:00** `validation/validate_tracing.py` · **15:01** `validation/results.json` ·
**15:02** `validation/RESULTS.md`

Validated against Bintu et al., *Science* 2018 chromatin tracing data (public:
`github.com/BogdanBintu/ChromatinImaging`) on three datasets × three random splits. Method: split
cells into halves A and B; feed only half A's contact frequencies (<150 nm) to the pipeline; score
against half B's measured median pairwise distances, which the model never saw.

Headline numbers, recorded as measured — including the unfavourable ones:

| Dataset | Model ρ | Baseline ρ | Model ρ, trend removed | Baseline, trend removed | Model / ceiling |
|---|---|---|---|---|---|
| IMR90 chr21:28–30 Mb (4,832 cells) | 0.86 ± 0.01 | 0.93 | **0.38 ± 0.02** | 0.01 | 39 % |
| A549 chr21:28–30 Mb (3,941 cells) | 0.83 ± 0.01 | 0.91 | **0.51 ± 0.03** | 0.00 | 54 % |
| IMR90 chr21:18–20 Mb (1,277 cells) | 0.51 ± 0.03 | 0.96 | **0.09 ± 0.02** | 0.00 | 36 % |

Stated limitations, from `validation/RESULTS.md`: on raw rank agreement the genomic-distance baseline
scores higher than the model; absolute distances are off by roughly 3× (the 72 nm bond length assumed
for 30 kb beads is too small for this data); and the input was imaging-derived contact frequencies
rather than sequencing Hi-C, so a direct Hi-C → imaging test remains to be done.

Contributors: _______________

**15:14 · Master README** — `README.md`

### 15:19:14 — Commit `af4d159` "v3.2: Master README, new modules (domains, genes, therapy, drug lab, PDF reports), validation suite, judges guide"

39 files — the state presented to the judges.

Contributors: _______________

---

## State at end of day

**Modules:** 20 in `chronocell/`, 10 in `ui/`.
**Tests:** 73 `def test_` functions across 5 files — `test_core.py` (16), `test_egnn.py` (5),
`test_v3.py` (14), `test_v31.py` (20), `test_v32.py` (18). The README quotes 83; the difference is
most likely parametrised cases, and the exact figure should be confirmed with
`python -m pytest --collect-only -q` before it is quoted anywhere.

**Six workspaces:** 3D structure · 4D dynamics · Compare · Drug lab · Genes · Guide.

**Public disclosure:** demonstrated to the HackBlitz judging panel on 27 September 2026, and pushed
to `github.com/WolframNU13/Team-NU-13`. If that repository is public, 27 September 2026 is the
worldwide publication date for everything in commit `af4d159`.

---

## Open items

- [ ] **Fill in every `Contributors:` line above, jointly.** This is the only part of the record that
      cannot be reconstructed from evidence later, and it is the part that matters most.
- [ ] **Add the co-originator as a git author.** Have them make commits under their own name and
      email from here on, so the record shows joint work going forward.
- [ ] Add `LICENSE`, `AUTHORS`, `CITATION.cff` — none currently exist, so the code is
      "all rights reserved" by default and authorship is recorded nowhere in the repo.
- [ ] Confirm whether `github.com/WolframNU13/Team-NU-13` is public or private.
- [ ] Check the HackBlitz registration terms for any IP assignment or licence clause.
- [ ] Confirm the real test count before quoting 83 anywhere.
- [ ] Direct Hi-C → imaging validation on IMR90 (named as the next step in `validation/RESULTS.md`).

---

## Maintaining this log

From here on, write entries **as work happens**, dated, with the contributor named. A
contemporaneous log is worth far more than a reconstructed one. Append to the bottom, never rewrite
history above, and commit each entry so its date is independently witnessed by git.

---

## 28 September 2026 — post-event cleanup (live log)

**Entries from here on are written as the work happens**, unlike the reconstructed sections above.
Times are the machine clock (IST, UTC+05:30); git commit times witness them independently.

Made by, for every entry in this section unless stated otherwise: **Claude (AI coding assistant,
Claude Code, model Claude Opus 5.5)**. Requested and approved by: **Shivoham Pandey**.

**15:02:34 · Commit `52c9174` on `main`: this log added to git.** It had existed only as an untracked
file. Committed at the requester's instruction before any cleanup began. Not pushed.

**15:02:35 · Branch `chore/hackathon-cleanup` created from `52c9174`.** All cleanup below happens on
this branch; `main` is untouched.

**15:05 · Repository location recorded.** The requester reports that the GitHub repository was
transferred and renamed from `github.com/WolframNU13/Team-NU-13` to
`github.com/Sh1voham/ChronoCell-5D`. It could not be independently verified from this machine: the
GitHub login here (WolframNU13) no longer has access, and the new repository is not public. The local
`origin` remote still points at the old name and was not changed.

**15:07 · Baseline before any change: 101 tests passed.**

**15:38:04 · Removed `legacy/app_v1.py`.** The original prototype app, imported by nothing. It
remains in git history: `git show 52c9174:legacy/app_v1.py`.

**15:38:04 · Stopped tracking `colab/chronocell_code.zip`.** A generated file (made by
`python colab/pack_code.py`) that `.gitignore` already listed. The local copy is kept.

**15:38:04 · `JUDGES_GUIDE.md` moved to `docs/OVERVIEW.md`,** then edited at **15:38:49**:
- the judge-specific wording and the 3-minute demo script were removed;
- "pitch" and "questions judges ask" were retitled "In 30 seconds" and "Frequently asked questions";
- the outdated accuracy paragraph ("imaging data weren't available") was replaced with the real
  microscopy validation results, including the unfavourable ones.

**15:38:05 — 15:38:17 · Stopped tracking `validation/data/*.csv` (16 MB) and added `validation/data/`
to `.gitignore`.** This is third-party data (Bintu et al., *Science* 2018) whose licence is
unconfirmed. `validation/validate_tracing.py` downloads it on first run. The local copies are kept.

**15:39:31 · Documentation and code references updated:**
- `APP_GUIDE.md` and `AUDIT.md`: references to `JUDGES_GUIDE.md` now point to `docs/OVERVIEW.md`,
  and the `legacy/app_v1.py` row is removed.
- `ARCHITECTURE.md`: folder name `Hack-a-thon/` changed to `ChronoCell-5D/`, and the legacy entries
  now point to git history.
- `chronocell/features.py` docstring: "Teammate 1 loop" changed to "original per-bin loop".
- `validation/validate_tracing.py` docstring corrected: it writes `results.json` only.
- `validation/RESULTS.md`: now notes the on-demand data download.
- `colab/ChronoCell5D_Colab.ipynb`: example Google Drive path changed to `MyDrive/ChronoCell-5D`.

**15:40 — 15:41 · Added `docs/images/fold.png` and `docs/images/compare.png`.** Rendered with
`chronocell.snapshot` from the synthetic reference model and the synthetic demo patients, and
captioned as synthetic in the README.

**15:42:20 · `README.md` rewritten.** The old file was wrapped in `<![CDATA[ ... ]]>`, which made
GitHub show it as one unformatted paragraph. Changes:
- new structure, with the two images and a pipeline diagram;
- a validation summary with its limits;
- test count corrected from 83 to 101;
- two claims the app does not support were removed: an "immunoglobulin cluster" preset and
  "persistence length";
- clone URL changed to `Sh1voham/ChronoCell-5D`;
- the "Team NU-13" section was removed at the requester's instruction.

**Licence: deliberately not added.** The old README claimed "MIT", but no `LICENSE` file exists, and
the open items above say to check the HackBlitz registration terms for an IP-assignment or licence
clause first. The claim was replaced with "no licence has been chosen yet". Choosing one remains an
open decision for the authors.

**Kept, by decision:**
- the two reference PDFs (at the requester's instruction);
- `UPDATES.md`;
- `.claude/launch.json` (local app-preview settings used by the assistant);
- the historical "Teammate" references in `AUDIT.md`, which describe the original specification.

**15:43 — 15:47 · Verification after the changes.** Every Python file compiles; the Colab notebook
and `validation/results.json` parse; **101 tests passed**, the same as the baseline. This includes
end-to-end runs through all six app pages.

**15:48 · Committed on `chore/hackathon-cleanup`** (the commit that adds this entry). Not pushed:
this machine's GitHub login has no access to the renamed repository.

Contributors to the ideas behind these changes: _______________ (the requester directed the cleanup;
fill in if others were involved).

## 28 September 2026 — reconstruction accuracy work (live log)

Made by Claude (Claude Code, Claude Opus 5.5); requested by Shivoham Pandey. Branch
`feat/microscopy-accuracy-v3.3`, created from `chore/hackathon-cleanup` at `124ed41` (which contains
`main` plus the cleanup commit).

**Goal set by the requester:** raise the accuracy measured against real microscopy (currently 36–54 %
of the experiment's own reproducibility, about 43 % on average) towards 80–90 %.

**Rules fixed before any tuning, so the final number is credible:**
- **Practice data (approved by the requester):** tuning uses only Bintu et al. datasets that were never part of the reported
  validation: K562 chr21:28–30 Mb, HCT116 chr21:28–30 Mb (untreated and 6 h auxin), and
  HCT116 chr21:34–37 Mb (untreated). The IMR90 cell-cycle set is excluded because it shares cell
  line and region with a test set.
- **Test data:** the three reported datasets (IMR90 chr21:28–30 Mb, IMR90 chr21:18–20 Mb,
  A549 chr21:28–30 Mb) stay untouched until the final run.
- **Headline number:** Σ model / Σ ceiling of the trend-removed Spearman ρ over the three test
  datasets, averaged over 3 random splits. This is the ceiling-weighted average of the per-dataset
  ratios, so the noisy 18–20 Mb set (ceiling 0.25) cannot dominate it. Per-dataset ratios, raw ρ,
  the genomic-distance baseline and Lin's CCC are reported alongside, favourable or not.
- **What the model sees:** only half A's contact frequencies (< 150 nm), as before.

**18:29 · Branch created.**

**18:35 · Branch renamed** from `feat/reconstruction-accuracy` to `feat/microscopy-accuracy-v3.3`, at the
requester's suggestion. It stays based on the cleanup commit so the restructured README is kept.

**Direction from the requester (18:35):**
- extend v3.2; keep all 101 tests passing;
- report two clearly separated scores:
  - *Contact map fit*: agreement with the input contact data;
  - *Independent microscopy accuracy*: against unseen imaging data;
- Phase 1: a microscopy validation pipeline and an ensemble Langevin population solver;
- Phase 2: physics fixes (ICE, positive fitted exponent, bending stiffness, nuclear confinement);
- Phase 3: UI additions (distance probe, slicing plane, timing table, loss charts);
- Phase 4: REST API and a JSON log without compliance claims.

**18:29 — 18:33 · Step 1 diagnostics, practice data only** (split 0; trend-removed Spearman ρ as a
% of each dataset's ceiling):

| | K562 28–30 Mb | HCT116 28–30 Mb |
|---|---|---|
| Ceiling (half A vs half B) | 0.983 | 0.930 |
| Contact frequencies alone, no 3D (f^-1/3) | 88 % | 82 % |
| Best single 3D structure, from half A's *true* medians | 63 % | 56 % |
| Current pipeline: shortest-path MDS start only | 74 % | 75 % |
| Current pipeline: after gradient fit | 37 % | 53 % |
| Current pipeline: full (with EGNN) | 37 % | 55 % |

Reading:
- the input carries enough information for more than 80 %;
- a single 3D structure is capped at roughly 60–75 %;
- the gradient stage of the current pipeline loses accuracy relative to its own starting point.

This motivates a population (ensemble) model.

**18:36 — 18:48 · First maximum-entropy Langevin ensemble (scratch prototype, practice data, split 0).**

How it works:
- 100 replica chains run overdamped Langevin dynamics as Gaussian chains, with soft excluded volume.
- Each pair of beads has a potential −ε_ij·φ(d), with φ a smooth step at the 150 nm contact radius.
- ε_ij is updated until the ensemble's contact frequency matches the input.
- The prediction is the ensemble's median distance.

This is a maximum-entropy inversion (the idea behind Zhang & Wolynes 2015), with lengths in units of
the contact radius. Results, untuned:

| | HCT116 28–30 Mb | K562 28–30 Mb |
|---|---|---|
| Trend-removed ρ vs half B (% of ceiling) | 0.694 (74.7 %) | 0.797 (81.1 %) |
| Raw ρ | 0.916 | 0.963 |
| Lin's CCC (nm) | 0.81 | 0.77 |
| Model / real scale | 1.01 | 1.16 |
| Contact-map fit ρ (ensemble vs input frequencies) | 0.794 | 0.754 |

Speed: about 5 ms per Langevin step on this CPU (float32, 100 replicas × 65 beads), about 35 s per
dataset.

**18:53 — 18:57 · Langevin prototype dropped.**
- On HCT116 + auxin it reached only 19 % of the ceiling, although contact frequencies alone give 82 %.
- Cause: the simulated time (about 60 units) is shorter than the chain's relaxation time (about 140),
  so long-range contacts never equilibrated. The model's contact frequencies ran about 50 % above the
  input.
- Overall on practice data: 63 %.

**18:57 — 19:10 · New population model: `chronocell/ensemble.py`.**
- **Model:** a maximum-entropy Gaussian polymer ensemble. This approach is prior art, not a novel
  method: HIPPS/DIMES, Shi & Thirumalai, PRX 2019 and Nat Commun 2023.
  - Each contact frequency is inverted to a pair spread through the Maxwell distribution.
  - A valid covariance is fitted by weighted least squares, weighting each pair by its binomial
    reliability.
  - Statistics are exact.
- **Trajectories:** 100 Langevin trajectories are sampled exactly (Ornstein–Uhlenbeck, mode by mode)
  from the fitted spring network.
- **Settings:** frozen after tuning on practice data only. Tested: 1,500 vs 4,000 iterations,
  binomial vs uniform weights, splits 0 and 1. The spread across settings was under 1 point.
- **Tests:** `tests/test_v33.py`, 7 tests on planted populations: ideal chain, partial loop, missing
  pairs, exact Langevin sampling, input checks. All pass.
- **`validation/validate_tracing.py`:**
  - now scores v3.2, the v3.3 ensemble (exact and 100 trajectories) and a no-3D reference;
  - reports contact-map fit and microscopy accuracy separately;
  - adds a `--practice` option.

**Practice-set results** (`validation/results_practice.json`; 4 datasets × 3 splits; trend-removed
Spearman ρ as a % of the ceiling, Σ model / Σ ceiling):

| | K562 | HCT116 | HCT116 + auxin | HCT116 34–37 Mb | **Overall** |
|---|---|---|---|---|---|
| v3.2 single structure | 31–39 % | 52–58 % | 47–50 % | 51–58 % | **48.3 %** |
| v3.3 ensemble (exact) | 92–93 % | 91–92 % | 87–90 % | 95–96 % | **92.0 %** |
| v3.3, 100 trajectories | 91–92 % | 92 % | 88–89 % | 95 % | **91.9 %** |
| No 3D (direct inversion) | | | | | **87.5 %** |

- Contact-map fit of the ensemble: 0.97–0.99.
- Raw ρ of the ensemble (0.92–0.98) now beats the genomic-distance baseline (0.77–0.93) on every
  practice set. v3.2 did not.

**Status:** the three held-out TEST datasets have NOT been run yet. The 80–90 % question is not
answered until they are. Next: run `python validation/validate_tracing.py` once, record the result
here, and update `validation/RESULTS.md` whatever it shows. The full 101-test suite was not rerun for
this commit; only new files and the validation script changed.

**19:53 — 19:55 · Held-out TEST run (once, settings frozen beforehand).** `python validation/validate_tracing.py`,
3 datasets × 3 splits → `validation/results.json`. Trend-removed Spearman ρ as a % of the ceiling:

| Test dataset | v3.2 single structure | **v3.3 ensemble** | 100 trajectories | No 3D (direct inversion) |
|---|---|---|---|---|
| IMR90 chr21:28–30 Mb | 39.1 % | **88.2 %** | 88.4 % | 81.4 % |
| A549 chr21:28–30 Mb | 54.2 % | **91.2 %** | 91.5 % | 87.9 % |
| IMR90 chr21:18–20 Mb (ceiling 0.25) | 35.2 % | **54.4 %** (±11 points across splits) | 53.2 % | 42.3 % |
| **Overall, Σ model / Σ ceiling (pre-registered)** | **45.2 %** | **85.6 %** | **85.7 %** | 79.7 % |

**Also recorded, favourable or not:**
- The unweighted mean of the per-dataset ratios is 77.9 %.
- Raw ρ of the v3.3 ensemble is 0.976 / 0.952 / 0.868. The genomic-distance baseline scores
  0.930 / 0.915 / **0.962**, so the model loses on raw ρ on the weak 18–20 Mb region.
- Lin's CCC (nm) is 0.97 / 0.93 / 0.46.
- Contact-map fit is 0.99 / 0.99 / 0.97.
- The ensemble's cell-to-cell spread (CV 0.42, fixed by the Gaussian model) is below the measured
  0.50–0.58.

**Conclusion:** the requester's target (overall 80–90 %) is met on held-out data under the rule fixed
before tuning: **85.6 %**. The weak-structure region remains far below it.

**19:58 · Documentation updated to match:**
- `validation/RESULTS.md` rewritten: protocol, the two scores, full tables and an honest reading.
- The accuracy sections of `README.md` and `docs/OVERVIEW.md` updated. Both say the population model
  runs on windows and is **not yet built into the app's pages**.
- Test count updated from 101 to 108.

**Full test suite: 108 passed** (the 101 existing tests plus 7 new).

**Found in passing, not a fix:** the requested "positive, data-fitted γ" is already how the code works
(`physics.contact_decay` / `domains.decay_exponent` fit P(s) ~ s^-γ from data, with γ > 0). The
negative-γ error was only in the pasted spec.

### Phase 2 — physics refinements (20:07)

Made by Claude (Claude Code, Claude Opus 5.5); requested by Shivoham Pandey.

- **ICE balancing, new `chronocell/normalize.py`.** ICE = Iterative Correction and Eigenvector
  decomposition (Imakaev et al., Nat Methods 2012).
  - Works on sparse contact lists, with cooler-style defaults: ignore 2 diagonals, MAD filter, and
    `min_nnz` re-filtering until stable.
  - The first version did not converge on sparse maps: row-sum CV was still 5 % after 500
    iterations. Two causes:
    - kept bins whose partners had been masked;
    - a stopping rule stricter than cooler's.
  - After fixing both, the synthetic chr22 map converges in 1,470 iterations (about 1 s). On a
    planted-bias test the bias is recovered with r = 0.999.
  - `balanced_contacts()` returns balanced (ci, cj, cm). `python -m chronocell.build_graph --balance`
    uses it.
- **Bending stiffness and nuclear-envelope confinement** added to `egnn.FitConfig` (`lambda_bend`,
  `bend_cos0`, `lambda_confine`, `confine_radius_nm`), with NumPy twins `physics.loss_bend` and
  `physics.loss_confinement`.
  - Measured on the practice sets: no accuracy change beyond ±1–2 points. Confinement is inactive at
    realistic radii.
  - **Both are off by default.** v3.2 behaviour is unchanged.
  - They are not added to the Gaussian population model, because the fitted couplings absorb any
    quadratic prior.
- **Positive, data-fitted γ:** already the case (see the Phase 1 note). No change needed.
- **Wording:** ChronoAgent's "Biophysical diagnosis" heading renamed to "Biophysical assessment"
  (agent, PDF output, test, APP_GUIDE, OVERVIEW), to avoid a medical-diagnosis reading. The exports
  already say "Research use only — not a clinical diagnostic" and contain no compliance claims.
- **`validation/TUNING.md` (new):** every practice-set experiment with its numbers, including the
  abandoned ones.
- **Tests: 112 passing** (4 new: ICE ×2, bend/confinement ×2). README count updated.

### Phase 3 — user interface (20:39)

Made by Claude (Claude Code, Claude Opus 5.5); requested by Shivoham Pandey.

- **Population model in the app** (3D structure → 03 Model & convergence → *Build population model*):
  - windows up to 400 beads; built from the window's sequencing counts (`ensemble.fit_from_counts`);
  - the adjacent-bead contact probability is an explicit user-visible assumption, and lengths are
    anchored to b₀;
  - *Population model* joins the structure switch; the view shows the representative member;
  - a note explains that members have no excluded volume;
  - multi-model PDB export of the 100 structures.
- **Speed-up:** the ensemble loss now uses the Gram matrix (A Aᵀ) instead of a (pairs × N) difference
  tensor. 300 beads went from more than 10 min to 22 s.
  - The held-out validation was re-run: the exact ensemble is unchanged (85.6 %).
  - The 100-trajectory sample moved from 85.7 % to 85.6 % (A549 from 91.5 % to 91.3 %) through
    sampling round-off. `validation/results.json` and `RESULTS.md` updated; the earlier numbers are
    left as logged above.
- **Two separate accuracy scores**, new `chronocell/accuracy.py`:
  - Contact-map fit, live on the window.
  - Microscopy accuracy, the method benchmark read from `validation/results.json` and labelled "not
    measured on this window".
  - Shown in 03, the PDF dossier (new section) and the JSON report (`accuracy` block).
- **Distance probe** (*Measure*): two beads by number. It shows the distance in the displayed
  structure, the separation along the DNA and the loci. With the population model it also gives the
  population median, the middle 50 % of cells and the contact probability. The pair is marked in 3D.
  - **Limitation, stated honestly:** beads are picked by number, not by clicking, because the app's
    3D chart widget does not return click events.
- **Slicing plane** (*Display*): x, y or z normal at a % of the fold's extent. Everything beyond it
  is hidden (tube faces, beads, line, context), and a translucent sheet marks the plane.
- **Execution telemetry table:** every reconstruction this session, per stage, with wall-clock
  time, ms per bead, device, final loss and contact-map fit. It uses real measurements.
- **Loss-convergence charts:** the existing v3.2 loss chart, plus a new chart for the population
  fit.
- **Docs:** APP_GUIDE §10 (real-microscopy accuracy) and new §20 (how to use v3.3); a Guide-page
  section "How accurate is it? Two scores, never mixed"; README and OVERVIEW updated. Version
  string changed to 3.3.
- **Verified in the running app** (in-app browser, chr22 30.5–32.0 Mb, 150 beads):
  - the population model builds;
  - contact-map fit 0.867; microscopy benchmark 85.6 % with its scope note;
  - the probe gives 196 nm in the displayed structure and a population median of 327 nm (IQR
    235–428 nm);
  - the slicing plane cuts the tube;
  - telemetry rows appear.
- **Tests: 114 passing** (2 new: end-to-end app test of the population model, scores, probe,
  slicing and dossier; and the PDF two-score section).

### Phase 4 — REST API and run log (21:41)

Made by Claude (Claude Code, Claude Opus 5.5); requested by Shivoham Pandey.

- **New `chronocell/api.py`:** handlers for `/api/v1/reconstruct` (population ≤ 400 beads or
  single ≤ 2,000), `/api/v1/metrics` and `/api/v1/benchmark`.
  - They are plain dict-in / dict-out functions; `create_app()` wraps them in FastAPI, and
    `python -m chronocell.api` serves them.
  - Responses keep the two accuracy scores separate.
  - Bad input, including model-level errors and non-object bodies, is rejected as a request error
    (HTTP 400), not a server error.
- **FastAPI is not installed and was NOT installed** (an install needs the requester's go-ahead).
  The wrapper's test is skipped until `pip install fastapi uvicorn`; the handlers are fully tested.
- **Run log:** JSON lines in `.chronocell_cache/api_run_log.jsonl` (git-ignored). Each line holds:
  - run id, UTC time, endpoint and software version;
  - parameters, as sizes only, never raw data;
  - the input's SHA-256;
  - run time, status and a result summary.
  It is documented as a reproducibility record, with no compliance claim.
- **Docs:** README (REST API section), APP_GUIDE, and an optional line in `requirements.txt`.
- **Tests: 117 passing, 1 skipped** (the FastAPI wrapper).

**21:41 · Pushed** branch `feat/microscopy-accuracy-v3.3` to github.com/Sh1voham/ChronoCell-5D
(private) at the requester's instruction. `origin` was updated from the old WolframNU13/Team-NU-13 URL.
No licence was added; `main` was not pushed.

---

## 4 October 2026 — v4 (reconstructed from git history)

**Reconstructed on 5 October 2026 from the commit messages of this branch, not written as the work
happened.** Commit times are the witnesses (IST). Every commit on `feat/v4-evidence` is authored
`DHRUV` (the machine's git identity), which does **not** record who contributed which idea; the
`Contributors:` lines are left blank for that reason. The measured results are in
`validation/RESULTS.md`; they are not repeated here.

**11:13 · `1e20719` Baseline:** ChronoCell-5D v3.3 as downloaded (118 tests passing).
**11:18 – 12:54 · Core of v4:**
- `e4feebc` scalable population model with exact per-pair uncertainty (Pillars 1, 2);
- `37cdab2` built-in `.hic` reader (local or remote by HTTP range) and the shared validation protocol;
- `35c1d8d` perturbation engine, SV files, chromosome names; cohesin parameters frozen (Pillar 4);
- `1217079` Gate 1 settings frozen before the Su chr21 test run;
- `141ea88` Self-Math PDB State Evaluator (Pillar 7);
- `5fa9472`, `ced6aaa` population model beyond 400 beads in the app; exact interval in the probe;
- `66023af` genome assemblies as configuration (hg38 + mm39) (Pillar 6);
- `f6516fe` reproducibility record, JSON + PDF (Pillar 8);
- `20493a0` held-out results of Gate 1 and Gate 4 (cohesin); Gate 2 recalibration frozen.

**16:12 – 16:32 · Structural variants and environment:** `50a1e27` SV test pre-registered; `38bfbbc`
Gate 4b measured, not validated; `7338fda`, `89fdaa2` variant files and the 04 Variant impact panel;
`1224f11` pinned versions, Dockerfile, CI.

**17:10 – 18:15 · Prediction without contact data:** `233fdaa` Gate 5 pre-registered; `780aba7`
predictor frozen on practice data; `9fef9e3`, `84effe1`, `c183b9b` tuning record and protocol;
`ebfe542` recalibrated interval in the probe; `9ba5061` RESULTS.md v4; `03f20c0` Gate 5 passes its
pre-registered rule, modestly; `e153602` the prediction in the app.

Contributors: _______________

---

## 4 – 5 October 2026 — remaining gates and features (live log)

Made by Claude (Claude Code, model Claude Opus 5.5), at the request of this machine's user (git
identity `DHRUV`). Every commit below is on `feat/v4-evidence`, local only: nothing was pushed and no
pull request was opened, as requested. Times are the machine clock (IST).

**Rules the work followed** (set by the requester): every number in the documents comes from a result
file through `validation/report.py`; pass rules are committed before the test that uses them; tuning
uses practice data only; failures are reported next to successes; nothing is pushed.

**Pre-registrations, each committed before its test:**
- `a372531` Gate 2b, intervals with Hi-C input (rule + practice-fitted Hi-C recalibration);
- `1af2b96` Gate 5m, the predictor on mouse (not run, see below);
- `b625717` Gate 2c, a per-pair reliability score (after a practice comparison of three candidates).

**Results:**
- `bea3fde` **Gate 2b fails**: with Hi-C input the stated intervals stay far below nominal even after
  a Hi-C recalibration fitted on practice data. The probe no longer applies the imaging recalibration
  to sequencing input; it shows the measured shortfall.
- `667d290` **Gate 2c fails**: the chosen per-pair score ranks error weakly with imaging-derived input and
  not at all with Hi-C input, below the pre-registered bar; no reliability score is shown in the app.
- `3297c52` **Gate 3 measured**: with imaging-derived input the population models are best or tied best
  on 27 of 28 test units; with Hi-C input on 18 of 26 (the no-3D inversion or PASTIS is ahead on seven
  genome-scale chromosomes). The full run was stopped by the machine's 2-hour job limit, so the plan
  was run in four parts and merged (`--merge`); PASTIS PM2 on the two 651-locus sets did not finish in
  an hour and is recorded as not finished.
- `21b7078` **Cost**: the scale benchmark is complete (v4 fits the whole synthetic chr22 in under nine
  minutes on this CPU). The per-chromosome runtime is **partial**: the 2-hour job limit stopped it after
  13 chromosomes; three were measured while the test suite ran and are kept but discarded from the
  table; 35 remain. `python validation/chromosome_runtime.py --resume` finishes it.

**Features** (`7263d77`):
- CTCF peaks by cell type: "ENCODE, by cell type" fetches the GRCh38 IDR peaks of IMR-90, A549, K562
  or HCT116 (listed in `chronocell/data/validation_sources.json`) once into the cache, checks the MD5
  the portal publishes, records the accession in the prediction's inputs and shows the citation.
  Upload and paste remain.
- Predicted and contact-built maps side by side, with their agreement; never blended.
- `python -m chronocell.predict`: the predicted map (`.npy`) plus a JSON record (inputs with SHA-256,
  model, validation reference, "predicted, not measured"); hg38 only.
- The cohesin-control caveat (Gate 5) in the predictor's banner and in the Guide, read from the result
  file.
- Downloads: one routine for the UCSC sequence and ENCODE peaks; resume after interruption or short
  reads; no file without a published MD5. Tests use a local HTTP server.

**Bugs found and fixed:**
- Building a second population model for a window that already had one (a prediction after a
  contact model, or the reverse) crashed: Streamlit forbids setting a widget's value after it is
  drawn. The choice is now applied before the widget on the next run.
- The variant-impact panel could take a *predicted* population as its base without saying so; it now
  uses populations built from contacts only. The PDB evaluator labels predicted populations.
- Download progress lagged the bytes received (it read the file size before the buffer was flushed).

**Checked by hand, not committed:** the real UCSC download of hg38 chr21 through the app's code
(12,709,705 bytes; MD5 184df2bd9b812b6e6b6da16c6021369e, equal to UCSC's), and the real ENCODE IMR-90
CTCF peaks (595,548 bytes; MD5 equal to the portal's; 374 peaks on chr21). The command line on those
real inputs reproduces the Gate 5 validation map for the IMR-90 test loci exactly (now a test that runs
when the validation data are present).

**Checked in the running app** (a separate instance on port 8502, driven in Chrome; reference chr22,
window 20.0–21.5 Mb, 150 beads): a population model built from contacts; Input "Sequence + CTCF" →
"ENCODE, by cell type" → IMR90 fetched the real ENCODE file (724 peaks on chr22); the chr22 sequence
downloaded from UCSC with a progress bar (12,255,678 bytes; MD5 equal to UCSC's md5sum.txt); the
prediction was built and labelled "predicted from sequence + CTCF (no contact data)"; the banner and
the cohesin caveat read Gate 5 from its result file; the side-by-side view showed both maps (rank
agreement 0.519, 0.132 beyond the trend, size ratio 2.23); the probe showed "coverage not tested" for
the predicted population and, after rebuilding from contacts (no crash), the Gate 2b shortfall for Hi-C
input; the Guide read its numbers from the result files. Seeing the side-by-side view on the synthetic
reference led to one change: it now says when the contacts are the synthetic reference, so the
agreement means nothing about real folding.

**Not done, with the reason:**
- **A second structural-variant test**: none could be pre-registered and run here (RESULTS.md,
  Gate 4b): the cleanest candidate (Firre deletion, GEO GSE98632) has raw reads only.
- **Mouse (mm39) prediction**: the pre-registered test needs a 4DN Data Portal access key (downloads
  answer HTTP 403 without one). Mouse stays without prediction.
- **PASTIS PM2 on the two 651-locus sets** (Gate 3): did not finish in 60 minutes; recorded as not
  finished.
- **35 of 45 per-chromosome runtimes**: stopped by the job limit (above); resumable.

**Housekeeping:** `truststore` listed as optional; `validation_sources.json` completed (PASTIS, UCSC
sequences, genome bundles, published checksums); README, ARCHITECTURE (new section 11), Guide updated;
references to a `validation/reproduce` script and a `paper/` folder, which never existed in this
repository, removed.

**Noticed, not touched:** the tracked file `AI-Powered Codon Optimization for Vaccines(1).pdf` was
deleted from the working tree during this session by something other than these commands. It is
left deleted-but-uncommitted for the owner to decide.

**Tests:** 226 passed, 1 skipped (the FastAPI wrapper), against 201 passed, 1 skipped before this work.

Contributors to the ideas behind these changes: _______________

## 5 October 2026 — Phase A (accuracy core) and Phase B (product features) (live log)

Made by Claude (Claude Code, model Claude Opus 5.5), at the request of this machine's user (git identity
`DHRUV`), who asked for the work to continue unattended. Every commit is on `feat/v4-evidence`, local only:
nothing was pushed and no pull request was opened. Times are the machine clock (IST). The measured results are
in `validation/RESULTS.md` (generated from the result files); this entry only lists what was done.

**Rules followed:** everything additive (existing pages, defaults, models, commands, endpoints and result files
unchanged); pass rules committed in `validation/frozen.py` before each test; tuning on practice data only; each
test run once; every number in the documents from a result file through `validation/report.py`; failures
reported next to successes.

**Phase A (accuracy core).**
- A0 `eb8887d`: PyTorch 2.11.0 with CUDA 12.8 on the RTX 5050 Laptop GPU; CPU fallback on GPU out-of-memory;
  a GPU / CPU agreement test.
- A6 `e07e340`: two untouched genome-scale sets registered with fixed roles; `f5a2295`, `fbc189d`: shared cached
  units (fits + held-out truth), ENCODE Hi-C sources, depth thinning.
- Gate 1c (sizes from Hi-C): pre-registered `c1726a9`, **fail** `41d46c2`.
- Gate 2d (conformal ranges): pre-registered `4dad69c`, **fail** `c3e466b`.
- Gate 2e (per-pair reliability): pre-registered `0a7998d`, **fail** (this entry's documentation commit).
- Gate 3b (learned correction): practice chose no correction, **not run** `8763ab1`.
- Gate 5b (prediction with cohesin peaks): pre-registered `cbebab1`, **fail** `4415433`.
- Gate 5m (the predictor on mouse): the 4DN files turned out to be public on 4DN Open Data, so the test that was
  pre-registered in `1af2b96` ran unchanged: runner `2f205fc`, **pass, modestly** `54075e2`; the predictor is now
  offered for mouse, labelled.

**Phase B (product features).**
- B8 `17c1578`: Research mode switch (default on), "SYNTHETIC ·" labels, warnings on synthetic results.
- B1–B7 library `817428f`: variant engine v2 (joins, two chromosomes, copy number, ranking), analysis suite
  (HiCCUPS-like loops, TopDom-like / Arrowhead-like domains), replicate-aware differential analysis, region-wise
  .hic / .cool / .mcool / .pairs reading with assembly checks, KR balancing, liftover, exports for IGV / Juicebox /
  HiGlass, offline HTML / PDF reports, the `chronocell` command, a local job queue, saved projects, three REST
  endpoints, `pyproject.toml`.
- App `80fb251`: new panels (01 → 06 Analysis suite, 03 → Differential analysis, 02 → 05 Variant engine v2,
  05 Genes → reference annotations), sidebar Projects and Jobs, a built-in bigWig reader, mouse prediction.
- Gate 4c (cohesin loss vs RAD21-degron Hi-C): pre-registered `a1d20fb`, **pass** `88e7436`.
- Gates 6 (loops) and 7 (differential FDR), and Gate 4d recorded as blocked: pre-registered `b22fe02`; Gate 6
  **fail**, Gate 7 **pass** (this entry's documentation commit).
- B9: chromosight and Mustache run from an isolated environment (`.chronocell_cache/tools-venv`); every other
  tool is recorded with the reason it could not run here (`validation/results_tools_b9.json`).

**Checked in the running app** (a separate instance on port 8502, driven in Chrome, real data): ENCODE HCT116
RAD21-mAC in situ Hi-C for chr21:28–30 Mb loaded as the contact map; the analysis suite with KR balancing gave the
same calls as the command line on the same map (1 loop, 2 · 6 boundaries, 14 domains); the differential tab read
four ENCODE maps by URL (two untreated, two 6 h auxin experiments) and tested 2,570 pixels (4 significant); a
population model was built on the GPU (17 s) and the engine applied a 300 kb deletion typed in genomic coordinates;
Research mode off hid the Drug lab; a project was saved. Three problems found this way were fixed: the analysis
panel called uploaded contacts "synthetic" when the structure was the reference model; the engine used the
reference model's synthetic H3K27ac as enhancers; saving a session crashed when an object could not be pickled
(after a module reload) and left a half-written folder. GTEx and ClinVar downloads were checked against their
published checksums.

**Interruptions.** This computer was suspended twice (about 07:35–12:00 and 12:45–15:40). Runs that stalled or
were killed (Gate 2d, the A3 practice, Gate 2e, the test suite) were completed by resumable or chunked reruns;
Gate 4c and Gate 5m crashed once (a dropped network read; an 8-column trace file) and were completed after the
fix. In each case the completed run reproduced the numbers the first run had printed.

**Not done, with the reason:** Gate 4d (two usable SV events found, three needed); a FASTQ → contacts pipeline
(the aligners need Linux / WSL, and no WSL distribution is installed); HiCCUPS, Arrowhead (Java), TopDom, dcHiC,
diffHic, multiHiCcompare (R), CHESS (pysam does not build), cooltools, hic_breakfinder, HiNT, Akita and Orca;
bigWig and beddb writing (compiled tools only; bedGraph and .mcool are written instead).

**Tests:** 262 tests in 20 files (226 passed, 1 skipped before Phase A). The last complete run passed every file (15:46–16:05, GPU). The final GPU run after the browser-check fixes passed 18 of 20 files (247 passed, 1 skipped) and was then stopped by Claude Code because the computer was low on memory; the two remaining files (test_v4_trust.py, test_v4_ui.py) passed in the previous complete run. The full CPU-only run (`CUDA_VISIBLE_DEVICES=` set to empty) was not repeated for the same reason; the GPU / CPU agreement test passes.

Contributors to the ideas behind these changes: _______________

## 6 October 2026 — Quantum lab (live log)

Made by Claude (Claude Code, model Claude Opus 5.5) at the request of this machine's user (git identity `DHRUV`),
who asked for a quantum-computing lab for a competition, with a quantum section in every workspace and one separate
tab. Local commits on `feat/v4-evidence` only; nothing pushed. The measured results are in `validation/RESULTS.md`
(Gate Q, generated from the result files); this entry only lists what was done.

**What was added (all additive, Research mode only, nothing runs until a Run button is pressed).**
- `chronocell/quantum/`: a statevector simulator written for ChronoCell (Qiskit's conventions; GPU through PyTorch;
  OpenQASM 2.0 export; QAOA; an approximate hardware-noise model), QUBO solvers (exact enumeration, exact dynamic
  programming for banded problems, simulated annealing, simulated quantum annealing), problem encodings (domain
  walls, variant set, drug combination, gene group, lattice folding), minimal-basis quantum chemistry from scratch
  with VQE (H2, HeH+), quantum-kernel SVM and the swap test, and quantum / classical walks on contact graphs.
- App: a new workspace **07 Quantum lab** and Quantum sections in 01 (domain walls, lattice fold), 02 (quantum walk,
  variant set), 03 (swap-test tab), 04 (drug combination, VQE) and 05 (gene classifier, gene group); a Guide section in
  plain words.
- `validation/quantum_crosscheck.py`: every circuit type checked against Qiskit 2.2.1 in an isolated environment
  (`.chronocell_cache/quantum-venv`; `requirements-quantum.txt`). The app does not need Qiskit.
- Gate Q (`validation/quantum_gateq.py`): practice on GM12878, rules committed in `validation/frozen.py`
  (`QUANTUM_GATEQ`) before the test, test run once on K562 / IMR-90 (Q1, Q2), the chemistry (Q3) and GM12878 → IMR-90
  genes (Q4). Practice choices: `validation/TUNING.md` §19.

**Course of the work.** The first practice grid for the domain QUBO ranked it below the classical insulation caller;
a per-boundary cost and a wider resolution range (practice only, with the classical callers' grids widened alongside)
brought it level on practice. A full search of every window's 2^19 states was too slow for the grid, so an exact
dynamic programme for banded QUBOs was written (equal to enumeration on every check) and used for the grid; that it
solves the problem in linear time is also reported as the honest classical comparison. The quantum-kernel grid was
widened once because its first best setting sat at the grid's edge. A VQE panel that would have computed on page load
was moved behind a button before the first commit.

**Gate Q.** Pre-registered in `eba6efa` (rules, practice results, the lab itself). Test run once (this entry's
commit): **Q1 pass** (QAOA reached the exact optimum of the held-out domain puzzles), **Q2 fail** (the domain QUBO's
calls agreed with ENCODE's Arrowhead calls less well than the classical TopDom-like and insulation callers, on both
cell lines), **Q3 pass** (UCCSD-VQE exact for H2 and HeH+), **Q4 fail** (the quantum-kernel gene classifier trailed
the classical RBF-SVM by slightly more than the margin; logistic regression beat both). Numbers: `validation/RESULTS.md`.
The test exposed one product problem, fixed afterwards and recorded there: the hardware-efficient VQE option let
HeH+ drift to three-electron states (energies below the true one); it now carries the standard electron-number penalty
and the panel shows the electron count.

## 7 October 2026 — Drug lab additions and quantum drug tabs (live log)

Made by Claude (Claude Code, model Claude Opus 5.5) at the request of this machine's user (git identity `DHRUV`), who
asked for more quantum tabs focused on drugs, new drugs with visualisations and information, accuracy tests for
everything, a slightly enhanced UI with the same layout, and for the work to go ahead without questions. Local commits
on `feat/v4-evidence` only; nothing pushed. Numbers: `validation/RESULTS.md` (generated); this entry lists what was done.

**Added (all additive; the Drug lab's default set and every earlier page and result are unchanged).**
- Drug lab: a *Drug set* switch (Core 4, the default, or Extended 12 with eight more chromatin drug classes) and a
  *Drug guide* (class cards with status and safety themes, a landscape chart, where each class acts on the fold, all
  twelve at full dose, pair synergy, PubChem molecule cards with 2D and 3D structures and drug-likeness rules).
- Three quantum drug tabs (in the Drug lab's Quantum section and in 07 Quantum lab): *Drug molecules* (STO-3G
  chemistry written from scratch, active-space UCCSD-VQE up to 12 qubits), *Heart safety* (a quantum-kernel hERG screen
  from SMILES, next to classical models) and *Docking* (QAOA maximum-weight clique on 20 qubits, PoseBusters examples).
- Accuracy tests, each pre-registered in `validation/frozen.py` before its test data were read: Gate 8 (the Drug lab
  against MINA chromatin tracing of IMR-90 cells treated with real drugs, 4DN), Q5 (hERG, TDC Wang -> TDC Karim), Q6
  (molecule energies against OpenFermion's independent data; stretched-molecule VQE), Q7 (re-docking, PoseBusters).
- Independent checks run only in `.chronocell_cache/quantum-venv`: RDKit (descriptors, de-duplication), OpenFermion
  (reference energies), the Basis Set Exchange (STO-3G exponents). New sources listed in
  `chronocell/data/validation_sources.json`.

**Course of the work.** The from-scratch chemistry matched seven textbook Hartree-Fock energies once its solver was given
a damped warm-up (DIIS alone had locked N2 onto an excited solution); VQE was made fast enough for 12 qubits with exact
closed-form rotations and an adjoint gradient. The first docking design (grid hot-spots) docked no practice complex;
hot-spots along protein H-bond vectors and a 20-vertex core subgraph reached 15 % on practice against 3 % for random
search. The hERG reader first split compound names on commas (637 of 655 rows) and was fixed before the settings were
chosen; the kernel grids were widened on practice when a best value sat at an edge. Gate 8's practice showed that the
simulator's agreement with alpha-amanitin came from a generic compaction pattern (shuffled targeting did as well), so its
rule asks that the real targeting beat shuffled targeting.

**Results** (pre-registered in `fa41e78`; each test run once; this entry's commit): **Gate 8 fail** (no drug met the
rule: where predicted and measured changes agreed, shuffled targeting agreed as well), **Q5 pass** (the quantum-kernel
hERG screen slightly ahead of the classical models on about 13,000 unseen compounds), **Q6 fail** (the chemistry matched
OpenFermion's independent LiH reference; VQE missed chemical accuracy on stretched N2 and HCN), **Q7 pass**, modestly
(QAOA found every best clique; the QAOA route docked more ligands than random search, but far fewer than on practice).

## 7 October 2026 — Round 2 of the quantum gates (live log)

Made by Claude (Claude Code, model Claude Opus 5.5) at the request of this machine's user (git identity `DHRUV`), who
asked to improve every quantum test and to try to turn every fail into a pass. Local commits on `feat/v4-evidence`
only; nothing pushed. A failed result is never re-run or overwritten: each failed gate gets a *new* method, developed on
the data its test had used (now practice), a rule committed before new test data are read, and one run. Numbers:
`validation/RESULTS.md` (generated).

- **Q6b (molecules), fail by one case.** ADAPT-VQE added (Drug molecules tab: *Circuit* switch; the fixed UCCSD circuit
  stays the default). On eleven new stretched cases it came within chemical accuracy in ten (stretched N2: 0.02 and
  0.29 mHa, where the fixed circuit is 126 mHa off); the eleventh, LiH at 3x, missed by 2.5 mHa. Found after the test:
  there the active space's lowest state is a triplet, and ADAPT-VQE had reached the lowest singlet exactly (the state a
  closed-shell VQE targets). The verdict stays fail; the app now reports the spin of the exact state and the lowest
  singlet.
- **Q4b (genes), pass, with a plain caveat.** Labels from ENCODE RNA-seq of each cell line, two more Hi-C features,
  training on three cell lines; on 489 genes of HMEC (a cell line never used) the quantum-kernel classifier was level
  with the classical ones, as the rule asks. All three are only modestly better than chance on a new cell line, and
  with a kernel estimated from 1,000 shots per entry (a real device's situation) the quantum classifier falls to chance.
- **Gate 6b (loop calls), pass.** Gate 6 (fail, unchanged) had chosen the loop caller's settings on three GM12878
  windows; re-chosen on all nine windows it had used, the same caller beat chromosight and Mustache on two cell lines
  never used (HMEC, HAP-1; ENCODE HiCCUPS loops as the reference), mostly through far fewer false calls.
- **Q2b (domains), pass.** The same domain QUBO and QAOA, settings re-chosen on all nine windows already seen: on two
  new cell lines (HMEC, HAP-1) the quantum route's domain calls beat both classical callers on HMEC and came within
  the margin on HAP-1. Plainly: at the new settings QAOA landed on the QUBO's exact optimum in only a minority of
  windows (12-19 %, against 98 % in Gate Q1); its best shots scored about as well as the optimum would have.
- **Q6c (molecules, reference corrected to the lowest singlet), fail.** The same ADAPT-VQE on twelve new cases: within
  chemical accuracy in nine, including N2 at 1.3x and 3.0x and LiH where a triplet lies lowest. It missed three bonds
  stretched to 2.2-2.8x (CO 168 mHa, HCN 16.5, CH2O 4.6): there the gradient-chosen circuit stops at a stationary point
  far from the ground state. A real limit of the method; molecules are not re-tested again.
- **App options from the passes** (defaults unchanged): the Analysis page's loop calling can use Gate 6b's setting
  (FDR 0.01), and the Quantum lab's domain-wall panel can load Q2b's settings.
- **Not retested, and why.** *Gate 8* (Drug lab vs drug-treated cells): a practice check of other chromatin marks
  (ENCODE IMR-90 H3K27me3, H3K9me3, H3K4me3, H3K36me3, H3K4me1, H3K9ac) as the drugs' targets on Gate 8's data
  (`results_gate8_marks_practice.json`) found more agreements than chance alone would give, but small ones, in
  directions that differ by drug and would have to be fitted on these very data; no fresh tracing after these drugs
  exists to test such a model (the only drug-treated tracing that matches a Drug lab class is one compound, A-485, a
  p300/CBP inhibitor, at the mouse Sox2 locus, 4DN). Gate 8 stays a fail and the Drug lab stays a mechanism simulator. *Gates 1c, 2d, 2e,
  5b* (Phase A) need new modelling (nanometre scale from Hi-C, per-band intervals, per-pair reliability, cohesin-peak
  prediction) and new held-out imaging data; they are left for later. *Gate 4d* stays blocked (no third event).
- **Q7b (docking), fail, as practice predicted.** With a Vina-like score and pose refinement the docking pipeline got
  much better: on the 85 complexes of the Astex Diverse set (never used) the QAOA route placed 35 % of ligands within
  2 A, against 5.5 % in Q7's test, and QAOA found the best clique in 95 % of graphs. But random search with the same
  score and the same refinement placed 58 %: once the score is reasonable, poses built from clique matches are worse
  starting points than many random placements around the site. Q7's earlier pass came from a crude score that also
  handicapped random search. The Docking tab's *Polish* button shows both side by side.

**Round 2 in one line.** Passed on new data: Q2b (domains), Q4b (genes, weakly), Gate 6b (loops). Failed: Q6b and
Q6c (molecules: 10 of 11 and 9 of 12; ADAPT-VQE stalls on nearly broken bonds), Q7b (docking: the quantum route loses
to random search). Not retested: Gate 8, Gates 1c / 2d / 2e / 5b, Gate 4d (reasons above). Every original verdict
stands.

## 8 October 2026 — New tabs and new methods (live log)

Made by Claude (Claude Code, model Claude Opus 5.5) at the request of this machine's user (git identity `DHRUV`), who
asked for new tabs and new methods to make the labs pass, over about eight hours, without questions. Local commits on
`feat/v4-evidence` only; nothing pushed. Same rules as round 2: every test pre-registered before its data are read,
run once, failures kept.

- **New tabs.** *08 Scoreboard*: every test in one place with charts. *ADMET profile* (07 Quantum lab and Drug lab →
  Quantum): 21 drug properties from a quantum-kernel model next to classical ones. *Noise & mitigation* (07 Quantum
  lab): a molecule's circuit on a simulated noisy chip, repaired by zero-noise extrapolation and symmetry
  verification.
- **New tests pre-registered.** Gate Q8 (ADMET, 21 official TDC scaffold test splits) and Gate Q9 (error mitigation,
  16 new molecules).
- **Gate Q9 (error mitigation), pass.** On 16 new molecules, symmetry verification + zero-noise extrapolation brought
  the energy of a 4-qubit VQE circuit on a simulated chip at today's best noise level within chemical accuracy in
  every case (median 0.41 mHa; about 29 times smaller than unmitigated).
- **Gate Q8 (ADMET), fail by one endpoint.** On the official held-out splits of 21 drug properties the quantum-kernel
  model was within the margin of the classical model on the same inputs for 16 (17 needed), ahead on several (P-gp,
  CYP3A4 substrate, intestinal absorption); it fell clearly behind on the two clearance endpoints and just past the
  margin on three more. The tab stays, with this standing shown.
- **Docking, once more (practice only).** Splitting random search's budget between random placements and the QAOA
  route's clique poses docked more than the clique route alone but still fewer than random search alone, so no new
  docking test was run: with a reasonable score, the quantum step does not help this docking task.
- **New in the app.** The Guide has a "New" section with a button to the Scoreboard; the Scoreboard downloads the whole
  record as a web page (*Download the evidence report*).
- **ADMET, once more (practice only).** A "projected" quantum kernel (each qubit's own state as the features) did worse
  than the plain quantum kernel on 4 of the 6 drug properties tried, so no new ADMET test was run; Gate Q8's fail
  stands.
- **Gate Q6d pre-registered (molecules, once more).** Round 2 had closed the molecule tests after Q6c; this session's
  request for new methods reopened them, with a new method on new cases. The circuit-growing method (ADAPT-VQE) now
  steps off a flat spot instead of stopping there, and starts from four arrangements of the electrons, keeping the
  best run that ends with the right spin. On all 38 molecule cases used so far it was within chemical accuracy every
  time, including Q6c's three misses (CO stretched 2.8x: 0.14 mHa, was 168). Test: fourteen new cases, all must pass.
