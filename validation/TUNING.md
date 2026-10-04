# Tuning log: practice datasets only

Every design choice for v3.3 was made on **practice datasets** that are not part of the reported
validation. The test datasets (IMR90 chr21:28–30 Mb, IMR90 chr21:18–20 Mb, A549 chr21:28–30 Mb) were
run once, after the settings below were frozen. The time-stamped record is in `UPDATES.md`
(28 September 2026).

- **Practice datasets:** Bintu et al. 2018:
  - K562 chr21:28–30 Mb;
  - HCT116 chr21:28–30 Mb, untreated and 6 h auxin (cohesin depleted);
  - HCT116 chr21:34–37 Mb, untreated.
- **Excluded:** the IMR90 cell-cycle set, because it shares cell line and region with a test set.
- **Score:** trend-removed Spearman ρ against half B as a % of the ceiling. Overall = Σ model / Σ
  ceiling.

## 1. Where the v3.2 pipeline loses accuracy (split 0)

| | K562 | HCT116 28–30 |
|---|---|---|
| Contact frequencies alone (f^-1/3, no 3D) | 88 % | 82 % |
| Best single 3D structure from half A's *true* medians | 63 % | 56 % |
| v3.2: shortest-path MDS start only | 74 % | 75 % |
| v3.2: after gradient fit | 37 % | 53 % |
| v3.2: full, with EGNN | 37 % | 55 % |

**Conclusion.** The input carries more than 80 % of the reproducible structure. One 3D structure
cannot hold it (cap about 60–75 %), and v3.2's gradient stage loses accuracy relative to its own MDS
start. This motivates a population model.

## 2. Langevin contact-potential ensemble (abandoned)

- **Setup:** 100 replica chains, overdamped Langevin, pair potentials −ε·φ(d) learned so the
  ensemble's contact frequencies match the input.
- **Result:** overall 63.2 % (K562 81 %, HCT116 75 %, HCT116 + auxin **19 %**, HCT116 34–37 Mb 78 %).
- **Cause:** about 60 time units were simulated against a chain relaxation time of about 140, so the
  ensemble never equilibrated. The model's long-range contact frequencies ran about 50 % above the
  input, and the contact-map fit was only 0.77.

## 3. Maximum-entropy Gaussian ensemble (adopted: `chronocell/ensemble.py`)

This is the HIPPS/DIMES approach (Shi & Thirumalai, PRX 2019; Nat Commun 2023).

| Variant (split 0) | Overall |
|---|---|
| Direct Maxwell inversion (no 3D) | 86.7 % |
| Gaussian ensemble, PSD projection only | 73.1 % |
| Gaussian ensemble, refined, binomial weights | **91.3 %** |
| Refined, 100 sampled structures (independent draws) | 78.6 % |

| Setting (splits 0 / 1) | Overall |
|---|---|
| Binomial weights, 1,500 iterations | 91.3 % / 92.2 % |
| Binomial weights, 4,000 iterations | 91.3 % / 92.1 % |
| Uniform weights, 1,500 iterations | 91.1 % / 91.9 % |

**Frozen:** binomial weights, 1,500 iterations, learning rate 0.01, 100 trajectories × 50 frames.

The final practice run (`python validation/validate_tracing.py --practice`, 3 splits) gave:

| Model | Overall |
|---|---|
| v3.2 | 48.3 % |
| v3.3 ensemble | **92.0 %** |
| v3.3, 100 trajectories | 91.9 % |
| No 3D | 87.5 % |

## 4. Phase 2 physics terms on the v3.2 single structure (split 0)

| `egnn.FitConfig` | Overall | K562 / HCT116 / auxin / 34–37 |
|---|---|---|
| Default | 47.7 % | 37 / 55 / 47 / 52 |
| Bending stiffness λ = 1, cos θ₀ = 0 | 46.4 % | 37 / 53 / 47 / 49 |
| Bending stiffness λ = 1, cos θ₀ = 0.3 | 48.2 % | 44 / 52 / 48 / 48 |
| Nuclear confinement λ = 10, R = 300 nm | 47.7 % | unchanged: every bead already inside R |
| Nuclear confinement λ = 10, R = 500 nm | 47.7 % | unchanged |

**Conclusion.** Neither term changes accuracy beyond noise (±1–2 points), so both stay **off by
default**. They remain available for simulation and what-if use (`lambda_bend`, `lambda_confine`).

They are not added to the population model. Any quadratic prior on a Gaussian ensemble is absorbed by
the fitted couplings, and a hard wall would break the exact Gaussian statistics.

---

# v4 (October 2026): practice-only choices

Every number below comes from a practice dataset and is logged in `validation/tuning_v4.json`
(written by `validation/tune_practice.py`, `validation/calibration.py --practice`,
`validation/perturbation.py --practice` and `validation/predictor.py --practice`). The frozen
values are in `validation/frozen.py`. The commit that froze each value comes before the test run
that uses it, so git history shows the order.

v4 practice data: the Bintu sets above, plus Su et al. 2020 (IMR-90): chr2, 935 loci at 250 kb
steps, and its p-arm replicate. Test data are Su chr21, its replicate, the genome-scale set, and
the Bintu test sets. The cohesin-depletion test set (HCT116 chr21:34–37 Mb, 6 h auxin) is test-only.

## 5. Whole-chromosome population model (Pillar 1), Su chr2, split 0

Input: half A's contact frequencies at r_c = half A's median adjacent distance. That rule was fixed
before any run, as the imaging analogue of the app's adjacent-probability anchor. Truth: half B's
medians. Pairs: "within" means inside the v3.3 model's ≤ 400-locus tiles; "cross" means across
tiles, which only a whole-chromosome model can predict.

| Variant | Within | Cross | All pairs | CCC (all) | Fit time | Best misfit |
|---|---|---|---|---|---|---|
| Full rank, float32 | 84.9 % | 78.2 % | 80.6 % | 0.956 | 32 s | 0.0079 |
| Full rank, float64 | 84.9 % | 78.2 % | 80.6 % | 0.956 | 73 s | 0.0079 |
| 3,000 iterations | 84.9 % | 78.2 % | 80.6 % | 0.956 | 61 s | 0.0079 |
| Learning rate 0.02 | 84.8 % | 78.2 % | 80.6 % | 0.956 | 32 s | 0.0079 |
| **Rank 256 + random walk** | **84.9 %** | **78.0 %** | **80.7 %** | **0.957** | **15 s** | 0.0091 |
| Rank 128 + random walk | 85.0 % | 78.0 % | 80.6 % | 0.957 | 12 s | 0.0091 |
| Coarse start (k = 3) | 84.9 % | 78.2 % | 80.6 % | 0.956 | 34 s | 0.0079 |
| Uniform weights | 82.8 % | 76.9 % | 79.0 % | 0.955 | 32 s | 0.0204 |

**Frozen** (`WHOLE_CHROMOSOME`): rank 256 plus the random-walk term, 1,500 iterations, learning
rate 0.01, binomial weights, float32, automatic start. Rank 256 keeps a margin over 128 at twice
its speed. None of the accuracy differences except uniform weights exceed 0.2 points.

## 6. Hi-C → imaging settings (Gate 1 Hi-C input), Su chr2, split 0

Input: Rao et al. 2014 IMR-90 in situ Hi-C summed onto the imaged loci by Su et al. Same truth as
§5. Columns give the whole-chromosome model's score unless marked.

| Counts | Zero-count pairs | p_adjacent | Within, windowed | Within, whole | Cross, whole | All pairs, whole | Raw ρ | CCC | Size ratio |
|---|---|---|---|---|---|---|---|---|---|
| raw | clipped | 0.3 | 38.7 % | 38.6 % | 22.3 % | 29.2 % | 0.839 | 0.180 | 1.84 |
| raw | clipped | 0.5 | 33.5 % | 28.1 % | −16.5 % | 3.9 % | 0.722 | 0.112 | 2.10 |
| raw | clipped | 0.7 | 35.3 % | 31.3 % | −10.5 % | 8.6 % | 0.764 | 0.092 | 2.32 |
| raw | clipped | 0.9 | 38.5 % | 37.8 % | −2.4 % | 16.5 % | 0.831 | 0.072 | 2.70 |
| **raw** | **unobserved** | **0.3** | 44.8 % | 45.3 % | 20.8 % | **32.0 %** | 0.836 | 0.226 | 1.65 |
| raw | unobserved | 0.5 | 40.9 % | 36.1 % | −8.2 % | 12.2 % | 0.788 | 0.159 | 1.83 |
| raw | unobserved | 0.7 | 43.7 % | 41.9 % | −6.3 % | 16.0 % | 0.816 | 0.128 | 2.02 |
| raw | unobserved | 0.9 | 53.4 % | 47.8 % | 6.4 % | 26.1 % | 0.851 | 0.091 | 2.37 |
| ICE | clipped | 0.3 | 29.4 % | 28.9 % | 3.7 % | 14.4 % | 0.782 | 0.165 | 1.85 |
| ICE | clipped | 0.5 | 28.1 % | 26.4 % | −19.1 % | 1.7 % | 0.716 | 0.115 | 2.07 |
| ICE | clipped | 0.7 | 32.1 % | 32.9 % | −17.8 % | 5.8 % | 0.756 | 0.092 | 2.29 |
| ICE | clipped | 0.9 | 35.0 % | 38.2 % | −6.9 % | 13.8 % | 0.803 | 0.069 | 2.67 |
| ICE | unobserved | 0.3 | 38.3 % | 37.7 % | 15.1 % | 25.3 % | 0.829 | 0.229 | 1.64 |
| ICE | unobserved | 0.5 | 36.9 % | 35.4 % | −5.5 % | 13.4 % | 0.794 | 0.168 | 1.80 |
| ICE | unobserved | 0.7 | 38.4 % | 37.7 % | −12.4 % | 11.1 % | 0.788 | 0.124 | 1.99 |
| ICE | unobserved | 0.9 | 48.0 % | 45.7 % | −0.8 % | 21.5 % | 0.828 | 0.089 | 2.33 |

**Frozen** (`GATE1`) when 13 of the 16 variants had finished: raw counts, zero-count pairs treated
as unobserved, p_adjacent 0.3. That was the best all-pairs score. The last three variants (ICE,
unobserved, 0.5 / 0.7 / 0.9) finished afterwards, all lower (13.4 / 11.1 / 21.5 %), so the
choice is unchanged.

Absolute size does not transfer. With the literature b₀ anchor the model is about 1.65× too large
here. A practice-calibrated anchor (measured adjacent distance / literature b₀ on chr2 =
461.0 / 146.2 nm = 3.153) is kept as a second, labelled anchor for the test; only CCC depends on it.

**Reading.** Hi-C ranks transfer much worse than imaging contacts on practice data: 32 % against
81 % of the ceiling on all pairs. p_adjacent mostly moves the cross-tile score.

## 7. Structural variants: spring-network surgery dropped

The first SV simulator rewired the fitted ensemble's spring network, with couplings
K = pinv(−½ J S J), and re-sampled. On synthetic planted ideal Gaussian chains, whose true K is
exactly zero beyond |i − j| = 1 (`validation/network_check.py` → `network_surgery_check.json`),
the fitted couplings beyond 5 beads have median |K| 0.5–0.9× the backbone median. Their
95th percentile is 1.7–4.1×. Rewiring that network moves noise, not structure.

It was replaced by the exact covariance-space construction now in `chronocell/perturb.py`. Pieces
keep their joint law, and only junction displacements are new bonds. It reduces to concatenation
for an ideal chain (tests: `tests/test_v4_perturb.py`). Its real-data test is Gate 4 (SV),
`validation/RESULTS.md`.

## 8. Cohesin-depletion parameters (Gate 4), practice pair HCT116 chr21:28–30 Mb untreated → 6 h auxin

`perturb.cohesin_loss_map`: log d′ = log d + shift(s) − λ · residual. Fitted by least squares of
log d_pred against the measured auxin medians, pooled over 3 splits
(`validation/perturbation.py --practice` → `chronocell/data/perturbation_params.json`):
λ = 0.448, shift(s) = −2.346 + 0.952 L − 0.0863 L², L = log₁₀ s, clamped to 30 kb – 1.92 Mb.

In-sample (practice; these are fits, not tests):

| | Change agreement (Spearman) | CCC vs auxin | Raw ρ vs auxin |
|---|---|---|---|
| Full model (shift + λ) | 0.112 | 0.926 | 0.909 |
| Trend only (λ = 0) | 0.312 | 0.915 | 0.895 |
| No change | — | 0.646 | 0.896 |

**Reading.** On its own training pair, the λ term does **not** improve agreement with the measured
change (0.112 against 0.312). The parameters were fitted to absolute auxin distances, not to the
change. The held-out test region, which has more domain structure, behaves differently
(`validation/RESULTS.md`, Gate 4). The effect is region-dependent and rests on one test region.

## 9. Per-bead reliability (Pillar 2): can the input alone say which beads will be wrong?

Spearman between each candidate score (higher = more reliable) and minus the bead's held-out
error, split 0:

| Dataset | a. misfit-based (`bead_reliability`) | b. sampling noise | c. product |
|---|---|---|---|
| K562 28–30 | +0.027 | +0.046 | +0.058 |
| HCT116 28–30 | +0.145 | −0.253 | −0.070 |
| HCT116 28–30 auxin | +0.034 | +0.240 | +0.092 |
| HCT116 34–37 | −0.364 | +0.251 | −0.039 |
| Su chr2 | +0.331 | −0.519 | +0.130 |
| Su chr2 p-arm rep. | +0.197 | +0.270 | +0.289 |

**Conclusion.** No candidate is consistently positive. Score a is kept in the API as
"fit consistency" (how well the bead's contacts were reproduced). It is **not** presented as a
reliability or accuracy estimate. The held-out test is in RESULTS.md, Gate 2.

## 10. Interval recalibration (Gate 2), fitted on practice

Quantile recalibration (Kuleshov, Fenner & Ermon, ICML 2018) on the PIT values of held-out
single-copy distances. Each practice dataset is weighted equally; the result is frozen in
`chronocell/data/calibration.json`. Coverage of the stated 50 / 80 / 90 % intervals, in-sample:

| Dataset | Raw | Recalibrated |
|---|---|---|
| K562 28–30 | 40 / 68 / 79 % | 47 / 78 / 88 % |
| HCT116 28–30 | 41 / 67 / 78 % | 48 / 77 / 88 % |
| HCT116 28–30 auxin | 39 / 66 / 77 % | 47 / 77 / 88 % |
| HCT116 34–37 | 46 / 74 / 85 % | 52 / 82 / 91 % |
| Su chr2 | 46 / 75 / 85 % | 53 / 83 / 93 % |
| Su chr2 p-arm rep. | 48 / 77 / 87 % | 54 / 84 / 93 % |

The raw intervals are too narrow everywhere. Cells vary more than the Gaussian model allows, as
v3.3 already showed with its fixed CV of 0.42.

## 11. No-Hi-C predictor (Gate 5), practice

`chronocell/predict.py`: CTCF peaks (ENCODE IDR) oriented by the JASPAR MA0139.1 motif, plus GC,
give 12 pair features. Ridge regression of the trend-removed log median distance. The CTCF motif
reached relative score ≥ 0.8 in 89 % (K562), 69 % (HCT116) and 94 % (IMR-90) of peaks.

Held-out trend-removed Spearman, leave-one-dataset-out over the four untreated practice datasets:

| Ridge λ | K562 28–30 | HCT116 28–30 | HCT116 34–37 | Su chr2 | Mean |
|---|---|---|---|---|---|
| 0.0001 | +0.251 | +0.130 | +0.414 | +0.049 | +0.211 |
| 0.001 | +0.251 | +0.131 | +0.410 | +0.054 | +0.212 |
| **0.01** | +0.254 | +0.140 | +0.391 | +0.077 | **+0.215** |
| 0.1 | +0.251 | +0.153 | +0.335 | +0.051 | +0.197 |
| 1 | +0.228 | +0.138 | +0.300 | +0.042 | +0.177 |

**Frozen** (`validation/predictor_model.json`): λ = 0.01, refitted on all four. Trend:
log d = 2.720 + 0.610 log₁₀ s (d in nm). The largest standardised coefficients are CTCF sites
between the loci (+0.084, falling with separation: −0.071 × log₁₀ s) and GC similarity (−0.034).
The orientation terms are small (|β| ≤ 0.024).
