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

## 12. Per-pair reliability (Gate 2c), practice

Gate 2's per-bead score did not rank held-out error (section 9). Three per-pair scores, computed from the
input and the fitted model only, were compared on the practice datasets (`validation/reliability.py --practice`
→ `results_reliability_practice.json`, split 0). Pair level: Spearman of the score vs minus the pair's
scale-free error within 10 separation strata, averaged; bead level: median pair score vs median pair error.

| Practice dataset | Input | input_se: pair / bead | misfit: pair / bead | combined: pair / bead |
|---|---|---|---|---|
| bintu_k562_28_30 | imaging | +0.191 / +0.331 | +0.124 / +0.253 | +0.178 / +0.348 |
| bintu_hct116_28_30 | imaging | +0.087 / +0.216 | +0.063 / +0.035 | +0.083 / +0.219 |
| bintu_hct116_28_30_auxin | imaging | +0.266 / +0.534 | +0.175 / +0.307 | +0.245 / +0.430 |
| bintu_hct116_34_37 | imaging | +0.009 / -0.281 | +0.067 / -0.147 | +0.066 / -0.208 |
| su_chr2 | imaging | -0.031 / -0.108 | +0.057 / +0.262 | +0.049 / +0.159 |
| su_chr2_parm_rep | imaging | +0.077 / +0.152 | +0.099 / +0.241 | +0.099 / +0.263 |
| bintu_k562_28_30 | hic | -0.079 / -0.040 | -0.080 / +0.194 | -0.095 / -0.015 |
| su_chr2 | hic | -0.048 / +0.002 | -0.020 / +0.019 | -0.024 / +0.009 |
| su_chr2_parm_rep | hic | -0.104 / -0.160 | -0.014 / +0.072 | -0.026 / +0.042 |

**Frozen** (`frozen.RELIABILITY`): the misfit score, which has the best worst case: on imaging input its pair-level minimum is +0.057 (combined +0.049, input_se −0.031) and its bead-level minimum −0.147 (combined −0.208, input_se −0.281), and on Hi-C practice it is the least negative on Su chr2 and its replicate (on K562, −0.080 against −0.079 for input_se). No candidate is positive on practice Hi-C input; the misfit score is tested there unchanged. Pass bar: ρ ≥ 0.20 on every test set with its 95 % interval above 0 (the usability bar report.py already applied to Gate 2's per-bead score), per input type and level.

**Reading.** Even the chosen score reaches the bar on none of the practice sets at pair level. A weak positive ranking on imaging-derived input is the best case the test can confirm.

## 13. Hi-C size calibration (Gate 1c, Phase A1), practice

`validation/hic_size_calibration.py --practice` → `results_hic_size_practice.json`. Model distances are multiplied
by exp(h), h from a least-squares fit of log(measured / model median) with equal weight per practice unit; six
forms of increasing complexity (one global factor; + separation; + separation²; + locus spacing; + depth (from
1/4 and 1/16 binomial thinning of the same maps); + protocol (in situ / intact)) were compared by
leave-one-dataset-out over five practice groups (Bintu K562, HCT116 28–30 and 34–37, Su chr2 with its
replicate, the genome-scale practice set). The practice table is in `RESULTS.md` (Gate 1c, generated).

**Frozen** (`chronocell/data/hic_size_calibration.json`): the single global factor (h = 0.895, ×2.45), the best
mean held-out CCC (0.688); every richer form extrapolated worse to the held-out group, most on Su chr2.

## 14. Distance ranges (Gate 2d, Phase A2), practice

`validation/intervals_v2.py --practice` → `results_intervals_v2_practice.json`. Split-conformal quantiles of
log(single-copy distance / model median) from half B's copies, compared four ways (per separation band or
pooled; with or without the Gate 1c size calibration) by leave-one-dataset-out; criterion: the smallest
worst-case |coverage at 90 % − 0.90| over held-out groups and bands, and within 0.005 of it the simpler variant.

| Variant | Imaging input: worst |cov90 − 0.90| | Hi-C input |
|---|---|---|
| per band | **0.066** (chosen) | 0.282 |
| pooled | 0.074 | **0.175** (chosen) |
| per band + size calibration | — | 0.281 |
| pooled + size calibration | — | 0.175 |

The size calibration is a global factor, so it cancels in the ratio (identical to the pooled variant within
binning); the tie rule keeps it out. Zero distances (identical rounded coordinates) count as below every range.

## 15. Per-pair reliability (Gate 2e, Phase A3), practice

`validation/reliability_v2.py --practice` (two chunks: imaging, then Hi-C; refits cached per unit) →
`results_reliability_v2_practice.json`; the table is generated in `RESULTS.md` (Gate 2e). Candidates: refit
spread (8 refits on resampled input), local evidence, their rank average, and Gate 2c's misfit.

**Frozen** (`frozen.RELIABILITY_V2`), by the best worst case per input: imaging input, misfit (the only candidate
positive on every practice group, +0.063 to +0.175); Hi-C input, the refit spread (worst −0.031; every
candidate is near zero on Hi-C input). No candidate reached the 0.30 bar on any practice group.

## 16. Learned correction (Gate 3b, Phase A4), practice

`validation/learned_correction.py --practice` → `results_learned_correction_practice.json` (generated in
`RESULTS.md`, Gate 3b). Ridge and a two-layer network (trained on the GPU) on pair features, leave-one-dataset-out.
Both raise Lin's CCC (0.42 → 0.68 / 0.69) and lower the held-out trend-removed ρ (0.709 → 0.654 / 0.682).

**Frozen**: no correction (the selection rule needs a better held-out pattern on every group). Gate 3b not run.

## 17. Prediction from more inputs (Gate 5b, Phase A5), practice

`validation/predictor_v2.py --practice` → `results_predictor_v2_practice.json` (generated in `RESULTS.md`, Gate
5b). Feature sets CTCF + GC, plus ATAC, H3K27ac, RAD21 or all three; ridge λ by leave-one-dataset-out over the
three non-IMR-90 untreated practice regions.

**Frozen** (`validation/predictor_v2_model.json`): CTCF + GC + RAD21 peaks, λ = 0.01 (held-out ρ +0.432 against
+0.294 for CTCF + GC). ATAC and H3K27ac did not help.

## 18. Phase B: settings chosen on practice data

- **Gate 4c (cohesin loss against Hi-C).** Nothing fitted: the Gate 4 parameters are used unchanged. The practice
  run (chr21:28–30 Mb, `results_cohesin_hic_practice.json`) only checked the pipeline; resolution (10 kb),
  window (2 Mb), read threshold (10) and the main pair (in situ) were fixed in the code before it ran.
- **Gate 6 (loops).** On three GM12878 windows (`results_gate6_practice.json`, generated in `RESULTS.md`): FDR 0.05 /
  0.1 / 0.2 per neighbourhood, raw or KR-balanced counts. **Frozen**: FDR 0.2 with KR balancing (precision 0.889,
  recall 0.533, F1 0.667); raw counts give more recall and much lower precision.
- **Gate 7 (differential).** On chr21:28–30 Mb (`results_gate7_practice.json`, generated in `RESULTS.md`): minimum
  mean count 5 or 10, distance normalisation on or off. **Frozen**: 5 reads with distance normalisation (mean FDP
  0.036 at nominal 0.05, recall 0.16 at ×2 and 0.94 at ×4; no discovery without a planted change).

## 19. Quantum lab (Gate Q): settings chosen on practice data

Everything "quantum" runs on ChronoCell's statevector simulator (`chronocell/quantum/sim.py`), which agrees with
Qiskit 2.2.1 on every circuit type the lab builds (`results_quantum_crosscheck.json`). Practice data only: GM12878
(Gate 6's practice map, chr1/2/3:100–110 Mb) and its ENCODE Arrowhead domains; the chemistry needs no data. The
tables are generated in `RESULTS.md` (block `gateq_practice`).

- **Domain QUBO (Q1/Q2).** 360 settings scored by the exact optimum of every 20-bin window (dynamic programming over
  the QUBO's band, equal to full enumeration): bin size 40 or 50 kb; minimum domain 2–4 bins; resolution γ; a cost
  per boundary (0–3, in units of the mean linear coefficient); weights as observed − γ·expected or as a log ratio.
  The first grid (no boundary cost) peaked at F1 0.279, below the insulation caller; adding the boundary cost and the
  wider γ range (still practice only) raised it. **Frozen**: 40 kb, minimum 3 bins, γ 1.5, boundary cost 1.5,
  difference weights (F1 0.400, precision 0.542, recall 0.317). The classical callers were given the same chance
  (grids widened alongside): insulation w 4, depth 0.25 (F1 0.395); TopDom-like window 3 (F1 0.303). Practice
  agreement with Arrowhead is low for every method (sparse reference: 6–14 domains per 10 Mb).
- **QAOA (Q1).** Depth 3 or 6, expectation or CVaR (α 0.1) objective, 4,096 shots, COBYLA 80 iterations per depth
  with INTERP growth, on the frozen QUBO's 33 practice windows (19 qubits each). Every setting sampled the optimum in
  33 of 33 windows. **Frozen**: p = 3, expectation (the lowest depth with the highest mean probability of the optimum,
  0.119; a uniform guess 1.9 × 10⁻⁶). Reported alongside: under the approximate noise model only 0.6 % of shots
  survive a depth-3 circuit (about 470 CX), and simulated annealing also found every optimum.
- **Gene classifier (Q4).** 106 GM12878 genes, 5-fold cross-validation, the same grid width for both kernels.
  **Frozen**: quantum kernel bandwidth 0.05, 1 repetition, C 100 (AUC 0.806); RBF-SVM γ 0.02, C 100 (AUC 0.796);
  logistic regression 0.762 for reference. The first, narrower grid put the quantum bandwidth at its edge, so both
  grids were widened before choosing.
- **Chemistry (Q3).** No fitting. The from-scratch STO-3G pipeline reproduced Szabo & Ostlund's H2 integrals and
  energies (E_HF −1.1167, E_FCI −1.1373 hartree at 1.4 bohr) and HeH+ (E_HF −2.86066), and the 15-term H2 qubit
  Hamiltonian.

## 20. Drug lab and the quantum drug tabs (Gate 8, Gates Q5-Q7): settings chosen on practice data

Tables generated in `RESULTS.md` (blocks `gate8` practice part and `qdrug_practice`).

- **Gate 8 (Drug lab vs tracing after real drugs).** Nothing in the simulator was fitted: the four core classes are the
  original ones and the eight extended classes were defined from their mechanisms in the literature before any
  drug-treated trace was read. Practice used the untreated traces split in two halves (no drug: the noise level of the
  comparison) and alpha-amanitin (a class not tested). The noise reached |rho| 0.14 with intervals spanning 0, and the
  alpha-amanitin agreement (rho 0.29) was matched by shuffled targeting (permutation p 0.86): the simulator's generic
  compaction pattern, not where it acts. The rule therefore asks for rho >= 0.20, an interval above 0 and a permutation
  p <= 0.05 per drug. The transcription-inhibitor class keeps its literature direction (compaction of active chromatin):
  the measured global size grew, but global scale is not trusted between imaging experiments.
- **Q5 (hERG).** The TDC file contains commas inside compound names; the first practice run split on commas and read 637
  of its 655 rows, so the reader was changed to a CSV parser before the settings were chosen (both runs are practice).
  Feature sets: 6 Lipinski-type descriptors, 8 hERG-relevant descriptors (with basic amines, halogens, sp3 fraction), or
  8 principal components of all 17. Grids for the quantum kernel (bandwidth, 1 or 2 repetitions, C) and the RBF kernel
  (gamma, C) were widened together three times when a best value sat at an edge. **Frozen**: 8 hERG descriptors; QSVM
  bandwidth 0.005, 2 repetitions, C 1000 (CV AUC 0.862; C at the top of its grid, where it trades off against the small
  bandwidth and the AUC moves in the third decimal); RBF gamma 0.002, C 100 (0.860); logistic regression 0.855.
- **Q6 (molecules).** No fitting. The STO-3G exponents equal the Basis Set Exchange values (basis_set_exchange 0.11, every element H-Ne,
  largest relative difference 4.5e-6: rounding of the published fit); the HF
  solver needed one practice fix (DIIS from the core guess locked N2 onto an excited solution: a damped warm-up and the
  lower of two routes are kept) and VQE an exact closed-form rotation and an adjoint gradient (speed only). Practice
  anchors: Szabo & Ostlund's seven HF energies and OpenFermion's H2 curve.
- **Q7 (docking).** On 39 usable practice complexes: grid-based hot-spots docked none correctly; hot-spots placed along
  protein N-H and C=O vectors, tau 1.0 A, 8 polar + 4 hydrophobic ligand features and a 20-vertex core subgraph (the
  vertices with the most consistent pairs) docked 15 % against 3 % for random search with the same score. QAOA: p 3 or 5,
  expectation or CVaR; **frozen** p 5, CVaR (best clique found in 97 % of practice graphs).

## 21. Round 2 of the quantum gates (Q2b, Q4b, Q6b, Q7b): settings chosen on practice data

Tables generated in `RESULTS.md` (block `round2_practice`). The data of the failed tests are practice here.

- **Q6b (molecules).** ADAPT-VQE with two operator pools on the six Q6 cases and nine more (equilibrium geometries, N2 at
  2x): the occupied-to-virtual pool stalls at a stationary point for stretched N2 (45 and 128 mEh); the generalized pool
  reached every case within 0.19 mEh. **Frozen**: generalized pool, gradient-norm threshold 1e-3; the cap of 100
  operators was reached twice on practice (still within 0.2 mEh) and was raised to 150 so that it binds less often.
- **Q4b (genes).** Labels: ENCODE CSHL long total RNA-seq of each cell line (one lab, one protocol), TPM >= 1. Features:
  Q4's five plus the TSS bin's local observed / expected contacts and the insulation slope across the TSS. Settings by
  leave-one-cell-line-out AUC over GM12878, K562 and IMR-90; both kernels' grids widened once together when the first
  run's best values sat at their edges. **Frozen**: seven features; QSVM bandwidth 0.005, 1 repetition, C 10,000;
  RBF-SVM gamma 0.0005, C 10,000 (both still at the edge where C and the kernel width trade off; the AUC moves in the
  third decimal there).
- **Q2b (domains).** Gate Q's domain settings were chosen on three GM12878 windows holding 41 reference boundaries
  and called too few boundaries on the test cell lines (recall 0.14-0.17). Here: every setting's exact optimum on all
  nine windows already seen (GM12878, K562, IMR-90), 672 settings (resolution 40-80 kb, minimum domain 2-4 bins, gamma,
  boundary cost, difference or log weights; gamma widened once when the best sat at 5.0, without moving the choice);
  the classical callers re-tuned on the same windows at each resolution (insulation: window and depth; TopDom-like:
  window), jointly over the three cell lines. Choice by the worst cell line's margin over the better classical caller.
  **Frozen**: 60 kb, minimum 3 bins, gamma 5.0, boundary cost 0.25, log weights; insulation (2, 0.15); TopDom-like 3.

## 22. Gate 6b (loop calls): settings chosen on practice data

Gate 6 chose FDR 0.2 with KR balancing on three GM12878 windows; on K562 that made 315 calls for 156 reference loops.
Here: the nine windows Gate 6 had used (GM12878, K562, IMR-90), a grid over FDR (0.001-0.2; the lower end added once
when the best sat at 0.01, the choice did not move), balancing (none, KR) and two optional post-filters (minimum donut
enrichment 2.5 or 3.0; minimum cluster size 2), chromosight and Mustache run on the same windows. Choice by the worst
cell line's margin over the better tool. **Frozen**: FDR 0.01, KR balancing, no post-filter.
