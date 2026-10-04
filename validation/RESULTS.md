# ChronoCell-5D: what was tested against real data, and how it came out

This file records every held-out test of ChronoCell-5D, including the ones it failed. Every number
in a table is generated from a result file by `python validation/report.py`. The text explains
them.

**Research software, not a diagnostic.** Nothing here is a clinical or regulatory validation.

## Rules every test follows

- **Two scores, never mixed.**
  - *Contact-map fit* compares the model with its own input. It shows the fit converged and is not
    evidence of accuracy.
  - *Microscopy accuracy* compares the model's distances with measurements it never saw.
  - The population's cell-to-cell spread is reported as *ensemble consistency*, never as accuracy.
- **Held out.** Imaged chromosome copies are split at random into halves A and B. The model sees
  only half A's contact frequencies (or sequencing Hi-C). Half B's measured median distances are
  the answer key.
- **Score.** Trend-removed Spearman ρ (agreement beyond "further along the DNA = further apart")
  as a percentage of the *ceiling*, which is how well half A's own medians agree with half B.
  Raw Spearman, Lin's concordance (CCC, absolute size in nm) and the median size ratio are also
  given.
- **Practice vs test.** Every setting was chosen on practice datasets
  (`validation/datasets.py`, role `practice`) and frozen in `validation/frozen.py` before the test
  datasets were run once. The tuning record is `validation/TUNING.md`.
- **Synthetic data** (the reference chr22 model, demo patients) are never used for an accuracy
  claim. They appear only in cost measurements and are labelled there.

## Data

| Data | Role | Cell line | Source and citation | Licence |
|---|---|---|---|---|
| Chromatin tracing, chr21 2–2.5 Mb regions, 30 kb steps | practice and test (see roles) | IMR-90, A549, K562, HCT116 (± auxin) | Bintu et al., *Science* 362, eaau1783 (2018); github.com/BogdanBintu/ChromatinImaging | no licence file (all rights reserved by default): downloaded on demand, not redistributed |
| Sequential-hybridisation tracing, chr21 (651 loci) and chr2 (935 loci); genome-scale DNA-MERFISH (1,041 loci) | chr2: practice; chr21 and genome: test | IMR-90 | Su et al., *Cell* 182, 1641 (2020); Zenodo 3928890 | CC-BY-4.0, downloaded on demand |
| In situ Hi-C (MAPQ ≥ 30) | input (with the tracing truth); GM12878 and K562 for the SV test | IMR-90, K562, GM12878 | Rao et al., *Cell* 159, 1665 (2014); GEO GSE63525 (read remotely by region) | public GEO data, citation required |
| CTCF ChIP-seq IDR peaks (GRCh38) | Gate 5 input | IMR-90, A549, K562, HCT116 | ENCODE (ENCFF670ULH, ENCFF624ZSR, ENCFF582SNT, ENCFF470EAN); ENCODE Project Consortium, *Nature* 583, 699 (2020) | public ENCODE data, citation requested |
| hg38 sequence (chr2, chr21) | Gate 5 input | — | UCSC hg38 (GRCh38) | freely available |
| CTCF motif MA0139.1 | Gate 5 input | — | JASPAR 2024, *Nucleic Acids Res* 52, D174 (2024) | CC BY 4.0 |

Every downloaded file is checked against the checksum its source publishes, and recorded with its
SHA-256 in `validation/data_manifest.json`. The data themselves are not committed.

## Summary

| Gate | Question | Result | Verdict |
|---|---|---|---|
| 1 | Does the whole-chromosome model at least match the windowed v3.3 model on the same regions? | Su chr21, 651 loci: 94.4 vs 94.4 % of the ceiling, and 96.2 vs 96.7 % on the replicate. It also predicts the cross-window pairs (89–93 %) that v3.3 cannot. | **Matches**, but slightly below: by 0.05–0.06 and 0.47–0.54 points, in every split |
| 1b | Do sequencing Hi-C contacts give the right distances? | Ranks: 80–85 % of the ceiling. Absolute sizes: no (CCC 0.22–0.33). | **Ranks yes, nanometres no** |
| 2 | Are the stated intervals honest? | Raw 90 % intervals cover 71–84 %. After recalibration fitted on practice data: 83–91 %. Short separations are worst. Per-bead reliability does not predict error. | **Recalibrated intervals usable for imaging-derived input; not for Hi-C input; no per-bead reliability** |
| 2b | Are they honest with sequencing Hi-C input? | Stated 90 % intervals hold 23–51 %; 43–73 % after a recalibration fitted on practice Hi-C | **Fail: with Hi-C input the range is the model's spread, not a 90 % range** |
| 2c | Can a per-pair score from the input say which distances are wrong? | Imaging input: weak (pairs +0.07 to +0.20, below the 0.20 bar); Hi-C input: none (−0.05 to +0.03) | **Fail: no reliability score in the app** |
| 3 | Benchmark against baselines and published tools | see below | see below |
| 4 | Does the cohesin-loss prediction match real RAD21 depletion? | Held-out region: change agreement 0.868 vs 0.336 for a trend-only shift | **Pass, on one region** |
| 4b | Does the SV simulator predict a real rearranged genome? | K562 chr9 deletions: 0.083 vs 0.149 (distance shift) vs 0.424 (no change) | **Not validated: mechanism simulator** |
| 5 | Can sequence + CTCF predict distances with no contact data? | 12–27 % of the ceiling on 4 of 5 test sets (0.5 % on the weak-structure set). Raw ρ gains over the separation baseline are ≤ 0.005. | **Pass (pre-registered rule), modest: a prior, not a substitute for contacts** |

## Gate 1 — whole chromosome at once (Pillar 1)

**Question.** v3.3 fits a population model to windows of at most 400 beads. Does the v4 model,
which fits a whole chromosome at once (`chronocell/population.py`), at least match it on the same
regions?

**Test** (`python validation/gate1.py --test`, run once after `GATE1` was frozen). Su et al. chr21,
651 loci at 50 kb, IMR-90, held out, plus its replicate; 3 random splits each. The windowed model
fits the fewest equal tiles of ≤ 400 loci. "Within tiles" compares the two models on the same
pairs. "Cross tiles" pairs can only be predicted by the whole-chromosome model.

<!-- BEGIN generated:gate1 -->
**su_chr21** (651 loci, 7541 copies; mean over 3 splits; test)

| Input | Pairs | Model | % of ceiling | Raw ρ | CCC (nm) | Size ratio |
|---|---|---|---|---|---|---|
| imaging | within tiles | windowed_v3_3 | 94.4 % | 0.967 | 0.934 | 0.93 |
| imaging | within tiles | whole_v4 | 94.4 % | 0.967 | 0.933 | 0.93 |
| imaging | within tiles | no_3d_direct_inversion | 92.2 % | 0.955 | 0.927 | 0.94 |
| imaging | within tiles | genomic_distance_only | -0.1 % | 0.731 | 0.751 | 1.00 |
| imaging | cross tiles | whole_v4 | 89.2 % | 0.906 | 0.833 | 0.94 |
| imaging | cross tiles | no_3d_direct_inversion | 82.8 % | 0.848 | 0.789 | 0.95 |
| imaging | cross tiles | genomic_distance_only | -0.0 % | 0.568 | 0.577 | 1.00 |
| imaging | all pairs | whole_v4 | 91.7 % | 0.954 | 0.921 | 0.93 |
| imaging | all pairs | no_3d_direct_inversion | 87.5 % | 0.929 | 0.901 | 0.94 |
| imaging | all pairs | genomic_distance_only | -0.1 % | 0.749 | 0.796 | 1.00 |
| hic | within tiles | windowed_v3_3[literature_b0] | 81.1 % | 0.909 | 0.287 | 0.53 |
| hic | within tiles | whole_v4[literature_b0] | 81.3 % | 0.909 | 0.286 | 0.53 |
| hic | within tiles | windowed_v3_3[practice_calibrated_b0] | 81.1 % | 0.909 | 0.306 | 1.68 |
| hic | within tiles | whole_v4[practice_calibrated_b0] | 81.3 % | 0.909 | 0.306 | 1.68 |
| hic | cross tiles | whole_v4[literature_b0] | 90.2 % | 0.928 | 0.218 | 0.61 |
| hic | cross tiles | whole_v4[practice_calibrated_b0] | 90.2 % | 0.928 | 0.121 | 1.93 |
| hic | all pairs | whole_v4[literature_b0] | 84.7 % | 0.930 | 0.328 | 0.57 |
| hic | all pairs | whole_v4[practice_calibrated_b0] | 84.7 % | 0.930 | 0.244 | 1.80 |

Whole − windowed, imaging input, within tiles, per split: -0.06, -0.05, -0.06 points.

**su_chr21_rep** (651 loci, 4539 copies; mean over 3 splits; test)

| Input | Pairs | Model | % of ceiling | Raw ρ | CCC (nm) | Size ratio |
|---|---|---|---|---|---|---|
| imaging | within tiles | windowed_v3_3 | 96.7 % | 0.977 | 0.924 | 0.93 |
| imaging | within tiles | whole_v4 | 96.2 % | 0.974 | 0.920 | 0.93 |
| imaging | within tiles | no_3d_direct_inversion | 95.6 % | 0.971 | 0.922 | 0.93 |
| imaging | within tiles | genomic_distance_only | -0.0 % | 0.682 | 0.696 | 1.02 |
| imaging | cross tiles | whole_v4 | 93.1 % | 0.920 | 0.781 | 0.91 |
| imaging | cross tiles | no_3d_direct_inversion | 89.5 % | 0.880 | 0.772 | 0.91 |
| imaging | cross tiles | genomic_distance_only | -0.0 % | 0.501 | 0.492 | 0.98 |
| imaging | all pairs | whole_v4 | 94.7 % | 0.960 | 0.892 | 0.92 |
| imaging | all pairs | no_3d_direct_inversion | 92.7 % | 0.945 | 0.886 | 0.92 |
| imaging | all pairs | genomic_distance_only | -0.0 % | 0.702 | 0.745 | 1.00 |
| hic | within tiles | windowed_v3_3[literature_b0] | 75.6 % | 0.875 | 0.259 | 0.53 |
| hic | within tiles | whole_v4[literature_b0] | 75.7 % | 0.875 | 0.258 | 0.53 |
| hic | within tiles | windowed_v3_3[practice_calibrated_b0] | 75.6 % | 0.875 | 0.283 | 1.67 |
| hic | within tiles | whole_v4[practice_calibrated_b0] | 75.7 % | 0.875 | 0.283 | 1.67 |
| hic | cross tiles | whole_v4[literature_b0] | 87.9 % | 0.906 | 0.221 | 0.62 |
| hic | cross tiles | whole_v4[practice_calibrated_b0] | 87.9 % | 0.906 | 0.112 | 1.95 |
| hic | all pairs | whole_v4[literature_b0] | 80.4 % | 0.907 | 0.312 | 0.58 |
| hic | all pairs | whole_v4[practice_calibrated_b0] | 80.4 % | 0.907 | 0.223 | 1.81 |

Whole − windowed, imaging input, within tiles, per split: -0.54, -0.47, -0.48 points.
<!-- END generated:gate1 -->

**Reading.**
- **Gate 1 passes in substance, not strictly.** Within tiles the whole-chromosome model is
  0.05–0.06 points below the windowed model on chr21 and 0.47–0.54 points below on the replicate.
  It is below in every split. The gap is small and consistent. We report it as "matches, slightly
  below", not "beats".
- **What the whole-chromosome fit adds** is the cross-window pairs: 89.2 % and 93.1 % of the
  ceiling, where v3.3 has no prediction. It also beats inverting each contact frequency on its own
  (no 3D) on every pair set.
- **Hi-C input (Gate 1b).** Sequencing Hi-C, binned on the imaged loci, gives the right ranks:
  80–85 % of the ceiling on all pairs, against 92–95 % from imaging-derived contacts. It does not
  give the right nanometres. With the literature b₀ the model is about 0.57× the measured size.
  With an anchor calibrated on practice data it is 1.8×, overcorrected. CCC is 0.22–0.33 either
  way. **Absolute distances from Hi-C input are not calibrated, and the app says so.**
- **Ensemble consistency** (cell-to-cell CV of the model, 0.43 on chr21) is reported next to these
  numbers and is not an accuracy.

## Gate 2 — are the stated intervals honest? (Pillar 2)

**What is stated.** For every pair of loci the population predicts a distribution of distances
across cells: a Maxwell law with per-axis σ_ij, exact from the fitted ensemble. Hence a mean, SD,
median and central 50 / 80 / 90 % interval (`population.pair_summary`), shown in the app's distance
probe.

**Test** (`python validation/calibration.py --test`, run once). Does the stated interval contain
that share of half B's *single-copy* distances? Recalibration (Kuleshov et al., ICML 2018, quantile
recalibration on the PIT histogram) was fitted on the six practice datasets, each weighted equally.
It was frozen in `chronocell/data/calibration.json` before this run. The last column is the Spearman
correlation between the per-bead reliability score (input-only) and minus the bead's held-out
error.

<!-- BEGIN generated:gate2 -->
test (run once; recalibration frozen); recalibration fitted on bintu_k562_28_30, bintu_hct116_28_30, bintu_hct116_28_30_auxin, bintu_hct116_34_37, su_chr2, su_chr2_parm_rep.

| Test dataset | Loci | Model | Stated 50 / 80 / 90 %: raw | Recalibrated | Reliability vs error (ρ) |
|---|---|---|---|---|---|
| bintu_imr90_28_30 | 65 | ensemble_v3_3 | 46 / 75 / 84 % | 52 / 82 / 91 % | +0.04 |
| bintu_imr90_18_20 | 65 | ensemble_v3_3 | 35 / 60 / 71 % | 41 / 71 / 83 % | -0.16 |
| bintu_a549_28_30 | 65 | ensemble_v3_3 | 42 / 69 / 79 % | 48 / 78 / 88 % | +0.14 |
| bintu_hct116_34_37_auxin | 83 | ensemble_v3_3 | 42 / 70 / 81 % | 49 / 79 / 90 % | +0.17 |
| su_chr21 | 651 | population_v4 | 43 / 71 / 82 % | 50 / 80 / 90 % | +0.02 |
| su_chr21_rep | 651 | population_v4 | 42 / 69 / 80 % | 49 / 79 / 90 % | +0.28 |

Coverage of the stated 90 % interval by genomic separation (raw → recalibrated):

| Test dataset | 0–0.1 Mb | 0.1–0.3 Mb | 0.3–1 Mb | 1–3 Mb | 3–10 Mb | 10–300 Mb |
|---|---|---|---|---|---|---|
| bintu_imr90_28_30 | 72 → 81 % | 80 → 88 % | 86 → 92 % | 90 → 94 % | — | — |
| bintu_imr90_18_20 | 64 → 76 % | 71 → 83 % | 72 → 84 % | 71 → 84 % | — | — |
| bintu_a549_28_30 | 72 → 80 % | 76 → 86 % | 81 → 90 % | 80 → 89 % | — | — |
| bintu_hct116_34_37_auxin | 75 → 84 % | 83 → 91 % | 83 → 91 % | 78 → 89 % | — | — |
| su_chr21 | 67 → 74 % | 71 → 79 % | 76 → 84 % | 81 → 89 % | 82 → 90 % | 82 → 91 % |
| su_chr21_rep | 65 → 74 % | 70 → 79 % | 75 → 84 % | 81 → 89 % | 81 → 90 % | 80 → 90 % |
<!-- END generated:gate2 -->

**Reading.**
- **The raw intervals are too narrow on every test set.** A stated 90 % interval holds 71–84 % of
  real single-cell distances. Real cells vary more than a Gaussian population allows.
- **Recalibration helps, and its parameters were fitted only on practice data.** The 90 % interval
  then holds 83–91 %, and the 50 % interval 41–52 %. It is applied in the app (probe interval)
  only to imaging-like contact input.
- **Where the intervals cannot be trusted:**
  - close pairs (< 100 kb): 74–84 % even after recalibration;
  - the low-structure IMR-90 18–20 Mb region: 83 % at the 90 % level;
  - **Hi-C input**: coverage is far below nominal, even after a Hi-C-specific recalibration
    (Gate 2b).
- **Per-bead reliability does not predict per-bead error** on test data (ρ from −0.16 to +0.28).
  The score is shown only as "fit consistency" and is not offered as a reliability. A second,
  per-pair attempt is Gate 2c below.

## Gate 2b — intervals with sequencing Hi-C input

**Question.** Gate 2 recalibrated the intervals for imaging-derived contacts. The app's usual input
is sequencing Hi-C. Does the same recalibration method, fitted on practice Hi-C input, make the
stated intervals honest for it?

**Test** (`python validation/calibration.py --test --input hic`, pre-registered in
`frozen.HIC_CALIBRATION` and committed with the practice-fitted recalibration before the run). Input:
Rao et al. 2014 Hi-C on the imaged loci, the Gate 1 Hi-C settings and the app's b₀ anchor. Truth:
half B's single-copy distances, split 0. Pass: on every test dataset the recalibrated 90 % interval
holds 83–97 % and the 50 % interval 40–60 %.

<!-- BEGIN generated:gate2b -->
test (run once; Hi-C recalibration frozen); Hi-C recalibration fitted on bintu_k562_28_30, su_chr2, su_chr2_parm_rep. Pass (pre-registered): recalibrated 90 % interval holds 83–97 % and 50 % interval 40–60 % on every test dataset.

| Test dataset | Loci | Model | Size ratio (model / measured) | Stated 50 / 80 / 90 %: raw | With the imaging recalibration | With the Hi-C recalibration | Within the rule |
|---|---|---|---|---|---|---|---|
| bintu_imr90_28_30 | 65 | ensemble_v3_3 | 0.49 | 16 / 30 / 39 % | 22 / 44 / 58 % | 40 / 58 / 66 % | no |
| bintu_imr90_18_20 | 65 | ensemble_v3_3 | 0.35 | 9 / 18 / 23 % | 13 / 26 / 36 % | 25 / 37 / 43 % | no |
| su_chr21 | 651 | population_v4 | 0.57 | 23 / 42 / 51 % | 30 / 56 / 70 % | 41 / 64 / 73 % | no |
| su_chr21_rep | 651 | population_v4 | 0.57 | 23 / 41 / 51 % | 29 / 55 / 69 % | 41 / 63 / 72 % | no |

Width of the stated 90 % interval (upper / lower bound): raw 4.71, Hi-C recalibrated 3.78. 0 of 4 test datasets within the rule. Verdict: **fail**.

Practice (in-sample, after the fit): bintu_k562_28_30 34 / 50 / 57 %; su_chr2 43 / 64 / 73 %; su_chr2_parm_rep 44 / 67 / 75 %.
<!-- END generated:gate2b -->

**Reading.**
- **Fail, on every test set.** With Hi-C input a stated 90 % interval holds 23–51 % of real
  single-cell distances. The recalibration fitted on practice Hi-C raises that to 43–73 %, still far
  below 90 %. The imaging recalibration, which the app applied to every input until now, gives
  36–70 %. The weak-structure IMR-90 18–20 Mb region is worst (model 0.35× the measured size).
- **Why, seen on practice data before the test.** The model built from Hi-C is about half the
  measured size (Gate 1b). 29–49 % of practice single-cell distances fell beyond the model's
  99.5th percentile. The recalibration works on a 200-bin histogram of those percentiles, so it
  cannot stretch the upper tail far enough. The method was kept as pre-registered, not changed
  after seeing this.
- **In the app.** For a population built from sequencing counts, the probe no longer shows the
  imaging recalibration. It shows the model's interval with the measured shortfall from this test,
  and says to read the range as the model's cell-to-cell spread, not as a 90 % range for real cells.
  Intervals of a predicted population (Gate 5) were not tested and are labelled so.

## Gate 2c — can the input say which distances will be wrong?

**Question.** Gate 2's per-bead score did not predict error. Can a per-*pair* score, computed from
the input and the fitted model only, rank which distances are off?

**Test** (`python validation/reliability.py`). Error of a pair: the absolute deviation of
log(model / measured median) from the median of that log ratio in its separation stratum (10
equal-count strata), so the known size and trend biases (Gate 1b) do not count as pair errors. Score: Spearman of the
candidate vs minus that error within each stratum, averaged; per bead, its median pair score vs its
median pair error. 95 % intervals from 200 resamples of the loci. Candidates, chosen on practice
data: the delta-method input noise of the pair's contact count, the model's misfit to the pair's
own input, and their combination. The test is pre-registered in `frozen.RELIABILITY`.

<!-- BEGIN generated:gate2c -->
Practice (all candidates; in-sample choice): pair-level stratified Spearman / bead-level Spearman.

| Practice dataset | Input | input_se | misfit | combined |
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

Test (run once): score **misfit**; pass on every test dataset: Spearman ≥ 0.20 and its 95 % interval above 0.

| Test dataset | Input | Pair level: stratified ρ [95 %] | Bead level: ρ [95 %] |
|---|---|---|---|
| bintu_imr90_28_30 | imaging | +0.112 [+0.018, +0.216] | +0.329 [+0.076, +0.581] |
| bintu_imr90_18_20 | imaging | +0.081 [-0.007, +0.135] | +0.248 [+0.013, +0.451] |
| bintu_a549_28_30 | imaging | +0.175 [+0.003, +0.273] | +0.373 [+0.097, +0.567] |
| bintu_hct116_34_37_auxin | imaging | +0.072 [+0.007, +0.131] | -0.033 [-0.204, +0.187] |
| su_chr21 | imaging | +0.130 [+0.102, +0.155] | +0.245 [+0.170, +0.312] |
| su_chr21_rep | imaging | +0.195 [+0.163, +0.221] | +0.452 [+0.383, +0.520] |
| bintu_imr90_28_30 | hic | +0.017 [-0.070, +0.097] | +0.107 [-0.143, +0.357] |
| bintu_imr90_18_20 | hic | +0.026 [-0.039, +0.106] | -0.016 [-0.285, +0.263] |
| su_chr21 | hic | -0.043 [-0.063, -0.022] | -0.099 [-0.186, -0.000] |
| su_chr21_rep | hic | -0.051 [-0.075, -0.027] | -0.049 [-0.137, +0.046] |

Verdicts: imaging pair **fail** (0 of 6); imaging bead **fail** (5 of 6); hic pair **fail** (0 of 4); hic bead **fail** (0 of 4).
<!-- END generated:gate2c -->

**Reading.**
- **Fail, by the pre-registered rule, for both input types at both levels.** No reliability score is
  shown in the app.
- **Imaging-derived input: a weak, consistent signal below the bar.** Pairs the model fits worst are
  somewhat more often wrong (pair level +0.07 to +0.20 within separation strata, every 95 % interval
  but one above 0). Per bead, 5 of 6 test sets reach the bar. The one that does not is the
  cohesin-depleted set (−0.03), so "every test set" is not met.
- **Sequencing Hi-C input, the app's usual case: no signal** (pair level −0.05 to +0.03). With Hi-C
  the model's misfit to its own input says nothing about where it is wrong, because the input itself
  departs from imaged distances (Gate 1b).
- **Multiplicity.** This is the second reliability attempt on Gate 2's test sets (the first was
  Gate 2's per-bead score). The score was chosen on practice data only, but these test sets were no
  longer untouched for this question. A pass would have needed confirming on new data; a fail does
  not.

## Gate 3 — benchmark (Pillar 3)

**Harness.** One command, `python -m validation.benchmark.run` (test) or `--practice`. Every
method gets the same input for a split and is scored against the same held-out truth, with 95 %
intervals from resampling half B's copies.

- **Baselines:** genomic distance only (a power law fitted on half A's *measured* medians, which
  the other methods never see); no 3D (each frequency inverted on its own).
- **ChronoCell:** v3.2 single structure; v3.3 windowed; v4 whole chromosome.
- **Published tools that run here:** PASTIS 0.4.0 (MDS and PM2; Varoquaux et al., *Bioinformatics*
  2014), its own code run unmodified from the official source distribution (SHA-256 checked).
- **Not run, with reasons** (no numbers are claimed for them):
  - ShRec3D: MATLAB only;
  - Chrom3D: C++ with Boost, needs lamina data;
  - 3DMax / LorDG: Java, no runtime here;
  - C.Origami: needs pyBigWig, which has no build on this Windows machine, and was trained on
    IMR-90 Hi-C including chr21, so an IMR-90 chr21 comparison would be in-sample.

<!-- BEGIN generated:gate3 -->
_validation/benchmark/results.json: not run (python -m validation.benchmark.run)._
<!-- END generated:gate3 -->

_Reading: pending (the held-out benchmark run had not finished when this was written)._

## Gate 4 — perturbations and variants (Pillar 4)

**Cohesin depletion** (`chronocell/perturb.py`, `validation/perturbation.py`). Parameters were
fitted on one practice pair (HCT116 chr21:28–30 Mb, untreated → 6 h auxin, RAD21 degraded;
`validation/TUNING.md` §8). They were tested once on the held-out region chr21:34–37 Mb. Input:
untreated half-A contacts. Truth: the auxin cells' measured medians.

**Structural variants** (`validation/sv_validation.py`, pre-registered in `frozen.SV_VALIDATION`).
K562 has two chr9 regions lost entirely, which include CDKN2A/CDKN2B (Zhou et al., *Genome Res*
2019). The test fits the population model to GM12878 (normal karyotype) Hi-C, applies the deletions,
and predicts K562's Hi-C counts for pairs that span the deletions. Nothing is fitted on K562.

<!-- BEGIN generated:gate4 -->
Cohesin depletion — bintu_hct116_34_37 -> bintu_hct116_34_37_auxin; test (held out; run once with frozen parameters).

| Prediction | Change agreement (Spearman) | Change RMSE (log) | CCC vs auxin (nm) | Raw ρ vs auxin | % of auxin ceiling |
|---|---|---|---|---|---|
| full model | 0.868 | 0.077 | 0.923 | 0.965 | 74.2 % |
| trend only | 0.336 | 0.145 | 0.827 | 0.899 | 74.5 % |
| no change | — | 0.298 | 0.724 | 0.908 | 74.5 % |
| data only reference | 0.336 | 0.145 | 0.888 | 0.910 | 77.6 % |

Structural variants — Zhou et al., Genome Res 29:472-484 (2019), K562: two chr9 regions lost entirely; window chr9:16,000,000-36,000,000 (800 bins of 25 kb); 6,475 pairs across the deletions (new separation 50 kb – 2 Mb).

| Prediction of K562 counts | Spearman [95 % block-bootstrap interval] |
|---|---|
| model | +0.083 [-0.068, +0.228] |
| distance shift | +0.149 [-0.018, +0.311] |
| no change | +0.424 [+0.329, +0.512] |
| model minus distance shift | -0.066 [-0.187, +0.052] |
| model trend removed | -0.068 [-0.233, +0.078] |

Check: K562 reads on the published deleted bins = 0.0004 of the kept bins (GM12878: 0.63). Verdict: **not validated (mechanism simulator)**.

Post hoc (after the verdict; cannot change it): K562 contacts across each junction as a fraction of a contiguous chain at the same separation (1 = joined on every copy):

| Junction (kept end – kept start) | s′ = 1 | 2 | 4 | 8 | 16 bins |
|---|---|---|---|---|---|
| 20.75 – 26.60 Mb | 0.000 | 0.000 | 0.000 | 0.000 | 0.003 |
| 28.55 – 31.62 Mb | 0.188 | 0.160 | 0.186 | 0.303 | 0.286 |
<!-- END generated:gate4 -->

**Reading.**
- **Cohesin loss: pass on the held-out region.** The predicted change agrees with the measured
  change far better than a separation-only shift (0.868 vs 0.336). The absolute auxin map is closer
  too (CCC 0.923 vs 0.724 with no change). The caveat is in TUNING.md §8: on its own practice pair
  the domain term did *not* help (0.112 vs 0.312). The effect is region-dependent, and one test
  region is not a general validation.
- **Structural variants: not validated.**
  - The published deletions are real in the data: K562 reads on the deleted bins are 0.04 % of the
    flanks.
  - The prediction is worse than a plain genomic-distance shift. Both lose to "no change".
  - The post hoc check shows why. Across the first junction, K562 has *no* contacts, so the kept
    pieces are not joined there. Across the second junction, contacts are 16–30 % of a contiguous
    chain, consistent with a join on only some copies.
  - The simulator assumes the variant list fully describes how the pieces are joined. A list of
    deleted intervals does not.
  - In the app the SV tools are labelled "mechanism simulator, not validated", with these numbers.
- **Drug-lab mechanisms** (loop-extrusion, compaction and similar what-ifs) have no matching
  perturbation data here and stay labelled "not validated".

## Gate 5 — prediction without contact data (Pillar 5)

**Model** (`chronocell/predict.py`). CTCF peaks of the same cell line (ENCODE), oriented by the
JASPAR CTCF motif, plus GC per locus give 12 pair features. A ridge regression on the trend-removed
log distance was fitted on the four untreated practice datasets (TUNING.md §11). It uses no
contact data from the target region. **Baseline:** the same training trend in genomic separation,
with no features. **Pass rule** (pre-registered): on ≥ 3 of 5 test datasets, % of ceiling > 0 with
its 95 % interval above 0, *and* raw Spearman above the baseline. The cohesin-depleted set is a
control (CTCF-anchored loops need cohesin) and is not in the rule.

<!-- BEGIN generated:gate5 -->
| Dataset | Role | Predictor: % of ceiling [95 %] | raw ρ | CCC | Baseline (trend, no data): raw ρ | Pass |
|---|---|---|---|---|---|---|
| bintu_imr90_28_30 | test | 26.7 [25.9, 28.7] | 0.933 | 0.856 | 0.930 | yes |
| bintu_imr90_18_20 | test | 0.5 [-13.4, 14.5] | 0.961 | 0.865 | 0.962 | no |
| bintu_a549_28_30 | test | 12.2 [9.0, 15.1] | 0.915 | 0.884 | 0.915 | yes |
| su_chr21 | test | 20.2 [18.7, 20.0] | 0.750 | 0.805 | 0.748 | yes |
| su_chr21_rep | test | 22.6 [21.1, 22.8] | 0.710 | 0.764 | 0.705 | yes |
| bintu_hct116_34_37_auxin | control | 26.8 [23.8, 28.6] | 0.973 | 0.769 | 0.973 | — |

4 of 5 test datasets pass both criteria (needed: 3). Verdict: **pass**.
<!-- END generated:gate5 -->

**Reading.**
- **Pass, by the pre-registered rule, 4 of 5, but modest.**
  - With no contact data at all, sequence + CTCF recover 12–27 % of the reproducible distance
    pattern beyond the separation trend on IMR-90 and A549 chr21:28–30 Mb, Su chr21 and its
    replicate. Contact data give 88–95 % (Gate 1, v3.3).
  - On the weak-structure IMR-90 18–20 Mb region it recovers nothing: 0.5 %, interval −13 to
    +14 %.
- **Overall rank agreement barely moves.** Raw ρ exceeds the separation-only baseline by
  0.000–0.005 (A549: 0.9150 vs 0.9148), so criterion (ii) is met only literally. The predictor's
  value is the pattern beyond the trend.
- **Absolute size** is reasonable (CCC 0.76–0.88), because the trend was learned from imaging.
- **The control did not behave as expected.** On cohesin-depleted cells the predictor scores
  26.8 %, as high as on untreated cells. Its signal is therefore not mainly loop extrusion: the
  CTCF-orientation coefficients are small (TUNING.md §11), and GC similarity and CTCF density
  between loci carry it. Read it as a compartment and insulation-level prior, not a loop predictor.
- **Intervals.** They are computed on a 20,000-pair subsample with resampled copies. For the
  651-locus sets the interval can sit just beside the all-pairs value (20.2 vs [18.7, 20.0]).
- **In the app** it is offered only where a window has no contact data, labelled "predicted from
  sequence + CTCF (no contact data)", with these numbers. When contacts are present the
  data-driven model is used. Blending prediction and data was not tested and is not offered.

## Cost (Pillar 1): runtime and peak memory against bead count

`python validation/scale_benchmark.py`: each run in its own process. Input: windows of the
**synthetic** reference chr22 contact map, up to the whole chromosome (5,082 beads at 10 kb). These
rows measure cost only, not accuracy.

<!-- BEGIN generated:scale -->
_validation/scale_benchmark.json: not run (python validation/scale_benchmark.py)._
<!-- END generated:scale -->

## Per-chromosome runtime (Pillar 6)

<!-- BEGIN generated:per_chromosome -->
_validation/chromosome_runtime.json: not run (python validation/chromosome_runtime.py)._
<!-- END generated:per_chromosome -->

## What is new here, and what is not

- **Not new:**
  - the maximum-entropy Gaussian polymer ensemble (HIPPS / DIMES: Shi & Thirumalai, *Phys Rev X*
    2019; *Nat Commun* 2023);
  - the tracing data (Bintu et al. 2018; Su et al. 2020) and the Hi-C (Rao et al. 2014);
  - quantile recalibration (Kuleshov et al. 2018);
  - CTCF-orientation features (the convergent-CTCF loop rule: Rao et al. 2014; Fudenberg et al.,
    *Cell Rep* 2016).
- **New in this repository:**
  - a low-rank-plus-random-walk parameterisation with exact block gradients, so the ensemble fits
    a whole chromosome (6,000 beads) on a CPU;
  - per-pair intervals with a held-out calibration test;
  - an exact covariance-space construction for rearranged chromosomes;
  - one harness that scores all of it, and published baselines, against held-out imaging, with
    the failures kept.

---

# v3.3 baseline (September 2026)

The v3.3 test results below were obtained with `python validation/validate_tracing.py`. They remain
valid for the windowed model (Bintu et al. 2018 regions, 150 nm contact radius, 3 splits).
Two limits listed at the end of this section are addressed above: the whole-chromosome
model (Gate 1) and the direct Hi-C → imaging test (Gate 1b).

## v3.3 results on the held-out test datasets (mean over 3 splits)

| Dataset (cells) | v3.2 single structure | **v3.3 population model** | 100 sampled trajectories | No 3D (direct inversion) | Ceiling ρ |
|---|---|---|---|---|---|
| IMR90 chr21:28–30 Mb (4,832) | 39 % | **88 %** | 88 % | 81 % | 0.98 |
| A549 chr21:28–30 Mb (3,941) | 54 % | **91 %** | 91 % | 88 % | 0.95 |
| IMR90 chr21:18–20 Mb (1,277) | 35 % | **54 %** (±11 points across splits) | 53 % | 42 % | 0.25 |
| **Overall (Σ model / Σ ceiling)** | **45 %** | **85.6 %** | **85.6 %** | 79.7 % | |

The unweighted mean of the three percentages is 78 %.

| Other measures, v3.3 | IMR90 28–30 | A549 | IMR90 18–20 |
|---|---|---|---|
| Contact-map fit (ρ vs input) | 0.99 | 0.99 | 0.97 |
| Raw ρ, model (trend kept) | **0.976** | **0.952** | 0.868 |
| Raw ρ, genomic-distance baseline | 0.930 | 0.915 | **0.962** |
| Lin's CCC in nm (1 = exact sizes) | 0.97 | 0.93 | 0.46 |
| Model / real size | 1.00 | 0.92 | 0.82 |
| Cell-to-cell spread (CV), model vs measured | 0.42 vs 0.50 | 0.42 vs 0.55 | 0.42 vs 0.58 |

Practice datasets, for comparison (Σ/Σ over 4 datasets × 3 splits): v3.2 48 %, v3.3 92 %.

## What the v3.3 model is

- A **population of structures**, not one fold. Every cell folds differently, and Hi-C and imaging
  average over thousands of cells.
- It is a maximum-entropy Gaussian polymer ensemble. The method follows HIPPS/DIMES (Shi &
  Thirumalai, PRX 2019; Nat Commun 2023) and is not a new method of ours.
- Each contact frequency is converted to a pair spread. A valid 3D ensemble is then fitted to all
  pairs at once, with each pair weighted by how reliably it was measured.
- 100 Langevin trajectories are drawn exactly from that ensemble.
- The prediction is the ensemble's median distance. It is computed exactly, and checked against the
  100 sampled trajectories, which give the same score.
- Code: `chronocell/ensemble.py`; tests: `tests/test_v33.py`.

## Honest reading (v3.3)

- **On regions with real structure the target is met.** The model recovers 88–91 % of the structure
  the experiment itself can reproduce, up from 39–54 % in v3.2.
  - Absolute sizes are now right: CCC 0.93–0.97, where v3.2 was about 3× too small.
  - Overall rank agreement now beats the genomic-distance baseline.
- **The 3D population adds real information.** It beats inverting each contact frequency on its own
  (no 3D) on every dataset, because fitting all pairs jointly corrects noisy ones.
- **Where the data are weak, the model is weak.** On IMR90 18–20 Mb, even the experiment's two
  halves barely agree once the trend is removed (ceiling 0.25). There the model reaches 54 %, and:
  - it swings ±11 points between splits;
  - its overall rank agreement (0.87) is **below** the genomic-distance baseline (0.96).
- **Cells vary more than the model.** Real cell-to-cell spread (CV 0.50–0.58) is larger than the
  Gaussian model's fixed 0.42. The model gets population medians right, not the full shape of each
  pair's distance distribution.
- **Limits:**
  - The input is imaging-derived contact frequencies at a known 150 nm radius, not sequencing Hi-C.
    For Hi-C the probability scale must be assumed (`ensemble.counts_to_probability`), and a direct
    Hi-C → imaging test is still to do.
  - The population model is intended for windows of up to a few hundred beads.
  - Sampled single structures have no excluded volume.
