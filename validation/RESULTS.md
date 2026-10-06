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
| 3 | Benchmark against baselines and published tools | Imaging input: best or tied best on 27 of 28 units (PASTIS PM2 ahead only on the weak IMR-90 18–20 Mb set). Hi-C input: best on 18 of 26; no-3D or PASTIS ahead on 7 genome-scale chromosomes and the 18–20 Mb set | **Ahead with imaging-derived input; mixed with Hi-C input**; PASTIS PM2 on the 651-locus sets not finished |
| 4 | Does the cohesin-loss prediction match real RAD21 depletion? | Held-out region: change agreement 0.868 vs 0.336 for a trend-only shift | **Pass, on one region** |
| 4b | Does the SV simulator predict a real rearranged genome? | K562 chr9 deletions: 0.083 vs 0.149 (distance shift) vs 0.424 (no change) | **Not validated: mechanism simulator** |
| 5 | Can sequence + CTCF predict distances with no contact data? | 12–27 % of the ceiling on 4 of 5 test sets (0.5 % on the weak-structure set). Raw ρ gains over the separation baseline are ≤ 0.005. | **Pass (pre-registered rule), modest: a prior, not a substitute for contacts** |

**Phase A and B gates** (October 2026; generated from the result files):

<!-- BEGIN generated:summary_ab -->
| Test (held-out, real data) | Measured | Verdict |
|---|---|---|
| Gate 1c: calibrated sizes from Hi-C (5 main sets) | CCC 0.43–0.72 after calibration (needed ≥ 0.8) | fail |
| Gate 2d: conformal distance ranges hold 85–95 % / 40–60 % in every band | imaging input 4 of 7 sets, Hi-C input 1 of 5 | imaging fail, hic fail |
| Gate 2e: per-pair reliability on untouched genome-scale sets (ρ ≥ 0.30) | su_genome · imaging +0.08; su_genome · hic +0.01; su_genome_amanitin · imaging +0.12; su_genome_amanitin · hic +0.01 | imaging fail, hic fail |
| Gate 3b: a learned correction on the population model | practice chose no correction | not run |
| Gate 4c: cohesin loss vs RAD21-degron Hi-C on 6 held-out regions | 5 of 6 regions beat trend only | pass |
| Gate 4d: SV effects on new events with Hi-C before and after | two usable events found, three needed | blocked; variant engine stays a mechanism simulator |
| Gate 5b: prediction with cohesin peaks | 0 of 5 test sets | fail |
| Gate 5m: the human predictor on mouse ES-cell tracing (4DN) | 17.0 % of the ceiling; 7.2 % of the ceiling | pass (modest) |
| Gate 6: loop calls vs ENCODE HiCCUPS calls, held-out cell lines | k562: F1 0.40 (chromosight 0.39, Mustache 0.49); imr90: F1 0.74 (chromosight 0.42, Mustache 0.47) | fail |
| Gate 7: false discoveries of the differential analysis (real replicates + planted changes) | mean FDP 0.002 at nominal 0.05; recall ×2 0.01, ×4 0.63 | pass |
<!-- END generated:summary_ab -->

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

## Gate 1c — calibrated sizes from Hi-C (Phase A1)

**Question.** Models built from sequencing Hi-C come out at 0.35–0.57× the imaged size (Gate 1b). Can a
size calibration learned on practice pairs of Hi-C and imaging bring them to the right nanometres?

**Test** (`python validation/hic_size_calibration.py --test`; rule in `frozen.HIC_SIZE`, committed with the
frozen calibration before any test unit was built). The calibration multiplies every model distance by
one factor, chosen by leave-one-dataset-out on practice units (Bintu K562 and HCT116 regions with Rao
2014, ENCODE in situ and intact Hi-C, each also thinned to 1/4 and 1/16 of its depth; Su chr2 and its
replicate; the genome-scale practice set). Main test sets: Bintu IMR-90 chr21:28–30 and 18–20 Mb, Su
chr21 and its replicate (all with Rao 2014 Hi-C), and the genome-scale set. Truth: half B's median
distances, split 0. Pass: on every main set, Lin's CCC ≥ 0.8, median size ratio 0.8–1.25, and the
pattern not lower than the uncalibrated model's.

<!-- BEGIN generated:gate1c -->
Practice (leave-one-dataset-out; CCC of the held-out group / its size ratio):

| Form | bintu_hct116_28_30 | bintu_hct116_34_37 | bintu_k562_28_30 | su_chr2 | su_genome_tx | Mean CCC |
|---|---|---|---|---|---|---|
| global (chosen) | 0.90 / 0.97 | 0.79 / 1.10 | 0.75 / 1.09 | 0.34 / 0.96 | 0.65 / 0.99 | 0.688 |
| sep | 0.89 / 1.02 | 0.71 / 1.19 | 0.69 / 1.17 | 0.31 / 0.94 | 0.64 / 1.00 | 0.649 |
| sep2 | 0.88 / 1.00 | 0.69 / 1.19 | 0.69 / 1.14 | 0.13 / 0.94 | 0.57 / 1.03 | 0.592 |
| sep2+step | 0.89 / 0.92 | 0.81 / 1.09 | 0.78 / 1.07 | -0.03 / 0.69 | 0.58 / 0.82 | 0.605 |
| sep2+step+depth | 0.88 / 0.93 | 0.84 / 1.06 | 0.83 / 1.09 | -0.04 / 0.71 | 0.44 / 0.68 | 0.591 |
| sep2+step+depth+protocol | 0.89 / 0.92 | 0.92 / 1.07 | 0.87 / 1.08 | -0.05 / 0.75 | 0.34 / 0.61 | 0.594 |

Test (run once): pass on every main set: CCC ≥ 0.8, size ratio 0.8–1.25, pattern unchanged.

| Set | Tier | CCC (nm): uncalibrated → calibrated | Size ratio | Trend-removed ρ | Within the rule |
|---|---|---|---|---|---|
| bintu_imr90_28_30 | main | 0.372 → 0.718 | 0.49 → 1.20 | 0.907 | no |
| bintu_imr90_18_20 | main | 0.040 → 0.624 | 0.35 → 0.85 | -0.049 | no |
| su_chr21 | main | 0.328 → 0.458 | 0.57 → 1.40 | 0.844 | no |
| su_chr21_rep | main | 0.311 → 0.428 | 0.57 → 1.41 | 0.798 | no |
| su_genome | main | 0.166 → 0.633 | 0.46 → 1.12 | 0.522 | no |
| bintu_imr90_28_30 · encode_imr90_intact | secondary | 0.629 → 0.446 | 0.61 → 1.49 | 0.912 | — |
| bintu_imr90_28_30 · encode_imr90_dilution | secondary | 0.126 → 0.809 | 0.34 → 0.83 | 0.671 | yes |
| bintu_imr90_18_20 · encode_imr90_intact | secondary | 0.069 → 0.825 | 0.42 → 1.04 | -0.000 | yes |
| bintu_imr90_18_20 · encode_imr90_dilution | secondary | 0.023 → 0.227 | 0.28 → 0.70 | 0.028 | — |
| bintu_a549_28_30 · encode_a549_insitu | secondary | 0.067 → 0.468 | 0.29 → 0.70 | 0.452 | — |
| su_genome_amanitin | secondary | 0.170 → 0.594 | 0.47 → 1.14 | 0.533 | — |

Verdict: **fail**.
<!-- END generated:gate1c -->

**Reading.**
- **Fail, on every main set.** The single factor (×2.45) raises CCC everywhere, from 0.04–0.37 to
  0.43–0.72, but no set reaches 0.8. Sizes land in range on three of the five (IMR-90 ×2 and the
  genome-scale set) and overshoot Su chr21 and its replicate (1.40×), whose uncalibrated models were
  already larger than the practice average.
- **One factor does not transfer between Hi-C datasets.** On the secondary sets each IMR-90 region
  reaches CCC 0.81–0.83 with one ENCODE file and only 0.23–0.45 with the other, and which file works
  swaps between the two regions (dilution Hi-C on 28–30 Mb, intact Hi-C on 18–20 Mb). The size error
  depends on the Hi-C dataset and the region in ways the practice units could not teach: every
  practice form with separation, locus-spacing, depth or protocol terms extrapolated worse than the
  single factor.
- **The pattern is untouched** (one factor for all pairs), as the rule required.
- **In the app** nothing changes: models built from Hi-C are still labelled as giving ranks, not
  nanometres (Gate 1b).

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

## Gate 2d — conformal distance ranges (Phase A2)

**Question.** Gates 2 and 2b showed the stated ranges are too narrow (Hi-C input) or uneven across separations
(imaging input). Do ranges learned on practice data from the observed spread of single-copy distances
around the model (split-conformal, in log-distance space) hold their stated coverage on held-out data?

**Test** (`python validation/intervals_v2.py --test`; rule in `frozen.INTERVALS_V2`, committed with the frozen
quantiles before the run). Practice chose, by leave-one-dataset-out, one set of quantiles per separation band
for imaging input and one pooled set for Hi-C input. Test sets: Gate 2's six (imaging input), Gate 2b's four
(Hi-C input) and the genome-scale set (both). Pass, per input type: in every set and every band with at least
200 pair-copies, the 90 % range holds 85–95 % and the 50 % range 40–60 %.

<!-- BEGIN generated:gate2d -->
Test (run once): variants imaging: bands, hic: pooled; pass on every set and band: 90 % ranges hold 85–95 %, 50 % ranges 40–60 %.

| Set · input | Coverage at 90 % / 50 % by separation band (width of the 90 % range, upper / lower) | Within the rule |
|---|---|---|
| bintu_imr90_28_30 · imaging | 0–0.1 Mb: 92 / 51 % (×9.5); 0.1–0.3 Mb: 92 / 53 % (×7.6); 0.3–1 Mb: 93 / 55 % (×6.5); 1–3 Mb: 96 / 63 % (×7.5) | no |
| bintu_imr90_18_20 · imaging | 0–0.1 Mb: 91 / 44 % (×9.5); 0.1–0.3 Mb: 90 / 45 % (×7.6); 0.3–1 Mb: 87 / 43 % (×6.5); 1–3 Mb: 89 / 46 % (×7.5) | yes |
| bintu_a549_28_30 · imaging | 0–0.1 Mb: 90 / 50 % (×9.5); 0.1–0.3 Mb: 91 / 50 % (×7.6); 0.3–1 Mb: 91 / 51 % (×6.5); 1–3 Mb: 93 / 54 % (×7.5) | yes |
| bintu_hct116_34_37_auxin · imaging | 0–0.1 Mb: 93 / 53 % (×9.5); 0.1–0.3 Mb: 95 / 56 % (×7.6); 0.3–1 Mb: 93 / 52 % (×6.5); 1–3 Mb: 93 / 52 % (×7.5) | yes |
| su_chr21 · imaging | 0–0.1 Mb: 84 / 42 % (×9.5); 0.1–0.3 Mb: 84 / 45 % (×7.6); 0.3–1 Mb: 86 / 47 % (×6.5); 1–3 Mb: 92 / 55 % (×7.5); 3–10 Mb: 94 / 55 % (×7.5); 10–300 Mb: 91 / 51 % (×6.2) | no |
| su_chr21_rep · imaging | 0–0.1 Mb: 85 / 40 % (×9.5); 0.1–0.3 Mb: 85 / 43 % (×7.6); 0.3–1 Mb: 86 / 45 % (×6.5); 1–3 Mb: 93 / 55 % (×7.5); 3–10 Mb: 93 / 54 % (×7.5); 10–300 Mb: 90 / 49 % (×6.2) | yes |
| bintu_imr90_28_30 · hic | 0–0.1 Mb: 89 / 49 % (×9.3); 0.1–0.3 Mb: 94 / 56 % (×9.3); 0.3–1 Mb: 94 / 58 % (×9.3); 1–3 Mb: 92 / 52 % (×9.3) | yes |
| bintu_imr90_18_20 · hic | 0–0.1 Mb: 84 / 38 % (×9.3); 0.1–0.3 Mb: 93 / 47 % (×9.3); 0.3–1 Mb: 94 / 50 % (×9.3); 1–3 Mb: 93 / 51 % (×9.3) | no |
| su_chr21 · hic | 0–0.1 Mb: 75 / 37 % (×9.3); 0.1–0.3 Mb: 85 / 46 % (×9.3); 0.3–1 Mb: 91 / 51 % (×9.3); 1–3 Mb: 92 / 52 % (×9.3); 3–10 Mb: 90 / 51 % (×9.3); 10–300 Mb: 89 / 47 % (×9.3) | no |
| su_chr21_rep · hic | 0–0.1 Mb: 68 / 29 % (×9.3); 0.1–0.3 Mb: 83 / 42 % (×9.3); 0.3–1 Mb: 91 / 50 % (×9.3); 1–3 Mb: 92 / 54 % (×9.3); 3–10 Mb: 90 / 51 % (×9.3); 10–300 Mb: 87 / 47 % (×9.3) | no |
| su_genome · imaging | 1–3 Mb: 83 / 44 % (×7.5); 3–10 Mb: 86 / 47 % (×7.5); 10–300 Mb: 90 / 51 % (×6.2) | no |
| su_genome · hic | 1–3 Mb: 83 / 43 % (×9.3); 3–10 Mb: 87 / 47 % (×9.3); 10–300 Mb: 93 / 54 % (×9.3) | no |

Verdicts: imaging input **fail**; hic input **fail**.
<!-- END generated:gate2d -->

**Reading.**
- **Fail for both input types.** With imaging-derived input 4 of 7 sets are inside the rule in every band; the
  misses are the IMR-90 28–30 Mb region (its widest band over-covers: 96 % / 63 %) and Su chr21 and the
  genome-scale set (short separations under-cover: 83–86 % at 90 %). With Hi-C input only 1 of 5 sets is
  inside: Su chr21 and its replicate hold 68–75 % of their 90 % ranges below 100 kb.
- **Much closer than before.** The same Hi-C sets held 23–51 % (raw) and 43–73 % (recalibrated) in Gate 2b;
  most bands now sit at 85–94 %. The remaining errors are at the shortest separations, where the model's
  spread and the cells' spread differ most.
- **Not trivially wide.** A 90 % range spans ×6–9.5 between its ends (reported per band), the real cell-to-cell
  spread of single-copy distances.
- The first run stalled when this computer was suspended and was killed before writing a result; the same
  frozen test was completed in a second, resumable run, which reproduced the stalled run's coverages exactly.
- **In the app** the probe is unchanged.

## Gate 2e — a per-pair reliability score on untouched data (Phase A3)

**Question.** Gate 2c's per-pair score ranked held-out errors weakly (imaging input) or not at all (Hi-C
input). Do the Phase A3 candidates do better on data no reliability test has touched?

**Test** (`python validation/reliability_v2.py --test`, run in chunks; rule in `frozen.RELIABILITY_V2`,
committed before the test). Scores chosen on practice per input (best worst case): imaging input, Gate 2c's
misfit; Hi-C input, the spread over 8 refits on resampled counts. Test sets: the genome-scale sets su_genome
and su_genome_amanitin. Pass, per input: ρ ≥ 0.30 with its 95 % interval above 0 on both.

<!-- BEGIN generated:gate2e -->
Practice (pair-weighted stratified ρ per group):

| Practice group · input | boot_sd | evidence | combined | misfit |
|---|---|---|---|---|
| bintu_hct116_28_30 · imaging | +0.138 | +0.092 | +0.129 | +0.063 |
| bintu_hct116_28_30_auxin · imaging | +0.232 | +0.213 | +0.241 | +0.175 |
| bintu_hct116_34_37 · imaging | -0.083 | +0.039 | -0.031 | +0.067 |
| bintu_k562_28_30 · imaging | +0.105 | +0.203 | +0.191 | +0.124 |
| su_chr2 · imaging | -0.029 | -0.019 | -0.037 | +0.063 |
| su_genome_tx · imaging | +0.124 | +0.097 | +0.126 | +0.081 |
| bintu_hct116_28_30 · hic | +0.003 | -0.044 | -0.032 | -0.027 |
| bintu_hct116_34_37 · hic | -0.017 | -0.065 | -0.062 | -0.042 |
| bintu_k562_28_30 · hic | +0.026 | -0.125 | -0.075 | -0.064 |
| su_chr2 · hic | +0.045 | -0.055 | -0.008 | -0.019 |
| su_genome_tx · hic | -0.031 | +0.017 | -0.011 | +0.104 |

Test (run once): score imaging input **misfit**, hic input **boot_sd**; pass per input type: ρ ≥ 0.30 with the 95 % interval above 0 on every new test set.

| New test set · input | Score | ρ [95 %] | Within the rule |
|---|---|---|---|
| su_genome · imaging | misfit | +0.083 [+0.062, +0.109] | no |
| su_genome · hic | boot_sd | +0.013 [-0.015, +0.050] | no |
| su_genome_amanitin · imaging | misfit | +0.116 [+0.097, +0.138] | no |
| su_genome_amanitin · hic | boot_sd | +0.007 [-0.022, +0.042] | no |

Verdicts: imaging input **fail**; hic input **fail**.
<!-- END generated:gate2e -->

**Reading.**
- **Fail for both input types, on both untouched sets.** With imaging-derived input the misfit score ranks pair
  errors at ρ +0.08 and +0.12 (intervals above 0, far below the 0.30 bar); with Hi-C input the refit spread is
  indistinguishable from 0. Every other candidate, reported alongside, also stays below 0.30 (the best, misfit
  with Hi-C input, +0.13 / +0.17, was not the practice choice).
- Together with Gates 2 and 2c this is the third pre-registered attempt: from the input alone, ChronoCell cannot
  say which of its distances will be wrong. **In the app** no reliability score is shown.
- The run was chunked (one part per set and input); the first job hit the 2-hour limit after its first part when
  this computer was suspended, and the remaining parts ran as separate chunks.

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
Run 2026-10-04T21:00:02+00:00 (2.5 h). All pairs, mean over splits; brackets: 95 % interval (resampling half B, split 0). Full tables: `validation/benchmark/RESULTS_TABLE.md`.

**Input: imaging-derived contacts** — % of ceiling (raw ρ)

| Dataset | genomic_distance_only | no_3d | v3_2_single | v3_3_windowed | v4_whole | pastis_mds | pastis_pm2 |
|---|---|---|---|---|---|---|---|
| bintu_imr90_28_30 | 0.6 [1.4, 1.8] (0.930) | 81.4 [82.5, 84.7] (0.962) | 44.3 [41.1, 43.7] (0.872) | 88.2 [89.1, 91.8] (0.976) | 88.2 [89.1, 91.8] (0.976) | 54.9 [55.2, 57.8] (0.911) | 76.7 [73.9, 76.9] (0.962) |
| bintu_imr90_18_20 | 0.4 [-3.4, 5.8] (0.962) | 43.0 [19.4, 65.5] (0.815) | 37.7 [13.5, 76.6] (0.523) | 55.4 [32.3, 88.1] (0.868) | 55.4 [32.3, 88.1] (0.868) | 51.0 [23.6, 69.1] (0.654) | 67.9 [28.6, 97.8] (0.848) |
| bintu_a549_28_30 | -1.0 [-2.9, -1.8] (0.915) | 87.9 [86.8, 91.0] (0.939) | 54.4 [48.4, 54.6] (0.827) | 91.2 [90.5, 95.0] (0.952) | 91.2 [90.5, 95.0] (0.952) | 61.9 [59.3, 66.1] (0.852) | 73.1 [70.6, 76.6] (0.906) |
| bintu_hct116_34_37_auxin | -0.1 [-0.2, 0.3] (0.973) | 66.5 [63.3, 65.9] (0.944) | 30.0 [28.6, 32.6] (0.752) | 79.7 [75.1, 78.8] (0.968) | 79.7 [75.1, 78.8] (0.968) | 34.2 [29.7, 34.3] (0.780) | 47.4 [37.9, 43.3] (0.862) |
| su_chr21 | -0.2 [-0.2, -0.2] (0.748) | 87.4 [86.9, 87.7] (0.929) | 37.2 [36.0, 37.2] (0.621) | — | 91.7 [91.3, 92.1] (0.954) | 65.6 [65.1, 66.2] (0.815) | — |
| su_chr21_rep | -0.0 [0.2, 0.2] (0.705) | 92.6 [92.2, 92.9] (0.943) | 37.9 [36.9, 37.8] (0.586) | — | 94.6 [94.4, 95.1] (0.959) | 64.4 [63.8, 64.6] (0.799) | — |
| su_genome:chr1 | 1.8 [1.1, 2.1] (0.775) | 92.9 [92.2, 93.6] (0.963) | 36.6 [35.0, 38.2] (0.668) | 94.1 [93.3, 94.9] (0.968) | 94.1 [93.3, 94.9] (0.968) | 57.5 [56.3, 59.1] (0.842) | 83.5 [81.9, 86.1] (0.938) |
| su_genome:chr3 | 2.1 [1.2, 2.6] (0.878) | 91.8 [90.9, 92.9] (0.957) | 38.0 [36.2, 40.3] (0.672) | 93.3 [92.4, 94.4] (0.961) | 93.3 [92.4, 94.4] (0.961) | 53.6 [52.1, 55.9] (0.810) | 82.1 [80.5, 84.1] (0.935) |
| su_genome:chr4 | -2.2 [-3.0, -1.5] (0.926) | 93.8 [92.9, 95.0] (0.976) | 28.8 [27.1, 30.1] (0.637) | 95.6 [94.7, 96.9] (0.980) | 95.6 [94.7, 96.9] (0.980) | 54.1 [52.2, 55.7] (0.827) | 84.0 [82.4, 86.0] (0.954) |
| su_genome:chr5 | 0.1 [-0.3, 0.7] (0.913) | 89.8 [88.5, 91.8] (0.974) | 32.1 [30.1, 34.2] (0.632) | 91.4 [90.1, 93.4] (0.976) | 91.4 [90.1, 93.4] (0.976) | 55.2 [53.1, 57.2] (0.832) | 84.4 [82.7, 85.8] (0.965) |
| su_genome:chr6 | 0.7 [0.4, 1.3] (0.889) | 90.9 [90.0, 91.7] (0.964) | 34.0 [32.4, 36.1] (0.640) | 91.7 [90.8, 92.6] (0.967) | 91.7 [90.8, 92.6] (0.967) | 53.1 [51.8, 54.9] (0.819) | 85.8 [84.4, 86.7] (0.957) |
| su_genome:chr7 | -0.7 [-1.2, -0.3] (0.906) | 94.9 [93.6, 96.2] (0.977) | 35.8 [33.6, 37.8] (0.643) | 95.6 [94.3, 96.9] (0.979) | 95.6 [94.3, 96.9] (0.979) | 54.5 [53.0, 56.3] (0.816) | 84.5 [83.1, 86.3] (0.961) |
| su_genome:chr8 | -0.4 [-1.2, 0.4] (0.928) | 94.2 [92.9, 95.7] (0.983) | 33.8 [32.0, 36.2] (0.635) | 95.1 [93.9, 96.7] (0.984) | 95.1 [93.9, 96.7] (0.984) | 46.2 [44.5, 48.1] (0.817) | 79.8 [77.8, 81.9] (0.959) |
| su_genome:chr9 | -0.1 [-0.7, 0.3] (0.874) | 96.0 [95.5, 96.8] (0.970) | 34.8 [32.9, 37.1] (0.630) | 96.7 [96.2, 97.5] (0.973) | 96.7 [96.2, 97.5] (0.973) | 48.3 [46.7, 50.4] (0.796) | 81.3 [80.2, 83.0] (0.944) |
| su_genome:chr10 | -0.1 [-0.9, 0.6] (0.910) | 92.9 [91.9, 93.8] (0.974) | 28.5 [26.5, 32.3] (0.604) | 93.4 [92.4, 94.4] (0.975) | 93.4 [92.4, 94.4] (0.975) | 55.2 [53.8, 57.5] (0.810) | 84.0 [82.4, 85.7] (0.961) |
| su_genome:chr11 | -2.0 [-2.8, -1.7] (0.882) | 94.7 [93.8, 95.7] (0.976) | 36.3 [34.6, 38.6] (0.667) | 95.3 [94.3, 96.2] (0.978) | 95.3 [94.3, 96.2] (0.978) | 58.8 [57.4, 60.8] (0.822) | 84.5 [83.3, 85.8] (0.955) |
| su_genome:chr12 | -0.4 [-0.9, 0.3] (0.902) | 93.8 [92.8, 94.8] (0.979) | 28.8 [26.2, 30.8] (0.661) | 94.3 [93.3, 95.2] (0.980) | 94.3 [93.3, 95.2] (0.980) | 53.9 [52.2, 56.1] (0.810) | 86.4 [84.6, 88.5] (0.961) |
| su_genome:chr13 | -0.2 [-1.6, 1.0] (0.894) | 89.9 [88.0, 93.9] (0.971) | 20.1 [17.0, 23.1] (0.594) | 89.9 [88.0, 93.9] (0.971) | 89.9 [88.0, 93.9] (0.971) | 58.5 [56.1, 61.9] (0.799) | 86.3 [83.7, 89.2] (0.952) |
| su_genome:chr14 | -1.2 [-3.0, -0.2] (0.837) | 95.5 [94.1, 96.8] (0.964) | 43.2 [41.9, 46.1] (0.581) | 95.6 [94.3, 96.9] (0.964) | 95.6 [94.3, 96.9] (0.964) | 63.3 [62.6, 67.0] (0.770) | 81.4 [80.3, 84.6] (0.917) |
| su_genome:chr15 | 0.2 [-1.2, 1.3] (0.720) | 96.2 [95.1, 97.6] (0.948) | 45.9 [43.6, 48.6] (0.642) | 96.2 [95.1, 97.6] (0.948) | 96.2 [95.1, 97.6] (0.948) | 68.2 [66.0, 70.3] (0.793) | 90.9 [88.8, 92.7] (0.930) |
| su_genome:chr16 | -2.6 [-3.5, -1.9] (0.782) | 95.2 [94.5, 96.6] (0.956) | 43.2 [42.1, 45.2] (0.683) | 95.3 [94.6, 96.7] (0.956) | 95.3 [94.6, 96.7] (0.956) | 67.4 [66.1, 68.5] (0.817) | 85.9 [84.5, 86.9] (0.934) |
| su_genome:chr17 | -1.1 [-1.7, -0.6] (0.644) | 96.2 [95.5, 96.8] (0.952) | 41.3 [39.5, 43.0] (0.600) | 96.2 [95.5, 96.8] (0.952) | 96.2 [95.5, 96.8] (0.952) | 65.6 [64.5, 66.4] (0.761) | 92.1 [91.1, 92.6] (0.931) |
| su_genome:chr18 | -5.9 [-7.1, -5.3] (0.698) | 95.0 [94.4, 95.9] (0.955) | 35.1 [33.7, 37.7] (0.598) | 95.4 [94.7, 96.3] (0.956) | 95.4 [94.7, 96.3] (0.956) | 55.4 [53.9, 56.5] (0.739) | 88.7 [87.2, 89.9] (0.939) |
| su_genome:chr19 | 0.0 [-1.2, 1.0] (0.523) | 96.0 [95.2, 96.5] (0.953) | 52.3 [50.4, 54.0] (0.611) | 96.1 [95.3, 96.7] (0.954) | 96.1 [95.3, 96.7] (0.954) | 74.4 [73.3, 75.8] (0.812) | 91.0 [90.0, 92.1] (0.931) |
| su_genome:chr20 | 1.6 [0.8, 2.5] (0.731) | 95.4 [94.3, 96.9] (0.969) | 45.8 [44.0, 48.0] (0.667) | 95.5 [94.3, 97.0] (0.969) | 95.5 [94.3, 97.0] (0.969) | 66.1 [64.9, 67.5] (0.829) | 86.0 [84.7, 87.9] (0.929) |
| su_genome:chr21 | -6.0 [-6.6, -5.3] (0.619) | 95.0 [94.5, 95.8] (0.967) | 47.7 [46.5, 49.0] (0.621) | 95.1 [94.7, 95.9] (0.967) | 95.1 [94.7, 95.9] (0.967) | 70.3 [69.5, 70.9] (0.793) | 89.6 [88.8, 90.3] (0.913) |
| su_genome:chr22 | 0.8 [-0.5, 1.0] (0.520) | 97.2 [96.6, 97.7] (0.964) | 63.8 [63.1, 65.1] (0.677) | 97.3 [96.7, 97.8] (0.964) | 97.3 [96.7, 97.8] (0.964) | 86.6 [86.0, 87.3] (0.887) | 94.2 [93.8, 95.0] (0.944) |
| su_genome:chrX | 0.5 [-0.4, 1.0] (0.921) | 90.8 [89.3, 93.1] (0.975) | 33.5 [31.5, 35.8] (0.640) | 90.9 [89.4, 93.1] (0.975) | 90.9 [89.4, 93.1] (0.975) | 41.7 [40.7, 45.1] (0.765) | 81.4 [79.6, 84.9] (0.950) |

**Input: sequencing Hi-C (Rao et al. 2014)** — % of ceiling (raw ρ)

| Dataset | no_3d | v3_2_single | v3_3_windowed | v4_whole | pastis_mds | pastis_pm2 |
|---|---|---|---|---|---|---|
| bintu_imr90_28_30 | 80.8 [79.8, 81.2] (0.968) | 67.8 [66.6, 68.4] (0.953) | 92.5 [91.2, 92.8] (0.987) | 92.5 [91.2, 92.8] (0.987) | 74.5 [72.5, 74.4] (0.963) | 85.1 [84.8, 86.8] (0.978) |
| bintu_imr90_18_20 | -18.5 [-54.2, 16.2] (0.839) | -1.6 [-44.5, 22.5] (0.564) | -24.2 [-64.7, 29.1] (0.882) | -24.2 [-64.8, 29.0] (0.882) | 0.5 [-21.2, 46.2] (0.688) | -10.5 [-51.8, 27.8] (0.806) |
| su_chr21 | 70.1 [69.3, 70.1] (0.852) | 56.0 [54.9, 56.1] (0.790) | — | 84.9 [84.0, 85.1] (0.931) | 72.8 [71.2, 72.2] (0.870) | — |
| su_chr21_rep | 66.9 [66.0, 67.0] (0.831) | 52.6 [51.5, 52.9] (0.766) | — | 80.5 [79.5, 80.9] (0.908) | 67.7 [66.3, 67.7] (0.841) | — |
| su_genome:chr1 | 45.4 [44.4, 48.2] (0.801) | 32.0 [30.0, 34.5] (0.789) | 53.1 [51.7, 56.2] (0.874) | 53.1 [51.7, 56.2] (0.874) | 39.3 [37.3, 42.0] (0.793) | 48.9 [47.6, 52.2] (0.854) |
| su_genome:chr3 | 50.3 [47.7, 52.4] (0.805) | 36.4 [33.5, 38.8] (0.804) | 61.9 [59.3, 64.3] (0.905) | 61.9 [59.3, 64.3] (0.905) | 42.6 [39.9, 44.9] (0.775) | 54.2 [50.7, 56.8] (0.845) |
| su_genome:chr4 | 41.0 [39.3, 43.3] (0.866) | 28.5 [26.6, 30.0] (0.875) | 57.2 [55.4, 59.7] (0.942) | 57.2 [55.4, 59.7] (0.942) | 36.0 [34.5, 37.9] (0.866) | 48.4 [46.4, 50.1] (0.929) |
| su_genome:chr5 | 43.3 [41.7, 45.5] (0.892) | 26.8 [24.7, 28.8] (0.854) | 55.7 [53.9, 58.1] (0.934) | 55.7 [53.9, 58.1] (0.934) | 34.5 [32.0, 36.8] (0.875) | 45.4 [43.9, 48.0] (0.938) |
| su_genome:chr6 | 39.0 [37.2, 40.6] (0.811) | 26.8 [24.3, 29.0] (0.780) | 45.8 [44.0, 47.5] (0.871) | 45.8 [44.0, 47.5] (0.871) | 33.3 [31.6, 35.3] (0.784) | 43.4 [41.3, 44.9] (0.867) |
| su_genome:chr7 | 36.5 [33.9, 39.3] (0.718) | 23.7 [20.7, 27.3] (0.750) | 42.7 [40.2, 45.4] (0.807) | 42.7 [40.2, 45.4] (0.807) | 29.6 [26.9, 32.5] (0.679) | 36.2 [33.2, 39.6] (0.761) |
| su_genome:chr8 | 22.6 [20.9, 25.2] (0.755) | 13.5 [11.1, 15.9] (0.815) | 40.7 [38.5, 43.1] (0.911) | 40.7 [38.5, 43.1] (0.911) | 14.3 [12.2, 17.1] (0.716) | 23.0 [21.1, 25.5] (0.777) |
| su_genome:chr9 | 19.0 [16.3, 21.0] (0.576) | 21.8 [19.1, 23.9] (0.736) | 28.5 [25.1, 30.3] (0.711) | 28.4 [25.1, 30.2] (0.711) | 18.4 [15.6, 20.4] (0.548) | 18.7 [15.8, 20.7] (0.565) |
| su_genome:chr10 | 37.0 [33.9, 40.0] (0.678) | 26.5 [24.0, 29.5] (0.805) | 47.1 [44.4, 49.5] (0.865) | 47.1 [44.4, 49.5] (0.865) | 34.2 [31.0, 37.3] (0.662) | 43.9 [41.1, 46.7] (0.736) |
| su_genome:chr11 | 36.7 [35.8, 38.7] (0.707) | 37.2 [35.5, 39.4] (0.766) | 49.2 [47.7, 51.5] (0.857) | 49.2 [47.7, 51.5] (0.857) | 32.9 [31.7, 34.9] (0.681) | 42.5 [40.3, 44.5] (0.770) |
| su_genome:chr12 | 36.3 [33.3, 39.1] (0.768) | 28.3 [25.7, 30.0] (0.802) | 50.5 [47.6, 52.9] (0.883) | 50.5 [47.6, 52.9] (0.884) | 33.1 [30.1, 35.7] (0.738) | 39.2 [36.3, 42.4] (0.793) |
| su_genome:chr13 | 36.9 [30.7, 43.5] (0.833) | 25.1 [18.8, 29.2] (0.788) | 36.2 [29.7, 41.2] (0.854) | 36.2 [29.8, 41.3] (0.854) | 33.2 [28.0, 39.3] (0.799) | 37.5 [31.6, 44.9] (0.843) |
| su_genome:chr14 | 62.0 [59.0, 67.7] (0.834) | 31.1 [28.1, 35.6] (0.734) | 53.2 [48.8, 58.7] (0.815) | 53.3 [48.9, 58.7] (0.816) | 58.1 [55.1, 63.1] (0.790) | 57.8 [54.3, 62.8] (0.816) |
| su_genome:chr15 | 54.5 [51.4, 57.9] (0.747) | 56.5 [53.9, 59.8] (0.781) | 56.3 [52.6, 59.3] (0.767) | 56.3 [52.6, 59.3] (0.766) | 60.1 [56.4, 63.8] (0.768) | 63.9 [60.3, 67.4] (0.805) |
| su_genome:chr16 | 59.8 [58.5, 61.8] (0.761) | 27.8 [26.5, 29.6] (0.608) | 64.4 [62.6, 66.3] (0.755) | 64.4 [62.6, 66.4] (0.755) | 58.5 [56.9, 60.5] (0.721) | 62.3 [60.4, 63.9] (0.747) |
| su_genome:chr17 | 42.7 [41.3, 43.9] (0.644) | 47.9 [46.3, 49.8] (0.710) | 57.8 [55.9, 59.2] (0.737) | 57.8 [55.9, 59.2] (0.737) | 40.8 [39.2, 42.1] (0.604) | 46.3 [44.7, 47.8] (0.650) |
| su_genome:chr18 | 47.6 [45.3, 50.0] (0.708) | 31.6 [29.1, 33.9] (0.619) | 56.1 [54.0, 58.1] (0.785) | 56.1 [54.1, 58.1] (0.785) | 43.0 [40.9, 45.3] (0.666) | 51.0 [48.5, 53.4] (0.747) |
| su_genome:chr19 | 80.5 [79.0, 82.9] (0.847) | 68.2 [66.4, 69.5] (0.745) | 71.1 [69.7, 72.2] (0.754) | 71.1 [69.7, 72.2] (0.754) | 80.1 [78.5, 82.3] (0.819) | 82.4 [80.9, 84.5] (0.838) |
| su_genome:chr20 | 56.8 [55.0, 60.5] (0.651) | 47.4 [44.8, 51.3] (0.693) | 59.8 [57.4, 62.5] (0.713) | 59.7 [57.3, 62.4] (0.713) | 57.1 [55.3, 60.9] (0.652) | 60.0 [58.5, 63.7] (0.704) |
| su_genome:chr21 | 84.3 [82.9, 85.4] (0.919) | 68.4 [67.5, 69.5] (0.812) | 84.0 [82.6, 85.3] (0.896) | 84.0 [82.6, 85.3] (0.896) | 81.4 [79.8, 82.7] (0.897) | 82.8 [81.0, 83.9] (0.906) |
| su_genome:chr22 | 83.0 [81.9, 84.3] (0.887) | 60.7 [59.5, 61.9] (0.706) | 68.5 [67.1, 69.8] (0.767) | 68.7 [67.3, 69.9] (0.768) | 84.3 [83.2, 85.0] (0.876) | 84.6 [83.3, 85.5] (0.881) |
| su_genome:chrX | 0.0 [-2.0, 1.5] (0.753) | 3.4 [-0.1, 6.7] (0.744) | 11.9 [8.1, 14.3] (0.879) | 11.9 [8.1, 14.3] (0.879) | -2.5 [-5.3, 0.3] (0.687) | 0.5 [-2.8, 2.9] (0.756) |

- Not finished: **su_chr21 · pastis_pm2**: started twice, never finished: stopped after 60 min of wall time on the imaging input (split 0), and earlier after about 53 min inside the first, interrupted full run (PASTIS 0.4.0 PM2, single-threaded, max_iter 5000, 651 loci). Not attempted on su_chr21_rep (the same 651 loci). PM2 ran on every other test unit.

'Where we lose' entries: 610 (listed in RESULTS_TABLE.md).
<!-- END generated:gate3 -->

**How it was run.** The single run of the whole plan was stopped at this machine's 2-hour job
limit, while PASTIS PM2 was still fitting Su chr21. That run writes its result file only at the end,
so nothing from it is used. The plan was then run as four parts (Bintu sets; Su chr21; its
replicate; genome-scale), each with the harness unchanged (`--datasets`, `--out`). The parts were
merged by `--merge`, which recomputes every summary from the parts' rows. Where the first run had
finished a unit, the parts reproduce its numbers exactly (fixed seeds). PM2 on the two 651-locus
sets is the one thing not measured (above).

**Reading.**
- **Imaging-derived input: the population models are best, or tied best, on 27 of 28 test units.**
  - The exception is the weak-structure IMR-90 18–20 Mb region: PASTIS PM2 67.9 % against 55.4 %.
    The 95 % intervals overlap almost entirely there (32–88 % and 29–98 %).
  - PASTIS ranges 47–94 % (PM2) and 34–87 % (MDS); the v3.2 single structure 20–64 %.
  - The population model's gain over inverting each frequency on its own (no 3D) is 2–13 points on
    the chr21 regions. On the genome-scale loci (3 Mb apart) it is 0.0–1.8 points: at that spacing
    the 3D model adds almost nothing to the input.
- **Sequencing Hi-C input: best on 18 of 26 units, beaten on 8.**
  - On the IMR-90 18–20 Mb region every method is near or below zero, and ours is lowest (−24 %,
    interval −65 to +29 %).
  - On seven genome-scale chromosomes (chr13, 14, 15, 19, 20, 21, 22) the no-3D inversion or PASTIS
    is higher. The largest gap is chr22: 68.7 % against 84.6 % for PASTIS PM2. The 3D fit loses
    most on the small, gene-dense chromosomes at genome scale.
  - Absolute size (CCC, nm) with Hi-C input is often better with PASTIS PM2 or the no-3D inversion
    (53 of the 610 "where we lose" entries), as expected from Gate 1b.
- **Intervals, and how this relates to Gate 2 / 2b.** Most of the 610 entries are interval coverage.
  - With imaging-derived input the raw 90 % intervals hold 70–85 %, and 83–91 % after the Gate 2
    recalibration.
  - With Hi-C input they hold 23–51 % raw and 36–70 % with that recalibration, on every Hi-C unit,
    the genome-scale ones included.
  - This is the same failure Gate 2b measured with its own pre-registered test. The benchmark adds
    that it holds at genome scale too.
- **Intervals in brackets** come from split 0 on a 20,000-pair subsample, while the point value is
  the mean over splits. On the Bintu sets the bracket can sit beside the point value; read them as
  the size of the truth noise, not as a range around the mean.

## Gate 3b — a learned correction on top of the population model (Phase A4)

**Question.** Can a correction learned on practice data (pair features: separation, model and input
residuals, evidence, depth, locus spacing) make the model best on more Hi-C benchmark units?

**Practice and decision.** A ridge and a small network (trained on the GPU) were compared with no correction
by leave-one-dataset-out. Both raised Lin's CCC (sizes) but lowered the held-out pattern, so the selection rule
fixed before the run (best mean held-out pattern, and better than no correction on every held-out group) chose
no correction. With nothing to test, the test was not run (`frozen.LEARNED_CORRECTION`).

<!-- BEGIN generated:gate3b -->
Practice (leave-one-dataset-out; mean held-out trend-removed ρ / CCC):

| Candidate | Pattern ρ | CCC |
|---|---|---|
| none (chosen) | +0.7093 | 0.423 |
| linear | +0.6540 | 0.678 |
| mlp | +0.6821 | 0.685 |

Test: **not run** — practice chose no correction (pre-specified selection rule). Pre-registered rule (frozen.LEARNED_CORRECTION): best on ≥ 90 % of Hi-C units and within 1 point of the current model on every imaging unit.
<!-- END generated:gate3b -->

**Reading.** A learned correction helps absolute sizes on practice (CCC 0.42 → 0.68) but not the pattern of
distances, which is what Gate 3 scores. The benchmark method `learned_correction` exists and returns the
base model with a note while no correction is frozen; the app is unchanged.

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
- **A second, independent SV test was looked for (October 2026) and none could be pre-registered
  and run here.** A useful test needs a clonal rearrangement with known joins, Hi-C of the same cell
  type with and without it, and contact maps (not raw reads) that can be read on this machine:
  - *Firre deletion* (Barutcu et al., *Nat Commun* 9:1444, 2018; GEO GSE98632): wild-type and
    knockout mouse fibroblasts, an 82 kb deletion (mm9 chrX:47,908,463–47,990,294), hemizygous in
    male cells. The design is right, but GEO holds raw reads only; building the maps means aligning
    hundreds of millions of reads, with no aligner or Hi-C pipeline on this machine. The deletion is
    also about two bins at the published 40 kb resolution, so the expected change is small.
  - *CTCF-motif edits in HAP1 cells* (Sanborn et al., *PNAS* 2015): tens of base pairs, far below
    the simulator's 5–25 kb beads; they test loop extrusion, not a bead-level rearrangement.
  - *Cancer lines with breakpoint-resolved SVs* (Dixon et al., *Nat Genet* 50:1388, 2018): mostly
    complex rearrangements with copy-number changes and no Hi-C of a matched normal of the same cell
    type, the same confound that weakens the K562-vs-GM12878 test above.
  - So the label stays: **mechanism simulator, not validated.**
- **Drug-lab mechanisms** (loop-extrusion, compaction and similar what-ifs) have no matching
  perturbation data here and stay labelled "not validated".

## Gate 4c — cohesin loss against sequencing Hi-C on held-out regions (Phase B1)

**Question.** Gate 4 checked the cohesin-loss prediction on one imaged region. Does the same model, with its
parameters unchanged, predict what RAD21 degradation does to Hi-C on regions it has never seen?

**Test** (`python validation/cohesin_hic.py --test`; rule in `frozen.COHESIN_HIC`, committed after the practice
run on the Gate 4 practice region and before any other auxin map was read). HCT116 RAD21-mAC cells, untreated
and 6 h auxin (Rao et al. 2017, ENCODE in situ Hi-C; intact Hi-C as a secondary pair). From the untreated map
only: v3.3 population of each 2 Mb window → the Gate 4 cohesin-loss map → predicted change of contact
probability. Score: Spearman of predicted and measured change (library-normalised auxin over untreated);
baseline: the same model without the domain term (trend only). Pass: at least 3 of 6 regions beat trend only,
with the 95 % interval of the difference above 0.

<!-- BEGIN generated:gate4c -->
Practice (the Gate 4 practice region; nothing fitted):

| Pair | Region | Model ρ | Trend-only ρ | Difference [95 %] |
|---|---|---|---|---|
| in situ | chr21:28-30 Mb | +0.219 | +0.227 | -0.009 [-0.115, +0.104] |
| intact | chr21:28-30 Mb | +0.410 | +0.053 | +0.357 [+0.250, +0.459] |

Test (run once): pass if ≥ 3 of 6 held-out regions beat trend only on the main pair (in situ), difference > 0 with its 95 % interval above 0.

| Pair | Region | Pairs | Model ρ | Trend-only ρ | Difference [95 %] | Beats trend only |
|---|---|---|---|---|---|---|
| in situ (main) | chr2:216-218 Mb | 17,394 | +0.486 | +0.092 | +0.394 [+0.325, +0.466] | yes |
| in situ (main) | chr5:140-142 Mb | 14,180 | +0.212 | +0.175 | +0.037 [-0.045, +0.132] | no |
| in situ (main) | chr7:130-132 Mb | 13,779 | +0.496 | +0.096 | +0.400 [+0.345, +0.444] | yes |
| in situ (main) | chr10:100-102 Mb | 18,378 | +0.416 | +0.255 | +0.161 [+0.082, +0.243] | yes |
| in situ (main) | chr12:52-54 Mb | 14,374 | +0.462 | +0.048 | +0.414 [+0.317, +0.511] | yes |
| in situ (main) | chr17:48-50 Mb | 18,424 | +0.416 | +0.301 | +0.115 [+0.047, +0.177] | yes |
| intact (secondary) | chr2:216-218 Mb | 9,175 | +0.601 | +0.063 | +0.537 [+0.466, +0.603] | yes |
| intact (secondary) | chr5:140-142 Mb | 6,583 | +0.494 | -0.043 | +0.537 [+0.455, +0.627] | yes |
| intact (secondary) | chr7:130-132 Mb | 6,962 | +0.708 | +0.106 | +0.602 [+0.552, +0.662] | yes |
| intact (secondary) | chr10:100-102 Mb | 9,305 | +0.504 | -0.119 | +0.624 [+0.556, +0.692] | yes |
| intact (secondary) | chr12:52-54 Mb | 6,199 | +0.548 | -0.138 | +0.686 [+0.625, +0.744] | yes |
| intact (secondary) | chr17:48-50 Mb | 8,842 | +0.503 | -0.158 | +0.660 [+0.580, +0.744] | yes |

5 of 6 main-pair regions beat trend only. Verdict: **pass**.
<!-- END generated:gate4c -->

**Reading.**
- **Pass: 5 of 6 held-out regions** (in situ, the main pair); the sixth (chr5:140–142 Mb, the protocadherin
  cluster) is positive but its interval includes 0. On the intact secondary pair every region beats trend only,
  by more.
- **What it shows.** Removing the domain pattern (the λ term fitted on imaging data in Gate 4) adds
  +0.04 to +0.62 of rank agreement with the measured Hi-C change over the separation shift alone: the model's
  cohesin-loss mechanism transfers from imaging to sequencing data and to new regions.
- **What it does not show.** One cell line and one perturbation. CTCF and WAPL loss (Nora 2017, Haarhuis
  2017) are different mechanisms the model does not describe and were not tested; Hafner 2023 has imaging,
  not Hi-C, after auxin.
- The run was interrupted once by a dropped network read in the secondary pair, after all six main-pair regions
  were computed; the completed run (resumable per region, retries on dropped reads) reproduced the printed
  values exactly.
- **In the app** the cohesin simulator keeps its label until Gate 4d passes as well (see below).

## Gate 4d — structural variants with Hi-C before and after (Phase B1): blocked

**Pre-registered rule** (`frozen.SV_GATE4D`): on at least 3 engineered or patient events with Hi-C of the same
cell type before and after, the engine's predicted change beats both the distance-shift and the no-change
baselines, with the 95 % interval of the difference above 0.

**Not run: not enough events.** Openly downloadable, checksum-published Hi-C before and after an SV with
stated breakpoints exists here for two events only: the DXZ4 macrosatellite deletion on the active and on the
inactive X of RPE-1 cells (ENCODE ENCSR624AJS and ENCSR566XIP, chrX:114,946,736–115,094,976 hg19). The 4DN
Treg CRISPR deletions (Ikzf2 enhancer, Foxp3) have matched controls but no breakpoints in their metadata, and
GEO maps publish no checksums. The variant engine therefore stays a **mechanism simulator, not validated**
everywhere in the app and in reports.

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

## Gate 5b — prediction without contact data, with more inputs (Phase A5)

**Question.** Does adding accessibility (ATAC-seq), active chromatin (H3K27ac) and cohesin (RAD21)
peaks of the same cell line lift the prediction without contact data towards half of the
reproducible pattern?

**Test** (`python validation/predictor_v2.py --test`; rule in `frozen.PREDICTOR_V2`, committed with
the frozen model before the run). Training: the three untreated non-IMR-90 practice regions (K562,
HCT116 ×2), so every test set is a held-out cell type. Feature set and ridge penalty by
leave-one-dataset-out on those three. Scoring as Gate 5. Pass: on at least 3 of the 5 test sets,
≥ 50 % of the pattern with its 95 % interval above 0 and above Gate 5's model.

<!-- BEGIN generated:gate5b -->
Practice (leave-one-dataset-out over K562 / HCT116 regions; best ridge per feature set):

| Feature set | Ridge | Held-out trend-removed ρ (mean) |
|---|---|---|
| ctcf | 0.001 | +0.294 |
| ctcf+atac | 0.001 | +0.256 |
| ctcf+h3k27ac | 0.001 | +0.270 |
| ctcf+rad21 (chosen) | 0.01 | +0.432 |
| ctcf+all | 0.01 | +0.414 |

Test (run once): pass on ≥ 3 of 5 sets: ≥ 50 % of the pattern, interval above 0, and above Gate 5's model.

| Test set | Predictor v2: % of ceiling [95 %] | Gate 5 model | Pass |
|---|---|---|---|
| bintu_imr90_28_30 | -10.5 [-12.6, -9.3] | 26.7 | no |
| bintu_imr90_18_20 | 0.2 [-13.8, 13.7] | 0.5 | no |
| bintu_a549_28_30 | 4.9 [0.8, 7.7] | 12.2 | no |
| su_chr21 | 9.4 [7.6, 9.0] | 20.2 | no |
| su_chr21_rep | 10.8 [9.2, 10.6] | 22.6 | no |
| bintu_hct116_34_37_auxin (control) | 45.1 | — | — |

0 of 5 test sets pass. Verdict: **fail**. Not run: Akita (Fudenberg et al., Nat Methods 2020); C.Origami (Tan et al., Nat Biotechnol 2023); Orca (Zhou, Nat Genet 2022) (reasons in validation/predictor_v2.py).
<!-- END generated:gate5b -->

**Reading.**
- **Fail, on every test set, and worse than Gate 5's simpler model on each.**
  - On practice, cohesin peaks helped the held-out region (leave-one-out ρ +0.43 against +0.29 for
    CTCF and GC alone). On new cell types they do not transfer: IMR-90 chr21:28–30 Mb even turns
    negative.
  - Gate 5's model, which also trained on IMR-90 (chr2), stays the better prior for IMR-90.
  - ATAC and H3K27ac did not help even on practice.
- **The control is not a held-out result.** The cohesin-depleted set is HCT116 chr21:34–37 Mb, the
  region the model trained on (untreated), so its 45 % mostly measures the training region again.
- **Engines not run:** Akita, C.Origami and Orca, for the reasons listed under the table. The frozen
  Gate 5 model remains the app's prediction; nothing from Gate 5b reaches the app.

## Gate 5m — the predictor on mouse (pre-registered, not run)

**Why.** The predictor was trained and tested on human data. The app offers it for hg38 only.

**Pre-registered test** (`validation/frozen.py`, `PREDICTOR_MOUSE`, committed before any mouse data
were read):
- **Truth:** ORCA chromatin tracing in mouse ES cells, two 3-Mb loci on chr6 and chr3 (Hafner et al.,
  *Mol Cell* 83:1377, 2023; 4DN sets 4DNESD28H8O7 and 4DNESWDXDZSE, untreated).
- **Inputs:** ENCODE CTCF IDR peaks of mouse ES-Bruce4 cells (ENCFF533APC, mm10; the only ENCODE mouse
  ES CTCF ChIP-seq, and a different ES line from the imaged cells) and the mm10 sequence.
- **Model:** the frozen human model, unchanged.
- **Rule:** Gate 5's criteria (i) and (ii), on both loci. The CTCF- and cohesin-degron lines are
  secondary sets and, with auxin, controls.

**Not run.** The 4DN Data Portal answers downloads with HTTP 403 unless the request carries a 4DN
account access key, which this machine does not have. An openly licensed alternative (Takei et al.,
*Nature* 590:344, 2021, DNA seqFISH+, Zenodo 3735329) gives spot positions, not traces. Turning spots
into single-copy traces is a separate analysis that would itself need validating. Until the test runs
and passes, mouse assemblies get no prediction, in the app or on the command line.

## Gate 5m — result (the predictor on mouse, run October 2026)

The section above was written before the run and is kept as it was. The 4DN files turned out to be public on
the 4DN AWS Open Data bucket (each file's `open_data_url`, with the MD5 the portal publishes), so the
pre-registered test ran unchanged, without an access key (`python validation/predictor_mouse.py --test`).

<!-- BEGIN generated:gate5m -->
Test (run once; the human model unchanged): pass if both test loci have % of ceiling > 0 with its 95 % interval above 0 and raw Spearman above the training-trend baseline.

| Set | Role | Traces | Loci | Predictor: % of ceiling (3 splits) | Trend baseline | Raw ρ: predictor / trend |
|---|---|---|---|---|---|---|
| 4DNESD28H8O7 (chr6 locus, untreated) | test | 877 | 76 | 17.0 [13.7, 24.1] | 1.2 | 0.869 / 0.866 |
| 4DNESWDXDZSE (chr3 locus, untreated) | test | 2,270 | 96 | 7.2 [5.4, 9.5] | 0.4 | 0.918 / 0.916 |
| 4DNESJ3TXVIR (chr6 locus, CTCF-AID untreated) | secondary | 5,171 | 76 | 17.6 | 2.4 | 0.903 / 0.899 |
| 4DNESQ49IXDU (chr6 locus, RAD21-AID untreated) | secondary | 4,233 | 76 | 13.7 | 2.7 | 0.890 / 0.886 |
| 4DNESTNG39BO (chr3 locus, CTCF-AID untreated) | secondary | 5,203 | 96 | 8.0 | 0.7 | 0.924 / 0.922 |
| 4DNESLRTOSQT (chr3 locus, RAD21-AID untreated) | secondary | 3,897 | 96 | 8.3 | 0.4 | 0.909 / 0.908 |
| 4DNES2KX6HQ5 (chr6 locus, CTCF-AID + auxin) | control | 6,088 | 76 | 15.3 | 1.8 | 0.917 / 0.913 |
| 4DNESMN7RCSB (chr6 locus, RAD21-AID + auxin) | control | 4,658 | 76 | -0.1 | 1.2 | 0.888 / 0.889 |
| 4DNESG62SAVA (chr3 locus, CTCF-AID + auxin) | control | 7,600 | 96 | 0.3 | 0.7 | 0.945 / 0.945 |
| 4DNESBH54BG2 (chr3 locus, RAD21-AID + auxin) | control | 7,735 | 96 | 4.4 | -0.1 | 0.937 / 0.938 |

Verdict: **pass** (95 % interval from split 0).
<!-- END generated:gate5m -->

**Reading.**
- **Pass, modestly.** The human model, unchanged, recovers 17 % (chr6 locus) and 7 % (chr3 locus) of the
  reproducible pattern in mouse ES cells, intervals above 0, and its raw Spearman is above the separation
  trend's on both, by 0.002–0.003: as in human (Gate 5), a prior, not a substitute for contacts.
- **Controls behave as expected where the signal is clear:** with cohesin (RAD21) degraded the chr6 locus
  drops to −0.1 %, and with CTCF degraded the chr3 locus to 0.3 %; the CTCF-degraded chr6 locus keeps 15 %,
  so the control is not uniform.
- The run crashed once on a secondary set whose tables have 8 columns (no Cell_ID) before writing a result;
  after the loader fix the test loci load identically (877 and 2,270 traces) and the numbers are those printed
  by the first run.
- **In the app** the predictor is now offered for mouse assemblies, labelled with these numbers; mouse CTCF
  peaks must be uploaded (the ENCODE list is human only).

## Gate 6 — loop calls against reference calls on held-out cell lines (Phase B3)

**Question.** Does ChronoCell's loop caller agree with the field's reference calls at least as well as the
comparable tools that run here?

**Test** (`python validation/loops_gate6.py --test`; rule in `frozen.LOOPS_GATE6`, committed after the practice
run on GM12878 and before any test window was read). Reference: the HiCCUPS loops ENCODE called on the same
maps. Held-out cell lines K562 and IMR-90, three 10 Mb windows each at 10 kb. Comparators: chromosight and
Mustache on the same windows. Pass: ChronoCell's F1 at least the best comparator's on both cell lines.

<!-- BEGIN generated:gate6 -->
Practice (GM12878, three 10 Mb windows; ChronoCell settings):

| Setting | Precision | Recall | F1 |
|---|---|---|---|
| fdr0.05_none | 0.417 | 0.750 | 0.536 |
| fdr0.05_kr | 0.900 | 0.450 | 0.600 |
| fdr0.1_none | 0.400 | 0.767 | 0.526 |
| fdr0.1_kr | 0.906 | 0.483 | 0.630 |
| fdr0.2_none | 0.331 | 0.833 | 0.474 |
| fdr0.2_kr (chosen) | 0.889 | 0.533 | 0.667 |

Test (run once): pass if ChronoCell's F1 is at least the best of chromosight and Mustache on both held-out cell lines (reference: ENCODE HiCCUPS calls on the same maps).

| Cell line | Method | Calls | Reference loops | Matched | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|
| k562 | chronocell | 315 | 156 | 94 | 0.298 | 0.603 | 0.399 |
| k562 | chromosight | 316 | 156 | 92 | 0.291 | 0.590 | 0.390 |
| k562 | mustache | 231 | 156 | 94 | 0.407 | 0.603 | 0.486 |
| imr90 | chronocell | 134 | 154 | 106 | 0.791 | 0.688 | 0.736 |
| imr90 | chromosight | 298 | 154 | 95 | 0.319 | 0.617 | 0.420 |
| imr90 | mustache | 255 | 154 | 97 | 0.380 | 0.630 | 0.474 |

Verdict: **fail**.
<!-- END generated:gate6 -->

**Reading.**
- **Fail by the rule (1 of 2 cell lines).** On IMR-90 ChronoCell's calls agree best with HiCCUPS (F1 0.74, precision
  0.79; chromosight 0.42, Mustache 0.47). On K562, Mustache is ahead (F1 0.49 against 0.40; chromosight 0.39):
  ChronoCell finds as many reference loops (recall 0.60) but makes more calls that HiCCUPS did not.
- **What the comparison measures.** Agreement with HiCCUPS on the same map; ChronoCell's caller follows the
  HiCCUPS recipe, so agreement is expected to favour it, and still it does not win on K562. The practice-chosen
  FDR (0.2) trades precision for recall; K562's karyotype (amplifications) may also inflate calls.
- **In the app** the analysis suite is labelled with this result; loop calls are offered as an analysis, not as a
  validated detector.

## Gate 7 — false-discovery control of the differential analysis (Phase B2)

**Question.** On real replicate maps with no biological difference, plus changes planted at known pixels,
does the differential analysis keep its false-discovery rate at the stated level, and how much does it find?

**Test** (`python validation/diff_gate7.py --test`; rule in `frozen.DIFF_GATE7`, committed after the practice
run and before any test region was read). Four independent untreated intact Hi-C experiments of the same
cells, two per condition; 100 pixels per region raised ×2 or ×4 in condition B (3 draws each); five 2 Mb test
regions. Pass: mean false-discovery proportion ≤ 0.05 at nominal 0.05.

<!-- BEGIN generated:gate7 -->
Practice (chr21:28-30 Mb):

| Setting | Mean FDP | Recall ×2 | Recall ×4 | No-change discoveries |
|---|---|---|---|---|
| min5_dist | 0.036 | 0.16 | 0.94 | 0 |
| min10_dist | 0.024 | 0.11 | 0.92 | 0 |
| min5_nodist | 0.054 | 0.09 | 0.88 | 0 |

Test (run once): pass if the mean false-discovery proportion over 30 spike-in runs is ≤ 0.05.

| Region | No-change discoveries | Spike-in FDP (×2 / ×4) | Recall ×2 / ×4 |
|---|---|---|---|
| chr3:180-182 Mb | 0 | 0.000 / 0.000 | 0.00 / 0.18 |
| chr8:126-128 Mb | 0 | 0.000 / 0.000 | 0.00 / 0.70 |
| chr11:65-67 Mb | 0 | 0.000 / 0.004 | 0.03 / 0.74 |
| chr14:90-92 Mb | 0 | 0.000 / 0.000 | 0.00 / 0.59 |
| chr19:12-14 Mb | 0 | 0.000 / 0.014 | 0.03 / 0.92 |

Mean FDP 0.002 [0.000, 0.004] at nominal 0.05; recall ×2 0.01, ×4 0.63. Not run: diffHic (R / Bioconductor not installed on this machine); multiHiCcompare (R / Bioconductor not installed on this machine); CHESS (pip install chess-hic fails: its dependency pysam (via FAN-C) does not build on Windows). Verdict: **pass**.
<!-- END generated:gate7 -->

**Reading.**
- **Pass: false discoveries are controlled.** Over 30 spike-in runs in five held-out regions the mean
  false-discovery proportion is 0.002 at a nominal 0.05, and no comparison without a planted change produced a
  discovery.
- **But the test is conservative and has little power with two replicates per side.** It finds 63 % of four-fold
  changes and about 1 % of two-fold changes. Read a short list of significant pixels as "strong changes only", not
  as the absence of smaller ones.
- On real data in the app (untreated vs 6 h auxin, two experiments each, chr21:28–30 Mb) it tested 2,570 pixels and
  called 4 (1 gained, 3 lost after auxin).
- diffHic, multiHiCcompare and CHESS could not run here (recorded in the result file).

## Benchmarks against other tools (Phase B9)

Every tool that runs here gets the same input and the same reference as ChronoCell (Gates 6 and 7); the others
are listed with the reason they could not run (`python validation/tools_b9.py`).

<!-- BEGIN generated:b9tools -->
| Task | Tool | Status | Where the comparison is |
|---|---|---|---|
| loops | HiCCUPS | not run: needs Java (Juicer tools); no Java runtime on this machine; ENCODE's HiCCUPS calls are Gate 6's reference | Gate 6 (validation/results_gate6.json) |
| loops | chromosight | run (1.6.3, isolated environment) | Gate 6 (validation/results_gate6.json) |
| loops | Mustache | run (1.3.3, isolated environment with NumPy 1.26; reads .mcool because hic-straw does not build here) | Gate 6 (validation/results_gate6.json) |
| TADs | TopDom | not run: R / Bioconductor package; R is not installed on this machine | ChronoCell's TopDom-like caller (analysis suite); no pre-registered boundary gate |
| TADs | Arrowhead | not run: needs Java (Juicer tools); no Java runtime on this machine | ChronoCell's Arrowhead-like caller; no gate |
| TADs | insulation (cooltools) | not run: cooltools has no Windows build (needs a C compiler) | ChronoCell's insulation score |
| compartments | dcHiC | not run: R / Bioconductor package; R is not installed on this machine | — |
| differential | diffHic | not run: R / Bioconductor package; R is not installed on this machine | Gate 7 (validation/results_gate7.json) |
| differential | multiHiCcompare | not run: R / Bioconductor package; R is not installed on this machine | Gate 7 |
| differential | CHESS | not run: pip install chess-hic fails building pysam (a FAN-C dependency) on Windows | Gate 7 |
| SV detection | HiNT | not run: not installable here (needs R and BWA / samtools) | — |
| SV detection | hic_breakfinder | not run: C++ (Eigen, BamTools) with no Windows binary; no compiler here | — |
| SV detection | EagleC | not run: no held-out SV truth set with Hi-C here (Gate 4d is blocked), so there is nothing to score it on | — |
| SV impact | Akita (basenji) | not run: TensorFlow model code not packaged for pip; no SV truth set (Gate 4d) | — |
| SV impact | Orca | not run: needs selene-sdk (no Windows build) and large weights; no SV truth set (Gate 4d) | — |
<!-- END generated:b9tools -->

## Gate Q — quantum lab: simulated quantum algorithms on real data (pre-registered)

The quantum lab (`chronocell/quantum`, the app's 07 Quantum lab and a Quantum section in each workspace) writes
ChronoCell problems in the forms a quantum computer takes and solves them on a statevector **simulator on this
computer** (not quantum hardware), next to the classical answer. Gate Q asks, on held-out real data, whether each
quantum route does what it claims; the rules are in `validation/frozen.py` (`QUANTUM_GATEQ`), committed before the
test, and the runner is `validation/quantum_gateq.py`.

- **Q1 (solver).** The domain-boundary QUBO of a 20-bin window (19 qubits) on held-out K562 and IMR-90 Hi-C: does
  QAOA's best shot reach the minimum energy (checked against all 2^19 states)?
- **Q2 (biology).** Do QAOA's domain boundaries agree with ENCODE's Arrowhead calls on the same maps at least as well
  as the app's classical callers (insulation, TopDom-like), within the pre-registered margin? Agreement with
  Arrowhead, the field's standard caller, is not biological truth.
- **Q3 (chemistry).** Does VQE reach the exact energy of H2 and HeH+ within chemical accuracy at every listed bond
  length?
- **Q4 (machine learning).** Does a quantum-kernel SVM predict which genes are expressed, from Hi-C features, as well
  as a classical RBF-kernel SVM, when trained on GM12878 and tested on IMR-90?

No speed-up is claimed: the sizes a simulator can hold are solved faster classically (for the banded domain QUBO a
dynamic programme finds the exact optimum in linear time), and the noise columns show what today's hardware errors
would do. Practice (settings; `validation/TUNING.md` §19):

<!-- BEGIN generated:gateq_practice -->
Domain QUBO, practice GM12878 (33 windows of 20 bins; 360 settings, exact optimum of each): chosen 40 kb bins, minimum domain 3 bins, gamma 1.5, boundary cost 1.5, difference weights: F1 0.400 (precision 0.542, recall 0.317). Classical callers on the same windows, best of their grids:

| Caller | Setting | Precision | Recall | F1 |
|---|---|---|---|---|
| Domain QUBO (exact optimum) | chosen | 0.542 | 0.317 | 0.400 |
| insulation | [4, 0.25] | 0.378 | 0.415 | 0.395 |
| topdom | [3] | 0.259 | 0.366 | 0.303 |

QAOA settings, practice (same windows; 4,096 shots; noise column: approximate hardware-noise model):

| Depth p | Objective | Optimum found | Mean P(optimum) | Uniform P(optimum) | With noise | Simulated annealing | SQA | Mean QAOA s | Mean CX |
|---|---|---|---|---|---|---|---|---|---|
| 3 | cvar | 100 % | 0.0764 | 1.9e-06 | 42 % | 100 % | 88 % | 1.9 | 473 |
| 3 | expectation | 100 % | 0.1185 | 1.9e-06 | 48 % | 100 % | 88 % | 2.5 | 473 |
| 6 | cvar | 100 % | 0.1094 | 1.9e-06 | 9 % | 100 % | 88 % | 7.0 | 947 |
| 6 | expectation | 100 % | 0.2322 | 1.9e-06 | 9 % | 100 % | 88 % | 8.2 | 947 |

Chosen: p = 3, expectation.

Gene classifier, practice: 106 GM12878 genes (57 % expressed), 5-fold cross-validated AUC: quantum kernel 0.806 ({'bandwidth': 0.05, 'reps': 1, 'C': 100.0}), RBF-SVM 0.796 ({'gamma': 0.02, 'C': 100.0}), logistic regression 0.762.

Chemistry anchors (Szabo & Ostlund): H2 R=1.4 bohr E_HF (Szabo & Ostlund: -1.1167): -1.11671; H2 R=1.4 bohr E_FCI (Szabo & Ostlund: -1.1373): -1.13728; HeH+ R=1.4632 bohr E_HF (Szabo & Ostlund: -2.86066): -2.86066.
<!-- END generated:gateq_practice -->

Test (run once):

<!-- BEGIN generated:gateq -->
Q1 (run once, held-out K562 and IMR-90 windows; 19 qubits each):

| Cell line | Windows | QAOA found the optimum | Mean P(optimum) | Random guessing found it | With noise | Simulated annealing | SQA | Mean QAOA s | Mean exact DP s |
|---|---|---|---|---|---|---|---|---|---|
| k562 | 33 | 100 % | 0.0768 | 3 % | 48 % | 100 % | 73 % | 2.9 | 0.0011 |
| imr90 | 33 | 97 % | 0.0520 | 0 % | 39 % | 100 % | 67 % | 2.6 | 0.0011 |

Pooled: 98.5 % of 66 windows (needed ≥ 90 %). Q1: **pass**.

Q2 (reference: ENCODE Arrowhead domains on the same maps; a call matches within one bin; pass if QAOA's F1 ≥ the better classical caller's − 0.05 on both cell lines):

| Cell line | Method | Precision | Recall | F1 |
|---|---|---|---|---|
| k562 | QAOA (simulated quantum) | 0.758 | 0.141 | 0.238 |
| k562 | exact optimum of the QUBO | 0.758 | 0.141 | 0.238 |
| k562 | insulation (classical) | 0.732 | 0.232 | 0.352 |
| k562 | topdom (classical) | 0.667 | 0.271 | 0.386 |
| imr90 | QAOA (simulated quantum) | 0.738 | 0.166 | 0.271 |
| imr90 | exact optimum of the QUBO | 0.744 | 0.171 | 0.278 |
| imr90 | insulation (classical) | 0.726 | 0.241 | 0.361 |
| imr90 | topdom (classical) | 0.693 | 0.278 | 0.397 |

Q2: **fail** (k562 fail, imr90 fail).

Q3 (UCCSD-VQE vs exact diagonalisation; noise column: hardware-efficient ansatz under the noise model):

| Molecule | Bond length (Å) | Ansatz | Hartree–Fock | FCI | VQE | Error (mHa) | With noise |
|---|---|---|---|---|---|---|---|
| H2 | 0.500 | UCCSD | -1.04300 | -1.05516 | -1.05516 | +1.4e-11 | -0.2798 |
| H2 | 0.741 | UCCSD | -1.11668 | -1.13727 | -1.13727 | +1.6e-11 | -0.5762 |
| H2 | 1.000 | UCCSD | -1.06611 | -1.10115 | -1.10115 | +1.3e-11 | -0.6832 |
| H2 | 1.500 | UCCSD | -0.91087 | -0.99815 | -0.99815 | +1.4e-11 | -0.7245 |
| H2 | 2.000 | UCCSD | -0.78379 | -0.94864 | -0.94864 | +1.2e-11 | -0.7246 |
| H2 | 2.500 | UCCSD | -0.70294 | -0.93605 | -0.93605 | +1.2e-11 | -0.7240 |
| H2 | 0.500 | hardware-efficient (2 layers) | -1.04300 | -1.05516 | -1.05516 | +2.7e-12 | -0.9524 |
| H2 | 0.741 | hardware-efficient (2 layers) | -1.11668 | -1.13727 | -1.13727 | +1.8e-12 | -1.0629 |
| H2 | 1.000 | hardware-efficient (2 layers) | -1.06611 | -1.10115 | -1.10115 | +2.1e-11 | -1.0458 |
| H2 | 1.500 | hardware-efficient (2 layers) | -0.91087 | -0.99815 | -0.99815 | +7.1e-12 | -0.9619 |
| H2 | 2.000 | hardware-efficient (2 layers) | -0.78379 | -0.94864 | -0.92454 | +24 | -0.8966 |
| H2 | 2.500 | hardware-efficient (2 layers) | -0.70294 | -0.93605 | -0.93164 | +4.4 | -0.9038 |
| HeH+ | 0.500 | UCCSD | -2.74612 | -2.76025 | -2.76025 | +4.1e-11 | -1.7161 |
| HeH+ | 0.774 | UCCSD | -2.86066 | -2.88071 | -2.88071 | +4e-11 | -2.1460 |
| HeH+ | 1.000 | UCCSD | -2.81376 | -2.83466 | -2.83466 | +3.7e-11 | -2.2457 |
| HeH+ | 1.500 | UCCSD | -2.70084 | -2.71027 | -2.71027 | +5.1e-11 | -2.2694 |
| HeH+ | 2.000 | UCCSD | -2.65539 | -2.65639 | -2.65639 | +3.8e-11 | -2.2574 |
| HeH+ | 0.500 | hardware-efficient (2 layers) | -2.74612 | -2.76025 | -2.76014 | +0.11 | -2.6218 |
| HeH+ | 0.774 | hardware-efficient (2 layers) | -2.86066 | -2.88071 | -2.92269 | -42 | -2.8223 |
| HeH+ | 1.000 | hardware-efficient (2 layers) | -2.81376 | -2.83466 | -3.03982 | -2.1e+02 | -2.9471 |
| HeH+ | 1.500 | hardware-efficient (2 layers) | -2.70084 | -2.71027 | -3.10132 | -3.9e+02 | -3.0149 |
| HeH+ | 2.000 | hardware-efficient (2 layers) | -2.65539 | -2.65639 | -3.10940 | -4.5e+02 | -3.0241 |

Q3: **pass** (worst UCCSD error 5.1e-11 mHa).

Q4 (trained on 106 GM12878 genes, tested on 380 IMR-90 genes, 70 % expressed; 95 % intervals from 1,000 resamples of genes):

| Classifier | Test AUC |
|---|---|
| Quantum-kernel SVM (simulated, exact kernel) | 0.623 [0.561, 0.683] |
| Quantum-kernel SVM, kernel from 1000 shots per entry | 0.628 |
| RBF-kernel SVM (classical) | 0.657 [0.596, 0.717] |
| Logistic regression (classical) | 0.749 [0.690, 0.803] |

Quantum minus RBF: -0.034 [-0.075, 0.001]. Q4: **fail** (needed: within 0.03 of the RBF-SVM and the lower 95 % bound above 0.5).
<!-- END generated:gateq -->

**Reading (written after the test run).**
- **Q1 passed.** QAOA on the simulator found the optimum of real-data domain puzzles (19 qubits) in all but one
  held-out window, where a random guess almost never does: the circuit concentrates probability on good answers, as
  the algorithm intends. It is no advantage: the same windows are solved exactly by a classical dynamic programme in
  about a millisecond, simulated annealing also found every optimum, and under the noise model most shots of a circuit
  this deep are noise (the "with noise" column).
- **Q2 failed** on both cell lines. QAOA's best shot was the QUBO's exact optimum in every window but one, so this is
  the formulation, not the quantum solver: its boundaries are precise but too few, the practice tie with the insulation
  caller did not carry over, and the TopDom-like caller agreed best with Arrowhead on both cell lines. The app keeps
  the quantum domain walls as a demonstration, labelled with this result.
- **Q3 passed.** UCCSD-VQE reached the exact energy at every geometry to numerical precision, as expected for two
  electrons in a minimal basis, where UCCSD can represent the exact state; the noise column shows how far the same
  circuit would drift on today's hardware. The hardware-efficient rows (reported, not gated) showed two things. On HeH+
  the circuit, which does not conserve the electron count, converged to three-electron states whose energies lie below
  the molecule's exact energy: physically wrong. After the test the app's hardware-efficient option was given the
  standard electron-number penalty (`chem.vqe`, default 2 hartree per electron squared; the rows above are the run as
  it happened). On stretched H2 it stopped in local minima, the known cost of generic circuits.
- **Q4 failed**, narrowly: the quantum-kernel SVM transferred from GM12878 to IMR-90 above chance (its interval is above
  0.5) but trailed the RBF-SVM by a little more than the margin, and plain logistic regression beat both. A kernel
  estimated from shots gave about the same result as the exact kernel. The training set is small (106 genes).
- **Overall.** The lab shows that ChronoCell's questions can be written as quantum algorithms and that those algorithms
  do what they are designed to do on a simulator; on these real-data tests they do not beat the classical methods, and
  the app says so on every quantum panel.

Summary:

<!-- BEGIN generated:summary_q -->
| Test (held-out, real data; simulated quantum) | Measured | Verdict |
|---|---|---|
| Gate Q1: QAOA (simulator) finds the optimum of the domain QUBO, 66 held-out windows | 98 % of windows (needed ≥ 90 %) | pass |
| Gate Q2: quantum domain calls vs classical callers (F1 vs ENCODE Arrowhead) | k562: QAOA 0.24 vs insulation 0.35, TopDom-like 0.39; imr90: QAOA 0.27 vs insulation 0.36, TopDom-like 0.40 | fail |
| Gate Q3: VQE within chemical accuracy (H2, HeH+) | worst error 5.1e-11 mHa (needed ≤ 1.6) | pass |
| Gate Q4: quantum-kernel gene classifier vs RBF-SVM, GM12878 → IMR-90 | AUC 0.623 vs 0.657 | fail |
<!-- END generated:summary_q -->

The simulator against Qiskit (a software check, not a gate; `validation/quantum_crosscheck.py`):

<!-- BEGIN generated:quantum_crosscheck -->
Qiskit 2.2.1 statevectors of the exported OpenQASM circuits vs ChronoCell's simulator:

| Circuit | Qubits | Gates | State fidelity |
|---|---|---|---|
| QAOA domain walls (11 qubits, p=2) | 11 | 135 | 1.000000000000 |
| UCCSD-VQE H2 | 4 | 198 | 1.000000000000 |
| Hardware-efficient VQE HeH+ | 4 | 20 | 1.000000000000 |
| ZZ feature map (5 qubits, 2 reps) | 5 | 80 | 1.000000000000 |
| Swap test with state preparation (5 qubits) | 5 | 14 | 1.000000000000 |
| random circuit, every gate | 5 | 60 | 1.000000000000 |

Verdict: agree.
<!-- END generated:quantum_crosscheck -->

## Gate 8 — the Drug lab against chromatin tracing after real drug treatment (pre-registered)

The Drug lab (`chronocell/therapy.py`) is a mechanism simulator: each drug class is reduced to where it acts and which
way it pushes. Gate 8 asks whether its predicted 3D change agrees with what real drugs do, using MINA chromatin tracing
of an 840 kb region of chrX in IMR-90 cells (Cheng et al., Genome Biol 2021; 4DN) before and after GSK126 (EZH2
inhibitor), Trichostatin A + sodium butyrate (HDAC inhibitors), 5-aza-2'-deoxycytidine (DNMT inhibitor) and DMOG
(demethylase blocker). The prediction uses only the untreated traces and IMR-90 H3K27ac; the comparison is the pattern
of pair-distance changes (Spearman rho, global scale removed), with the trace-bootstrap interval and a permutation check
that shuffles where the drug acts. The rule (`validation/frozen.py`, `DRUG_GATE8`) was committed before any treated
test data were read; practice used only untreated split halves and alpha-amanitin. Runner: `validation/drug_gate8.py`.

<!-- BEGIN generated:gate8 -->
Practice (no treated test data): untreated traces split in two halves (no drug, so any agreement is noise) and alpha-amanitin (transcription-inhibitor class):

| Comparison | Allele | Drug-lab class | Spearman rho | 95 % interval | Shuffled-target mean | Permutation p |
|---|---|---|---|---|---|---|
| untreated half vs half | Xa | ezh2 | -0.137 |  [-0.265, 0.087] | -0.125 | 0.935 |
| untreated half vs half | Xa | hdac | -0.119 |  [-0.246, 0.126] | -0.112 | 0.667 |
| untreated half vs half | Xa | dnmt | -0.120 |  [-0.256, 0.097] | -0.120 | 1.000 |
| untreated half vs half | Xa | kdm | +0.021 |  [-0.145, 0.179] | +0.013 | 0.189 |
| untreated half vs half | Xa | txn | -0.048 |  [-0.161, 0.131] | -0.081 | 0.179 |
| untreated half vs half | Xi | ezh2 | +0.120 |  [-0.144, 0.308] | +0.139 | 0.965 |
| untreated half vs half | Xi | hdac | +0.135 |  [-0.153, 0.332] | +0.161 | 0.806 |
| untreated half vs half | Xi | dnmt | +0.103 |  [-0.122, 0.293] | +0.103 | 1.000 |
| untreated half vs half | Xi | kdm | +0.012 |  [-0.181, 0.207] | +0.002 | 0.249 |
| untreated half vs half | Xi | txn | +0.069 |  [-0.124, 0.178] | +0.105 | 0.826 |
| alpha-amanitin vs untreated | Xa | txn | +0.293 |  [0.165, 0.343] | +0.308 | 0.856 |
| alpha-amanitin vs untreated | Xi | txn | +0.171 |  [0.065, 0.246] | +0.219 | 0.955 |

Test (run once; pass per drug on Xa: rho >= 0.2, 95 % interval above 0, permutation p <= 0.05; gate: at least 3 of 4 drugs):

| Drug (class) | Allele | Spearman rho | 95 % interval | Shuffled-target mean | Permutation p | Measured global log change | Traces (untreated / treated) | Pass |
|---|---|---|---|---|---|---|---|---|
| GSK126 (EZH2 inhibitor) | Xa | +0.033 |  [-0.122, 0.263] | -0.003 | 0.030 | -0.022 | 508 / 191 | no |
| GSK126 (EZH2 inhibitor) | Xi | +0.361 |  [0.104, 0.467] | +0.349 | 0.159 | -0.077 | 512 / 191 | — |
| TSA + NaBu (HDAC inhibitors) | Xa | +0.293 |  [0.089, 0.359] | +0.274 | 0.109 | +0.221 | 508 / 219 | no |
| TSA + NaBu (HDAC inhibitors) | Xi | -0.052 |  [-0.200, 0.174] | -0.084 | 0.189 | +0.170 | 512 / 201 | — |
| 5-aza-dC (DNMT inhibitor) | Xa | +0.045 |  [-0.119, 0.200] | +0.045 | 1.000 | +0.067 | 508 / 118 | no |
| 5-aza-dC (DNMT inhibitor) | Xi | -0.012 |  [-0.173, 0.186] | -0.012 | 1.000 | +0.155 | 512 / 118 | — |
| DMOG (demethylase blocker) | Xa | -0.102 |  [-0.172, 0.007] | -0.092 | 0.657 | +0.129 | 508 / 383 | no |
| DMOG (demethylase blocker) | Xi | -0.057 |  [-0.162, 0.113] | -0.086 | 0.055 | +0.065 | 512 / 336 | — |

0 of 4 drugs pass. Verdict: **fail**.
<!-- END generated:gate8 -->

**Reading (written after the test run).** Gate 8 **failed**: no drug met the rule on the active X. Where the
simulator's predicted pattern did agree with the measured one (the HDAC inhibitors on the active X; GSK126 on the
inactive X, reported only), shuffling where the drug acts agreed about as well, so the agreement comes from the
simulator's generic way of opening or compacting a fold, not from targeting the right loci; the same was seen with
alpha-amanitin in practice. The overall size of the region grew after the HDAC inhibitors and after 5-aza-dC, the
direction these opening classes predict, and also after DMOG, which the compacting demethylase-blocker class does not
predict; overall size is not part of the rule because it is not comparable between imaging experiments. The Drug lab
therefore stays labelled "mechanism simulator, not validated", now with this measurement behind the label, and the
Drug guide shows the result. The test covered one 840 kb region in one cell line; drugs that act mostly elsewhere in the
genome would not show here.


## Gates Q5–Q7 — quantum drug tabs: heart safety, molecules, docking (pre-registered)

Three drug-focused quantum tabs, each tested on held-out data against classical methods on the same input
(`validation/quantum_drug_gates.py`; rules in `validation/frozen.py`, `QUANTUM_DRUG_GATES`):
- **Q5, heart safety.** A quantum-kernel SVM screens molecules for hERG blocking from SMILES descriptors
  (`chronocell/quantum/molfeat.py`, checked against RDKit). Trained on TDC hERG, tested on the independent
  TDC hERG_Karim set, against an RBF-SVM and logistic regression.
- **Q6, molecules.** The from-scratch STO-3G chemistry (`chronocell/quantum/molecules.py`) against OpenFermion's
  independently computed reference energies (test: LiH), and active-space UCCSD-VQE on stretched molecules against
  the exact energy in the same space.
- **Q7, docking.** Rigid re-docking of PoseBusters ligands as a maximum-weight clique solved by QAOA on 20 qubits,
  against the same pipeline with classical clique enumeration and against random search with the same score.

Practice (settings; `validation/TUNING.md` §20):

<!-- BEGIN generated:qdrug_practice -->
Q5 practice: TDC hERG, 655 compounds read (0 unreadable), 69 % blockers. SMILES descriptors against RDKit (share exact or within tolerance; Pearson r): mw 100 % (r 1.000); hbd 100 % (r 1.000); hba 100 % (r 1.000); tpsa 69 % (r 0.991); rot_bonds 79 % (r 0.986); rings 98 % (r 0.994); aromatic_rings 100 % (r 1.000); heavy_atoms 100 % (r 1.000); fsp3 100 % (r 1.000).

| Features | Quantum-kernel SVM (CV AUC) | RBF-SVM | Logistic regression |
|---|---|---|---|
| lipinski6 | 0.853 | 0.849 | 0.849 |
| herg8 | 0.862 | 0.860 | 0.855 |
| pca8 | 0.858 | 0.858 | 0.852 |

Chosen: {'features': 'herg8', 'qsvm': {'bandwidth': 0.005, 'reps': 2, 'C': 1000.0}, 'rbf': {'gamma': 0.002, 'C': 100.0}}.

Q6 practice: STO-3G Hartree-Fock against Szabo & Ostlund (Table 3.13, printed to 1 mEh): H2 +0.29 mEh; CO +0.42 mEh; N2 +0.16 mEh; CH4 +0.15 mEh; NH3 -0.07 mEh; H2O +0.06 mEh; HF +0.21 mEh. Against OpenFermion's H2 data (HF / FCI, mEh): 0.5 A +5.8e-05 / +5.1e-05; 0.7414 A +5.0e-05 / +3.9e-05; 1.0 A +7.0e-05 / +4.9e-05; 1.5 A +1.2e-04 / +4.2e-05; 2.0 A +2.0e-04 / +1.9e-05; 2.5 A +2.5e-04 / -3.1e-06.

Q7 practice (39 usable complexes; settings {'tau': 1.0, 'max_polar': 8, 'max_hydrophobic': 4, 'per_type': 8, 'qubits': 20, 'cliques': 30, 'random_poses': 1000}):

| QAOA | Max clique found | Mean P(best clique) | Uniform | Docked (QAOA route) | Docked (classical cliques) | Docked (random search) | Mean QAOA s |
|---|---|---|---|---|---|---|---|
| p3_expectation | 82 % | 0.0011 | 1.4e-06 | 15 % | 15 % | 3 % | 3.5 |
| p3_cvar | 79 % | 0.0011 | 1.4e-06 | 15 % | 15 % | 3 % | 3.6 |
| p5_expectation | 90 % | 0.0013 | 1.4e-06 | 15 % | 15 % | 3 % | 8.5 |
| p5_cvar | 97 % | 0.0027 | 1.4e-06 | 15 % | 15 % | 3 % | 8.8 |

Chosen: {'p': 5, 'objective': 'cvar'}.
<!-- END generated:qdrug_practice -->

Test (run once):

<!-- BEGIN generated:qdrug -->
Q5 (trained on 655 TDC hERG compounds, tested on 12989 hERG_Karim compounds not in the training set (456 removed), 50 % blockers):

| Classifier | Test AUC (95 % interval) |
|---|---|
| Quantum-kernel SVM (simulated) | 0.710 [0.701, 0.718] |
| RBF-kernel SVM (classical) | 0.695 [0.686, 0.703] |
| Logistic regression (classical) | 0.686 [0.676, 0.694] |

Quantum minus RBF +0.015 [0.013, 0.017]. Q5: **pass** (needed: within 0.03 of the RBF-SVM, lower bound above 0.5).

Q6 (independent reference: OpenFermion's stored data; tolerance 0.1 mEh):

| Reference | HF (ours) | HF (reference) | Difference (mEh) | FCI (ours) | FCI (reference) | Difference (mEh) |
|---|---|---|---|---|---|---|
| H1-Li1_sto-3g_singlet_1.45.hdf5 | -7.862568 | -7.862568 | +2.2e-04 | -7.880982 | -7.880982 | +2.2e-04 |

Active-space UCCSD-VQE on stretched molecules (chemical accuracy 1.6 mHa):

| Molecule | Bond × equilibrium | Active space | Qubits | Parameters | HF | Active-space FCI | VQE | Error (mHa) |
|---|---|---|---|---|---|---|---|---|
| H2O | 1.5 | 4e, 4o | 8 | 26 | -74.74717 | -74.78666 | -74.78666 | +0.003 |
| NH3 | 1.5 | 6e, 5o | 10 | 54 | -55.12443 | -55.21169 | -55.21114 | +0.544 |
| N2 | 1.5 | 6e, 6o | 12 | 117 | -106.94677 | -107.50580 | -107.42744 | +78.353 |
| HF | 2.0 | 2e, 2o | 4 | 3 | -98.30081 | -98.47485 | -98.47485 | -0.000 |
| CH2O | 1.3 | 4e, 4o | 8 | 26 | -112.21697 | -112.32501 | -112.32499 | +0.021 |
| HCN | 1.3 | 4e, 4o | 8 | 26 | -91.42908 | -91.60457 | -91.60188 | +2.698 |

Q6: **fail** (reference within tolerance; VQE outside chemical accuracy).

Q7 (200 usable test complexes of 213; 20-qubit interaction graphs):

| Measure | Value |
|---|---|
| QAOA found the maximum-weight clique | 100 % (needed ≥ 80 %) |
| Mean probability of the best clique (uniform guess) | 0.0069 (2.5e-03) |
| Docked within 2 A: QAOA route | 5.5 % |
| Docked within 2 A: classical cliques, same graph | 6.5 % |
| Docked within 2 A: random search, same score, 1,000 poses | 2.0 % |
| Mean QAOA time per complex | 8.6 s |

Q7: **pass** (solver pass, docking pass: QAOA route at least as good as random search).
<!-- END generated:qdrug -->

**Reading (written after the test run).**
- **Q5 passed.** On almost 13,000 compounds the screen had never seen, from a different laboratory's measurements, the
  quantum-kernel SVM ranked hERG blockers better than the RBF-SVM and logistic regression trained on the same 655
  compounds; the margin is small but its interval excludes zero. Transfer between the two datasets is modest for every
  model (AUC about 0.7, against 0.86 in cross-validation on the training set). The kernel is computed exactly on a
  simulator, so this is a well-chosen kernel, not a quantum speed-up.
- **Q6 failed** on its VQE half. The from-scratch chemistry agreed with OpenFermion's independently computed LiH energies
  far inside the tolerance, and UCCSD-VQE reached chemical accuracy on four of the six stretched molecules; it missed on
  nitrogen at 1.5 times its bond length (breaking a triple bond, where a single Hartree-Fock reference describes the
  electrons badly and a disentangled UCCSD started from it stalls) and, by a little, on hydrogen cyanide at 1.3 times.
  This is the known limit of the method, now measured; the panel shows the exact active-space energy next to VQE.
- **Q7 passed** its pre-registered rule, modestly: QAOA found the best clique in every test graph, and the QAOA route
  docked more test ligands correctly than random search with the same score, but far fewer than on practice, about as
  many as the classical clique route on the same graph. The practice settings were tuned on 39 complexes and did not
  carry over; with a 20-qubit graph and a simple score, this docking is a demonstration of the formulation, not a
  docking tool.

Summary:

<!-- BEGIN generated:summary_qd -->
| Test (held-out, real data) | Measured | Verdict |
|---|---|---|
| Gate 8: Drug lab vs chromatin tracing after real drug treatment (IMR-90 chrX, 4 drugs) | 0 of 4 drugs met the rule | fail |
| Gate Q5: quantum-kernel hERG screen vs RBF-SVM (TDC hERG → hERG_Karim) | AUC 0.710 vs 0.695 | pass |
| Gate Q6: molecule energies vs OpenFermion; stretched-molecule VQE | worst reference difference 2.2e-04 mHa; worst VQE error 78.35 mHa | fail |
| Gate Q7: QAOA max-clique docking (PoseBusters) | docked 6 % vs random search 2 %; clique found 100 % | pass |
<!-- END generated:summary_qd -->

## Cost (Pillar 1): runtime and peak memory against bead count

`python validation/scale_benchmark.py`: each run in its own process. Input: windows of the
**synthetic** reference chr22 contact map, up to the whole chromosome (5,082 beads at 10 kb). These
rows measure cost only, not accuracy.

<!-- BEGIN generated:scale -->
Machine: Windows-10-10.0.26200-SP0, AMD64 Family 25 Model 117 Stepping 2, AuthenticAMD, 16 logical CPUs, 23.1 GB RAM, no GPU (CPU only); torch 2.11.0+cpu, Python 3.10.11. Synthetic input: cost only, not accuracy. Each row ran in its own process, one at a time; the last column is the whole-machine CPU load measured just before the row (other work slows timings).

| Beads | Model | Fit (s) | Sampling (s) | Total (s) | Peak memory (MB) | Device | Rank | CPU load before |
|---|---|---|---|---|---|---|---|---|
| 200 | v3.3 | 7.5 | 7.7 | 15.2 | 1,156 | cpu | 200 | 11 % |
| 200 | v4 | 1.6 | 4.4 | 6.0 | 1,934 | cpu | 200 | 19 % |
| 400 | v3.3 | 15.8 | 17.2 | 33.2 | 1,162 | cpu | 400 | 18 % |
| 400 | v4 | 3.8 | 5.6 | 9.6 | 1,973 | cpu | 256 | 12 % |
| 800 | v3.3 | 75.0 | 49.8 | 125.1 | 1,525 | cpu | 800 | 22 % |
| 800 | v4 | 11.7 | 6.6 | 18.6 | 2,064 | cpu | 256 | 23 % |
| 1,600 | v3.3 | 439.5 | 218.2 | 659.2 | 4,699 | cpu | 1600 | 18 % |
| 1,600 | v4 | 40.5 | 8.9 | 50.9 | 2,372 | cpu | 256 | 15 % |
| 3,200 | v3.3 | — | — | — | — | — | dense v3.3 fit not attempted above 1600 beads | |
| 3,200 | v4 | 158.6 | 12.9 | 175.2 | 3,049 | cpu | 256 | 22 % |
| 5,082 | v3.3 | — | — | — | — | — | dense v3.3 fit not attempted above 1600 beads | |
| 5,082 | v4 | 485.1 | 22.6 | 513.7 | 4,321 | cpu | 256 | 15 % |
<!-- END generated:scale -->

**Reading.** The v4 model fits the whole synthetic chr22 (5,082 beads) in under 9 minutes on this
CPU, at about 4,300 MB peak memory. The dense v3.3 fit grows much faster with size: at 1,600 beads it takes
about 13 times as long as v4 and twice the memory, which is why the app switches to v4 above 400
beads. These are cost figures on synthetic input; accuracy is Gates 1–3.

## Per-chromosome runtime (Pillar 6)

<!-- BEGIN generated:per_chromosome -->
Machine: Windows-10-10.0.26200-SP0, AMD64 Family 25 Model 117 Stepping 2, AuthenticAMD, 16 logical CPUs, 23.1 GB RAM, no GPU (CPU only); torch 2.11.0+cpu, Python 3.10.11. Synthetic input (one planted reference per chromosome at the app's default resolution, at most 6,000 beads), fitted with the app's settings; cost only, not accuracy. Each row ran in its own process, one at a time; the last column is the whole-machine CPU load measured just before the row (other work slows timings).

| Assembly | Chromosome | Resolution | Bins | Assembled beads | Model | Fit (s) | Total (s) | Peak memory (MB) | CPU load before |
|---|---|---|---|---|---|---|---|---|---|
| hg38 | chr4 | 40 kb | 4,756 | 4,745 | population_v4 | 407.8 | 435.6 | 4,088 | 22 % |
| hg38 | chr5 | 40 kb | 4,539 | 4,532 | population_v4 | 353.9 | 381.3 | 3,914 | 15 % |
| hg38 | chr6 | 40 kb | 4,271 | 4,253 | population_v4 | 299.9 | 324.7 | 3,704 | 15 % |
| hg38 | chr7 | 40 kb | 3,984 | 3,974 | population_v4 | 302.9 | 329.8 | 3,555 | 25 % |
| hg38 | chr8 | 25 kb | 5,806 | 5,791 | population_v4 | 615.9 | 651.3 | 4,488 | 14 % |
| hg38 | chr9 | 25 kb | 5,536 | 4,873 | population_v4 | 548.9 | 581.9 | 4,679 | 7 % |
| hg38 | chr10 | 25 kb | 5,352 | 5,332 | population_v4 | 537.6 | 571.9 | 4,578 | 18 % |
| hg38 | chr11 | 25 kb | 5,404 | 5,382 | population_v4 | 541.1 | 573.8 | 4,606 | 20 % |
| hg38 | chr12 | 25 kb | 5,332 | 5,326 | population_v4 | 528.8 | 559.8 | 4,549 | 14 % |
| hg38 | chr13 | 20 kb | 5,719 | 4,900 | population_v4 | 607.7 | 644.1 | 4,680 | 16 % |

10 chromosomes; total 84.2 min of fitting; slowest 651 s; largest peak memory 4,680 MB.

**Status: partial: the run was stopped at this machine's 2-hour job limit after 13 chromosomes; continue with python validation/chromosome_runtime.py --resume.** Not measured yet (35 of 45): hg38 chr1, hg38 chr2, hg38 chr3, hg38 chr14, hg38 chr15, hg38 chr16, hg38 chr17, hg38 chr18, hg38 chr19, hg38 chr20, hg38 chr21, hg38 chr22, hg38 chrX, hg38 chrY, mm39 chr1, mm39 chr2, mm39 chr3, mm39 chr4, mm39 chr5, mm39 chr6, mm39 chr7, mm39 chr8, mm39 chr9, mm39 chr10, mm39 chr11, mm39 chr12, mm39 chr13, mm39 chr14, mm39 chr15, mm39 chr16, mm39 chr17, mm39 chr18, mm39 chr19, mm39 chrX, mm39 chrY.
- Discarded: hg38 chr1 (472 s): measured while the full test suite ran on the same machine (02:33-02:50 IST, 5 October 2026); discarded so that --resume re-measures it.
- Discarded: hg38 chr2 (544 s): measured while the full test suite ran on the same machine (02:33-02:50 IST, 5 October 2026); discarded so that --resume re-measures it.
- Discarded: hg38 chr3 (488 s): measured while the full test suite ran on the same machine (02:33-02:50 IST, 5 October 2026); discarded so that --resume re-measures it.
<!-- END generated:per_chromosome -->

**Reading.** Partial. On the hg38 chromosomes measured so far, the app's default resolution (about
4,000–5,800 beads) means 5–11 minutes and 3,600–4,700 MB of peak memory per chromosome on this CPU. The run was
stopped by this machine's 2-hour job limit; `python validation/chromosome_runtime.py --resume`
measures the rest (and the three discarded rows) and `python validation/report.py` updates this
table.

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
