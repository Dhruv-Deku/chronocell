"""
Settings frozen after practice-only tuning (validation/TUNING.md, sections 5-8). Test scripts import
these and nothing else; changing a value here after a test run invalidates that run, and git history
records when each value was set.
"""

# Pillar 1: whole-chromosome population model (chronocell.population.PopulationConfig fields).
# Practice evidence (Su chr2, 935 loci, split 0; validation/tuning_v4.json): full rank, rank 256 and
# rank 128 (+ random walk) gave 84.9 / 84.9 / 85.0 % of the ceiling within windows and 80.6 / 80.7 /
# 80.6 % on all pairs; 3,000 iterations, lr 0.02, float64 and the coarse start changed nothing;
# uniform weights were worse (82.8 / 79.0 %). Rank 256 is kept as a margin over 128 at 2x its speed.
WHOLE_CHROMOSOME = {"rank_cap": 256, "local_term": True, "iterations": 1500, "learning_rate": 0.01,
                    "weighting": "binomial", "dtype": "float32", "init": "auto"}

# Gate 1 (validation/gate1.py --test). Imaging input: contact radius = half A's median adjacent
# distance (pre-registered rule), whole-chromosome model = WHOLE_CHROMOSOME, windowed model = v3.3
# defaults. Hi-C input (Rao 2014 binned on the imaged loci by Su et al.): practice Su chr2, split 0;
# 13 of 16 variants had finished when this was frozen (validation/tuning_v4.json): raw counts beat ICE
# in all 13; zero-count pairs as
# unobserved beat clipping; p_adjacent 0.3 gave the best all-pairs score (32.0 % of the ceiling,
# whole model). The nm anchor is the literature b0 for the locus spacing (physics.bond_length_for)
# and, separately, that b0 times the practice ratio of measured adjacent distance to literature b0
# on chr2 (461.0 nm / 146.2 nm = 3.153); only CCC depends on it.
GATE1 = {"population": WHOLE_CHROMOSOME, "p_adjacent": 0.3, "hic_zeros": "unobserved",
         "hic_anchors": {"literature_b0": 1.0, "practice_calibrated_b0": 3.153}}

# Gate 4, structural variants (validation/sv_validation.py). Pre-registered before that script read any
# K562 counts on chr9 (the only K562 data read before were chr21 counts for the practice benchmark).
# Event: the two chr9 regions "lost entirely" in K562 (Zhou et al., Genome Res 29:472-484, 2019; hg19),
# which remove CDKN2A / CDKN2B / MTAP. Reference = GM12878 (normal karyotype), variant = K562, both Rao
# et al. 2014 in situ Hi-C (GSE63525, MAPQ >= 30 maps), raw counts. Deleted bins: midpoint inside a
# published interval; the kept runs are assumed joined in order (simple interstitial deletions).
# Model: v4 population fitted on the GM12878 window with the Gate 1 Hi-C settings, then the deletions
# applied by perturb.deletion_pieces / derive (no parameter is fitted on K562).
# Predictions of K562 counts for pairs on different kept runs ("spanning pairs") with new separation
# min..max bins: ours = p_after; baselines = distance shift (GM12878's mean count at the new
# separation) and no change (GM12878's count for the pair).
# Metrics: Spearman of prediction vs K562 counts (primary); trend-removed Spearman (each map divided by
# its own mean over within-run pairs at the same new separation; ours only, the baselines carry no
# pattern beyond separation); 95 % intervals from a block bootstrap (blocks of block_bins x block_bins).
# Pass ("validated on this event"): (i) ours beats the distance shift, 95 % interval of the Spearman
# difference above 0, and (ii) the trend-removed Spearman's 95 % interval is above 0. Anything else:
# the SV simulator stays "mechanism simulator, not validated" in the app.
SV_VALIDATION = {"reference": "rao2014_gm12878", "variant": "rao2014_k562", "chrom": "chr9", "assembly": "hg19",
                 "start": 16_000_000, "end": 36_000_000, "binsize": 25_000,
                 "deleted": [(20_750_000, 26_590_000), (28_560_000, 31_620_000)],
                 "deleted_source": "Zhou et al., Genome Res 29:472-484 (2019), K562: two chr9 regions lost entirely",
                 "population": WHOLE_CHROMOSOME, "p_adjacent": GATE1["p_adjacent"], "hic_zeros": GATE1["hic_zeros"],
                 "min_new_separation_bins": 2, "max_new_separation_bins": 80, "block_bins": 8, "bootstrap": 1000,
                 "seed": 0}

# Gate 5 (Pillar 5, no-Hi-C prediction; validation/predictor.py). Pre-registered before any test run.
# Model: chronocell.predict (CTCF peaks oriented by the JASPAR MA0139.1 motif + GC per locus -> ridge
# regression of the trend-removed log median distance). Inputs per dataset: ENCODE CTCF IDR peaks of
# the same cell line (validation/datasets.py CTCF_PEAKS) and the hg38 sequence. Training: the untreated
# practice datasets, each weighted equally; the ridge penalty is chosen on practice data only, by
# leave-one-dataset-out cross-validation over ridge_grid, then the model is refitted on all of them and
# frozen in validation/predictor_model.json. Test (run once): the untreated test datasets; the cohesin-
# depleted test set is reported as a control (CTCF-anchored loops need cohesin) and is not in the rule.
# Scoring: validation/benchmark/run.py protocol (half B medians as truth, all pairs, mean over splits,
# 95 % interval from split 0). Baseline: the same training trend log d = a + b log10 s with no features
# (genomic distance only, no data from the target region).
# Pass (keep in the UI): on at least 3 of the 5 test datasets, (i) % of ceiling > 0 with the 95 %
# interval above 0, and (ii) raw Spearman above the baseline's. Otherwise: documented, not in the UI.
PREDICTOR = {"train": ["bintu_k562_28_30", "bintu_hct116_28_30", "bintu_hct116_34_37", "su_chr2"],
             "test": ["bintu_imr90_28_30", "bintu_imr90_18_20", "bintu_a549_28_30", "su_chr21", "su_chr21_rep"],
             "control": ["bintu_hct116_34_37_auxin"],
             "motif": "JASPAR MA0139.1", "motif_min_relative_score": 0.8, "pseudocount": 0.25,
             "summit_half_width_bp": 100, "ridge_grid": [1e-4, 1e-3, 1e-2, 1e-1, 1.0], "pass_min_datasets": 3}

# Gate 2b (Pillar 2, sequencing Hi-C input; validation/calibration.py --input hic). Pre-registered before
# the practice fit and before any test run. Gate 2 recalibrated the intervals for imaging-derived
# contacts only; with Hi-C input the practice benchmark (K562 chr21:28-30 Mb) covered 33 % raw and 50 %
# with the imaging recalibration at the stated 90 % level. Question: does a recalibration fitted on
# practice Hi-C input make the stated intervals honest for Hi-C input?
# Input: the Gate 1 Hi-C settings (raw counts, zero-count pairs unobserved, p_adjacent 0.3) and the
# app's length anchor (adjacent beads at the literature b0 for the locus step); Rao et al. 2014 Hi-C
# of the same cell line on the imaged loci. Model: v3.3 up to 400 loci, WHOLE_CHROMOSOME above (as
# Gate 2). Truth: half B's single-copy distances, split 0. Recalibration: the Gate 2 method unchanged
# (quantile recalibration of the pooled PIT histogram, every practice dataset weighted equally),
# frozen in chronocell/data/calibration_hic.json before the test run.
# Pass ("usable for sequencing input"): on EVERY test dataset the recalibrated stated 90 % interval
# holds 83-97 % of half B's single-copy distances AND the stated 50 % interval holds 40-60 %.
# Otherwise: the app shows the interval for sequencing input with the measured coverage and the
# warning that it is far below nominal. Reported alongside, not in the rule: raw coverage, coverage
# with the imaging recalibration, and the width of the recalibrated 90 % interval (upper / lower).
HIC_CALIBRATION = {"practice": ["bintu_k562_28_30", "su_chr2", "su_chr2_parm_rep"],
                   "test": ["bintu_imr90_28_30", "bintu_imr90_18_20", "su_chr21", "su_chr21_rep"],
                   "split": 0, "p_adjacent": GATE1["p_adjacent"], "hic_zeros": GATE1["hic_zeros"],
                   "anchor": "literature_b0", "pass_90": (0.83, 0.97), "pass_50": (0.40, 0.60)}

# Gate 5m (mouse transfer of the Gate 5 predictor). Pre-registered before any mouse data were read, and
# NOT RUN: the traces are on the 4DN Data Portal, whose downloads need a (free) 4DN account access key that
# this machine does not have. The app keeps the predictor hg38-only until this test has run and passed.
# Truth: ORCA chromatin tracing in mouse ES cells (E14TG2a), Hafner et al., Mol Cell 83:1377 (2023), GRCm38
# (mm10), two 3-Mb loci: 4DN experiment sets 4DNESD28H8O7 (chr6, untreated) and 4DNESWDXDZSE (chr3,
# untreated) = the test; the same loci in the CTCF-AID and RAD21-AID lines without auxin are reported as
# secondary; with auxin (CTCF or cohesin degraded) as controls, not in the rule.
# Inputs: ENCODE CTCF IDR thresholded peaks of mouse ES-Bruce4 (ENCFF533APC, experiment ENCSR000CCB, mm10,
# MD5 d2c0c0c2ed1d8e11989362580ca85af7; a different ES line from the imaged cells, the only ENCODE mouse ES
# CTCF ChIP-seq) and the UCSC mm10 sequence (md5sum.txt). Model: validation/predictor_model.json unchanged
# (trained on human data only; nothing is refitted on mouse). Scoring and pass rule: Gate 5's, unchanged
# (half-B medians as truth, all pairs, 3 splits, 95 % interval from split 0): pass if BOTH test loci have
# (i) % of ceiling > 0 with its 95 % interval above 0 and (ii) raw Spearman above the training-trend
# baseline's. Pass -> the predictor may be offered for mouse assemblies (mm39), labelled with this result.
PREDICTOR_MOUSE = {"status": "pre-registered, not run (4DN access key needed)",
                   "test": ["4DNESD28H8O7", "4DNESWDXDZSE"],
                   "secondary": ["4DNESJ3TXVIR", "4DNESQ49IXDU", "4DNESTNG39BO", "4DNESLRTOSQT"],
                   "control": ["4DNES2KX6HQ5", "4DNESMN7RCSB", "4DNESG62SAVA", "4DNESBH54BG2"],
                   "peaks": ("ENCFF533APC", "ENCSR000CCB", "d2c0c0c2ed1d8e11989362580ca85af7"), "assembly": "mm10",
                   "model": "validation/predictor_model.json", "splits": 3, "pass": "both test loci meet Gate 5 (i) and (ii)"}

# Gate 2c (Pillar 2, per-pair reliability; validation/reliability.py). Pre-registered after the practice
# comparison (validation/results_reliability_practice.json) and before any test run.
# Score: "misfit" = minus |log sigma_model^2 - log sigma_target^2| of the pair, from the input and the fitted
# model only. It has the best worst case on practice data: imaging input, pair level misfit +0.057 to +0.175,
# combined +0.049 to +0.245, input_se -0.031 to +0.266; bead-level minimum -0.147 vs -0.208 and -0.281. On
# practice Hi-C input every candidate was negative (misfit the least negative on Su chr2 and its replicate,
# -0.080 vs input_se -0.079 on K562); misfit is tested there unchanged.
# Error of a pair: |r - median of r in its separation stratum| with r = log(model / half-B median), 10
# equal-count strata. Pair level: Spearman(score, -error) within strata, averaged; bead level: median pair
# score vs median pair error. 95 % intervals: 200 resamples of the loci. Test sets: Gate 2's six (imaging
# input) and Gate 2b's four (Hi-C input), split 0.
# Pass, separately for each input type and level: on EVERY test set Spearman >= min_spearman (0.20, the
# usability bar validation/report.py already applied to Gate 2's per-bead score) with the 95 % interval
# above 0. Only a passing input type and level may show a reliability in the app (probe: pair level; 3D
# view: bead level), labelled with its measured standing; otherwise none is shown.
RELIABILITY = {"score": "misfit", "min_spearman": 0.20, "strata": 10, "bootstrap": 200, "split": 0,
               "test_imaging": ["bintu_imr90_28_30", "bintu_imr90_18_20", "bintu_a549_28_30",
                                "bintu_hct116_34_37_auxin", "su_chr21", "su_chr21_rep"],
               "test_hic": ["bintu_imr90_28_30", "bintu_imr90_18_20", "su_chr21", "su_chr21_rep"]}

# Gate 5b (Phase A5, prediction without contact data from more inputs; validation/predictor_v2.py).
# Pre-registered after the practice comparison (results_predictor_v2_practice.json) and before any test run.
# Model frozen in validation/predictor_v2_model.json: Gate 5's 12 CTCF + GC features plus ENCODE RAD21 (cohesin)
# peak features (per-locus signal similarity and anchors, RAD21 sites between the loci, each also times
# log10 separation), ridge 0.01, chosen by leave-one-dataset-out over the three non-IMR-90 untreated practice
# regions it was trained on (K562 chr21:28-30, HCT116 chr21:28-30 and 34-37), so every test set is a held-out
# cell type. ATAC and H3K27ac did not help on practice and are not used. The frozen Gate 5 model is unchanged
# and stays the app's default.
# Test sets and scoring: Gate 5's (PREDICTOR["test"], the cohesin-depleted set as a control, half-B medians,
# all pairs, 3 splits for the Bintu sets, 95 % interval from split 0).
# Pass ("may be offered as an option"): on at least 3 of the 5 test sets, % of the reproducible pattern >= 50
# with its 95 % interval above 0 AND higher than the frozen Gate 5 model's % on that set
# (validation/results_predictor.json).
PREDICTOR_V2 = {"train": ["bintu_k562_28_30", "bintu_hct116_28_30", "bintu_hct116_34_37"], "test": PREDICTOR["test"],
                "control": PREDICTOR["control"], "marks": ["rad21"], "min_percent": 50.0, "min_sets": 3,
                "beat_gate5": True}

# Gate 1c (Phase A1, calibrated sizes from Hi-C; validation/hic_size_calibration.py). Pre-registered after the
# practice comparison (results_hic_size_practice.json) and before any test unit was built. Calibration frozen in
# chronocell/data/hic_size_calibration.json: model distance x exp(h); leave-one-dataset-out over the practice
# units (Bintu K562 / HCT116 regions with Rao 2014, ENCODE in situ and intact Hi-C at 1, 1/4, 1/16 depth; Su chr2
# and its replicate; the genome-scale practice set) chose the single global factor (h = 0.8953, x2.45): every
# form with separation, locus-spacing, depth or protocol terms extrapolated worse to the held-out dataset.
# Test sets (split 0, Hi-C input with the Gate 1 settings and the app's b0 anchor): main = Bintu IMR-90
# chr21:28-30 and 18-20 Mb (Rao 2014), Su chr21 and its replicate (Rao 2014 binned by Su et al.), the genome-scale
# set su_genome (chromosome units pooled); secondary, reported only = the Bintu IMR-90 regions with ENCODE intact
# and dilution Hi-C, Bintu A549 with ENCODE in situ Hi-C, and su_genome_amanitin (treated cells, untreated Hi-C).
# Pass ("calibrated Hi-C sizes" may be offered as an option): on EVERY main set, Lin's CCC >= 0.8 AND median size
# ratio (calibrated / measured) in 0.8-1.25 AND the trend-removed Spearman not lower than the uncalibrated model's
# (identical by construction: one factor for all pairs; checked to 1e-9).
HIC_SIZE = {"min_ccc": 0.8, "size_ratio": (0.8, 1.25), "pattern_tolerance": 1e-9, "split": 0,
            "main": ["bintu_imr90_28_30", "bintu_imr90_18_20", "su_chr21", "su_chr21_rep", "su_genome"],
            "secondary": ["bintu_imr90_28_30 intact", "bintu_imr90_28_30 dilution", "bintu_imr90_18_20 intact",
                          "bintu_imr90_18_20 dilution", "bintu_a549_28_30 in situ", "su_genome_amanitin"]}

# Gate 3b (Phase A4, a learned correction on top of the population model; validation/learned_correction.py and the
# benchmark method "learned_correction"). Rule, from the Phase A brief and fixed before the practice run: best % of
# the reproducible pattern on >= 90 % of the Hi-C test units of Gate 3's plan, AND no loss of more than 1 point
# against the current model on any imaging-input unit. Practice (results_learned_correction_practice.json;
# leave-one-dataset-out over the practice units) compared no correction, a ridge correction and a small MLP
# (trained on the GPU) on pair features; the pre-specified selection (highest mean held-out pattern rho, and
# beating "none" on every held-out group) chose NO correction: both learned corrections raised Lin's CCC (sizes)
# but lowered the held-out pattern. There is nothing to test, so Gate 3b is not run and its test units stay unused.
LEARNED_CORRECTION = {"hic_best_fraction": 0.90, "imaging_max_loss_points": 1.0, "status": "not run",
                      "reason": "practice chose no correction (pre-specified selection rule)"}

# Gate 2d (Phase A2, honest uncertainty ranges; validation/intervals_v2.py). Pre-registered after the practice
# comparison (results_intervals_v2_practice.json) and before any test unit's truth was read for this question.
# Intervals frozen in chronocell/data/intervals_v2.json: split-conformal quantiles of log(single-copy distance /
# model median distance), learned on practice units only (Bintu K562 / HCT116 regions, Su chr2 and its replicate,
# the genome-scale practice set; split 0), separately per input type. Leave-one-dataset-out on practice chose, by
# the smallest worst-case |coverage at 90 % - 0.90| (ties within 0.005 to the simplest variant): imaging input =
# one set of quantiles per separation band (0-0.1, 0.1-0.3, 0.3-1, 1-3, 3-10, > 10 Mb); Hi-C input = one set
# pooled over bands. The A1 size calibration does not change these intervals (a global factor cancels in the
# ratio) and is not used. A zero distance (identical rounded coordinates) counts as below every interval.
# Test sets (split 0, truth = half B's single copies): imaging input = Gate 2's six test sets and su_genome
# (chromosome units pooled); Hi-C input = Gate 2b's four (Rao 2014 Hi-C, Gate 1 settings, the app's b0 anchor)
# and su_genome. Pass, separately per input type: on EVERY test set and in every separation band with at least
# 200 pair-copies, the 90 % range holds 85-95 % AND the 50 % range holds 40-60 %. The median width of the 90 %
# range (upper / lower bound) is reported next to it, so a trivially wide interval cannot pass silently.
# A passing input type may be offered in the distance probe as an option, its note reading the coverage from
# results_intervals_v2.json; otherwise the probe is unchanged.
INTERVALS_V2 = {"cov90": (0.85, 0.95), "cov50": (0.40, 0.60), "min_pair_copies": 200, "split": 0,
                "test_imaging": ["bintu_imr90_28_30", "bintu_imr90_18_20", "bintu_a549_28_30",
                                 "bintu_hct116_34_37_auxin", "su_chr21", "su_chr21_rep", "su_genome"],
                "test_hic": ["bintu_imr90_28_30", "bintu_imr90_18_20", "su_chr21", "su_chr21_rep", "su_genome"]}

# Gate 4c (Phase B1, cohesin loss against sequencing Hi-C on more regions; validation/cohesin_hic.py).
# Pre-registered after the practice run on the Gate 4 practice region (chr21:28-30 Mb,
# results_cohesin_hic_practice.json) and before any auxin Hi-C of another region was read. Nothing is fitted:
# the Gate 4 cohesin parameters (chronocell/data/perturbation_params.json, fitted on imaging) are used unchanged.
# Data: HCT116 RAD21-mAC, untreated vs 6 h auxin (Rao et al. 2017), ENCODE GRCh38 maps read by region. Main pair:
# in situ Hi-C ENCSR123UVP (ENCFF750AOC) -> ENCSR637QCS (ENCFF301BWY), the data the brief names; secondary,
# reported only: intact Hi-C ENCSR958BEA (ENCFF528XGK) -> ENCSR087JOM (ENCFF317OIA, 5-Ph-IAA).
# Regions: six 2 Mb windows at 10 kb on chromosomes the method has not seen (hg38 starts below).
# Prediction from the untreated map only (Gate 1 Hi-C settings, v3.3 population, r_c 150 nm); measured change =
# log2 of library-normalised (count + 1), auxin over untreated, on pairs >= 2 bins apart with >= 10 reads in the
# two maps together. Score: Spearman of predicted vs measured change; baseline: the same model with lam = 0
# (trend only). A region "beats trend only" if the difference is > 0 with its 95 % interval (200 resamples of
# loci) above 0. Pass: at least 3 of the 6 regions beat trend only on the main pair. Until Gate 4c and Gate 4d
# pass, the app keeps the cohesin simulator's current label.
COHESIN_HIC = {"pairs": ["in situ", "intact"], "main_pair": "in situ", "min_regions": 3,
               "regions": [("chr2", 216_000_000), ("chr5", 140_000_000), ("chr7", 130_000_000),
                           ("chr10", 100_000_000), ("chr12", 52_000_000), ("chr17", 48_000_000)],
               "res": 10_000, "window": 2_000_000, "min_reads": 10, "min_sep_bins": 2, "boot": 200, "r_c_nm": 150.0}

# Gate 2e (Phase A3, a per-pair reliability score on data no reliability test has touched;
# validation/reliability_v2.py). Pre-registered after the practice comparison (results_reliability_v2_practice.json,
# run in two chunks, imaging then Hi-C) and before any test unit was scored. Candidates (input only): spread of
# the log median over 8 refits on resampled input (boot_sd), contacts around the pair (evidence), their rank
# average (combined), and Gate 2c's misfit. Chosen per input by the best worst case over the practice groups:
# imaging input = misfit (the only candidate positive on every practice group, +0.063 to +0.175; boot_sd -0.083,
# evidence -0.019, combined -0.037 at worst); Hi-C input = boot_sd (worst -0.031; evidence -0.125, combined
# -0.075, misfit -0.064). No candidate reached 0.30 on any practice group, so a pass is not expected.
# Test sets: the genome-scale su_genome and su_genome_amanitin (chromosome units, split 0, pooled pair-weighted;
# 95 % interval by resampling the units), each with imaging and with Hi-C input. Gate 2's six sets are not reused.
# Error and scoring as Gate 2c (scale-free pair error; Spearman of score vs minus error within 10 separation strata).
# Pass, separately per input type: on BOTH test sets rho >= 0.30 with the 95 % interval above 0. Only a passing input
# type may show a reliability in the probe and the 3D view, labelled with its measured standing; otherwise nothing.
RELIABILITY_V2 = {"score": {"imaging": "misfit", "hic": "boot_sd"}, "min_rho": 0.30, "boot_refits": 8, "strata": 10,
                  "split": 0, "test": ["su_genome", "su_genome_amanitin"]}

# Gate 7 (Phase B2, false-discovery control of the differential analysis; validation/diff_gate7.py). Pre-registered
# after the practice run on chr21:28-30 Mb (results_gate7_practice.json) and before any test region was read.
# Data: four independent untreated intact Hi-C experiments of HCT116 RAD21-AID cells (ENCODE, Aiden lab), split into
# condition A = ENCSR401WPL + ENCSR002OIN and condition B = ENCSR579TBL + ENCSR697MNL: no biological change between
# A and B. Spike-ins: in each region 100 pixels (3-100 bins apart, mean >= 20 reads) get their counts in both B
# replicates multiplied by 2 or by 4 (3 random draws each); every other significant pixel is a false discovery.
# Settings (chosen on practice: the most recall with mean FDP <= 0.05): min_count 5, distance normalisation on,
# FDR 0.05, pairs up to 2 Mb at 10 kb. Practice: mean FDP 0.036, recall 0.16 (x2) and 0.94 (x4), no discovery in
# the no-change comparison.
# Pass: mean false-discovery proportion over all spike-in runs of the five test regions <= 0.05 (the stated level),
# reported with its 95 % interval (regions resampled), recall at x2 and x4, and the discoveries in each no-change
# comparison. diffHic and multiHiCcompare (R / Bioconductor) and CHESS are not run: R is not installed here, and
# installing CHESS (pip install chess-hic) fails building pysam (a FAN-C dependency) on Windows; recorded in the
# result file.
DIFF_GATE7 = {"fdr": 0.05, "setting": {"min_count": 5.0, "distance_normalise": True}, "spike": 100, "folds": [2.0, 4.0],
              "seeds": [0, 1, 2], "res": 10_000, "window": 2_000_000,
              "regions": [("chr3", 180_000_000), ("chr8", 126_000_000), ("chr11", 65_000_000), ("chr14", 90_000_000),
                          ("chr19", 12_000_000)],
              "not_run": {"diffHic": "R / Bioconductor not installed on this machine",
                          "multiHiCcompare": "R / Bioconductor not installed on this machine",
                          "CHESS": "pip install chess-hic fails: its dependency pysam (via FAN-C) does not build on Windows"}}

# Gate 6 (Phase B3, loop calls against reference calls on held-out cell lines; validation/loops_gate6.py). Pre-registered
# after the practice run on GM12878 (results_gate6_practice.json) and before any test window was read.
# Reference: the HiCCUPS loops ENCODE called on the same maps ("hic-loop-calling-step", GRCh38, portal MD5). So the
# gate measures agreement with HiCCUPS, the field's standard caller, not biological truth.
# ChronoCell settings chosen on practice (best F1 of 6): FDR 0.2 per neighbourhood with KR balancing (practice
# precision 0.889, recall 0.533, F1 0.667 over three 10 Mb GM12878 windows). On one practice window chromosight
# (F1 0.49) and Mustache (F1 0.44) were run only to check that they work here (Mustache needs NumPy < 2, pinned in
# the isolated tools environment, and reads the .mcool because hic-straw does not build).
# Test: K562 (ENCSR545YBD: map ENCFF616PUW, loops ENCFF693XIL) and IMR-90 (ENCSR852KQC: ENCFF188SSH, ENCFF527JOL),
# windows chr4:100-110, chr7:100-110, chr11:60-70 Mb at 10 kb; loops 30 kb-2 Mb apart inside the windows; a call
# matches a reference loop when both anchors are within 25 kb (one-to-one). Comparators on the same windows:
# chromosight 1.6.3 and Mustache 1.3.3 with their defaults (ICE-balanced .mcool). Not run: HiCCUPS (Java, Juicer
# tools; it is the reference), TopDom (R), Arrowhead (Java), dcHiC (R), cooltools (no Windows build).
# Pass: on BOTH test cell lines, ChronoCell's F1 is at least the higher F1 of chromosight and Mustache (precision and
# recall reported alongside).
LOOPS_GATE6 = {"choice": {"fdr": 0.2, "balance": "kr"}, "tolerance_bp": 25_000, "min_sep": 30_000, "max_sep": 2_000_000,
               "res": 10_000, "window": 10_000_000,
               "test": {"k562": [("chr4", 100_000_000), ("chr7", 100_000_000), ("chr11", 60_000_000)],
                        "imr90": [("chr4", 100_000_000), ("chr7", 100_000_000), ("chr11", 60_000_000)]},
               "tools": {"chromosight": "1.6.3", "mustache": "1.3.3"}}

# Gate 4d (Phase B1, structural-variant effects on new events with Hi-C before AND after, same cell type).
# Pre-registered before any variant-state Hi-C was read. Rule: on at least 3 events, the variant engine's predicted
# change (sv_engine, the event's joins, a population fitted to the "before" map) agrees with the measured change
# (after / before) better than BOTH the distance-shift and the no-change baselines, with the 95 % interval of each
# difference above 0 (block bootstrap over loci, as Gate 4b). Status: BLOCKED - openly downloadable,
# checksum-published before / after Hi-C with stated breakpoints was found for two events only (ENCODE RPE-1 DXZ4
# deletion on Xa, ENCSR624AJS, and on Xi, ENCSR566XIP; chrX:114,946,736-115,094,976 hg19); the 4DN Treg CRISPR
# deletions (4DNESKIGI9SP, 4DNESHQUSH4D, control 4DNES4ZVGG33) state no breakpoints, and GEO maps publish no
# checksums. Until a third event is available and the test passes, the engine stays "mechanism simulator, not
# validated".
SV_GATE4D = {"min_events": 3, "baselines": ["distance_shift", "no_change"], "status": "blocked",
             "candidates": {"ENCSR624AJS": "RPE-1 DXZ4 deletion, active X", "ENCSR566XIP": "RPE-1 DXZ4 deletion, inactive X"},
             "reason": "two events with stated breakpoints and checksum-published before / after Hi-C; three are needed"}


# Gate Q (quantum lab, chronocell/quantum; validation/quantum_gateq.py). Everything "quantum" runs on the statevector
# SIMULATOR of this computer (checked against Qiskit 2.2.1 on every circuit type: results_quantum_crosscheck.json).
# Pre-registered after the practice runs (results_gateq_practice_{tad,qaoa,qsvm,chem}.json) and before any test
# window, test reference file or test gene was read.
# Q1/Q2 data: the ENCODE GRCh38 maps of Gate 6 read by region (10 Mb at 10 kb, coarsened); practice GM12878
# chr1/2/3:100-110 Mb; test K562 and IMR-90 chr4:100-110, chr7:100-110, chr11:60-70 Mb, tiled into 20-bin windows
# (19 qubits each) with flanks. Reference: ENCODE Arrowhead "contact domains" on the same maps (preferred default;
# GM12878 ENCFF531LSJ, K562 ENCFF271SAF, IMR-90 ENCFF166QGX; portal MD5); boundaries = domain starts and ends;
# a call matches within one bin (one-to-one); only boundaries strictly inside each window are scored, for every method.
# Domain QUBO chosen on practice (exact optimum, 360 settings): 40 kb bins, minimum domain 3 bins, gamma 1.5, boundary
# cost 1.5, difference weights (practice F1 0.400). Classical callers given the same chance (best of their grids on the
# same practice windows): insulation w 4, depth 0.25 (F1 0.395); TopDom-like window 3 (F1 0.303).
# QAOA chosen on practice (4 settings, 33 windows, 4,096 shots, COBYLA 80 iterations per depth, INTERP growth): every
# setting found the optimum in 33/33 windows; the cheapest with the highest mean P(optimum) at the lowest depth:
# p = 3, expectation objective (mean P(optimum) 0.119 vs 2e-6 for a uniform guess).
# Q1 pass: on the test windows (pooled, both cell lines) the best of the 4,096 noiseless shots has the minimum energy
# (exact enumeration of all 2^19 states) in >= 90 % of windows.
# Q2 pass: on BOTH test cell lines, F1 of QAOA's best-shot boundaries >= max(F1 insulation, F1 TopDom-like) - 0.05.
# Q3 pass: UCCSD-VQE energy within 1.6 mHa of FCI (exact diagonalisation, same Hamiltonian) at every listed geometry of
# H2 and HeH+ (STO-3G from scratch; practice reproduced Szabo & Ostlund's E_HF -1.1167 / E_FCI -1.1373 for H2 at
# 1.4 bohr and E_HF -2.86066 for HeH+).
# Q4 data: protein-coding genes with TSS > 500 kb inside the windows; features from the map (compartment eigenvector,
# insulation, coverage, distance to an insulation boundary, gene density); label GTEx v10 median TPM >= 1 in the
# matching cell type. Train on all 106 practice GM12878 genes (EBV lymphocytes), test on IMR-90 genes (cultured
# fibroblasts); K562 has no GTEx match. Settings chosen by 5-fold CV on practice (QSVM 0.806, RBF 0.796, logistic 0.762).
# Q4 pass: test AUC of the quantum-kernel SVM >= RBF-SVM AUC - 0.03, and its 95 % interval (1,000 gene resamples)
# above 0.5. Reported, not gated: the noise model's effect, simulated / quantum-inspired annealing, the QSVM with a
# kernel estimated from 1,000 shots per entry, timings (no other heavy job runs during the test).
QUANTUM_GATEQ = {"tad": {"res": 40_000, "min_size": 3, "gamma": 1.5, "boundary_cost": 1.5, "weight": "difference"},
                 "qaoa": {"p": 3, "objective": "expectation"}, "qaoa_fixed": {"shots": 4096, "maxiter": 80},
                 "noise": {"p1": 1e-3, "p2": 1e-2, "readout": 2e-2},
                 "classical": {"insulation": [4, 0.25], "topdom": [3]},
                 "q1_min_hit_rate": 0.9, "q2_margin": 0.05,
                 "geometries": {"H2": [0.5, 0.7414, 1.0, 1.5, 2.0, 2.5], "HeH+": [0.5, 0.7743, 1.0, 1.5, 2.0]},
                 "qsvm": {"bandwidth": 0.05, "reps": 1, "C": 100.0}, "rbf": {"gamma": 0.02, "C": 100.0},
                 "q4_margin": 0.03, "qsvm_shots": 1000}

# Gate 8 (the Drug lab's mechanism simulator against chromatin tracing after real drug treatment;
# validation/drug_gate8.py). Pre-registered after the practice run (results_gate8_practice.json: untreated traces split
# in halves, and alpha-amanitin) and before any test (drug-treated) trace was read. Data: MINA tracing of chrX:77.66-78.50
# Mb (28 loci of 30 kb) in IMR-90 (Cheng et al., Genome Biol 2021; 4DN, MD5). Prediction from the untreated traces only
# (MDS of the median distance matrix, IMR-90 H3K27ac peaks as the signal, the class at full dose, efficacy 0.8,
# mechanism-only mode). Score: Spearman rho between predicted and measured centred log pair-distance changes; 95 %
# interval from 200 trace bootstraps; permutation p from 200 shuffles of the signal track (does WHERE the drug acts
# matter?). Practice: untreated half-vs-half gave |rho| up to 0.14 with intervals spanning 0 (noise); alpha-amanitin gave
# rho 0.29 on Xa, but shuffled targeting did as well (permutation p 0.86): a generic compaction pattern, not targeting.
# Hence: per drug, on the active X (Xa), pass if rho >= 0.20, its 95 % interval is above 0 AND permutation p <= 0.05.
# Gate pass: at least 3 of the 4 test drugs (GSK126 -> ezh2, TSA + NaBu -> hdac, 5-aza-dC -> dnmt, DMOG -> kdm). The
# classes' targets and directions were fixed from the literature before any treated data were read (the core four are
# the original ones). Xi is reported, not gated. Until it passes, the Drug lab keeps "mechanism simulator, not validated".
DRUG_GATE8 = {"allele": "Xa", "min_rho": 0.20, "max_perm_p": 0.05, "min_drugs": 3, "efficacy": 0.8, "min_detected": 0.7,
              "bootstrap": 200, "permutations": 200,
              "test_drugs": ["GSK126 (EZH2 inhibitor)", "TSA + NaBu (HDAC inhibitors)", "5-aza-dC (DNMT inhibitor)",
                             "DMOG (demethylase blocker)"]}

# Gates Q5-Q7 (quantum drug tabs; validation/quantum_drug_gates.py). Pre-registered after the practice runs
# (results_qdrug_practice_{q5,q6,q7}.json) and before any test compound, reference file or test complex was read.
# Q5 (hERG): practice TDC hERG (655 rows), 5-fold CV over three feature sets and matched grids for both kernels (the grids
# were widened three times on practice when a best setting sat at an edge; the quantum kernel's C stays at the top of its
# grid, where at its small bandwidth C and bandwidth trade off and the AUC changes in the third decimal). Chosen: the
# 8 hERG-relevant descriptors; QSVM bandwidth 0.005, 2 repetitions, C 1000 (CV AUC 0.862); RBF-SVM gamma 0.002, C 100
# (0.860); logistic regression 0.855. Test: train on all of TDC hERG, score TDC hERG_Karim minus compounds whose RDKit
# canonical SMILES occur in training. Pass: QSVM AUC >= RBF AUC - 0.03 and its 95 % interval (1,000 resamples) above 0.5.
# Q6 (molecules): practice reproduced Szabo & Ostlund's seven STO-3G HF energies within their printed precision and
# OpenFermion's H2 HF / FCI energies within 3e-4 mEh. Test: OpenFermion's LiH (STO-3G, 1.45 A) HF and FCI within 0.1 mEh,
# AND active-space UCCSD-VQE within chemical accuracy (1.6 mHa) of the active-space FCI for six stretched cases never
# run before.
# Q7 (docking): settings chosen on 39 usable practice complexes (1 in 5 by a hash of the PDB id): directional hot-spots,
# tau 1.0 A, 8 polar + 4 hydrophobic features, 8 hot-spots per type, 20-vertex core subgraph, 30 cliques, best score;
# practice docked 15 % (QAOA route) vs 3 % (random search, 1,000 placements). QAOA p 5 with CVaR found the maximum-weight
# clique in 97 % of practice graphs. Pass: on the test complexes QAOA finds the maximum-weight clique in >= 80 % of
# graphs AND the QAOA route docks (RMSD <= 2 A) at least as many as random search.
QUANTUM_DRUG_GATES = {
    "q5": {"features": "herg8", "qsvm": {"bandwidth": 0.005, "reps": 2, "C": 1000.0}, "rbf": {"gamma": 0.002, "C": 100.0},
           "margin": 0.03},
    "q6": {"openfermion_test": ["H1-Li1_sto-3g_singlet_1.45.hdf5"], "reference_tol_mEh": 0.1,
           "vqe_cases": [["H2O", 1.5, [4, 4]], ["NH3", 1.5, [6, 5]], ["N2", 1.5, [6, 6]], ["HF", 2.0, [2, 2]],
                         ["CH2O", 1.3, [4, 4]], ["HCN", 1.3, [4, 4]]]},
    "q7": {"settings": {"tau": 1.0, "max_polar": 8, "max_hydrophobic": 4, "per_type": 8, "qubits": 20, "cliques": 30,
                        "random_poses": 1000},
           "qaoa": {"p": 5, "objective": "cvar"}, "min_hit_rate": 0.8},
}

# Round 2 of the quantum gates (validation/quantum_round2.py; October 2026): new methods for the gates that failed,
# each developed on data already seen (the failed tests' data are practice now), pre-registered here part by part
# BEFORE its test data are read, and run once. The original results (results_gateq.json, results_qdrug.json) stay.
# Q6b (molecules): ADAPT-VQE (molecules.adapt_vqe) in place of the fixed-order UCCSD circuit. Practice
# (results_round2_practice_q6b.json): the six Q6 cases and nine more; the pool of generalized singles and doubles
# reached chemical accuracy in 15 of 15 (largest error 0.19 mEh; N2 at 1.5x 0.04, HCN at 1.3x 0.00), the occupied ->
# virtual pool in 13 of 15 (it stalls at a stationary point for stretched N2, 45 and 128 mEh), fixed-order UCCSD in 11
# of 15. Chosen: generalized pool, stop when the pool gradient norm < 1e-3; two practice cases reached the 100-operator
# cap (still within 0.2 mEh), so the cap is raised to 150 so that it binds less often. Test: eleven cases never run
# (molecule, bond length x equilibrium, active electrons and orbitals). Pass: ADAPT-VQE within chemical accuracy
# (1.6 mHa) of the active-space FCI in every case. Fixed-order UCCSD on the same cases is reported, not gated.
QUANTUM_ROUND2 = {
    "q6b": {"settings": {"pool": "gsd", "grad_tol": 1e-3, "max_operators": 150}, "tolerance_mEh": 1.6,
            "cases": [["N2", 1.8, [6, 6]], ["N2", 2.5, [6, 6]], ["CO", 1.6, [6, 6]], ["CO", 2.2, [6, 6]],
                      ["HCN", 1.7, [6, 6]], ["H2O", 2.2, [8, 6]], ["NH3", 1.8, [6, 6]], ["CH2O", 1.7, [6, 6]],
                      ["HF", 2.6, [6, 4]], ["LiH", 3.0, [2, 4]], ["CH4", 1.4, [4, 4]]]},
}
# Q4b (genes; pre-registered after results_round2_practice_q4b.json, before any HMEC file was read): labels from
# ENCODE CSHL whole-cell long total RNA-seq of each cell line (GENCODE V29, mean TPM over replicates >= 1; they agree
# with Q4's GTEx labels for 89 % of GM12878 and 61 % of IMR-90 genes); Q4's five Hi-C features plus local_oe and
# insulation_slope; training on every gene of the seen windows (GM12878 106, K562 381, IMR-90 380). Leave-one-cell-
# line-out AUC on practice: QSVM 0.701 (bandwidth 0.005, 1 repetition, C 10,000), RBF-SVM 0.694 (gamma 0.0005,
# C 10,000), logistic 0.688; with Q4's five features 0.699 / 0.696 / 0.691. The grids were widened once (both kernels
# alike); the best settings stay at the edge where C and the kernel width trade off and the AUC moves in the third
# decimal. Test: HMEC (ENCODE in situ Hi-C ENCFF943JRY; RNA-seq ENCFF798WGM), genes of six 10 Mb regions (chr1, chr2,
# chr3:100-110 Mb, chr4, chr7:100-110 Mb, chr11:60-70 Mb). Pass: QSVM AUC >= RBF AUC - 0.03 and its 95 % interval
# (1,000 gene resamples) above 0.5 (Q4's rule). Reported: logistic regression, the QSVM with 1,000-shot kernel entries.
QUANTUM_ROUND2["q4b"] = {"features": "q4b", "qsvm": {"bandwidth": 0.005, "reps": 1, "C": 10000.0},
                         "rbf": {"gamma": 0.0005, "C": 10000.0}, "margin": 0.03, "qsvm_shots": 1000, "test_cell": "hmec",
                         "test_regions": [["chr1", 100_000_000], ["chr2", 100_000_000], ["chr3", 100_000_000],
                                          ["chr4", 100_000_000], ["chr7", 100_000_000], ["chr11", 60_000_000]]}
# Q6c (molecules, reference corrected; pre-registered after the Q6b test and its post-hoc check,
# results_round2_q6b_posthoc.json). Q6b scored ADAPT-VQE against the lowest state of the electron-number / Sz = 0
# sector; for LiH at 3x that state is a triplet, while a VQE started from a closed-shell determinant with
# spin-conserving excitations targets the lowest singlet (which it reached to 1e-4 mEh). Q6b stays a fail. Q6c: the
# same method and settings, scored against the lowest singlet of the active space (molecules.fci_singlet), on twelve
# cases never run. Pass: within chemical accuracy (1.6 mHa) in every case. Also reported: the lowest state of the sector
# and its spin, and fixed-order UCCSD.
QUANTUM_ROUND2["q6c"] = {"settings": {"pool": "gsd", "grad_tol": 1e-3, "max_operators": 150}, "tolerance_mEh": 1.6,
                         "reference": "lowest singlet",
                         "cases": [["N2", 1.3, [6, 6]], ["N2", 3.0, [6, 6]], ["CO", 1.3, [6, 6]], ["CO", 2.8, [6, 6]],
                                   ["HCN", 2.2, [6, 6]], ["H2O", 1.8, [8, 6]], ["NH3", 2.5, [6, 6]],
                                   ["CH2O", 2.2, [6, 6]], ["HF", 3.5, [6, 4]], ["LiH", 2.5, [2, 4]], ["LiH", 4.0, [2, 4]],
                                   ["CH4", 2.0, [4, 4]]]}

# Gate 6b (loop calls; validation/loops_gate6b.py). Pre-registered after results_gate6b_practice.json and before any
# HMEC or HAP-1 loop file or window of the test regions was read. Gate 6 (fail) stays as it is. Practice: the nine
# windows already used (GM12878 practice, K562 and IMR-90 test windows of Gate 6); grid FDR x balancing x minimum donut
# enrichment x minimum cluster size (FDR widened once when the best sat at 0.01, the edge); chromosight and Mustache run
# on the same windows. Choice (largest worst-cell margin over the better tool): FDR 0.01, KR balancing, no extra
# filter: F1 0.575 / 0.506 / 0.765 (GM12878 / K562 / IMR-90) against the better tool's 0.539 / 0.486 / 0.474.
# Test: HMEC (map ENCFF943JRY, HiCCUPS loops ENCFF999UXN) and HAP-1 (ENCFF898HRO, ENCFF557WIU), three 10 Mb windows
# each that no other round-2 test reads. Pass (Gate 6's rule): on BOTH, ChronoCell's F1 >= the higher of chromosight's
# and Mustache's.
LOOPS_GATE6B = {"choice": {"fdr": 0.01, "balance": "kr", "min_oe": 0.0, "min_cluster": 1}, "tolerance_bp": 25_000,
                "min_sep": 30_000, "max_sep": 2_000_000, "res": 10_000, "window": 10_000_000,
                "test": {"hmec": [("chr5", 100_000_000), ("chr12", 60_000_000), ("chr17", 40_000_000)],
                         "hap1": [("chr5", 100_000_000), ("chr12", 60_000_000), ("chr17", 40_000_000)]},
                "tools": {"chromosight": "1.6.3", "mustache": "1.3.3"}}
# Q2b (domains; pre-registered after results_round2_practice_q2b.json, before any HMEC or HAP-1 domain file was read;
# the HMEC Hi-C counts of the same six regions are read by the Q4b test, which runs first and uses no domain call).
# Practice: the domain QUBO's exact optimum on the nine seen windows (GM12878, K562, IMR-90) for 672 settings (gamma
# widened once to 6 and 8 when the first run's best sat at 5.0, the edge; the choice did not move), the classical
# callers' grids on the same windows at the same resolution (each tuned jointly over the three cell lines). Choice
# (largest worst-cell margin): 60 kb bins, minimum domain 3 bins, gamma 5.0, boundary cost 0.25, log weights; F1
# 0.358 / 0.464 / 0.481 against the better classical caller's 0.311 / 0.393 / 0.387 (insulation window 2, depth 0.15;
# TopDom-like window 3). QAOA as Gate Q (p 3, expectation objective, 4,096 shots, 80 COBYLA iterations per depth).
# Test: HMEC and HAP-1, six 10 Mb regions each (chr1, chr2, chr3, chr4, chr7:100-110 Mb, chr11:60-70 Mb), reference
# ENCODE Arrowhead domains called on each map. Pass (Q2's rule): on BOTH, F1 of QAOA's best-shot boundaries >= the
# higher F1 of the two classical callers - 0.05. Reported: QAOA's hit rate on the windows, the exact optimum's F1.
QUANTUM_ROUND2["q2b"] = {"tad": {"res": 60_000, "min_size": 3, "gamma": 5.0, "boundary_cost": 0.25, "weight": "log"},
                         "qaoa": {"p": 3, "objective": "expectation"}, "noise": {"p1": 1e-3, "p2": 1e-2, "readout": 2e-2},
                         "classical": {"insulation": [2, 0.15], "topdom": [3]}, "margin": 0.05,
                         "test_cells": ["hmec", "hap1"],
                         "test_regions": [["chr1", 100_000_000], ["chr2", 100_000_000], ["chr3", 100_000_000],
                                          ["chr4", 100_000_000], ["chr7", 100_000_000], ["chr11", 60_000_000]]}
