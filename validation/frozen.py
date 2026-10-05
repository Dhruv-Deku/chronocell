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
