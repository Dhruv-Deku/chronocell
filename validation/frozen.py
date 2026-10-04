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
