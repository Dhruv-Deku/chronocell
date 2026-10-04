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
