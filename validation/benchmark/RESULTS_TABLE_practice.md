# Benchmark: practice datasets (tuning allowed; in-sample)

Generated 2026-10-04T07:44:08+00:00 by `python -m validation.benchmark.run --practice` in 879 s. Every number is computed by the harness from the data; nothing is typed in.

Columns: % of ceiling = trend-removed Spearman / ceiling (distance-pattern recovery beyond the genomic-separation trend); raw ρ = Spearman; CCC = Lin's concordance in nm (absolute size); coverage = share of held-out single-copy distances inside the stated 50 / 80 / 90 % intervals (raw → recalibrated). Brackets: 95 % interval from resampling half B's copies (split 0).

## bintu_k562_28_30 · input: imaging · all pairs

| Method | % of ceiling | raw ρ | CCC (nm) | size ratio | coverage 50/80/90 % |
|---|---|---|---|---|---|
| v3_3_windowed | 92.3 [91.7, 93.0] | 0.980 [0.977, 0.981] | 0.924 [0.914, 0.930] | 0.90 | 40/68/79 → 47/77/88 |
| v4_whole | 92.3 [91.7, 93.0] | 0.980 [0.977, 0.981] | 0.924 [0.914, 0.930] | 0.90 | 40/68/79 → 47/77/88 |
| no_3d | 88.3 [87.8, 88.9] | 0.970 [0.967, 0.970] | 0.915 [0.907, 0.920] | 0.90 | 40/68/79 → 47/77/88 |
| pastis_pm2 | 59.2 [57.8, 60.7] | 0.887 [0.882, 0.891] | 0.248 [0.242, 0.255] | 1.77 |  |
| pastis_mds | 45.3 [43.8, 46.7] | 0.834 [0.829, 0.839] | 0.564 [0.552, 0.576] | 1.22 |  |
| v3_2_single | 35.2 [33.7, 36.4] | 0.798 [0.792, 0.801] | 0.111 [0.110, 0.113] | 0.32 |  |
| genomic_distance_only | 1.1 [0.8, 1.5] | 0.928 [0.923, 0.930] | 0.916 [0.910, 0.918] | 0.98 |  |

## bintu_k562_28_30 · input: hic · all pairs

| Method | % of ceiling | raw ρ | CCC (nm) | size ratio | coverage 50/80/90 % |
|---|---|---|---|---|---|
| v3_3_windowed | 53.7 [51.7, 55.7] | 0.940 [0.936, 0.942] | 0.244 [0.242, 0.247] | 0.45 | 14/26/33 → 19/37/50 |
| v4_whole | 53.7 [51.7, 55.7] | 0.940 [0.936, 0.942] | 0.244 [0.242, 0.247] | 0.45 | 14/26/33 → 19/37/50 |
| no_3d | 51.7 [49.8, 53.4] | 0.925 [0.921, 0.926] | 0.279 [0.277, 0.282] | 0.44 | 15/28/35 → 20/39/51 |
| pastis_mds | 46.5 [45.1, 48.2] | 0.911 [0.907, 0.913] | 0.392 [0.387, 0.396] | 0.53 |  |
| pastis_pm2 | 44.5 [42.5, 46.5] | 0.918 [0.914, 0.920] | 0.724 [0.720, 0.729] | 0.72 |  |
| v3_2_single | 43.6 [42.4, 45.0] | 0.898 [0.895, 0.900] | 0.215 [0.213, 0.217] | 0.39 |  |

## Where we lose

- bintu_k562_28_30 · imaging · all_pairs · v3_3_windowed coverage at 50 %: 40.3 %
- bintu_k562_28_30 · imaging · all_pairs · v3_3_windowed coverage at 80 %: 67.6 %
- bintu_k562_28_30 · imaging · all_pairs · v3_3_windowed coverage at 90 %: 78.6 %
- bintu_k562_28_30 · imaging · all_pairs · v4_whole coverage at 50 %: 40.3 %
- bintu_k562_28_30 · imaging · all_pairs · v4_whole coverage at 80 %: 67.6 %
- bintu_k562_28_30 · imaging · all_pairs · v4_whole coverage at 90 %: 78.6 %
- bintu_k562_28_30 · hic · all_pairs · lin_ccc_nm: no_3d 0.279 > v3_3_windowed 0.244
- bintu_k562_28_30 · hic · all_pairs · lin_ccc_nm: pastis_mds 0.392 > v3_3_windowed 0.244
- bintu_k562_28_30 · hic · all_pairs · lin_ccc_nm: pastis_pm2 0.724 > v3_3_windowed 0.244
- bintu_k562_28_30 · hic · all_pairs · v3_3_windowed coverage at 50 %: 13.8 %
- bintu_k562_28_30 · hic · all_pairs · v3_3_windowed coverage at 80 %: 26.0 %
- bintu_k562_28_30 · hic · all_pairs · v3_3_windowed coverage at 90 %: 33.2 %
- bintu_k562_28_30 · hic · all_pairs · v3_3_windowed coverage_recalibrated at 50 %: 18.5 %
- bintu_k562_28_30 · hic · all_pairs · v3_3_windowed coverage_recalibrated at 80 %: 37.4 %
- bintu_k562_28_30 · hic · all_pairs · v3_3_windowed coverage_recalibrated at 90 %: 49.9 %
- bintu_k562_28_30 · hic · all_pairs · v4_whole coverage at 50 %: 13.8 %
- bintu_k562_28_30 · hic · all_pairs · v4_whole coverage at 80 %: 26.0 %
- bintu_k562_28_30 · hic · all_pairs · v4_whole coverage at 90 %: 33.2 %
- bintu_k562_28_30 · hic · all_pairs · v4_whole coverage_recalibrated at 50 %: 18.5 %
- bintu_k562_28_30 · hic · all_pairs · v4_whole coverage_recalibrated at 80 %: 37.4 %
- bintu_k562_28_30 · hic · all_pairs · v4_whole coverage_recalibrated at 90 %: 49.9 %

## Methods not run

- **ShRec3D**: Published as MATLAB code (Lesne et al., Nat Methods 2014); no installable implementation here. ChronoCell's v3.2 shortest-path MDS start is ShRec3D-style but is not the published tool.
- **Chrom3D**: C++ program (Paulsen et al., Genome Biol 2017) needing Boost to build and lamina-association data; no compiler toolchain or LAD data here.
- **3DMax / LorDG**: Java programs (Oluwadare et al.); no Java runtime on this machine.
