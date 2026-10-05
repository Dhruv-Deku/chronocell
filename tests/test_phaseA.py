"""Phase A (accuracy core): data loading, GPU / CPU agreement, and the new calibration and interval pieces.
Every test runs without the downloaded validation data unless it says otherwise."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

VALIDATION = Path(__file__).resolve().parent.parent / "validation"
sys.path.insert(0, str(VALIDATION))


def test_genome_scale_columns_are_found_by_name(tmp_path):
    import datasets as D
    p = tmp_path / "g.tsv"
    # the extra-column files put other columns between the ones the loader needs
    p.write_text("x(nm)\tZ(nm)\ty(nm)\tgenomic coordinate\ttranscription\thomolog number\tcell number\t"
                 "experiment number\tdistance to speckle (nm)\n"
                 "10\t1\t20\tchr1:0-100000\t1\t1\t5\t2\t3\n"
                 "11\t2\t21\tchr1:3000000-3100000\t0\t2\t5\t2\t4\n")
    df = D._read_genome_columns(p)
    assert list(df.columns) == ["z", "x", "y", "locus", "homolog", "cell", "exp"]
    assert df["z"].tolist() == [1, 2] and df["x"].tolist() == [10, 11] and df["homolog"].tolist() == [1, 2]
    assert df["exp"].tolist() == [2, 2]


def test_phase_a_datasets_have_fixed_roles():
    import datasets as D
    assert D.REGISTRY["su_genome_tx"].role == "practice"
    assert D.REGISTRY["su_genome_amanitin"].role == "test"
    for k in ("su_genome_tx", "su_genome_amanitin"):
        assert D.REGISTRY[k].paired_hic == "su_hic_genome" and D.REGISTRY[k].kind == "su_genome"


def _planted(n: int = 60, seed: int = 0):
    from chronocell import ensemble as E
    rng = np.random.default_rng(seed)
    s = np.abs(np.subtract.outer(np.arange(n), np.arange(n))).astype(float)
    sig = 0.6 * np.maximum(s, 1) ** 0.35 * np.exp(rng.normal(0, 0.08, (n, n)))
    f = E.contact_probability_from_sigma((sig + sig.T) / 2)
    np.fill_diagonal(f, np.nan)
    return f


@pytest.mark.skipif(not __import__("torch").cuda.is_available(), reason="no CUDA GPU")
def test_gpu_and_cpu_fits_agree():
    """Stated tolerance (fixed before use in any gate): per-pair median distance within 1 %, the median
    pair within 1e-4, best misfit within 1e-3 (relative), for the v3.3 and v4 models."""
    from chronocell import ensemble as E, population as P
    f = _planted()
    iu = np.triu_indices(len(f), 1)
    fits = {
        "v3.3": [E.fit_ensemble(f, 2000.0, r_c_nm=150.0, cfg=E.EnsembleConfig(seed=0, device=d)) for d in ("cpu", "cuda")],
        "v4": [P.fit_population(f, 2000.0, r_c_nm=150.0, cfg=P.PopulationConfig(seed=0, device=d, rank_cap=32))
               for d in ("cpu", "cuda")],
    }
    for name, (cpu, gpu) in fits.items():
        assert cpu.config["device_used"] == "cpu" and gpu.config["device_used"].startswith("cuda"), name
        rel = np.abs(cpu.median_distance_nm[iu] / gpu.median_distance_nm[iu] - 1)
        assert rel.max() <= 1e-2 and np.median(rel) <= 1e-4, (name, rel.max(), np.median(rel))
        a, b = cpu.history["best_loss"][0], gpu.history["best_loss"][0]
        assert abs(a - b) <= 1e-3 * abs(a), (name, a, b)


def test_auto_device_uses_the_cpu_without_cuda(monkeypatch):
    import torch
    from chronocell import ensemble as E, population as P
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    f = _planted(30)
    assert E.fit_ensemble(f, 2000.0, cfg=E.EnsembleConfig(iterations=50)).config["device_used"] == "cpu"
    assert P.fit_population(f, 2000.0, cfg=P.PopulationConfig(iterations=50, rank_cap=16)).config["device_used"] == "cpu"


def test_cuda_out_of_memory_falls_back_to_the_cpu_only_on_auto():
    import torch
    from types import SimpleNamespace
    from chronocell import ensemble as E
    calls = []

    def fit(cfg):
        calls.append(cfg.device)
        if cfg.device != "cpu":
            raise torch.cuda.OutOfMemoryError("CUDA out of memory (test)")
        return SimpleNamespace(config={"device_used": "cpu"})

    res = E.with_cpu_fallback(fit, E.EnsembleConfig(device="auto"))
    assert calls == ["auto", "cpu"] and res.config["device_fallback"].startswith("CUDA out of memory")
    calls.clear()
    with pytest.raises(torch.cuda.OutOfMemoryError):
        E.with_cpu_fallback(fit, E.EnsembleConfig(device="cuda"))          # an explicit device is never overridden
    assert calls == ["cuda"]


def test_binomial_thinning_keeps_the_expected_share_of_reads():
    import phase_a as A
    rng = np.random.default_rng(3)
    c = rng.poisson(50, (40, 40)).astype(float)
    c = np.triu(c, 1) + np.triu(c, 1).T + np.diag(np.diag(c))
    t = A.thin_counts(c, 0.25, seed=0)
    iu = np.triu_indices(40, 1)
    assert np.allclose(t, t.T) and (t <= c).all() and np.array_equal(np.diag(t), np.diag(c))
    assert abs(t[iu].sum() / c[iu].sum() - 0.25) < 0.01
    assert np.array_equal(A.thin_counts(c, 0.25, seed=0), t) and A.thin_counts(c, 1.0) is c


def test_phase_a_units_respect_dataset_roles():
    import datasets as D
    import phase_a_units as U
    for s in U.PRACTICE_HIC + U.PRACTICE_IMAGING:
        assert D.REGISTRY[s.dataset.split(":")[0]].role == "practice", s
    for s in U.TEST_HIC_MAIN + U.TEST_HIC_SECONDARY + U.TEST_IMAGING:
        assert D.REGISTRY[s.dataset.split(":")[0]].role == "test", s


def test_learned_correction_features_and_identity():
    from chronocell import learned_correction as LC
    rng = np.random.default_rng(5)
    n = 30
    s = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) * 30_000.0
    med = 100.0 * np.maximum(s / 30_000, 1) ** 0.4 * np.exp(rng.normal(0, 0.05, (n, n)))
    med = (med + med.T) / 2
    np.fill_diagonal(med, 0.0)
    f = np.clip(rng.uniform(0.01, 0.6, (n, n)), 0, 1)
    f[3, 7] = np.nan                                                    # an unobserved pair
    feats = LC.features(med, f, f * 1000, s, 150.0, True, 2000.0, 30_000.0)
    iu = np.triu_indices(n, 1)
    assert feats.shape == (n, n, len(LC.FEATURES)) and np.isnan(feats[5, 2]).all()     # lower triangle empty
    assert np.isfinite(feats[iu][:, [0, 1, 4, 5, 6, 7, 8]]).all()
    zero = {"kind": "linear", "mu": [0.0] * len(LC.FEATURES), "sd": [1.0] * len(LC.FEATURES), "intercept": 0.0,
            "beta": [0.0] * len(LC.FEATURES)}
    assert np.allclose(LC.correct(med, feats, zero)[iu], med[iu])                 # zero correction = the model
    shift = dict(zero, intercept=np.log(2.0))
    assert np.allclose(LC.correct(med, feats, shift)[iu], 2 * med[iu])


def test_learned_correction_method_without_a_frozen_correction_is_the_base_model(monkeypatch):
    sys.path.insert(0, str(VALIDATION.parent))
    from chronocell import learned_correction as LC
    from validation.benchmark import methods as M
    monkeypatch.setattr(LC, "load", lambda path=None: {"chosen": "none"})
    f = _planted(40)
    n = len(f)
    sep = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) * 30_000
    inp = M.Input("imaging", n, f, np.full((n, n), 500.0), f * 500, 150.0, np.zeros((n, n)), sep, [(0, n)], 0)
    pred = M.learned_correction(inp)
    assert pred.median_nm.shape == (n, n) and "no correction" in pred.notes[0]


def test_conformal_intervals_quantiles_coverage_and_tie_rule():
    import intervals_v2 as I
    rng = np.random.default_rng(2)
    r = rng.normal(0.0, 0.5, 200_000)
    h = np.histogram(np.clip(r, I.GRID[0], I.GRID[-1]), bins=I.GRID)[0].astype(float)
    lo, hi = I.quantiles(h, [0.05, 0.95])
    assert abs(lo + 0.5 * 1.645) < 0.01 and abs(hi - 0.5 * 1.645) < 0.01
    assert abs(I.coverage(h, lo, hi) - 0.90) < 0.002
    h2 = h.copy()
    h2[0] += h.sum() / 9                                    # zero distances (-inf, clipped to the grid floor)
    assert I.coverage(h2, lo, hi) < 0.82                    # count as misses, never dropped
    table = {"bands|hic": {"worst_abs_err_90": 0.20}, "pooled|hic": {"worst_abs_err_90": 0.175},
             "bands+A1|hic": {"worst_abs_err_90": 0.20}, "pooled+A1|hic": {"worst_abs_err_90": 0.1749},
             "bands|imaging": {"worst_abs_err_90": 0.066}, "pooled|imaging": {"worst_abs_err_90": 0.074}}
    assert I.choose(table, "hic") == "pooled|hic"           # within the tie margin: the simpler variant
    assert I.choose(table, "imaging") == "bands|imaging"    # better by more than the margin


def test_fofct_reader_pools_files_and_keeps_nm(tmp_path):
    import predictor_mouse as PM
    head = ("##FOF-CT_version=v0.1,,,,,,,,\n##XYZ_unit=micron,,,,,,,,\n\"#Software_Authors: x, y\",,,,,,,,\n"
            "##columns=(Spot_ID, Trace_ID, X, Y, Z,Chrom, Chrom_Start, Chrom_End, Cell_ID)\n")
    rows = []
    for t in (1, 2):
        for k in range(4):
            if t == 2 and k == 3:
                continue                                    # one locus missing in trace 2
            rows.append(f"{len(rows)},{t},{100.0 * k + t},{0.0},{0.0},chr6,{1001 + 30000 * k},{31001 + 30000 * k},{t}")
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    a.write_text(head + "\n".join(rows) + "\n")
    b.write_text(head + "\n".join(rows[:4]) + "\n")
    tr = PM.read_fofct([a, b])
    assert tr.xyz.shape == (3, 4, 3)                        # trace ids made unique per file
    assert tr.meta["scale_to_nm"] == 1.0 and tr.meta["locus_bp"] == 30000
    assert np.allclose(np.diff(tr.starts), 30000) and tr.chrom[0] == "chr6"
    assert np.isnan(tr.xyz[1, 3]).all() and np.isclose(tr.xyz[0, 2, 0], 201.0)
