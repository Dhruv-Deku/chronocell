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
