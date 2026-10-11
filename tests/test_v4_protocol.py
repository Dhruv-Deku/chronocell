"""Validation protocol helpers (validation/protocol.py) that every held-out number depends on."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
import protocol as PR  # noqa: E402


def test_fast_nanmedian_is_exactly_numpys():
    rng = np.random.default_rng(3)
    for shape, frac in (((301, 257), 0.3), ((40, 6, 9), 0.6), ((7, 11), 1.0), ((5, 3), 0.0)):
        a = rng.normal(size=shape).astype(np.float32)
        a[rng.random(shape) < frac] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            ref = np.nanmedian(a, axis=0)
        assert np.array_equal(PR.nanmedian0(a), ref, equal_nan=True)


def test_half_stats_and_scores_on_a_planted_population():
    rng = np.random.default_rng(0)
    steps = rng.normal(size=(800, 29, 3)) * 100.0
    xyz = np.concatenate([np.zeros((800, 1, 3)), np.cumsum(steps, axis=1)], axis=1)
    xyz[rng.random(xyz.shape[:2]) < 0.1] = np.nan                     # 10 % of loci undetected
    a, b = PR.split(len(xyz), 0)
    ha, hb = PR.half_stats(xyz[a], 150.0), PR.half_stats(xyz[b], 150.0)
    assert ha.median.shape == (30, 30) and np.allclose(np.diag(ha.median), 0)
    assert np.nanmin(ha.freq) >= 0 and np.nanmax(ha.freq) <= 1
    sep = PR.separation(30)
    mask = np.triu(np.ones((30, 30), bool), 1)
    s = PR.scores(hb.median, hb.median, sep, mask)
    assert s["spearman"] > 0.999 and s["lin_ccc_nm"] > 0.999                  # a map scored against itself
