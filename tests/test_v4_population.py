"""v4 Pillar 1 and 2: the scalable population model (chronocell.population) and its per-pair uncertainty."""

from __future__ import annotations

import numpy as np
import pytest
import torch

from chronocell import ensemble as E, population as P
from tests.test_v33 import _population


def _random_problem(n=37, r=9, seed=0):
    rng = np.random.default_rng(seed)
    A = torch.tensor(rng.normal(size=(n, r)), dtype=torch.float64)
    phi = torch.tensor(rng.normal(size=n - 1) * 0.3 - 1.0, dtype=torch.float64)
    T = rng.normal(size=(n, n)) * 0.5 + 1.0
    T = torch.tensor((T + T.T) / 2, dtype=torch.float64)
    Wn = rng.random((n, n)) * (rng.random((n, n)) < 0.8)
    Wn = (Wn + Wn.T) / 2
    np.fill_diagonal(Wn, 0.0)
    return A, phi, T, torch.tensor(Wn, dtype=torch.float64), float((Wn > 0).sum()) / 2


def _autograd_loss(A, phi, T, W, Pn):
    A = A.clone().requires_grad_(True)
    phi = phi.clone().requires_grad_(True)
    n = A.shape[0]
    cum = torch.cat([torch.zeros(1, dtype=A.dtype), torch.cumsum(torch.exp(phi), 0)])
    s = ((A[:, None] - A[None]) ** 2).sum(-1) + (cum[None] - cum[:, None]).abs()
    iu = torch.triu_indices(n, n, 1)
    s_u = s[iu[0], iu[1]].clamp_min(1e-9)
    loss = (W[iu[0], iu[1]] * (torch.log(s_u) - T[iu[0], iu[1]]) ** 2).sum() / Pn
    loss.backward()
    return loss.item(), A.grad, phi.grad


@pytest.mark.parametrize("block", [5, 16, 64])
def test_blockwise_gradients_match_autograd(block):
    A, phi, T, W, Pn = _random_problem()
    val, gA, gphi = P._loss_grad(A, phi, T, W, Pn, block)
    ref_val, ref_A, ref_phi = _autograd_loss(A, phi, T, W, Pn)
    assert val == pytest.approx(ref_val, rel=1e-10)
    assert torch.allclose(gA, ref_A, rtol=1e-8, atol=1e-10)
    assert torch.allclose(gphi, ref_phi, rtol=1e-8, atol=1e-10)


def test_full_rank_population_matches_v33_on_a_planted_ideal_chain():
    d = _population(n=14)
    freq, truth = (d < 1.0).mean(0), np.median(d, axis=0) * 150.0
    iu = np.triu_indices(14, 1)
    cfg = P.PopulationConfig(iterations=600, replicas=40, frames=20, device="cpu")
    res = P.fit_population(freq, 6000, r_c_nm=150.0, cfg=cfg)
    v33 = E.fit_ensemble(freq, 6000, r_c_nm=150.0, cfg=E.EnsembleConfig(iterations=600, replicas=40, frames=20))
    assert res.config["rank"] == 14 and not res.config["local_term"]
    assert np.corrcoef(np.log(res.median_distance_nm[iu]), np.log(truth[iu]))[0, 1] > 0.98
    assert np.allclose(res.median_distance_nm[iu], v33.median_distance_nm[iu], rtol=0.05)
    assert res.trajectories_nm.shape == (20, 40, 14, 3) and res.representative_nm.shape == (14, 3)
    assert res.spread_cv == pytest.approx(0.42, abs=0.05)


def test_truncated_rank_with_random_walk_fits_a_loop_population():
    d = _population(n=60, loop=(10, 45), seed=4)
    freq, truth = (d < 1.0).mean(0), np.median(d, axis=0) * 150.0
    iu = np.triu_indices(60, 1)
    res = P.fit_population(freq, 6000, 150.0, P.PopulationConfig(iterations=800, rank_cap=12, replicas=30, frames=10))
    assert res.config["rank"] == 12 and res.config["local_term"]
    assert res.history["best_loss"][0] < 0.1 * res.history["loss"][0]
    assert E._spearman(res.median_distance_nm[iu], truth[iu]) > 0.95
    assert res.median_distance_nm[10, 45] < res.median_distance_nm[10, 30]


def test_coarse_to_fine_start_converges():
    d = _population(n=90, seed=5)
    freq = (d < 1.0).mean(0)
    cfg = P.PopulationConfig(iterations=300, init="coarse", coarse_beads=30, coarse_iterations=200, rank_cap=40,
                             replicas=20, frames=5)
    res = P.fit_population(freq, 6000, 150.0, cfg)
    assert res.config["init"] == "coarse" and "coarse_loss" in res.history
    assert res.contact_fit > 0.97


def test_counts_entry_point_anchors_adjacent_pairs_to_b0():
    d = _population(n=40, cells=4000, seed=6)
    f = (d < 1.0).mean(0)
    i, j = np.triu_indices(40, 1)
    cnt = np.random.default_rng(0).poisson(f[i, j] * 400).astype(float)
    keep = cnt > 0
    res = P.fit_population_from_counts(i[keep], j[keep], cnt[keep], 40, b0_nm=50.0,
                                       cfg=P.PopulationConfig(iterations=400, replicas=20, frames=5))
    adj = np.diag(res.median_distance_nm, 1)
    assert np.median(adj) == pytest.approx(50.0, rel=0.15)
    assert res.config["p_adjacent_assumed"] == 0.5 and res.config["input"] == "sequencing counts"


# ---------------------------------------------------------------- uncertainty (Pillar 2)
def test_maxwell_moments_and_quantiles_match_sampling():
    rng = np.random.default_rng(0)
    r = np.linalg.norm(rng.normal(size=(400_000, 3)), axis=1)
    assert r.mean() == pytest.approx(P.MAXWELL_MEAN, rel=5e-3)
    assert r.std() == pytest.approx(P.MAXWELL_SD, rel=5e-3)
    assert np.median(r) == pytest.approx(P.MAXWELL_MEDIAN, rel=5e-3)
    for level in (0.5, 0.8, 0.9):
        lo, hi = P.central_interval(level)
        assert np.mean((r > lo) & (r < hi)) == pytest.approx(level, abs=3e-3)
    assert P.maxwell_cdf(np.array([P.MAXWELL_MEDIAN]))[0] == pytest.approx(0.5, abs=1e-5)


def test_pair_summary_is_consistent_for_v33_and_v4_results():
    d = _population(n=14, seed=7)
    freq = (d < 1.0).mean(0)
    v33 = E.fit_ensemble(freq, 6000, cfg=E.EnsembleConfig(iterations=300, replicas=20, frames=5))
    v4 = P.fit_population(freq, 6000, cfg=P.PopulationConfig(iterations=300, replicas=20, frames=5))
    for res in (v33, v4):
        s = P.pair_summary(res, 2, 9, level=0.9)
        assert s.lower < s.median < s.upper and s.sd > 0
        assert s.median == pytest.approx(float(res.median_distance_nm[2, 9]), rel=1e-4)
        assert s.mean == pytest.approx(s.median * P.MAXWELL_MEAN / P.MAXWELL_MEDIAN, rel=1e-9)
    # members drawn from the v4 model reproduce the stated 90 % interval
    rng = np.random.default_rng(1)
    x = v4.model.sample(20_000, rng)
    dd = np.linalg.norm(x[:, 2] - x[:, 9], axis=1)
    s = P.pair_summary(v4, 2, 9, level=0.9)
    assert np.mean((dd > s.lower) & (dd < s.upper)) == pytest.approx(0.9, abs=0.01)


def test_rmsf_and_reliability_are_well_defined():
    rng = np.random.default_rng(2)
    base = np.cumsum(rng.normal(size=(30, 3)), axis=0) * 50
    members = base[None] + rng.normal(size=(50, 30, 3)) * np.linspace(1, 20, 30)[None, :, None]
    rot = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    rot *= np.sign(np.linalg.det(rot))
    members[::2] = members[::2] @ rot.T                      # rigid rotations must not count as spread
    f = P.rmsf_nm(members, base)
    assert f.shape == (30,) and np.corrcoef(f, np.linspace(1, 20, 30))[0, 1] > 0.9
    d = _population(n=14, seed=8)
    freq = (d < 1.0).mean(0)
    freq[3, :] = freq[:, 3] = np.nan                         # bead 3 never observed
    res = P.fit_population(freq, 6000, cfg=P.PopulationConfig(iterations=300, replicas=10, frames=3))
    rel = P.bead_reliability(res, freq, 6000)
    assert rel.shape == (14,) and np.all((rel >= 0) & (rel <= 1)) and rel[3] == 0.0 and rel.max() > 0.5


def test_excluded_volume_members_report_the_change_they_make():
    d = _population(n=40, seed=9)
    res = P.fit_population((d < 1.0).mean(0), 6000, 150.0, P.PopulationConfig(iterations=300, replicas=8, frames=2))
    relaxed, info = P.excluded_volume_members(res.frames_nm, b0_nm=60.0)
    assert relaxed.shape == res.frames_nm.shape
    assert info["overlaps_after"] <= info["overlaps_before"]
    assert 0.0 <= info["median_relative_pair_distance_change"] < 5.0
    assert 0.85 < info["median_bond_over_b0_after"] < 1.15


def test_app_defaults_match_the_frozen_validation_settings():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import frozen
    assert P.WHOLE_CHROMOSOME_DEFAULTS == frozen.WHOLE_CHROMOSOME
    assert P.config_for(5000).rank_cap == frozen.WHOLE_CHROMOSOME["rank_cap"]


def test_recalibrated_intervals_match_the_validation_protocol():
    import json
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import protocol as PR
    cal = json.loads(P.CALIBRATION_PATH.read_text(encoding="utf-8"))
    pits = PR.isotonic_recalibration(np.asarray(cal["pit_histogram"]))
    for level in (0.5, 0.8, 0.9):
        mine = P.recalibrated_interval(level)
        ref = tuple(float(P.maxwell_quantile(q)) for q in pits(level))
        assert np.allclose(mine, ref)
        raw = P.central_interval(level)
        assert mine[1] > raw[1]                    # the practice data have a heavier upper tail than Maxwell
    s = P.summarise_sigma(np.array([100.0]), 150.0, 0.9, recalibrated=True)
    assert s["upper"][0] == pytest.approx(100.0 * P.recalibrated_interval(0.9)[1])
