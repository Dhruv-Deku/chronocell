"""v4 Pillar 7: Self-Math PDB State Evaluator (chronocell.analytics.pdb_evaluator): every metric checked
against geometry whose answer is known, PDB unit handling, and the rule-based state on structures built
to have (or lack) the stated features."""

from __future__ import annotations

import numpy as np
import pytest

from chronocell import formats, genome, physics, synthetic
from chronocell.analytics import pdb_evaluator as PE


@pytest.fixture(scope="module")
def globule():
    return synthetic.fractal_globule(800, b0=50.0, seed=3)


def test_gyration_tensor_asphericity_and_acylindricity():
    rod = np.c_[np.arange(200.0), np.zeros(200), np.zeros(200)]
    rg, lam = PE.gyration(rod)
    assert rg == pytest.approx(physics.radius_of_gyration(rod))
    assert PE.asphericity(lam) == pytest.approx(1.0) and lam[1] - lam[2] == pytest.approx(0.0, abs=1e-9)
    ball = np.random.default_rng(0).normal(size=(50_000, 3))
    assert PE.asphericity(PE.gyration(ball)[1]) < 0.01
    disc = np.random.default_rng(1).normal(size=(50_000, 3)) * [1.0, 1.0, 0.0]
    _, lam_d = PE.gyration(disc)
    assert lam_d[1] - lam_d[2] == pytest.approx(lam_d[1], rel=1e-9)          # flat: acylindricity = l2
    # the specified Delta is the relative shape anisotropy kappa^2 of physics.gyration_shape
    x = np.random.default_rng(2).normal(size=(300, 3)) * [3.0, 1.5, 0.7]
    assert PE.asphericity(PE.gyration(x)[1]) == pytest.approx(physics.gyration_shape(x).anisotropy, rel=1e-9)


def test_hull_volume_and_packing_density():
    cube = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float) * 1000.0   # 1 um cube (nm)
    assert PE.hull_volume(cube) == pytest.approx(1e9)
    pts = np.vstack([cube, np.random.default_rng(3).random((92, 3)) * 1000.0])                   # 100 points inside
    e = PE.evaluate(pts)
    assert e.hull_volume_um3 == pytest.approx(1.0) and e.packing_density_per_um3 == pytest.approx(100.0)
    assert np.isnan(PE.hull_volume(np.zeros((10, 3))))                                           # degenerate


def test_distance_decay_exponent_recovers_known_scaling(globule):
    e = PE.evaluate(globule)
    assert e.gamma_d == pytest.approx(1 / 3, abs=0.05)                  # fractal globule: nu = 1/3
    assert e.gamma_c == pytest.approx(PE.ALPHA * e.gamma_d)
    rng = np.random.default_rng(4)
    g = [PE.evaluate(np.cumsum(rng.normal(size=(1500, 3)), axis=0) * 40.0).gamma_d for _ in range(8)]
    assert np.mean(g) == pytest.approx(0.5, abs=0.06)                   # ideal random walk: nu = 1/2
    straight = np.c_[np.arange(300.0) * 50, np.zeros(300), np.zeros(300)]
    assert PE.evaluate(straight).gamma_d == pytest.approx(1.0, abs=1e-6)


def test_states_follow_the_rules_on_constructed_structures(globule):
    normal = PE.evaluate(globule)
    assert normal.state == "Normal" and 0.75 <= normal.gamma_c <= 1.10
    compact = synthetic.relax(globule.mean(0) + (globule - globule.mean(0)) * 0.75, 50.0, 40.0, iters=60)
    sen = PE.evaluate(compact)
    assert sen.state == "Senescent"
    flags = {c.name: c.met for c in sen.criteria if c.state == "Senescent"}
    assert sum(flags.values()) >= 2
    # every criterion reports its value, threshold and where the threshold comes from
    for c in sen.criteria:
        assert c.source in ("specification", "assumption") and c.threshold
    assert PE.evaluate(compact).summary() == sen.summary()             # deterministic: same input, same output


def test_thresholds_are_parameters_not_hidden_constants(globule):
    strict = PE.Thresholds(gamma_normal=(1.2, 1.3))
    assert PE.evaluate(globule, strict).state != "Normal"


def test_pdb_units_and_models(globule):
    ch = genome.chrom("chr22")
    gc = np.full(ch.n_bins, 0.4)
    text, _ = formats.write_pdb(globule[:300], 100, gc, gc, 1.0, source="t", method="t", chrom=ch)
    models, note = PE.coords_from_pdb_text(text)
    assert "declared" in note and len(models) == 1 and np.allclose(models[0], globule[:300], atol=2e-3)
    angstrom = "\n".join(f"ATOM  {k + 1:5d}  CA  ALA A{k + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
                         for k, (x, y, z) in enumerate(globule[:50] * 10.0 / 10.0))
    models, note = PE.coords_from_pdb_text(angstrom)
    assert "Angstrom" in note and np.allclose(models[0], globule[:50] / 10.0, atol=1e-3)
    frames = np.stack([globule[:120], globule[:120] * 1.1])
    traj = formats.write_pdb_trajectory(frames, np.arange(1, 121), np.full(120, "CH22"), ch, "t", "t")
    models, _ = PE.coords_from_pdb_text(traj)
    assert len(models) == 2 and np.allclose(models[1], frames[1], atol=2e-3)
    with pytest.raises(ValueError):
        PE.coords_from_pdb_text("HEADER nothing\nEND\n")


def test_population_evaluation_summarises_members(globule):
    rng = np.random.default_rng(5)
    members = np.stack([globule[:200] + rng.normal(scale=5.0, size=(200, 3)) for _ in range(6)])
    out = PE.evaluate_population(members)
    assert out["members"] == 6 and sum(out["states"].values()) == 6
    assert out["rg_nm"]["p10"] <= out["rg_nm"]["median"] <= out["rg_nm"]["p90"]


def test_bad_input_is_rejected():
    with pytest.raises(ValueError):
        PE.evaluate(np.zeros((5, 3)))
    with pytest.raises(ValueError):
        PE.evaluate(np.full((20, 3), np.nan))
