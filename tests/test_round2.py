"""Round 2 of the quantum gates (October 2026): the docking score and refinement, virtual hydrogens, the analysis
suite's loop-FDR option and the generated report blocks. Each piece is checked against exact expectations; the
defaults of every earlier function are checked unchanged."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from chronocell import analysis as AN
from chronocell.quantum import docking as DK

ROOT = Path(__file__).resolve().parent.parent


def _pocket():
    """A small synthetic 'protein': a ring of backbone-like N/C/O atoms around the origin."""
    el, xyz, nm, rs = [], [], [], []
    for k in range(12):
        t = 2 * np.pi * k / 12
        for e, r, name in (("N", 5.0, "N"), ("C", 5.6, "CA"), ("O", 6.2, "O")):
            el.append(e)
            xyz.append((r * np.cos(t), r * np.sin(t), 0.6 * (k % 3) - 0.6))
            nm.append(name)
            rs.append("GLY")
    return DK.Protein(np.array(el), np.array(xyz, float), np.array(nm), np.array(rs))


def _ligand():
    mol = """lig
  test

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7000   -1.2000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
M  END
"""
    return DK.read_sdf(mol)


def test_vina_score_terms_and_refinement_never_worsen():
    prot = DK.add_virtual_hydrogens(_pocket())
    assert (prot.el == "H").sum() == 12                       # one virtual H per backbone N (glycine)
    h = prot.xyz[prot.el == "H"]
    n = prot.xyz[(prot.el == "N")]
    assert np.allclose(np.linalg.norm(h - n, axis=1), 1.0)
    lig = _ligand()
    lt = DK.ligand_types(lig)
    assert lt.don.tolist() == [False, False, True, True] and lt.acc.tolist() == [False, False, True, False]
    assert lt.hyd.tolist() == [False, False, False, False]      # both carbons touch N or O
    pt = DK.protein_types(prot, np.zeros(3), 20.0)
    centre = lt.xyz - lt.xyz.mean(0)
    s_mid = DK.vina_score(centre, lt, pt)
    s_clash = DK.vina_score(centre + np.array([5.0, 0.0, 0.0]), lt, pt)   # on top of the ring atoms
    s_far = DK.vina_score(centre + np.array([0.0, 0.0, 30.0]), lt, pt)    # beyond the 8 A cutoff
    assert s_far == 0.0 and s_clash > s_mid
    pose, s = DK.refine(centre + np.array([1.0, 0.5, 0.0]), lt, pt, maxfev=200)
    assert s <= DK.vina_score(centre + np.array([1.0, 0.5, 0.0]), lt, pt) + 1e-9
    assert np.allclose(np.linalg.norm(pose - pose.mean(0), axis=1), np.linalg.norm(centre, axis=1))   # rigid
    p0, s0 = DK.random_search_vina(centre, np.zeros(3), lt, pt, 50, 0, np.random.default_rng(0))
    assert np.isfinite(s0) and p0.shape == centre.shape
    assert DK.strip_hydrogens(prot).el.tolist() == _pocket().el.tolist()
    assert DK.add_virtual_hydrogens(prot) is prot             # structures with hydrogens are left as they are


def test_loop_fdr_option_keeps_the_default():
    rng = np.random.default_rng(0)
    n = 120
    i, j = np.triu_indices(n, 1)
    d = j - i
    lam = 200.0 * (d + 1.0) ** -1.0
    for a, b in ((20, 45), (60, 90), (30, 100)):
        near = (np.abs(i - a) <= 1) & (np.abs(j - b) <= 1)
        lam[near] *= 6.0
    v = rng.poisson(lam).astype(float)
    ok = v > 0
    ci, cj, cm = i[ok], j[ok], v[ok]
    default = AN.run_suite(ci, cj, cm, n, 10_000)
    same = AN.run_suite(ci, cj, cm, n, 10_000, loop_fdr=0.1)
    strict = AN.run_suite(ci, cj, cm, n, 10_000, loop_fdr=0.01)
    assert [(c.i, c.j) for c in default.loops] == [(c.i, c.j) for c in same.loops]
    assert len(strict.loops) <= len(default.loops)
    assert all(c.q["donut"] <= 0.01 for c in strict.loops)


def test_round2_report_blocks_render():
    sys.path.insert(0, str(ROOT / "validation"))
    import report as R
    for fn in (R.round2_practice, R.round2, R.summary_r2, R.gate6b):
        text = fn()
        assert isinstance(text, str) and text
    assert "round2" in R.BLOCKS and "gate6b" in R.BLOCKS
