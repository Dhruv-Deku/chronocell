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
