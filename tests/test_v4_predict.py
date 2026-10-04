"""No-contact-data predictor (chronocell/predict.py, Gate 5): motif orientation, features, ridge, the
frozen model in the app, and the population fitted to a predicted distance map."""

from __future__ import annotations

import gzip

import numpy as np
import pytest

from chronocell import predict as PD


def _random_seq(n: int, seed: int = 0) -> bytes:
    return np.random.default_rng(seed).choice(np.frombuffer(b"ACGT", np.uint8), size=n).tobytes()


def _plant(seq: bytes, pos: int, site: str) -> bytes:
    return seq[:pos] + site.encode() + seq[pos + len(site):]


def test_motif_orientation_is_read_from_the_sequence():
    pwm = PD.load_pwm()
    cons = "".join("ACGT"[k] for k in PD.read_jaspar(PD.JASPAR_PATH.read_text()).argmax(0))
    rc = cons[::-1].translate(str.maketrans("ACGT", "TGCA"))
    seq = _random_seq(20_000, 1)
    seq = _plant(seq, 5_000, cons)
    seq = _plant(seq, 12_000, rc)
    starts, ends = np.array([4_900, 11_900, 17_000]), np.array([5_200, 12_200, 17_300])
    strand, score = PD.orient_peaks(seq, starts, ends, np.array([-1, -1, -1]), pwm, 0.8, 100)
    assert strand[0] == 1 and strand[1] == -1 and score[0] > 0.99 and score[1] > 0.99
    assert strand[2] == 0 or score[2] >= 0.8               # random sequence: usually no match above 0.8


def test_peaks_are_read_from_narrowpeak_bed_and_gzip():
    text = ("track name=ctcf\n"
            "chr21\t1000\t1400\tp1\t0\t.\t5.0\t4\t3\t200\n"
            "21\t5000\t5300\tp2\t0\t.\t5.0\t4\t3\t-1\n"
            "chr22\t10\t20\tother\n"
            "chr21\t9000\t9100\n")
    for data in (text, text.encode(), gzip.compress(text.encode())):
        pk = PD.read_peaks(data, "chr21", lambda c: c if c.startswith("chr") else "chr" + c)
        assert pk["starts"].tolist() == [1000, 5000, 9000] and pk["summits"].tolist() == [200, -1, -1]


def test_features_and_ridge_recover_a_planted_effect():
    rng = np.random.default_rng(2)
    starts = np.arange(40) * 30_000
    loci = PD.Loci(starts, starts + 30_000, rng.integers(0, 2, 40), rng.integers(0, 2, 40), np.zeros(40, int),
                   rng.normal(0.42, 0.03, 40))
    mid = np.sort(rng.integers(0, 40 * 30_000, 30))
    i, j = np.triu_indices(40, 1)
    X = PD.pair_features(loci, mid, i, j)
    assert X.shape == (len(i), 12) and np.isfinite(X).all()
    between = np.searchsorted(mid, loci.starts[j]) - np.searchsorted(mid, loci.ends[i])
    assert np.allclose(X[:, 4], np.log1p(between))
    beta_true = np.zeros(12)
    beta_true[4] = 0.3
    y = ((X - X.mean(0)) / X.std(0)) @ beta_true + rng.normal(0, 0.01, len(i))
    beta, mu, sd = PD.fit_ridge([X], [y], 1e-4)
    assert abs(beta[4] - 0.3) < 0.05 and np.max(np.abs(np.delete(beta, 4))) < 0.06


def test_frozen_model_predicts_a_valid_map_on_any_sequence():
    model, meta = PD.load_model()
    assert meta["assembly"] == "hg38" and meta["validation"]["result_file"] == "validation/results_predictor.json"
    seq = _random_seq(3_000_000, 3)
    starts = 1_000_000 + np.arange(50) * 30_000
    pk = {"starts": np.array([1_100_000, 1_600_000, 2_200_000]), "ends": np.array([1_100_400, 1_600_400, 2_200_400]),
          "summits": np.array([200, 200, -1])}
    d, info = PD.predict_window(starts, starts + 30_000, seq, pk, model, meta["settings"])
    assert d.shape == (50, 50) and np.allclose(d, d.T) and np.all(np.diag(d) == 0)
    off = d[np.triu_indices(50, 1)]
    assert np.isfinite(off).all() and (off > 0).all()
    assert info["loci"] == 50 and info["peaks"] == 3
    # the separation trend dominates: farther along the DNA is farther apart on average
    assert np.mean(np.diag(d, 40)) > np.mean(np.diag(d, 1))


def test_population_follows_a_predicted_map():
    torch = pytest.importorskip("torch")  # noqa: F841
    from chronocell import population as P
    rng = np.random.default_rng(4)
    n = 40
    s = np.abs(np.subtract.outer(np.arange(n), np.arange(n))).astype(float)
    d = 120.0 * np.maximum(s, 1) ** 0.4 * np.exp(rng.normal(0, 0.05, (n, n)))
    d = (d + d.T) / 2
    np.fill_diagonal(d, 0.0)
    res = P.fit_population_from_medians(d)
    iu = np.triu_indices(n, 1)
    from scipy.stats import spearmanr
    assert spearmanr(res.median_distance_nm[iu], d[iu]).correlation > 0.95
    assert res.config["input"] == "distance map (no contact data)"
    with pytest.raises(ValueError):
        P.fit_population_from_medians(np.zeros((2, 2)))
