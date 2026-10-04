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


def _cli_inputs(tmp_path, n: int = 3_000_000, records=("chr21",)):
    fa = tmp_path / "seq.fa.gz"
    with gzip.open(fa, "wb") as fh:
        for k, name in enumerate(records):
            seq = _random_seq(n, 10 + k)
            fh.write(b">" + name.encode() + b" test\n" + b"\n".join(seq[i:i + 60] for i in range(0, n, 60)) + b"\n")
    peaks = tmp_path / "ctcf.narrowPeak"
    peaks.write_text("".join(f"chr21\t{1_000_000 + k * 150_000}\t{1_000_000 + k * 150_000 + 400}\tp{k}\t0\t.\t5\t5\t5\t200\n"
                             for k in range(10)) + "chr22\t5\t500\tother\t0\t.\t5\t5\t5\t100\n")
    return fa, peaks


def test_cli_writes_the_map_and_a_record_with_hashes(tmp_path, capsys):
    import hashlib
    import json
    fa, peaks = _cli_inputs(tmp_path, records=("chr20", "chr21"))
    out = tmp_path / "map.npy"
    PD.main(["--chrom", "21", "--start", "1000000", "--end", "2500000", "--peaks", str(peaks), "--fasta", str(fa),
             "--out", str(out)])
    d = np.load(out)
    assert d.shape == (50, 50) and np.allclose(d, d.T) and (d[np.triu_indices(50, 1)] > 0).all()
    rec = json.loads(out.with_suffix(".json").read_text())
    assert rec["status"] == "predicted, not measured"
    assert rec["region"] == {"assembly": "hg38", "chrom": "chr21", "start": 1_000_000, "end": 2_500_000, "bin_bp": 30_000,
                             "loci": 50}
    assert rec["inputs"]["peaks"]["sha256"] == hashlib.sha256(peaks.read_bytes()).hexdigest()
    assert rec["inputs"]["peaks"]["peaks"] == 10                     # the chr22 line is not on this chromosome
    assert rec["inputs"]["sequence"]["sha256"] == hashlib.sha256(fa.read_bytes()).hexdigest()
    assert rec["output"]["sha256"] == hashlib.sha256(out.read_bytes()).hexdigest()
    assert rec["model"]["sha256"] == hashlib.sha256(PD.MODEL_PATH.read_bytes()).hexdigest()
    assert rec["validation"]["result_file"] == "validation/results_predictor.json" and "warning" not in rec["validation"]
    # the same inputs as the app's prediction give the same map
    model, meta = PD.load_model()
    seq = PD.read_fasta_record(fa, "chr21")
    starts = 1_000_000 + np.arange(50) * 30_000
    d2, _ = PD.predict_window(starts, starts + 30_000, seq, PD.read_peaks(peaks.read_bytes(), "chr21"), model,
                              meta["settings"])
    assert np.allclose(d, d2)
    assert "PREDICTED, not measured" in capsys.readouterr().out


def test_cli_refuses_what_was_not_validated(tmp_path):
    import json
    fa, peaks = _cli_inputs(tmp_path)
    with pytest.raises(SystemExit, match="hg38 only"):
        PD.main(["--chrom", "chr1", "--start", "0", "--end", "300000", "--peaks", str(peaks), "--fasta", str(fa),
                 "--assembly", "mm39", "--out", str(tmp_path / "m.npy")])
    with pytest.raises(SystemExit, match="outside the chromosome"):
        PD.main(["--chrom", "chr21", "--start", "0", "--end", "999000000", "--peaks", str(peaks), "--fasta", str(fa),
                 "--out", str(tmp_path / "m.npy")])
    with pytest.raises(SystemExit, match="shorter than the region"):
        PD.main(["--chrom", "chr21", "--start", "2000000", "--end", "4000000", "--peaks", str(peaks), "--fasta", str(fa),
                 "--out", str(tmp_path / "m.npy")])
    # an untested bin size runs, with a warning in the record
    PD.main(["--chrom", "chr21", "--start", "1000000", "--end", "2500000", "--bin", "100000", "--peaks", str(peaks),
             "--fasta", str(fa), "--out", str(tmp_path / "m.npy")])
    assert "was not tested" in json.loads((tmp_path / "m.json").read_text())["validation"]["warning"]


def test_cli_reproduces_the_gate5_map_on_real_inputs(tmp_path):
    """The app / CLI path (read_peaks, offset GC index) gives exactly the map scored in Gate 5 for the IMR-90
    test loci. Needs the downloaded validation inputs (validation/data, git-ignored); skipped without them."""
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parent.parent / "validation"
    peaks, fasta = root / "data" / "encode" / "ENCFF670ULH.bed.gz", root / "data" / "hg38" / "chr21.fa.gz"
    if not (peaks.exists() and fasta.exists() and (root / "data" / "IMR90_chr21-28-30Mb.csv").exists()):
        pytest.skip("validation inputs not downloaded (python validation/predictor.py --test fetches them)")
    sys.path.insert(0, str(root))
    import datasets as D
    import predictor as VP
    tr = D.load("bintu_imr90_28_30")
    d = VP.design("bintu_imr90_28_30", tr)
    model = PD.Predictor.from_json((root / "predictor_model.json").read_text())
    n = d["loci"].n
    want = np.exp(VP._matrix(n, d["i"], d["j"], model.log_distance(d["X"], d["sep"].astype(float))))
    start, step = int(tr.starts[0]), int(tr.starts[1] - tr.starts[0])
    PD.main(["--chrom", "chr21", "--start", str(start), "--end", str(start + n * step), "--bin", str(step),
             "--peaks", str(peaks), "--fasta", str(fasta), "--out", str(tmp_path / "m.npy")])
    got = np.load(tmp_path / "m.npy")
    iu = np.triu_indices(n, 1)
    assert np.allclose(got[iu], want[iu], rtol=1e-12)


def test_cohesin_control_caveat_is_read_from_the_result_file():
    from chronocell import accuracy as ACC
    import ui.predict_view as PV
    pred = ACC.load_benchmark()["models"]["predicted_sequence_ctcf"]
    r = __import__("json").loads((ACC.VALIDATION / "results_predictor.json").read_text(encoding="utf-8"))
    ctl = r["settings"]["control"][0]
    want = r["summary"][ctl]["sequence + CTCF"]["all_pairs"]["predictor"]["percent_of_ceiling"]
    assert pred["control_percent_of_ceiling"][ctl] == round(want, 1)
    txt = PV.control_caveat(pred)
    assert f"{want:.1f} %" in txt and "not CTCF loops" in txt
    # a control that scored low would be read the other way, with no claim about compartments
    low = dict(pred, control_percent_of_ceiling={ctl: 1.0})
    assert "not CTCF loops" not in PV.control_caveat(low)
    assert PV.control_caveat(dict(pred, control_percent_of_ceiling={})) is None
