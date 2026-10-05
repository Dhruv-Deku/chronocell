"""Phase B: product features. Each test checks that the default app is what it was before Phase B and that
the new feature works when switched on."""

from __future__ import annotations

from pathlib import Path

import pytest

APP = str(Path(__file__).resolve().parent.parent / "app.py")


@pytest.fixture()
def app(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    import ui.common as C
    import ui.settings as SET
    import ui.states_panel as SP
    monkeypatch.setattr(C, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "DEMO_ROOT", tmp_path / "demo")
    monkeypatch.setattr(SET, "PATH", tmp_path / "settings.json")
    at = AppTest.from_file(APP, default_timeout=900)
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at


def _text(at) -> str:
    return "".join(m.value for m in at.markdown)


def test_settings_default_and_save(tmp_path, monkeypatch):
    import ui.settings as SET
    monkeypatch.setattr(SET, "PATH", tmp_path / "s.json")
    assert SET.load() == {"research_mode_default": True}               # unchanged app by default
    SET.save(research_mode_default=False, unknown_key=1)
    assert SET.load() == {"research_mode_default": False}
    assert "unknown_key" not in (tmp_path / "s.json").read_text()


def test_is_synthetic():
    from types import SimpleNamespace as NS
    from ui.common import SYNTHETIC_PREFIX, is_synthetic
    ref = NS(structure_label="x", inputs=(), is_reference=True)
    demo = NS(structure_label="Healthy · demo/synthetic_demo_healthy.pdb", inputs=(), is_reference=False)
    real = NS(structure_label="patient.pdb", inputs=(("structure", "patient.pdb", "0" * 64),), is_reference=False)
    assert is_synthetic(ref) and is_synthetic(demo) and not is_synthetic(real)
    assert is_synthetic(None, SYNTHETIC_PREFIX + "anything") and not is_synthetic(None, "anything")


def test_research_mode_switch_and_synthetic_labels(app):
    at = app
    assert at.toggle(key="research_mode").value is True                # default: everything as before
    assert any("Drug lab" in o for o in at.segmented_control(key="workspace").options)
    at.toggle(key="research_mode").set_value(False).run()
    assert not at.exception, [e.value for e in at.exception]
    assert not any("Drug lab" in o for o in at.segmented_control(key="workspace").options)
    at.segmented_control(key="workspace").set_value("Compare").run()
    assert any(str(o).startswith("SYNTHETIC · Reference model") for o in at.selectbox(key="pe_source").options)
    at.button(key="pe_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    txt = _text(at)
    assert "Research mode is off" in txt and "computed on SYNTHETIC data" in txt
    assert "Senescent: at least 2 of its 3 criteria" not in txt          # the rule table is hidden
    at.toggle(key="research_mode").set_value(True).run()
    assert "Senescent: at least 2 of its 3 criteria" in _text(at)


def test_gate4c_measured_change_and_scoring():
    import sys
    sys.path.insert(0, str(Path(APP).parent / "validation"))
    import numpy as np
    import cohesin_hic as CH
    rng = np.random.default_rng(0)
    n = 40
    sep = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) + 1.0
    u = rng.poisson(2000.0 / sep).astype(float)
    a = 3 * u                                                    # deeper library, same shape: no change
    meas, mask = CH.measured_change(u, a)
    assert mask.any() and not mask[np.tril_indices(n, 1)].any()
    assert np.allclose(meas[mask], np.median(meas[mask]), atol=0.35)
    planted = rng.normal(0, 1, (n, n))
    planted = planted + planted.T
    trend = -np.log2(sep)
    s = CH.score_region({"model": planted, "trend_only": trend}, planted, mask)
    assert s["model"] > 0.99 and s["difference"] > 0.5 and s["difference_ci95"][0] > 0


def _chain_S(n, b2=1.0, seed=0):
    """Pair variances of a random Gaussian chain with some structure (nm² up to a scale)."""
    import numpy as np
    rng = np.random.default_rng(seed)
    steps = rng.normal(0, 1, (n, 3)) * np.sqrt(b2 / 3) * (1 + 0.5 * np.sin(np.arange(n) / 5))[:, None]
    x = np.cumsum(steps, axis=0)
    X = np.stack([np.cumsum(rng.normal(0, 1, (n, 3)) * 0.6, axis=0) for _ in range(200)])
    d2 = ((X[:, :, None, :] - X[:, None, :, :]) ** 2).sum(-1).mean(0) / 3
    return d2 + 0.0 * x[0, 0]


def test_sv_engine_walks_joins_into_derivatives():
    from chronocell import sv_engine as SV
    sizes = {"A": 100, "B": 80}
    d, lost, _ = SV.walk(sizes, [SV.Join("A", 30, "L", "A", 50, "R")])                    # deletion
    assert [s.label() for s in d[0].segments] == ["A:0-30", "A:50-100"] and [x.label() for x in lost] == ["A:30-50"]
    d, lost, _ = SV.walk(sizes, [SV.Join("A", 30, "L", "A", 50, "L"), SV.Join("A", 30, "R", "A", 50, "R")])   # inversion
    assert [s.label() for s in d[0].segments] == ["A:0-30", "A:30-50(-)", "A:50-100"] and not lost
    d, lost, _ = SV.walk(sizes, [SV.Join("A", 40, "L", "B", 20, "R"), SV.Join("A", 40, "R", "B", 20, "L")])   # balanced t
    labels = sorted(" + ".join(s.label() for s in x.segments) for x in d)
    assert labels == ["A:0-40 + B:20-80", "B:0-20 + A:40-100"] and all(x.junctions for x in d)
    with pytest.raises(ValueError):
        SV.walk(sizes, [SV.Join("A", 40, "L", "B", 20, "R"), SV.Join("A", 40, "L", "B", 30, "R")])
    assert SV.breakend_sides("BND", "N[chr2:5[") == ("L", "R") and SV.breakend_sides("BND", "]chr2:5]N") == ("R", "L")
    assert SV.breakend_sides("BND", "+-") == ("L", "R")
    segs = SV.parse_segments("A:0-30 + B:5-20(-)")
    assert segs[1].reverse and segs[1].beads[0] == 19


def test_sv_engine_matches_perturb_on_simple_variants():
    import numpy as np
    from chronocell import perturb as PT, sv_engine as SV
    n, r_c = 60, 2.0
    S = _chain_S(n)
    for op, params, joins in (("deletion", {"a": 20, "b": 30}, [SV.Join("A", 20, "L", "A", 30, "R")]),
                              ("inversion", {"a": 20, "b": 35}, [SV.Join("A", 20, "L", "A", 35, "L"),
                                                                SV.Join("A", 20, "R", "A", 35, "R")])):
        old = PT.variant_impact_from_variance(S, r_c, op, params)
        src = SV.Source("A", S * r_c ** 2, r_c, resolution=10_000)
        kt = SV.karyotype_from_joins({"A": n}, joins, zygosity="homozygous")
        new = SV.impact([src], kt)
        ok = np.isfinite(old.log2_fc) & np.isfinite(new.log2_fc)
        assert ok.sum() > 100 and np.allclose(old.log2_fc[ok], new.log2_fc[ok], atol=1e-8), op
        assert np.array_equal(np.isfinite(old.log2_fc), np.isfinite(new.log2_fc))
    # heterozygous deletion: half the copies change
    src = SV.Source("A", S * r_c ** 2, r_c, resolution=10_000)
    het = SV.impact([src], SV.karyotype_from_joins({"A": n}, [SV.Join("A", 20, "L", "A", 30, "R")]))
    hom = SV.impact([src], SV.karyotype_from_joins({"A": n}, [SV.Join("A", 20, "L", "A", 30, "R")], "homozygous"))
    i, j = 10, 45
    assert np.isclose(het.p_after[i, j], 0.5 * (hom.p_after[i, j] + hom.p_before[i, j]))
    assert het.copy_number[25] == 1 and hom.copy_number[25] == 0 and het.copy_number[5] == 2


def test_sv_engine_translocation_outputs_and_ranking():
    import numpy as np
    from chronocell import sv_engine as SV
    SA, SB = _chain_S(60, seed=1), _chain_S(50, seed=2)
    A = SV.Source("A", SA * 4.0, 2.0, chrom="chr9", resolution=10_000)
    B = SV.Source("B", SB * 4.0, 2.0, chrom="chr22", resolution=10_000)
    kt = SV.karyotype_from_joins({"A": 60, "B": 50}, [SV.Join("A", 30, "L", "B", 25, "R"), SV.Join("A", 30, "R", "B", 25, "L")])
    imp = SV.impact([A, B], kt)
    new = imp.new_contacts()
    assert np.nanmax(new[imp.index("A", 29), 60:]) > 0.01                  # beads joined across the junction touch
    assert np.isnan(imp.log2_fc[0, 70])                                      # no fold change between chromosomes
    genes = SV.affected_genes(imp, {"A": [("GENE1", 29, 25, 33), ("FAR", 2, 0, 4)], "B": [("GENE2", 26, 26, 28)]})
    names = {g["gene"]: g["effects"] for g in genes}
    assert "breakpoint inside the gene" in names["GENE1"] and "next to a new junction" in names["GENE2"]
    ep = SV.enhancer_promoter(imp, [("A", "GENE1", 29)], [("B", 26), ("A", 10)])
    assert any(r["change"] == "new across a junction" for r in ep)
    sc = SV.ranking_score(genes, ep, SV.boundary_changes(imp), interval=False)
    assert sc["score"] >= 3 and sc["genes_deleted_or_broken"] >= 1


def test_sv_engine_copy_number_paths():
    import numpy as np
    from chronocell import sv_engine as SV
    n = 200
    i, j = np.triu_indices(n, 1)
    keep = (j - i) <= 20
    i, j = i[keep], j[keep]
    w = np.ones(n)
    w[80:120] = 0.5                                                        # one copy of 80-120 lost
    cm = 100.0 / (j - i) * np.sqrt(w[i] * w[j]) * 2
    est, segs = SV.estimate_copy_number(i, j, cm, n)
    assert any(a <= 85 and b >= 115 and cn == 1 for a, b, cn in segs), segs
    kt = SV.karyotype_from_copy_number(60, [(20, 30, 1.0), (40, 50, 4.0)])
    S = _chain_S(60)
    imp = SV.impact([SV.Source("A", S, 1.0, resolution=10_000)], kt)
    assert imp.copy_number[25] == 1 and imp.copy_number[45] == 4 and imp.copy_number[5] == 2
    assert np.nanmean(imp.log2_fc[40:50, 40:50]) > 0.5                     # extra copies add internal contacts
    g = SV.vcf_genotypes("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                         "chr1\t100\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=500\tGT:CN\t1/1:0\n")
    assert g[3] == {"gt": "1/1", "cn": 0.0} and SV.zygosity("1/1") == "homozygous" and SV.zygosity("0|1") == "heterozygous"
    assert SV.read_cnv_bed("chr1\t1000\t5000\t3\nchr2\t0\t10\t1\n", "chr1", 0, 10, 1000) == [(1, 5, 3.0)]


def _planted_map(n=300, seed=0, loops=((40, 80), (120, 170), (200, 260)), domains=(0, 60, 140, 220, 300)):
    import numpy as np
    rng = np.random.default_rng(seed)
    i, j = np.triu_indices(n, 1)
    s = (j - i).astype(float)
    lam = 4000.0 / s ** 1.0
    dom = np.searchsorted(domains, i, side="right") == np.searchsorted(domains, j, side="right")
    lam = lam * np.where(dom, 2.5, 1.0)
    for a, b in loops:
        near = (np.abs(i - a) <= 1) & (np.abs(j - b) <= 1)
        lam = lam * np.where((i == a) & (j == b), 12.0, np.where(near, 4.0, 1.0))
    c = rng.poisson(lam).astype(float)
    keep = c > 0
    return i[keep], j[keep], c[keep]


def test_analysis_suite_finds_planted_loops_and_domains():
    import numpy as np
    from chronocell import analysis as AN
    n = 300
    ci, cj, cm = _planted_map(n)
    calls = AN.call_loops(ci, cj, cm, n, 10_000, max_sep_bp=1_500_000)
    found = {(c.i, c.j) for c in calls}
    hits = sum(any(abs(a - x) <= 1 and abs(b - y) <= 1 for x, y in found) for a, b in ((40, 80), (120, 170), (200, 260)))
    assert hits == 3 and len(calls) <= 5, [(c.i, c.j, round(c.oe_donut, 1)) for c in calls]
    td, sig = AN.topdom(ci, cj, cm, n, window=5)
    assert sum(any(abs(b - x) <= 3 for x in td) for b in (60, 140, 220)) >= 2, td
    ps = AN.p_of_s(ci, cj, cm, n, 10_000)
    assert -1.6 < ps["slope"] < -0.6
    rep = AN.run_suite(ci, cj, cm, n, 10_000)
    assert rep.summary()["loops"] >= 3 and len(rep.insulation_boundaries) >= 2
    k = AN.kernels(2, 5)
    assert k["donut"][5, 5] == 0 and k["lower_left"][8, 2] == 1 and k["lower_left"][2, 8] == 0
    assert np.allclose(AN._bh(np.array([0.01, 0.04, 0.03, 0.5])), [0.04, 0.16 / 3, 0.16 / 3, 0.5])


def test_differential_controls_false_discoveries_and_finds_planted_changes():
    import numpy as np
    from chronocell import differential as DF
    rng = np.random.default_rng(3)
    n = 120
    sep = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) + 1.0
    base = 3000.0 / sep * np.exp(rng.normal(0, 0.3, (n, n)))
    base = (base + base.T) / 2

    def rep(scale):
        m = rng.poisson(base * scale * np.exp(rng.normal(0, 0.1, (n, n)))).astype(float)
        m = np.triu(m, 1)
        return m + m.T
    A = [rep(1.0), rep(1.3)]
    B = [rep(0.8), rep(1.1)]
    null = DF.differential_pixels(A, B, 10_000, fdr=0.05, max_sep=40)
    assert null.statistics and null.summary()["significant"] <= 0.01 * len(null.pixels) + 2
    iu = [(i, i + k) for i in range(5, 100, 7) for k in (5, 12)]
    Bp = [DF.plant_changes(m, iu, 4.0, rng) for m in B]
    res = DF.differential_pixels(A, Bp, 10_000, fdr=0.05, max_sep=40)
    sig = {(int(r.i), int(r.j)) for r in res.significant.itertuples()}
    tp = len(sig & set(iu))
    assert tp >= 0.5 * len(iu) and (len(sig) - tp) <= 0.1 * max(len(sig), 1) + 2
    one = DF.differential_pixels(A[:1], B[:1], 10_000)
    assert not one.statistics and "two replicates" in one.notes[-1]
    assert DF.aggregate(np.ones((10, 10)), 2).shape == (5, 5) and DF.aggregate(np.ones((10, 10)), 2)[0, 0] == 4
    full = DF.compare(A, Bp, 10_000, max_sep_bp=400_000)
    assert full.tracks["insulation"] is not None and full.boundaries is not None and full.loops is not None
    bed = DF.to_bedpe(res.significant.head(3), "chr1", 1_000_000, 10_000)
    assert bed.count("\n") == min(3, len(res.significant)) and bed.startswith("chr1\t")


def test_liftover_chain_parsing_and_mapping():
    from chronocell import liftover as LO
    chain = (b"chain 1000 chr1 1000 + 100 400 chrA 2000 + 500 800 1\n100 50 50\n150\n\n"
             b"chain 500 chr2 1000 + 0 100 chrB 1000 - 0 100 2\n100\n")
    ch = LO.Chain.parse(chain)
    assert ch.position("chr1", 100) == ("chrA", 500, "+")
    assert ch.position("chr1", 210) is None                         # in the gap
    assert ch.position("chr1", 260) == ("chrA", 660, "+")           # second block: 250 -> 650
    assert ch.position("chr2", 10) == ("chrB", 1000 - 10 - 1, "-")
    assert ch.interval("chr1", 120, 180) == ("chrA", 520, 580)
    txt, lost = ch.bed("chr1\t120\t180\tx\nchr1\t205\t215\ty\n")
    assert txt.startswith("chrA\t520\t580\tx") and lost == 1


def test_contact_reader_cool_pairs_and_errors(tmp_path):
    import gzip
    import h5py
    import numpy as np
    from chronocell import contacts_io as IO
    # a minimal cooler (schema v3) with chr1 of 10 bins at 1 kb
    path = tmp_path / "t.cool"
    n = 10
    b1, b2 = np.triu_indices(n)
    cnt = (10 + b1 + b2).astype(np.int32)
    with h5py.File(path, "w") as f:
        f.attrs["bin-size"] = 1000
        g = f.create_group("chroms")
        g["name"] = np.array([b"chr1"])
        g["length"] = np.array([248_956_422])
        f.create_group("bins")["start"] = np.arange(n) * 1000
        px = f.create_group("pixels")
        px["bin1_id"], px["bin2_id"], px["count"] = b1, b2, cnt
        ix = f.create_group("indexes")
        ix["chrom_offset"] = np.array([0, n])
        ix["bin1_offset"] = np.searchsorted(b1, np.arange(n + 1))
    c = IO.read_region(str(path), "1", 2000, 6000, 1000)
    assert c.n == 4 and "hg38" in c.notes[0]
    d = c.dense()
    assert d[0, 1] == 10 + 2 + 3 and d[1, 0] == d[0, 1]
    c2 = IO.read_region(str(path), "chr1", 0, 10000, 2000)
    assert c2.n == 5 and c2.dense()[0, 0] == 10 + 0 + 0 + (10 + 0 + 1) + (10 + 1 + 1)
    with pytest.raises(IO.ContactFileError, match="hg38.*session uses hg19"):
        IO.read_region(str(path), "chr1", 0, 5000, 1000, expected_assembly="hg19")
    with pytest.raises(IO.ContactFileError, match="not in this file"):
        IO.read_region(str(path), "chr7", 0, 5000, 1000)
    with pytest.raises(IO.ContactFileError, match="does not divide"):
        IO.read_region(str(path), "chr1", 0, 5000, 1500)
    bad = tmp_path / "bad.cool"
    bad.write_bytes(b"not hdf5")
    with pytest.raises(IO.ContactFileError, match="corrupt"):
        IO.read_region(str(bad), "chr1", 0, 5000, 1000)
    pairs = tmp_path / "p.pairs.gz"
    with gzip.open(pairs, "wt") as fh:
        fh.write("## pairs format v1.0\n#chromsize: chr1 248956422\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n")
        for k in range(50):
            fh.write(f"r{k}\tchr1\t{1000 + 37 * k}\tchr1\t{3000 + 41 * k}\t+\t-\n")
        fh.write("x\tchr2\t5\tchr2\t9\t+\t+\n")
    p = IO.read_region(str(pairs), "chr1", 0, 10_000, 1000)
    assert p.cm.sum() == 50 and "streamed" in p.notes[-1]
    k = IO.read_region(str(path), "chr1", 0, 10_000, 1000, normalization="kr")
    assert k.weights is not None and np.isfinite(k.weights).sum() >= 5
    with pytest.raises(IO.ContactFileError, match="not a contact file"):
        IO.kind_of("x.bam")


def _write_fixture_mcool(path, seed=0, extra_loop=None):
    from chronocell import exports as EX
    ci, cj, cm = _planted_map(300, seed=seed, loops=((40, 80), (120, 170), (200, 260)) + ((extra_loop,) if extra_loop else ()))
    L = 248_956_422
    EX.write_mcool(path, [("chr1", L)], {"chr1": (ci + 10_000, cj + 10_000, cm)}, 10_000, factors=(1, 2))
    return path


def test_mcool_export_roundtrip_cli_analyze_diff_report(tmp_path):
    import json
    import numpy as np
    from chronocell import cli, contacts_io as IO
    f1 = _write_fixture_mcool(tmp_path / "a1.mcool", 1)
    assert IO.resolutions(str(f1)) == [10_000, 20_000] and IO.detect_assembly(IO.chromosomes(str(f1))) == "hg38"
    region = "chr1:100000000-103000000"
    c = IO.read_region(str(f1), "chr1", 100_000_000, 103_000_000, 10_000)
    assert c.n == 300 and c.cm.sum() > 1000
    out = tmp_path / "an"
    cli.main(["analyze", str(f1), "--region", region, "--res", "10000", "--norm", "kr", "--out", str(out)])
    summ = json.loads((out / "summary.json").read_text())["summary"]
    assert summ["loops"] >= 3 and (out / "loops.bedpe").read_text().count("\n") == summ["loops"]
    assert (out / "run.json").exists() and json.loads((out / "run.json").read_text())["inputs"][0]["sha256"]
    files = [_write_fixture_mcool(tmp_path / f"{k}.mcool", s, (60, 110) if k.startswith("b") else None)
             for k, s in (("a2", 2), ("b1", 3), ("b2", 4))]
    dout = tmp_path / "df"
    cli.main(["diff", "--a", str(f1), str(files[0]), "--b", str(files[1]), str(files[2]), "--region", region,
              "--res", "10000", "--max-sep", "800000", "--out", str(dout)])
    ds = json.loads((dout / "summary.json").read_text())["summary"]
    assert ds["statistics"] and ds["replicates_a"] == 2 and (dout / "differential_pixels_all.tsv").exists()
    loops = __import__("pandas").read_csv(dout / "differential_loops.csv")
    assert ((loops["change"] == "gained in B") & (abs(loops["i"] - 60) <= 2) & (abs(loops["j"] - 110) <= 2)).any()
    for folder in (out, dout):
        cli.main(["report", str(folder), "--pdf"])
        html = (folder / "report.html").read_text(encoding="utf-8")
        assert "Reproducibility record" in html and "<script" not in html and (folder / "report.pdf").stat().st_size > 1000
    assert "Gate 7" in (dout / "report.html").read_text(encoding="utf-8")
    sheet = tmp_path / "s.csv"
    sheet.write_text("sample,path,region,resolution,condition\n" + "\n".join(
        f"{p.stem},{p},{region},10000,{'ctrl' if p.stem.startswith('a') else 'treat'}" for p in [f1] + files) + "\n")
    bout = tmp_path / "batch"
    cli.main(["batch", str(sheet), "--out", str(bout)])
    assert (bout / "a1" / "summary.json").exists() and (bout / "diff_ctrl_vs_treat" / "summary.json").exists()
    cli.main(["batch", str(sheet), "--out", str(bout)])          # resumes: nothing redone, no error


def test_cli_impact_with_joins_and_api_endpoints(tmp_path):
    import json
    import numpy as np
    from chronocell import api, cli
    f = _write_fixture_mcool(tmp_path / "v.mcool", 5)
    out = tmp_path / "imp"
    cli.main(["impact", str(f), "--region", "chr1:100000000-101500000", "--res", "10000", "--joins", "A:40:L-A:60:R",
              "--zygosity", "homozygous", "--assembly", "hg38", "--out", str(out)])
    rank = __import__("pandas").read_csv(out / "variant_ranking.csv")
    assert len(rank) == 1 and rank["score"].iloc[0] >= 0
    s = json.loads((out / "summary.json").read_text())
    assert "not validated" in s["standing"] and s["variants"][0]["lost_pieces"] == ["A:40-60"]
    cli.main(["report", str(out)])
    assert "Mechanism simulator, not validated" in (out / "report.html").read_text(encoding="utf-8")
    ci, cj, cm = _planted_map(120, seed=7, loops=((20, 60),), domains=(0, 50, 120))
    cont = {"i": ci.tolist(), "j": cj.tolist(), "count": cm.tolist(), "n": 120, "resolution": 10_000, "chrom": "chr1"}
    log = api.AuditLog(tmp_path / "log.jsonl")
    r = api.analyze({"contacts": cont}, log=log)
    assert r["summary"]["loops"] >= 1 and r["run_id"]
    r = api.diff({"condition_a": [cont, cont], "condition_b": [cont, cont], "max_sep_bp": 300_000}, log=log)
    assert r["summary"]["statistics"] and r["summary"]["significant"] == 0
    r = api.impact({"sources": [{"name": "A", "contacts": cont}], "joins": [{"source1": "A", "cut1": 30, "side1": "L",
                    "source2": "A", "cut2": 45, "side2": "R"}]}, log=log)
    assert r["summary"]["derivatives"] == ["A:0-30 + A:45-120"] and r["top_changes"]
    with pytest.raises(api.RequestError):
        api.impact({"sources": [{"name": "A", "contacts": cont}]}, log=log)
    assert [x["endpoint"] for x in log.read()][-1] == "/api/v1/impact" and log.read()[-1]["status"] == "rejected"


def test_job_queue_runs_stops_and_recovers(tmp_path):
    import json
    import time
    from chronocell import jobs as JQ
    q = JQ.Queue(tmp_path / "jobs", max_cpu=1)
    folder = tmp_path / "rep"
    folder.mkdir()
    (folder / "summary.json").write_text(json.dumps({"summary": {"loops": 1}}))
    j1 = q.submit(["report", str(folder)], "cpu", "report")
    j2 = q.submit(["report", str(folder)], "gpu", "report gpu")
    j3 = q.submit(["report", str(folder)], "cpu", "to stop")
    q.stop(j3["id"])
    q.work(poll=0.2, until_empty=True)
    states = {j["id"]: j["state"] for j in q.jobs()}
    assert states[j1["id"]] == "done" and states[j2["id"]] == "done" and states[j3["id"]] == "stopped"
    assert (folder / "report.html").exists() and "attempt 1" in q.log(j1["id"])
    j4 = q.submit(["report", str(folder)], "cpu")
    job = q.get(j4["id"])
    job.update(state="running", pid=999_999_9)                 # a worker that died mid-job
    q._put(job)
    assert q.recover() == [j4["id"]] and q.get(j4["id"])["attempt"] == 1


def test_projects_save_and_reopen(tmp_path):
    import numpy as np
    from chronocell import projects as PRJ
    d = PRJ.save("My study", {"chrom": "chr21", "res": 10_000, "obj": object(), "pair": [1, 2]},
                 {"c:/data/x.mcool": b"abc"}, {"telemetry": [1, 2], "arr": np.arange(3)}, root=tmp_path)
    assert d.name == "My study" and PRJ.list_projects(tmp_path)[0]["files"] == 1
    p = PRJ.load("My study", root=tmp_path)
    assert p["settings"] == {"chrom": "chr21", "res": 10_000, "pair": [1, 2]} and p["files"] == {"x.mcool": b"abc"}
    assert np.array_equal(p["results"]["arr"], np.arange(3))
    (d / "data" / "x.mcool").write_bytes(b"abd")
    with pytest.raises(ValueError, match="changed"):
        PRJ.load("My study", root=tmp_path)
    with pytest.raises(ValueError):
        PRJ.save("///", {}, root=tmp_path)


@pytest.fixture()
def app_b(tmp_path, monkeypatch):
    """The app with the Phase B platform folders in tmp_path and no worker process started."""
    from streamlit.testing.v1 import AppTest
    import ui.common as C
    import ui.settings as SET
    import ui.states_panel as SP
    from chronocell import jobs as JQ, projects as PRJ
    monkeypatch.setattr(C, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "DEMO_ROOT", tmp_path / "demo")
    monkeypatch.setattr(SET, "PATH", tmp_path / "settings.json")
    monkeypatch.setattr(PRJ, "ROOT", tmp_path / "projects")
    monkeypatch.setattr(JQ, "ROOT", tmp_path / "jobs")
    monkeypatch.setattr(JQ, "start_worker", lambda root=None: 0)
    at = AppTest.from_file(APP, default_timeout=900)
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at


def test_app_differential_tab_on_replicate_maps(app_b, tmp_path):
    at = app_b
    files = [_write_fixture_mcool(tmp_path / f"{k}.mcool", s, (60, 110) if k.startswith("b") else None)
             for k, s in (("a1", 1), ("a2", 2), ("b1", 3), ("b2", 4))]
    at.segmented_control(key="workspace").set_value("Compare").run()
    assert [t.label for t in at.tabs][:3] == ["Side by side", "Self-Math PDB State Evaluator", "Differential analysis"]
    at.text_area(key="df_a_paths").set_value(f"{files[0]}\n{files[1]}")
    at.text_area(key="df_b_paths").set_value(f"{files[2]}\n{files[3]}")
    at.text_input(key="df_chrom").set_value("chr1")
    at.number_input(key="df_start").set_value(100_000_000)
    at.number_input(key="df_end").set_value(102_000_000).run()
    at.button(key="df_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    s = at.session_state["df_last"]["result"]["summary"]
    assert s["statistics"] and s["replicates_a"] == 2 and s["tested_pixels"] > 1000
    assert "Gate 7" in _text(at)


def test_app_analysis_suite_and_engine_v2(app_b):
    at = app_b
    at.button(key="an_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    r = at.session_state["an_last"]["result"]
    assert r["summary"]["bins"] > 50 and "Gate 6" in _text(at)
    at.session_state["custom_window"] = (2000, 2150)
    at.segmented_control(key="region_choice").set_value("custom").run()
    next(b for b in at.button if (b.label or "").startswith("Build population model")).click().run()
    assert not at.exception, [e.value for e in at.exception]
    at.segmented_control(key="workspace").set_value("4D dynamics").run()
    assert not at.exception, [e.value for e in at.exception]
    at.text_area(key="ve_joins").set_value("A:40:L-A:60:R")
    at.button(key="ve_go").click().run()
    assert not at.exception, [e.value for e in at.exception]
    res = at.session_state["ve_last"]["results"]
    assert res[0][1]["summary"]["lost_pieces"] == ["A:40-60"]
    assert "Mechanism simulator, not validated" in _text(at)


def test_app_projects_and_jobs_panels(app_b, tmp_path):
    from chronocell import jobs as JQ, projects as PRJ
    at = app_b
    at.text_input(key="prj_name").set_value("Study one").run()
    at.button(key="prj_save").click().run()
    assert not at.exception, [e.value for e in at.exception]
    assert [p["name"] for p in PRJ.list_projects()] == ["Study one"]
    at.selectbox(key="prj_pick").set_value("Study one")
    at.button(key="prj_open").click().run()
    assert not at.exception and "Opened Study one" in _text(at)
    at.text_input(key="job_cmd").set_value(f"report {tmp_path}").run()
    at.button(key="job_submit").click().run()
    assert not at.exception, [e.value for e in at.exception]
    jobs = JQ.Queue().jobs()
    assert len(jobs) == 1 and jobs[0]["state"] == "queued" and jobs[0]["args"][0] == "report"


def test_bigwig_reader_and_track_ingest(tmp_path):
    import numpy as np
    from chronocell import bigwig as BW, genome, ingest
    p = tmp_path / "t.bw"
    BW.write_bigwig_bedgraph(p, [("chr22", 50_818_468)], {"chr22": [(20_000_000, 20_005_000, 2.0), (20_005_000, 20_010_000, 4.0),
                                                                     (30_000_000, 30_001_000, 8.0)]})
    b = BW.BigWig(str(p))
    s, e, v = b.intervals("22", 19_000_000, 21_000_000)
    assert list(v) == [2.0, 4.0]
    assert np.allclose(b.bin_means("chr22", 20_000_000, 20_020_000, 10_000), [3.0, np.nan], equal_nan=True)
    ch = genome.chrom("chr22", 10_000)
    vals, note = ingest.read_track(p.read_bytes(), "t.bw", ch)
    assert len(vals) == ch.n_bins and vals[2000] == 3.0 and vals[3000] == 8.0
    assert np.isnan(vals[100]) and ("built-in" in note or "bigWig" in note)
    with pytest.raises(ValueError):
        BW.BigWig(str(tmp_path / "t.bw").replace("t.bw", "missing.bw")) if False else BW.BigWig.__init__(object.__new__(BW.BigWig), str(_bad(tmp_path)))


def _bad(tmp_path):
    q = tmp_path / "bad.bw"
    q.write_bytes(b"\x00" * 64)
    return q


def test_annotations_from_local_files(tmp_path, monkeypatch):
    import gzip
    from chronocell import annotations as AO, pipelines as PL
    monkeypatch.setattr(AO, "CACHE", tmp_path)
    AO.clinvar_genes.cache_clear()
    AO.gtex_median_tpm.cache_clear()
    assert AO.clinvar_genes(False) is None and AO.gtex_median_tpm(False) is None
    (tmp_path / "gene_specific_summary.txt").write_text(
        "Overview of ClinVar\n#Symbol\tGeneID\tTotal_submissions\tAlleles_reported_Pathogenic_Likely_pathogenic\n"
        "TP53\t7157\t5000\t1200\nFOO\t1\t3\t0\n")
    with gzip.open(tmp_path / AO.Path(AO.GTEX_URL).name, "wt") as fh:
        fh.write("#1.2\n2\t2\nName\tDescription\tLiver\tLung\nENSG1\tTP53\t10.5\t20.0\nENSG2\tFOO\t0.1\t0.2\n")
    assert AO.clinvar_genes(False) == {"TP53": 1200, "FOO": 0}
    g = AO.gtex_for(["TP53", "BAR"], "Lung")
    assert list(g["gene"]) == ["TP53"] and float(g["Lung"].iloc[0]) == 20.0
    assert AO.cosmic_from_file(b"Gene Symbol,Tier\nTP53,1\nMYC,1\n", "cgc.csv") == {"TP53", "MYC"}
    ann = PL.annotate_genes(["TP53", "FOO"])
    assert ann["TP53"]["clinvar_pathogenic_alleles"] == 1200 and ann["TP53"]["curated_list"] == "cancer"
    AO.clinvar_genes.cache_clear()
    AO.gtex_median_tpm.cache_clear()


def test_gate6_matching_and_gate7_summary():
    import sys
    import pandas as pd
    sys.path.insert(0, str(Path(APP).parent / "validation"))
    import diff_gate7 as G7
    import loops_gate6 as G6
    ref = pd.DataFrame({"chrom": ["chr4"] * 3, "a": [1_000_000, 2_000_000, 3_000_000], "b": [1_500_000, 2_400_000, 3_900_000]})
    calls = pd.DataFrame({"chrom": ["chr4"] * 4, "a": [1_010_000, 1_020_000, 2_000_000, 5_000_000],
                          "b": [1_490_000, 1_500_000, 2_430_000, 5_500_000]})
    m = G6.match(calls, ref)
    assert m["matched"] == 1                  # one-to-one: two calls near loop 1 count once; loop 2 is 30 kb off
    assert m["precision"] == 0.25 and abs(m["recall"] - 1 / 3) < 1e-9
    rows = [{"region": "r1", "fold": 1.0, "discoveries": 0, "false": 0, "planted": 0, "found": 0},
            {"region": "r1", "fold": 2.0, "discoveries": 10, "false": 1, "planted": 100, "found": 9},
            {"region": "r2", "fold": 4.0, "discoveries": 50, "false": 0, "planted": 100, "found": 50},
            {"region": "r2", "fold": 2.0, "discoveries": 0, "false": 0, "planted": 100, "found": 0}]
    s = G7.summarise(rows)
    assert abs(s["mean_fdp"] - 0.1 / 3) < 1e-9 and s["null_discoveries"] == {"r1": 0} and s["recall"]["4.0"] == 0.5
