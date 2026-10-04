"""v4 Pillar 4: ensemble perturbations (cohesin depletion, structural variants on the spring network),
SV file parsing and chromosome-name normalisation."""

from __future__ import annotations

import numpy as np
import pytest

from chronocell import ensemble as E, genome, perturb as PT, svio
from tests.test_v33 import _population


@pytest.fixture(scope="module")
def chain_result():
    d = _population(n=40, cells=6000, seed=11)                          # planted ideal Gaussian chains
    return E.fit_ensemble((d < 1.0).mean(0), 6000, r_c_nm=150.0,
                          cfg=E.EnsembleConfig(iterations=800, replicas=20, frames=5))


# ---------------------------------------------------------------- chromosome names and SV files
@pytest.mark.parametrize("raw,norm", [("chr9", "chr9"), ("9", "chr9"), ("Chr9", "chr9"), ("CHRX", "chrX"), ("x", "chrX"),
                                      ("MT", "chrM"), ("chrMT", "chrM"), ("NC_000009.12", "chr9"),
                                      ("NC_000023.11", "chrX"), ("chr1_KI270706v1_random", "chr1_KI270706v1_random")])
def test_chromosome_names_normalise(raw, norm):
    assert genome.normalize_chrom(raw) == norm


def test_unknown_chromosome_name_raises():
    with pytest.raises(KeyError):
        genome.normalize_chrom("scaffold_banana")


def test_vcf_and_bedpe_parsing():
    vcf = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
           "9\t21900000\tdel1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=22100000\n"
           "chr22\t23180000\tbnd1\tN\tN]chr9:130854064]\t.\tPASS\tSVTYPE=BND\n"
           "chr4\t89700000\tdup1\tN\t<DUP:TANDEM>\t.\tPASS\tSVTYPE=DUP;SVLEN=1100000\n"
           "chr2\t5000\tinv1\tN\t<INV>\t.\tPASS\tEND=9000\n"
           "chr7\t100\tcnv\tN\t<CNV>\t.\tPASS\tSVTYPE=CNV;END=200\n"
           "chr7\t100\tshort\n")
    vf = svio.read_any(vcf.encode(), "calls.vcf")
    kinds = [(v.kind, v.chrom, v.start, v.end) for v in vf.variants]
    assert kinds[0] == ("DEL", "chr9", 21_900_000, 22_100_000)
    assert vf.variants[1].kind == "BND" and vf.variants[1].mate_chrom == "chr9" and vf.variants[1].mate_pos == 130_854_063
    assert kinds[2] == ("DUP", "chr4", 89_700_000, 90_800_000) and kinds[3] == ("INV", "chr2", 5000, 9000)
    assert len(vf.skipped) == 2 and "CNV" in vf.skipped[0]
    bedpe = ("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\n"
             "chr9\t21899990\t21900010\tchr9\t22099990\t22100010\tsv1\t0\t+\t-\n"
             "chr9\t1000\t1010\tchr9\t5000\t5010\tsv2\t0\t-\t+\n"
             "9\t1\t2\t22\t5\t6\tsv3\t0\t+\t+\n"
             "chr9\t1\t2\tchr9\t5\t6\tsv4\n")
    bf = svio.read_any(bedpe.encode(), "calls.bedpe")
    assert [v.kind for v in bf.variants] == ["DEL", "DUP", "BND"] and bf.variants[2].mate_chrom == "chr22"
    assert len(bf.skipped) == 1


# ---------------------------------------------------------------- structural variants (piecewise construction)
def _ideal(n=40, c=0.12):
    return c * np.abs(np.subtract.outer(np.arange(n), np.arange(n))).astype(float)


def test_identity_rearrangement_reproduces_the_reference():
    S = _ideal()
    Sd, origin, kind = PT.derive(S, [PT.Piece(np.arange(40), "native", PT.NATIVE)])
    assert np.allclose(Sd, S) and (origin == np.arange(40)).all()


def test_deletion_is_concatenation_for_an_ideal_chain_and_leaves_pieces_unchanged():
    S = _ideal()
    a, b = 15, 25
    pieces, _ = PT.pieces_for("deletion", 40, {"a": a, "b": b})
    Sd, origin, _ = PT.derive(S, pieces)
    ref = dict((int(o), k) for k, o in enumerate(origin))
    sb = np.median(np.diag(S, 1))
    for i, j in [(10, 30), (0, 39), (14, 25)]:
        assert Sd[ref[i], ref[j]] == pytest.approx(S[i, a - 1] + sb + S[b, j], rel=1e-9)
    assert Sd[ref[2], ref[8]] == pytest.approx(S[2, 8]) and Sd[ref[27], ref[35]] == pytest.approx(S[27, 35])


def test_several_deletions_at_once_concatenate_the_kept_runs():
    S = _ideal()
    pieces, desc = PT.pieces_for("deletion", 40, {"segments": [(5, 10), (15, 18)]})
    assert [p.beads.tolist()[:1] + p.beads.tolist()[-1:] for p in pieces] == [[0, 4], [10, 14], [18, 39]]
    assert desc == "deletion of beads 5-9, 15-17"
    Sd, origin, _ = PT.derive(S, pieces)
    kept = np.r_[0:5, 10:15, 18:40]
    assert (origin == kept).all()
    rank = np.abs(np.subtract.outer(np.arange(len(kept)), np.arange(len(kept))))
    assert np.allclose(Sd, 0.12 * rank)                       # an ideal chain simply loses the deleted beads
    single = PT.pieces_for("deletion", 40, {"segments": [(15, 25)]})[0]
    assert [p.beads.tolist() for p in single] == [p.beads.tolist() for p in PT.pieces_for("deletion", 40, {"a": 15, "b": 25})[0]]
    with pytest.raises(ValueError):
        PT.pieces_for("deletion", 40, {"segments": [(0, 40)]})
    with pytest.raises(ValueError):
        PT.pieces_for("deletion", 40, {"segments": [(12, 12)]})


def test_deletion_on_a_fitted_population(chain_result):
    imp = PT.variant_impact(chain_result, "deletion", {"a": 15, "b": 25})
    assert np.isnan(imp.p_after[16, 30]) and np.isfinite(imp.p_after[10, 30])     # deleted beads have no contacts
    assert imp.log2_fc[12, 28] > 1.0                                              # flanks gain contact
    assert abs(imp.log2_fc[2, 8]) < 1e-9                                          # inside one piece: unchanged
    assert imp.top_changes(5)[0]["log2_fc"] > 1.0


def test_inversion_creates_junction_contacts_and_breaks_old_bonds(chain_result):
    imp = PT.variant_impact(chain_result, "inversion", {"a": 10, "b": 30})
    assert imp.log2_fc[9, 29] > 1.5 and imp.log2_fc[10, 30] > 1.5                # new neighbours at both junctions
    assert imp.log2_fc[9, 10] < -1.0 and imp.log2_fc[29, 30] < -1.0              # old bonds broken
    assert abs(imp.log2_fc[12, 20]) < 1e-9                                        # inside the inverted piece


def test_duplication_and_translocation_stay_physical(chain_result):
    dup = PT.variant_impact(chain_result, "duplication", {"a": 10, "b": 20})
    assert np.nanmax(dup.p_after) <= 4.0 + 1e-9                                   # <= 2 x 2 copies summed
    assert np.nanmin(dup.p_after) >= 0.0
    assert dup.p_after[12, 15] > dup.p_before[12, 15]                             # the copy adds contacts
    tra = PT.variant_impact(chain_result, "translocation", {"breakpoint": 25, "partner_beads": 30})
    assert np.isnan(tra.p_after[30, 35]) and np.allclose(tra.p_after[5, 20], tra.p_before[5, 20])
    with pytest.raises(ValueError):
        PT.variant_impact(chain_result, "teleport", {})


def test_enhancer_promoter_pairs_and_bootstrap(chain_result):
    imp = PT.variant_impact(chain_result, "deletion", {"a": 15, "b": 25})
    pairs = PT.enhancer_promoter_pairs(imp, [("GENE1", 12), ("GENE2", 3)], [28, 18, 5])
    genes = {(p["gene"], p["enhancer_bead"], p["change"]) for p in pairs}
    assert ("GENE1", 28, "gained") in genes and ("GENE1", 18, "removed") in genes
    d = _population(n=30, cells=3000, seed=12)
    f = (d < 1.0).mean(0)
    i, j = np.triu_indices(30, 1)
    cnt = np.random.default_rng(0).poisson(f[i, j] * 300).astype(float)
    keep = cnt > 0

    def fit(ci, cj, cm):
        return E.fit_from_counts(ci, cj, cm, 30, b0_nm=50.0, cfg=E.EnsembleConfig(iterations=300, replicas=5, frames=2))

    band = PT.bootstrap_impact(fit, (i[keep], j[keep], cnt[keep]), "deletion", {"a": 10, "b": 15}, reps=3)
    assert band["lo"].shape == (30, 30) and np.nanmean(band["hi"] - band["lo"]) >= 0


# ---------------------------------------------------------------- cohesin depletion map
def test_cohesin_loss_map_shifts_trend_and_shrinks_the_domain_pattern():
    rng = np.random.default_rng(0)
    n = 30
    sep = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) * 30_000
    base = 150.0 * np.maximum(sep / 30_000, 1) ** 0.4
    pattern = np.exp(rng.normal(0, 0.2, (n, n)))
    pattern = np.sqrt(pattern * pattern.T)
    d = base * pattern
    np.fill_diagonal(d, 0.0)
    p = PT.CohesinParams(lam=0.5, c=(0.1, 0.0, 0.0), log_sep_range=(4.0, 7.0))
    out = PT.cohesin_loss_map(d, sep, p)
    iu = np.triu_indices(n, 1)
    _, r0 = PT.trend_residual(np.log(np.where(d > 0, d, 1.0)), sep)
    _, r1 = PT.trend_residual(np.log(np.where(out > 0, out, 1.0)), sep)
    assert np.allclose(r1[iu], 0.5 * r0[iu], atol=1e-9)                        # half the domain pattern lost
    assert np.mean(np.log(out[iu] / d[iu])) == pytest.approx(0.1, abs=0.02)


# ---------------------------------------------------------------- windowed impact and variant files -> scenarios
@pytest.mark.parametrize("op,params", [("deletion", {"a": 15, "b": 25}), ("inversion", {"a": 12, "b": 22}),
                                       ("duplication", {"a": 14, "b": 20}),
                                       ("deletion", {"segments": [(10, 13), (20, 26)]})])
def test_windowed_impact_matches_the_full_window(chain_result, op, params):
    full = PT.variant_impact(chain_result, op, params)
    win = PT.variant_impact(chain_result, op, params, window=(4, 36))
    a, b = full.log2_fc[4:36, 4:36], win.log2_fc
    assert (np.isfinite(a) == np.isfinite(b)).all()
    ok = np.isfinite(a)
    assert np.allclose(a[ok], b[ok], atol=1e-9)


def test_windowed_impact_on_a_v4_model():
    from chronocell import population as P
    d = _population(n=60, cells=4000, seed=5)
    res = P.fit_population((d < 1.0).mean(0), 4000, r_c_nm=150.0,
                           cfg=P.PopulationConfig(iterations=300, replicas=4, frames=2, rank_cap=32))
    full = PT.variant_impact(res, "deletion", {"a": 25, "b": 35})
    win = PT.variant_impact(res, "deletion", {"a": 25, "b": 35}, window=(10, 50))
    a, b = full.log2_fc[10:50, 10:50], win.log2_fc
    ok = np.isfinite(a)
    assert (np.isfinite(b) == ok).all() and np.allclose(a[ok], b[ok], atol=1e-6)
    assert np.nanmax(win.log2_fc) > 1.0                     # flanks of a 10-bead deletion gain contact


def test_variants_map_onto_the_loaded_window():
    text = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr22\t20000000\td1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20500000\n"
            "chr22\t30000000\tv1\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=30200000\n"
            "chr22\t23290000\tb1\tN\tN[chr9:130700000[\t.\tPASS\tSVTYPE=BND\n"
            "chr9\t100\tx\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000\n"
            "chr22\t49000000\tfar\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=49100000\n")
    vf = svio.read_any(text, "pasted")
    assert vf.format == "VCF" and len(vf.variants) == 5
    res, bin0, n = 10_000, 1_500, 2_000                    # window chr22:15.0-35.0 Mb
    got = [svio.to_scenario(v, "chr22", res, bin0, n, partners=genome.main_chromosomes("hg38")) for v in vf.variants]
    assert got[0] == ("deletion", {"a": 500, "b": 550})
    assert got[1] == ("inversion", {"a": 1500, "b": 1520})
    assert got[2][0] == "translocation" and got[2][1]["partner"] == "chr9" and got[2][1]["breakpoint"] == 828      # VCF POS 23,290,000 is 0-based 23,289,999
    assert got[3] == "on chr9, not chr22" and got[4] == "outside the loaded window"
    assert svio.deletion_segments(vf.variants, "chr22", res, bin0, n) == [(500, 550)]
    assert "DEL chr22:20,000,001-20,500,000" in svio.label(vf.variants[0])
    bed = svio.read_any("chr22\t20000000\t20000100\tchr22\t20400000\t20400100\tdelA\t0\t+\t-\n", "pasted")
    assert bed.format == "BEDPE" and svio.to_scenario(bed.variants[0], "chr22", res, bin0, n)[0] == "deletion"
