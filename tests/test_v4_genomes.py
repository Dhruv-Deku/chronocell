"""v4 Pillar 6: genomes as configuration (chronocell/data/genomes/<assembly>/manifest.json), every main
chromosome of hg38 and mm39, chromosome-name aliases, and the app switching assembly."""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pytest

from chronocell import formats, genes as G, genome, ingest, physics, synthetic, viz


def test_assemblies_are_discovered_from_manifests():
    a = genome.assemblies()
    assert {"hg38", "mm39"} <= set(a)
    for asm in a.values():
        assert asm.licence and asm.citation and asm.main_chromosomes and asm.path("chromosomes").exists()
    assert genome.assembly().name == "hg38" and genome.assembly("GRCm39").name == "mm39"
    assert genome.main_chromosomes("hg38") == genome.MAIN_CHROMOSOMES
    with pytest.raises(KeyError):
        genome.assembly("xenopus99")


@pytest.mark.parametrize("asm", ["hg38", "mm39"])
def test_every_main_chromosome_loads_with_a_bounded_default_resolution(asm):
    for name in genome.main_chromosomes(asm):
        ch = genome.chrom(name, None, asm)
        assert ch.assembly == asm and ch.size > 0 and 0 < ch.n_bins <= genome.MAX_DEFAULT_BEADS
        assert ch.assembled_mask().any() and len(ch.bands) >= 1


@pytest.mark.parametrize("asm,raw,ucsc", [
    ("hg38", "NC_000022.11", "chr22"), ("hg38", "CM000684.2", "chr22"), ("hg38", "22", "chr22"),
    ("hg38", "chrx", "chrX"), ("hg38", "MT", "chrM"),
    ("mm39", "CM001012.3", "chr19"), ("mm39", "NC_000085.7", "chr19"), ("mm39", "19", "chr19"), ("mm39", "X", "chrX"),
])
def test_chromosome_aliases_per_assembly(asm, raw, ucsc):
    assert genome.normalize_chrom(raw, asm) == ucsc
    if ucsc in genome.main_chromosomes(asm):                 # chrM is named but not modelled
        assert genome.chrom(raw, None, asm).name == ucsc


def test_mouse_genes_structures_and_files():
    ch = genome.chrom("chr19", 40_000, "mm39")
    assert len(G.table("mm39")) > 15_000 and G.find("Sox17", "mm39")["chrom"] == "chr1"
    ref = synthetic.build(ch, seed=3)
    tab = G.accessibility_table(ref.coords, ref.epi, ref.valid, physics.bond_length_for(ch.resolution), ch, 0)
    assert len(tab) > 100 and set(tab["chrom"] if "chrom" in tab else ["chr19"]) == {"chr19"}
    pdb, _ = formats.write_pdb(ref.coords[:200], 0, ref.gc, ref.epi, 10.0, source="t", method="t", chrom=ch)
    title = next(line for line in pdb.splitlines() if line.startswith("TITLE"))
    assert "MOUSE CHR19 (GRCM39)" in title and len(title) <= 80 and formats.validate_pdb(pdb).ok
    buf = io.BytesIO()
    formats.write_bundle(buf, formats.StructureBundle("chr19", ch.resolution, ref.coords[None], np.zeros(1), ["t"],
                                                      assembly="mm39"))
    b, _ = formats.read_bundle(buf.getvalue(), "x.npz")
    assert b.assembly == "mm39" and b.chrom == "chr19"
    with pytest.raises(ValueError):
        formats.StructureBundle("chr22", 10_000, ref.coords[None, :100], np.zeros(1), ["t"], assembly="mm39").validate()


def test_tracks_with_any_chromosome_spelling_land_on_the_right_bins():
    ch = genome.chrom("chr19", 40_000, "mm39")
    bg = "\n".join(f"NC_000085.7\t{s}\t{s + 40_000}\t3.0" for s in range(0, 400_000, 40_000)).encode()
    arr, _ = ingest.read_track(bg, "t.bedGraph", ch)
    assert np.allclose(arr[:10], 3.0) and np.isnan(arr[11])
    bg38 = "\n".join(f"22\t{s}\t{s + 10_000}\t2.0" for s in range(0, 50_000, 10_000)).encode()
    arr38, _ = ingest.read_track(bg38, "t.bedGraph", genome.chrom("chr22"))
    assert np.allclose(arr38[:5], 2.0)


def test_track_chart_names_the_chromosome_shown():
    ch = genome.chrom("chr7")
    fig = viz.tracks_chart(np.arange(10), np.full(10, 0.4), np.ones(10), [], chrom=ch)
    titles = [ax.title.text for ax in (fig.layout.xaxis, fig.layout.xaxis2) if ax.title.text]
    assert any("chr7" in t for t in titles) and not any("chr22" in t for t in titles)


def test_app_switches_to_the_mouse_assembly(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    import ui.common as C
    import ui.states_panel as SP
    monkeypatch.setattr(C, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "SLOT_ROOT", tmp_path / "empty")
    monkeypatch.setattr(SP, "DEMO_ROOT", tmp_path / "demo")
    at = AppTest.from_file(str(Path(__file__).resolve().parent.parent / "app.py"), default_timeout=900)
    at.run()
    assert not at.exception
    text = lambda: "".join(m.value for m in at.markdown)  # noqa: E731
    assert "GRCh38 · chr22" in text()                         # unchanged default
    at.selectbox(key="assembly_choice").set_value("mm39").run()
    assert not at.exception, [e.value for e in at.exception]
    assert at.session_state["chrom_choice"] == "chr19" and "GRCm39 · chr19" in text() and "mouse chr19" in text()
    at.segmented_control(key="workspace").set_value("Genes").run()
    assert not at.exception, [e.value for e in at.exception]
