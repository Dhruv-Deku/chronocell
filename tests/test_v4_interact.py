"""Selection plumbing (ui/interact.py): clicked 3D beads, map pixels and table rows become beads,
each selection is applied once, and beads are labelled with the genes they overlap."""

from __future__ import annotations

from chronocell import genes as G, genome
from ui import interact as I


def test_each_selection_is_applied_once():
    ss = {"viewport": {"selection": {"points": [{"customdata": [2005.0]}]}}}
    assert I.fresh(ss, "viewport") == [{"customdata": [2005.0]}]
    assert I.fresh(ss, "viewport") == []                        # kept by Streamlit between reruns: not re-applied
    ss["viewport"] = {"selection": {"points": [{"customdata": [2007.0]}]}}
    assert len(I.fresh(ss, "viewport")) == 1
    assert I.fresh({}, "missing") == []
    t = {"genes_table": {"selection": {"rows": [3], "columns": []}}}
    assert I.fresh(t, "genes_table", "rows") == [3]
    assert I.fresh(t, "genes_table", "rows") == []


def test_points_become_beads_and_pairs():
    assert I.bead_from_point({"customdata": [2005.0, 1.0, 2.0]}, 2000, 50) == 5
    assert I.bead_from_point({"customdata": 2049}, 2000, 50) == 49
    assert I.bead_from_point({"customdata": [2050]}, 2000, 50) is None        # outside the view
    assert I.bead_from_point({"customdata": [float("nan")]}, 2000, 50) is None
    assert I.bead_from_point({"x": 1.0}, 2000, 50) is None                     # a trace without bins
    res = 10_000
    px = lambda g: (g + 0.5) * res / 1e6                                       # noqa: E731  pixel centre, Mb
    assert I.pair_from_map_point({"x": px(2010), "y": px(2003)}, 2000, 50, res) == (3, 10)
    assert I.pair_from_map_point({"x": px(2010), "y": px(2010)}, 2000, 50, res) is None   # diagonal
    assert I.pair_from_map_point({"x": px(1990), "y": px(2003)}, 2000, 50, res) is None   # outside
    assert I.pair_from_map_point({"z": 1.0}, 2000, 50, res) is None


def test_beads_are_labelled_with_their_genes():
    ch = genome.chrom("chr22", 10_000)
    lab = I.gene_labels(ch, 0, ch.n_bins)
    g = G.find("BCR", ch.assembly)
    assert "BCR" in lab[int(g["tss"]) // ch.resolution]
    assert lab[0] == ""                                                        # 22p: no annotated genes
    some = [x for x in lab if "+" in x]
    assert all(x.count(",") <= 1 for x in some)                                 # at most two names, then '+k'
