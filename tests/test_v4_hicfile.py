"""chronocell.hicfile: the dependency-free .hic reader, on a small file written to the published layout
(versions 8 and 9, list and dense blocks), plus an optional live check against a public Rao et al. 2014
map (set CHRONOCELL_NETWORK_TESTS=1)."""

from __future__ import annotations

import io
import os
import struct
import zlib

import numpy as np
import pytest

from chronocell import hicfile as H


def _cstr(s: str) -> bytes:
    return s.encode() + b"\0"


def _block_list(xs, ys, cs, version: int) -> bytes:
    """Type-1 (list of rows) block, int16 positions and counts."""
    rows: dict[int, list[tuple[int, int]]] = {}
    for x, y, c in zip(xs, ys, cs):
        rows.setdefault(int(y), []).append((int(x), int(c)))
    out = struct.pack("<iii", len(xs), 0, 0)
    out += struct.pack("<bbbb", 0, 0, 0, 1) if version >= 9 else struct.pack("<bb", 0, 1)
    out += struct.pack("<h", len(rows))
    for y, items in sorted(rows.items()):
        out += struct.pack("<hh", y, len(items))
        for x, c in items:
            out += struct.pack("<hh", x, c)
    return zlib.compress(out)


def _block_dense(x0, y0, grid: np.ndarray, version: int) -> bytes:
    """Type-2 (dense) block; -32768 marks an empty cell."""
    h, w = grid.shape
    vals = np.where(grid > 0, grid, -32768).astype("<i2").ravel()
    out = struct.pack("<iii", int((grid > 0).sum()), x0, y0)
    out += struct.pack("<bbbb", 0, 0, 0, 2) if version >= 9 else struct.pack("<bb", 0, 2)
    out += struct.pack("<ih", vals.size, w) + vals.tobytes()
    return zlib.compress(out)


def write_hic(version: int, binsize: int, xs, ys, cs, dense: bool = False) -> bytes:
    chroms = [("ALL", 1000), ("chr21", 2_000_000)]
    head = bytearray(b"HIC\0" + struct.pack("<i", version) + struct.pack("<q", 0) + _cstr("hg38"))
    if version >= 9:
        head += struct.pack("<qq", 0, 0)
    head += struct.pack("<i", 0)
    head += struct.pack("<i", len(chroms))
    for name, length in chroms:
        head += _cstr(name) + (struct.pack("<q", length) if version >= 9 else struct.pack("<i", length))
    head += struct.pack("<i", 1) + struct.pack("<i", binsize) + struct.pack("<i", 0)
    xs, ys, cs = np.asarray(xs), np.asarray(ys), np.asarray(cs)
    if dense:
        n = int(max(xs.max(), ys.max())) + 1
        grid = np.zeros((n, n), int)
        grid[ys, xs] = cs
        block = _block_dense(0, 0, grid, version)
    else:
        block = _block_list(xs, ys, cs, version)
    bbc, bcc = 10_000, 1                       # one block holds every bin
    matrix_pos = len(head)
    zoom = (_cstr("BP") + struct.pack("<i", 0) + struct.pack("<ffff", 0, 0, 0, 0)
            + struct.pack("<iiii", binsize, bbc, bcc, 1))
    block_pos = matrix_pos + 12 + len(zoom) + 16
    zoom += struct.pack("<iqi", 0, block_pos, len(block))
    matrix = struct.pack("<iii", 1, 1, 1) + zoom
    body = bytes(head) + matrix + block
    master = len(body)
    entry = _cstr("1_1") + struct.pack("<qi", matrix_pos, len(matrix))
    footer = (struct.pack("<q", 4 + len(entry)) if version >= 9 else struct.pack("<i", 4 + len(entry)))
    footer += struct.pack("<i", 1) + entry
    data = bytearray(body + footer)
    data[8:16] = struct.pack("<q", master)
    return bytes(data)


@pytest.mark.parametrize("version,dense", [(8, False), (8, True), (9, False), (9, True)])
def test_reader_round_trips_a_spec_file(tmp_path, version, dense):
    rng = np.random.default_rng(version)
    n = 40
    i, j = np.triu_indices(n)                               # upper triangle incl. diagonal, as Juicer stores
    keep = rng.random(len(i)) < 0.6
    xs, ys = j[keep], i[keep]                               # x >= y
    cs = rng.integers(1, 300, keep.sum())
    path = tmp_path / f"t{version}{int(dense)}.hic"
    path.write_bytes(write_hic(version, 5000, xs, ys, cs, dense))
    hf = H.HicFile(str(path))
    assert hf.header.version == version and hf.header.genome == "hg38" and hf.header.bp_resolutions == [5000]
    m = hf.dense("21", 5000, 0, n * 5000)                    # "21" resolves to "chr21"
    truth = np.zeros((n, n))
    truth[ys, xs] = cs
    truth = truth + np.triu(truth, 1).T
    assert np.array_equal(m, truth)
    sub = hf.dense("chr21", 5000, 50_000, 120_000)          # bins 10..23
    assert np.array_equal(sub, truth[10:24, 10:24])
    with pytest.raises(KeyError):
        hf.matrix("chr21", "chr21", 10_000)
    with pytest.raises(KeyError):
        hf.chrom_index("chr7")
    hf.f.close()


def test_read_counts_accepts_bytes_and_drops_the_diagonal():
    data = write_hic(8, 10_000, [0, 1, 3, 3], [0, 0, 1, 3], [5, 7, 2, 9])
    ci, cj, cm = H.read_counts(data, "chr21", 10_000)
    assert sorted(zip(ci.tolist(), cj.tolist(), cm.tolist())) == [(0, 1, 7.0), (1, 3, 2.0)]


def test_not_a_hic_file_is_rejected(tmp_path):
    p = tmp_path / "x.hic"
    p.write_bytes(b"NOTHIC" + b"\0" * 100)
    with pytest.raises(ValueError):
        H.HicFile(str(p))


@pytest.mark.skipif(not os.environ.get("CHRONOCELL_NETWORK_TESTS"), reason="set CHRONOCELL_NETWORK_TESTS=1 for live checks")
def test_live_rao2014_imr90_region():
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_IMR90_combined_30.hic"
    hf = H.HicFile(url)
    assert hf.header.genome == "hg19" and 5000 in hf.header.bp_resolutions
    m = hf.dense("21", 5000, 29_370_000, 31_325_000)
    assert m.shape == (391, 391) and np.allclose(m, m.T) and np.median(np.diag(m, 1)) > 10


def test_ingest_reads_hic_without_hic_straw():
    from chronocell import genome, ingest
    data = write_hic(8, 5000, [1, 2, 7], [0, 1, 2], [4, 6, 3])
    ci, cj, cm, notes = ingest.read_contacts(data, "patient.hic", genome.chrom("chr21", 5000))
    assert sorted(zip(ci.tolist(), cj.tolist(), cm.tolist())) == [(0, 1, 4.0), (1, 2, 6.0), (2, 7, 3.0)]
    assert "built-in" in notes[0] or "hic" in notes[0]
