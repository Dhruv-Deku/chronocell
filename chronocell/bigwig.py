"""
A small pure-Python bigWig reader (Phase B5), for ChIP-seq / ATAC-seq signal tracks when the compiled pyBigWig
package is not available (it does not install on Windows without a compiler). Local files or URLs (HTTP range
requests through chronocell.hicfile's range reader): only the index and the data blocks of the requested region
are read.

Format (Kent et al., Bioinformatics 26:2204, 2010): a header, a B+ tree of chromosome names, an R-tree index of
data blocks, and zlib-compressed blocks of bedGraph, variableStep or fixedStep records. Zoom levels are not used:
values come from the full-resolution data.
"""

from __future__ import annotations

import struct
import zlib

import numpy as np

BIGWIG_MAGIC = 0x888FFC26
CHROM_TREE_MAGIC = 0x78CA8C91
CIR_TREE_MAGIC = 0x2468ACE0


class BigWig:
    def __init__(self, src: str):
        from .hicfile import _RangeFile
        self.f = _RangeFile(src)
        head = self._read(0, 64)
        magic = struct.unpack("<I", head[:4])[0]
        if magic != BIGWIG_MAGIC:
            if struct.unpack(">I", head[:4])[0] == BIGWIG_MAGIC:
                raise ValueError("Big-endian bigWig files are not supported.")
            raise ValueError("Not a bigWig file (bad magic number).")
        (_, self.version, self.zoom_levels, self.chrom_tree, self.full_data, self.full_index, _, _, _, _,
         self.uncompress_buf) = struct.unpack("<IHHQQQHHQQI", head[:56])
        self.chroms = self._chrom_tree()

    def _read(self, pos: int, n: int) -> bytes:
        self.f.seek(pos)
        return self.f.read(n)

    def _chrom_tree(self) -> dict[str, tuple[int, int]]:
        h = self._read(self.chrom_tree, 32)
        magic, block, key_size, val_size, items = struct.unpack("<IIIIQ", h[:24])
        if magic != CHROM_TREE_MAGIC:
            raise ValueError("Corrupt bigWig: bad chromosome tree.")
        out: dict[str, tuple[int, int]] = {}

        def node(pos: int) -> None:
            is_leaf, _, count = struct.unpack("<BBH", self._read(pos, 4))
            size = key_size + (8 if is_leaf else 8)
            data = self._read(pos + 4, count * size)
            for k in range(count):
                rec = data[k * size:(k + 1) * size]
                key = rec[:key_size].rstrip(b"\x00").decode()
                if is_leaf:
                    cid, csize = struct.unpack("<II", rec[key_size:key_size + 8])
                    out[key] = (cid, csize)
                else:
                    node(struct.unpack("<Q", rec[key_size:key_size + 8])[0])
        node(self.chrom_tree + 32)
        return out

    def _chrom_id(self, chrom: str) -> int:
        for cand in (chrom, chrom.removeprefix("chr"), "chr" + chrom.removeprefix("chr")):
            if cand in self.chroms:
                return self.chroms[cand][0]
        raise KeyError(f"{chrom} is not in this bigWig ({', '.join(list(self.chroms)[:8])}...).")

    def _blocks(self, cid: int, start: int, end: int) -> list[tuple[int, int]]:
        h = self._read(self.full_index, 48)
        if struct.unpack("<I", h[:4])[0] != CIR_TREE_MAGIC:
            raise ValueError("Corrupt bigWig: bad data index.")
        out = []

        def overlaps(sc, sb, ec, eb) -> bool:
            return (sc, sb) < (cid, end) and (ec, eb) > (cid, start)

        def node(pos: int) -> None:
            is_leaf, _, count = struct.unpack("<BBH", self._read(pos, 4))
            size = 32 if is_leaf else 24
            data = self._read(pos + 4, count * size)
            for k in range(count):
                rec = data[k * size:(k + 1) * size]
                sc, sb, ec, eb = struct.unpack("<IIII", rec[:16])
                if not overlaps(sc, sb, ec, eb):
                    continue
                if is_leaf:
                    out.append(struct.unpack("<QQ", rec[16:32]))
                else:
                    node(struct.unpack("<Q", rec[16:24])[0])
        node(self.full_index + 48)
        return out

    def intervals(self, chrom: str, start: int, end: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(starts, ends, values) of every record overlapping [start, end)."""
        cid = self._chrom_id(chrom)
        S, E, V = [], [], []
        for off, size in self._blocks(cid, start, end):
            raw = self._read(off, size)
            if self.uncompress_buf:
                raw = zlib.decompress(raw)
            chrom_id, s0, _, step, span, typ, _, n = struct.unpack("<IIIIIBBH", raw[:24])
            if chrom_id != cid:
                continue
            body = raw[24:]
            if typ == 1:
                rec = np.frombuffer(body[:12 * n], dtype=np.dtype([("s", "<u4"), ("e", "<u4"), ("v", "<f4")]))
                s, e, v = rec["s"].astype(np.int64), rec["e"].astype(np.int64), rec["v"].astype(float)
            elif typ == 2:
                rec = np.frombuffer(body[:8 * n], dtype=np.dtype([("s", "<u4"), ("v", "<f4")]))
                s = rec["s"].astype(np.int64)
                e, v = s + span, rec["v"].astype(float)
            elif typ == 3:
                v = np.frombuffer(body[:4 * n], dtype="<f4").astype(float)
                s = s0 + step * np.arange(n, dtype=np.int64)
                e = s + span
            else:
                raise ValueError(f"Unknown bigWig section type {typ}.")
            keep = (e > start) & (s < end)
            S.append(s[keep])
            E.append(e[keep])
            V.append(v[keep])
        if not S:
            return np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0)
        return np.concatenate(S), np.concatenate(E), np.concatenate(V)

    def bin_means(self, chrom: str, start: int, end: int, binsize: int) -> np.ndarray:
        """Mean value per bin over the bases that have data (NaN where none), as pyBigWig stats(type='mean')."""
        nb = -(-(end - start) // binsize)
        s, e, v = self.intervals(chrom, start, end)
        tot = np.zeros(nb)
        cov = np.zeros(nb)
        for a, b, x in zip(np.maximum(s, start), np.minimum(e, end), v):
            i0, i1 = (a - start) // binsize, (b - 1 - start) // binsize
            for i in range(int(i0), int(i1) + 1):
                lo, hi = max(a, start + i * binsize), min(b, start + (i + 1) * binsize)
                if hi > lo:
                    tot[i] += x * (hi - lo)
                    cov[i] += hi - lo
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(cov > 0, tot / cov, np.nan)


def write_bigwig_bedgraph(path, chroms: list[tuple[str, int]], records: dict[str, list[tuple[int, int, float]]]) -> None:
    """A minimal bigWig writer (one zlib-compressed bedGraph block per chromosome, no zoom levels),
    used to test the reader and to export small tracks."""
    import io
    buf = io.BytesIO()
    names = [c for c, _ in chroms]
    key_size = max(len(c) for c in names)
    header_size, chrom_tree_pos = 64, 64
    # chromosome B+ tree: one leaf node
    tree = struct.pack("<IIIIQQ", CHROM_TREE_MAGIC, len(names), key_size, 8, len(names), 0)
    tree += struct.pack("<BBH", 1, 0, len(names))
    for k, (c, L) in enumerate(chroms):
        tree += c.encode().ljust(key_size, b"\x00") + struct.pack("<II", k, L)
    data_pos = chrom_tree_pos + len(tree)
    blocks, body = [], b""
    pos = data_pos + 8
    for k, c in enumerate(names):
        recs = sorted(records.get(c, []))
        if not recs:
            continue
        sec = struct.pack("<IIIIIBBH", k, recs[0][0], recs[-1][1], 0, 0, 1, 0, len(recs))
        sec += b"".join(struct.pack("<IIf", s, e, v) for s, e, v in recs)
        comp = zlib.compress(sec)
        blocks.append((k, recs[0][0], k, recs[-1][1], pos, len(comp), len(sec)))
        body += comp
        pos += len(comp)
    index_pos = pos
    index = struct.pack("<IIQIIIIQII", CIR_TREE_MAGIC, max(1, len(blocks)), len(blocks), blocks[0][0], blocks[0][1],
                        blocks[-1][2], blocks[-1][3], index_pos, 1, 0)
    index += struct.pack("<BBH", 1, 0, len(blocks))
    for b in blocks:
        index += struct.pack("<IIIIQQ", *b[:6])
    unc = max(b[6] for b in blocks)
    header = struct.pack("<IHHQQQHHQQIQ", BIGWIG_MAGIC, 4, 0, chrom_tree_pos, data_pos, index_pos, 0, 0, 0, 0, unc, 0)
    buf.write(header.ljust(header_size, b"\x00"))
    buf.write(tree)
    buf.write(struct.pack("<Q", len(blocks)))
    buf.write(body)
    buf.write(index)
    with open(path, "wb") as fh:
        fh.write(buf.getvalue())
