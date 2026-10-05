"""
A small, dependency-free reader for Juicer .hic files (format versions 6-9), local or remote.

Only what ChronoCell needs: the header (genome, chromosomes, resolutions), the master index, and the
observed (raw) counts of one chromosome pair at one base-pair resolution inside a region. Remote
files (http/https) are read with HTTP range requests in cached chunks, so a region of a 10-40 GB
public map downloads only the blocks that cover it.

The layout follows the published specification (Durand et al., Cell Systems 2016; "HiC Format
Specification" in github.com/aidenlab/hic-format) and the reference reader `straw` (Aiden lab):
  header  : "HIC\\0", version, master-index pointer, genome id, [v9: NVI pointer + length],
            attributes, chromosomes (name, length), BP resolutions, FRAG resolutions
  footer  : at the master pointer: byte count, entries "c1_c2" -> (file position, size)
  matrix  : c1, c2, resolutions; per resolution: unit, sums, bin size, block bin/column counts,
            block index (number, file position, byte size)
  block   : zlib-compressed records (v6: (x, y, count) triples; v7+: offsets and either a list of
            rows or a dense grid; v9 adds integer-width flags)
Normalisation vectors are not read: ChronoCell uses raw counts (and ICE in chronocell.normalize).
"""

from __future__ import annotations

import io
import struct
import urllib.request
import zlib
from dataclasses import dataclass, field

import numpy as np

UA = {"User-Agent": "ChronoCell-5D (research; hic reader)"}


class _RangeFile:
    """Read-only, seekable view of a local file or an http(s) URL (range requests, cached chunks)."""

    def __init__(self, src: str, chunk: int = 1 << 18):
        self.src = src
        self.remote = src.startswith(("http://", "https://"))
        self.pos = 0
        self.chunk = chunk
        self._cache: dict[int, bytes] = {}
        self.bytes_fetched = 0
        self._fh = None if self.remote else open(src, "rb")

    def seek(self, pos: int) -> None:
        self.pos = int(pos)

    def tell(self) -> int:
        return self.pos

    def _fetch(self, start: int, end: int) -> bytes:
        import http.client
        import time
        req = urllib.request.Request(self.src, headers={**UA, "Range": f"bytes={start}-{end - 1}"})
        last = None
        for attempt in range(6):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = r.read()
                self.bytes_fetched += len(data)
                return data
            except (OSError, http.client.HTTPException) as exc:     # transient network / server errors (incl. a
                last = exc                                         # connection dropped mid-read): retry with backoff
                time.sleep(min(30.0, 2.0 * 2 ** attempt))
        raise OSError(f"range request failed for {self.src} [{start}, {end}): {last}")

    def read(self, n: int) -> bytes:
        if not self.remote:
            self._fh.seek(self.pos)
            data = self._fh.read(n)
            self.pos += len(data)
            return data
        if n > 4 * self.chunk:                  # big reads (blocks) go straight through
            data = self._fetch(self.pos, self.pos + n)
            self.pos += len(data)
            return data
        out = bytearray()
        while n > 0:
            k = self.pos // self.chunk
            if k not in self._cache:
                self._cache[k] = self._fetch(k * self.chunk, (k + 1) * self.chunk)
            buf = self._cache[k]
            off = self.pos - k * self.chunk
            take = buf[off:off + n]
            if not take:
                break
            out += take
            self.pos += len(take)
            n -= len(take)
        return bytes(out)

    def close(self) -> None:
        if self._fh:
            self._fh.close()


def _i32(f) -> int:
    return struct.unpack("<i", f.read(4))[0]


def _i64(f) -> int:
    return struct.unpack("<q", f.read(8))[0]


def _f32(f) -> float:
    return struct.unpack("<f", f.read(4))[0]


def _cstr(f) -> str:
    out = bytearray()
    while True:
        c = f.read(1)
        if not c or c == b"\0":
            return out.decode("utf-8", "replace")
        out += c


@dataclass
class HicHeader:
    version: int
    genome: str
    master: int
    chromosomes: list[tuple[str, int]]
    bp_resolutions: list[int]
    attributes: dict[str, str] = field(default_factory=dict)


class HicFile:
    """Raw-count access to a .hic map. `HicFile(path_or_url).matrix("chr21", "chr21", 5000, region)`."""

    def __init__(self, src: str):
        self.f = _RangeFile(src)
        self.header = self._read_header()
        self._index = None

    # ---- header and master index -----------------------------------------------------------
    def _read_header(self) -> HicHeader:
        f = self.f
        f.seek(0)
        if f.read(4)[:3] != b"HIC":
            raise ValueError("Not a .hic file (missing HIC magic).")
        version = _i32(f)
        if version < 6:
            raise ValueError(f".hic version {version} is not supported (6-9 are).")
        master = _i64(f)
        genome = _cstr(f)
        if version >= 9:
            _i64(f)                             # normalisation-vector index position
            _i64(f)                             # and length
        attrs = {}
        for _ in range(_i32(f)):
            k = _cstr(f)
            attrs[k] = _cstr(f)
        chroms = []
        for _ in range(_i32(f)):
            name = _cstr(f)
            length = _i64(f) if version >= 9 else _i32(f)
            chroms.append((name, length))
        res = [_i32(f) for _ in range(_i32(f))]
        return HicHeader(version, genome, master, chroms, res, {k: v[:200] for k, v in attrs.items()})

    def _master_index(self) -> dict[str, tuple[int, int]]:
        if self._index is None:
            f = self.f
            f.seek(self.header.master)
            _i64(f) if self.header.version >= 9 else _i32(f)        # byte count of the footer
            idx = {}
            for _ in range(_i32(f)):
                key = _cstr(f)
                idx[key] = (_i64(f), _i32(f))
            self._index = idx
        return self._index

    def chrom_index(self, name: str) -> int:
        names = [c for c, _ in self.header.chromosomes]
        for cand in (name, name.removeprefix("chr"), "chr" + name.removeprefix("chr"),
                     "MT" if name in ("chrM", "M") else name):
            if cand in names:
                return names.index(cand)
        raise KeyError(f"{name} is not in this .hic file ({', '.join(names[:8])}...).")

    # ---- one matrix, one resolution -------------------------------------------------------
    def _zoom(self, i1: int, i2: int, binsize: int):
        key = f"{i1}_{i2}"
        idx = self._master_index()
        if key not in idx:
            raise KeyError(f"No matrix {key} in this .hic file.")
        f = self.f
        f.seek(idx[key][0])
        _i32(f), _i32(f)
        for _ in range(_i32(f)):                                     # resolutions of this matrix
            unit = _cstr(f)
            _i32(f)                                                  # zoom index
            for _k in range(4):
                _f32(f)                                              # sumCounts, occupied cells, sd, p95
            bs, bbc, bcc = _i32(f), _i32(f), _i32(f)
            nb = _i32(f)
            raw = f.read(16 * nb)
            if unit == "BP" and bs == binsize:
                arr = np.frombuffer(raw, dtype=np.dtype([("n", "<i4"), ("pos", "<i8"), ("size", "<i4")]))
                return bbc, bcc, {int(r["n"]): (int(r["pos"]), int(r["size"])) for r in arr}
        raise KeyError(f"Resolution {binsize:,} bp is not in this .hic file ({self.header.bp_resolutions}).")

    def _block(self, pos: int, size: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        self.f.seek(pos)
        b = zlib.decompress(self.f.read(size))
        v = self.header.version
        n_rec = struct.unpack_from("<i", b, 0)[0]
        if v < 7:
            rec = np.frombuffer(b, dtype=np.dtype([("x", "<i4"), ("y", "<i4"), ("c", "<f4")]), count=n_rec, offset=4)
            return rec["x"].astype(np.int64), rec["y"].astype(np.int64), rec["c"].astype(np.float64)
        x0, y0 = struct.unpack_from("<ii", b, 4)
        if v < 9:
            use_float = struct.unpack_from("<b", b, 12)[0] != 0          # 0: counts stored as int16
            int_x = int_y = False
            kind = struct.unpack_from("<b", b, 13)[0]
            p = 14
        else:
            use_float = struct.unpack_from("<b", b, 12)[0] != 0
            int_x = struct.unpack_from("<b", b, 13)[0] != 0
            int_y = struct.unpack_from("<b", b, 14)[0] != 0
            kind = struct.unpack_from("<b", b, 15)[0]
            p = 16
        xs, ys, cs = [], [], []
        if kind == 1:                                                 # list of rows
            fy, sy = ("<i", 4) if int_y else ("<h", 2)
            fx, sx = ("<i", 4) if int_x else ("<h", 2)
            fc, sc = ("<f", 4) if use_float else ("<h", 2)
            rows = struct.unpack_from(fy, b, p)[0]
            p += sy
            for _ in range(rows):
                y = struct.unpack_from(fy, b, p)[0]
                p += sy
                cols = struct.unpack_from(fx, b, p)[0]
                p += sx
                dt = np.dtype([("x", fx), ("c", fc)])
                rec = np.frombuffer(b, dtype=dt, count=cols, offset=p)
                p += cols * dt.itemsize
                xs.append(x0 + rec["x"].astype(np.int64))
                ys.append(np.full(cols, y0 + y, dtype=np.int64))
                cs.append(rec["c"].astype(np.float64))
        elif kind == 2:                                               # dense grid
            n_pts = struct.unpack_from("<i", b, p)[0]
            w = struct.unpack_from("<h", b, p + 4)[0]
            p += 6
            vals = np.frombuffer(b, dtype="<f4" if use_float else "<i2", count=n_pts, offset=p).astype(np.float64)
            k = np.arange(n_pts)
            good = ~np.isnan(vals) if use_float else vals != -32768
            xs.append(x0 + (k % w)[good])
            ys.append(y0 + (k // w)[good])
            cs.append(vals[good])
        else:
            raise ValueError(f"Unknown .hic block type {kind}.")
        if not xs:
            return np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0)
        return np.concatenate(xs), np.concatenate(ys), np.concatenate(cs)

    def matrix(self, chrom1: str, chrom2: str, binsize: int, start1: int = 0, end1: int | None = None,
               start2: int | None = None, end2: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Observed counts (bin1, bin2, count) with bin = position // binsize, inside the region.
        For intra-chromosomal maps each pair is returned once (bin1 <= bin2)."""
        i1, i2 = self.chrom_index(chrom1), self.chrom_index(chrom2)
        flip = i1 > i2
        if flip:
            i1, i2 = i2, i1
            start1, end1, start2, end2 = start2, end2, start1, end1
        len1 = self.header.chromosomes[i1][1]
        len2 = self.header.chromosomes[i2][1]
        start2 = start1 if start2 is None and i1 == i2 else (start2 or 0)
        end1 = len1 if end1 is None else end1
        end2 = (end1 if i1 == i2 else len2) if end2 is None else end2
        b1lo, b1hi = start1 // binsize, (end1 - 1) // binsize
        b2lo, b2hi = start2 // binsize, (end2 - 1) // binsize
        bbc, bcc, blocks = self._zoom(i1, i2, binsize)
        intra = i1 == i2
        if self.header.version >= 9 and intra:
            # v9 intra maps index blocks by (depth = log2 distance from the diagonal, position along it),
            # as in straw's getBlockNumbersForRegionFromBinPosition
            lo_pad = (b1lo + b2lo) // 2 // bbc
            hi_pad = (b1hi + b2hi) // 2 // bbc + 1
            near = int(np.log2(1 + abs(b1lo - b2hi) / np.sqrt(2) / bbc))
            far = int(np.log2(1 + abs(b1hi - b2lo) / np.sqrt(2) / bbc))
            d_lo = min(near, far)
            if (b1lo > b2hi and b1hi < b2lo) or (b1hi > b2lo and b1lo < b2hi):
                d_lo = 0
            d_hi = max(near, far) + 1
            want = {d * bcc + p for d in range(d_lo, d_hi + 1) for p in range(lo_pad, hi_pad + 1)}
        else:
            def numbers(alo, ahi, blo, bhi):
                return {r * bcc + c for r in range(blo // bbc, bhi // bbc + 1) for c in range(alo // bbc, ahi // bbc + 1)}
            want = numbers(b1lo, b1hi, b2lo, b2hi)
            if intra:
                want |= numbers(b2lo, b2hi, b1lo, b1hi)
        X, Y, C = [], [], []
        for num in sorted(want):
            if num not in blocks:
                continue
            x, y, c = self._block(*blocks[num])
            X.append(x)
            Y.append(y)
            C.append(c)
        if not X:
            return np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0)
        x, y, c = np.concatenate(X), np.concatenate(Y), np.concatenate(C)
        if intra:
            a, b = np.minimum(x, y), np.maximum(x, y)
            keep = (b >= b1lo) & (a <= b1hi) & (a >= min(b1lo, b2lo)) & (b <= max(b1hi, b2hi))
            keep &= ((a >= b1lo) & (a <= b1hi) & (b >= b2lo) & (b <= b2hi)) | ((a >= b2lo) & (a <= b2hi) & (b >= b1lo) & (b <= b1hi))
            x, y, c = a[keep], b[keep], c[keep]
        else:
            keep = (x >= b1lo) & (x <= b1hi) & (y >= b2lo) & (y <= b2hi)
            x, y, c = x[keep], y[keep], c[keep]
        if flip:
            x, y = y, x
        return x, y, c

    def dense(self, chrom: str, binsize: int, start: int, end: int) -> np.ndarray:
        """Symmetric dense matrix of raw counts for one region of one chromosome."""
        x, y, c = self.matrix(chrom, chrom, binsize, start, end)
        lo = start // binsize
        n = (end - 1) // binsize - lo + 1
        m = np.zeros((n, n))
        np.add.at(m, (x - lo, y - lo), c)
        off = x != y
        np.add.at(m, (y[off] - lo, x[off] - lo), c[off])
        return m


def read_counts(src: str | bytes, chrom: str, binsize: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Whole-chromosome raw counts (i < j) at `binsize` from a path, URL or the file's bytes."""
    if isinstance(src, (bytes, bytearray)):
        import os
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".hic", delete=False) as fh:
            fh.write(src)
            path = fh.name
        try:
            hf = HicFile(path)
            x, y, c = hf.matrix(chrom, chrom, binsize)
            hf.f.close()
        finally:
            os.unlink(path)
    else:
        hf = HicFile(src)
        x, y, c = hf.matrix(chrom, chrom, binsize)
    keep = x != y
    return x[keep], y[keep], c[keep]
