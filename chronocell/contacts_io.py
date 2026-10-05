"""
Region-wise contact reading from files on disk or URLs (Phase B4), without loading whole matrices.

    .hic                 chronocell.hicfile (versions 6-9; local or HTTP range requests): only the blocks
                         covering the region are read
    .cool / .mcool       h5py, path-based: only the pixel rows of the region (bin1_offset index) are read;
                         an .mcool uses the requested resolution (or the finest one that divides it, summed)
    .pairs / .pairs.gz   4DN pairs format (header "#columns: readID chrom1 pos1 chrom2 pos2 ..."), streamed in
                         chunks; pairs are binned at the requested resolution on the fly

Every reader returns `Contacts` (region-local bins ci <= cj, counts, n, notes). Checks with plain-language
errors: the file type, the chromosome name (spellings are normalised; a missing chromosome lists the ones
present), the resolution (lists the ones available), the assembly (chromosome lengths are compared with known
assemblies: hg19, hg38, mm10, mm39) and corrupt files.
Normalisation: "none" (raw counts, the default where validated settings use them), "ice" or "kr"
(chronocell.normalize), returned as balancing weights so raw counts stay available.
"""

from __future__ import annotations

import gzip
import io
import os
from dataclasses import dataclass, field

import numpy as np

# chr1 and chrX lengths of the assemblies we can recognise (UCSC chrom.sizes)
KNOWN_LENGTHS = {"hg19": {"chr1": 249_250_621, "chrX": 155_270_560}, "hg38": {"chr1": 248_956_422, "chrX": 156_040_895},
                 "mm10": {"chr1": 195_471_971, "chrX": 171_031_299}, "mm39": {"chr1": 195_154_279, "chrX": 169_476_592}}
CHUNK_LINES = 1_000_000


class ContactFileError(ValueError):
    """A contact file could not be used; the message says why and what to do."""


@dataclass
class Contacts:
    ci: np.ndarray
    cj: np.ndarray
    cm: np.ndarray
    n: int
    chrom: str
    start: int
    resolution: int
    weights: np.ndarray | None = None
    notes: list[str] = field(default_factory=list)

    def dense(self) -> np.ndarray:
        m = np.zeros((self.n, self.n))
        np.add.at(m, (self.ci, self.cj), self.cm)
        off = self.ci != self.cj
        np.add.at(m, (self.cj[off], self.ci[off]), self.cm[off])
        return m


def kind_of(path: str) -> str:
    low = path.lower().split("?")[0]
    if low.endswith(".hic"):
        return "hic"
    if low.endswith(".mcool"):
        return "mcool"
    if low.endswith(".cool"):
        return "cool"
    if low.endswith((".pairs", ".pairs.gz", ".pairsam", ".pairsam.gz")):
        return "pairs"
    raise ContactFileError(f"'{os.path.basename(path)}' is not a contact file this reader knows "
                           "(.hic, .cool, .mcool, .pairs, .pairs.gz).")


def _norm(name: str) -> str:
    n = str(name).strip()
    if n.lower().startswith("chr"):
        return "chr" + n[3:]
    if n in ("MT", "M"):
        return "chrM"
    return "chr" + n


def detect_assembly(lengths: dict[str, int]) -> str | None:
    norm = {_norm(k): int(v) for k, v in lengths.items()}
    for asm, ref in KNOWN_LENGTHS.items():
        if any(c in norm for c in ref) and all(norm.get(c) == L for c, L in ref.items() if c in norm):
            return asm
    return None


def check_assembly(lengths: dict[str, int], expected: str | None) -> list[str]:
    """Notes about the file's assembly; raises if it is a known assembly different from `expected`."""
    if not lengths:
        return ["No chromosome lengths in the file: assembly not checked."]
    found = detect_assembly(lengths)
    if expected and found and found != expected:
        raise ContactFileError(f"This file is {found} (its chromosome lengths match {found}), but the session uses "
                               f"{expected}. Lift the coordinates over or load the matching assembly.")
    if found is None:
        return ["Assembly not recognised from chromosome lengths (not hg19, hg38, mm10 or mm39)."]
    return [f"Assembly {found} (from chromosome lengths)."]


def _find_chrom(names: list[str], chrom: str) -> int:
    norm = [_norm(x) for x in names]
    want = _norm(chrom)
    if want not in norm:
        raise ContactFileError(f"{chrom} is not in this file. It has: {', '.join(names[:12])}"
                               + ("..." if len(names) > 12 else "") + ".")
    return norm.index(want)


def _h5(path: str):
    try:
        import h5py
    except ImportError as exc:
        raise ContactFileError("Reading .cool / .mcool needs h5py (pip install h5py).") from exc
    try:
        return h5py.File(path, "r")
    except OSError as exc:
        raise ContactFileError(f"'{os.path.basename(path)}' is not a readable HDF5 cooler file "
                               f"(corrupt or incomplete download?): {exc}") from exc


def _cool_group(f, resolution: int | None):
    if "resolutions" in f:
        avail = sorted(int(r) for r in f["resolutions"].keys())
        if resolution is None:
            return f["resolutions"][str(avail[0])]
        if resolution in avail:
            return f["resolutions"][str(resolution)]
        finer = [r for r in avail if resolution % r == 0]
        if not finer:
            raise ContactFileError(f"The .mcool has resolutions {', '.join(f'{r:,}' for r in avail)} bp; none equals or "
                                   f"divides {resolution:,} bp.")
        return f["resolutions"][str(max(finer))]
    if "pixels" not in f or "bins" not in f:
        raise ContactFileError("Not a cooler file (no 'pixels' / 'bins' groups).")
    if resolution is not None:
        binsize = int(f.attrs.get("bin-size", 0) or 0)
        if binsize and resolution % binsize:
            raise ContactFileError(f"The .cool is at {binsize:,} bp, which does not divide {resolution:,} bp.")
    return f


def _open_text(path: str):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _pairs_header(path: str) -> list[str]:
    out = []
    with _open_text(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            out.append(line.rstrip("\n"))
    return out


def chromosomes(path: str) -> dict[str, int]:
    k = kind_of(path)
    if k == "hic":
        from .hicfile import HicFile
        try:
            return {c: int(L) for c, L in HicFile(path).header.chromosomes if c.upper() != "ALL"}
        except (ValueError, OSError) as exc:
            raise ContactFileError(f"Could not read the .hic header: {exc}") from exc
    if k in ("cool", "mcool"):
        with _h5(path) as f:
            g = _cool_group(f, None)
            names = [x.decode() if isinstance(x, bytes) else str(x) for x in g["chroms"]["name"][:]]
            return dict(zip(names, (int(v) for v in g["chroms"]["length"][:])))
    out = {}
    for line in _pairs_header(path):
        if line.startswith("#chromsize:"):
            parts = line.split()
            if len(parts) >= 3:
                out[parts[1]] = int(parts[2])
    return out


def resolutions(path: str) -> list[int]:
    k = kind_of(path)
    if k == "hic":
        from .hicfile import HicFile
        return sorted(HicFile(path).header.bp_resolutions)
    if k == "mcool":
        with _h5(path) as f:
            return sorted(int(r) for r in f["resolutions"].keys())
    if k == "cool":
        with _h5(path) as f:
            return [int(f.attrs.get("bin-size", 0))]
    return []           # pairs: any resolution


def read_region(path: str, chrom: str, start: int, end: int, resolution: int, normalization: str = "none",
                expected_assembly: str | None = None) -> Contacts:
    """Raw counts of chrom:start-end at `resolution` (region-local bins), with optional balancing weights."""
    if end <= start:
        raise ContactFileError("The region is empty (end <= start).")
    k = kind_of(path)
    notes = check_assembly(chromosomes(path), expected_assembly)
    lo = start // resolution
    n = (end - 1) // resolution - lo + 1
    if k == "hic":
        from .hicfile import HicFile
        hf = HicFile(path)
        if resolution not in hf.header.bp_resolutions:
            raise ContactFileError(f"The .hic has resolutions {', '.join(f'{r:,}' for r in sorted(hf.header.bp_resolutions))} "
                                   f"bp, not {resolution:,} bp.")
        names = [c for c, _ in hf.header.chromosomes]
        name = names[_find_chrom(names, chrom)]
        try:
            x, y, c = hf.matrix(name, name, resolution, start, end)
        except (ValueError, KeyError, OSError) as exc:
            raise ContactFileError(f"Could not read {chrom} from the .hic file: {exc}") from exc
        notes.append(f".hic version {hf.header.version}, {hf.f.bytes_fetched:,} bytes read")
        i, j, v = x - lo, y - lo, c
    elif k in ("cool", "mcool"):
        i, j, v, note = _cool_region(path, chrom, start, end, resolution)
        notes.append(note)
    else:
        i, j, v, note = _pairs_region(path, chrom, start, end, resolution)
        notes.append(note)
    a, b = np.minimum(i, j), np.maximum(i, j)
    v = np.asarray(v, float)
    ok = (a >= 0) & (b < n) & (v > 0)
    key = a[ok] * n + b[ok]
    uniq, inv = np.unique(key, return_inverse=True)
    sums = np.bincount(inv, weights=v[ok])
    out = Contacts((uniq // n).astype(np.int64), (uniq % n).astype(np.int64), sums, n, _norm(chrom), lo * resolution,
                   resolution, None, notes)
    if normalization != "none":
        from . import normalize as NZ
        off = out.ci != out.cj
        if normalization == "ice":
            r = NZ.ice_balance(out.ci[off], out.cj[off], out.cm[off], n)
            out.weights = np.where(np.isfinite(r.bias), 1.0 / r.bias, np.nan)
            out.notes.append(f"ICE weights ({'converged' if r.converged else 'not converged'}, {int(r.masked.sum())} bins masked)")
        elif normalization == "kr":
            r = NZ.kr_balance(out.ci[off], out.cj[off], out.cm[off], n)
            out.weights = r.weights
            out.notes.append(f"KR weights ({'converged' if r.converged else 'not converged'}, {int(r.masked.sum())} bins masked)")
        else:
            raise ContactFileError(f"Unknown normalisation '{normalization}' (none, ice or kr).")
    return out


def _cool_region(path: str, chrom: str, start: int, end: int, resolution: int):
    with _h5(path) as f:
        g = _cool_group(f, resolution)
        binsize = int(g.attrs.get("bin-size", 0)) or int(np.diff(g["bins"]["start"][:2])[0])
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in g["chroms"]["name"][:]]
        cid = _find_chrom(names, chrom)
        offsets = g["indexes"]["chrom_offset"][:]
        c_lo = int(offsets[cid])
        b_lo = c_lo + start // binsize
        b_hi = min(int(offsets[cid + 1]), c_lo + (end - 1) // binsize + 1)
        bin1_offset = g["indexes"]["bin1_offset"]
        p0, p1 = int(bin1_offset[b_lo]), int(bin1_offset[b_hi])
        b1 = g["pixels"]["bin1_id"][p0:p1].astype(np.int64)
        b2 = g["pixels"]["bin2_id"][p0:p1].astype(np.int64)
        cnt = g["pixels"]["count"][p0:p1].astype(np.float64)
    keep = (b2 >= b_lo) & (b2 < b_hi)
    pos1 = (b1[keep] - c_lo) * binsize
    pos2 = (b2[keep] - c_lo) * binsize
    lo = start // resolution
    factor = resolution // binsize
    note = (f"cooler at {binsize:,} bp" + (f", summed x{factor} to {resolution:,} bp" if factor > 1 else "")
            + f"; {p1 - p0:,} pixels read")
    return pos1 // resolution - lo, pos2 // resolution - lo, cnt[keep], note


def _pairs_region(path: str, chrom: str, start: int, end: int, resolution: int):
    import pandas as pd
    head = _pairs_header(path)
    cols = next((h.split()[1:] for h in head if h.startswith("#columns:")), None)
    if cols is None:
        cols = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]
    try:
        ic1, ip1, ic2, ip2 = cols.index("chrom1"), cols.index("pos1"), cols.index("chrom2"), cols.index("pos2")
    except ValueError as exc:
        raise ContactFileError("The .pairs header has no chrom1 / pos1 / chrom2 / pos2 columns.") from exc
    want = {chrom, _norm(chrom), chrom.removeprefix("chr")}
    lo = start // resolution
    I, J, total, seen = [], [], 0, set()
    with _open_text(path) as fh:
        reader = pd.read_csv(fh, sep="\t", comment="#", header=None, usecols=[ic1, ip1, ic2, ip2],
                             chunksize=CHUNK_LINES, dtype={ic1: str, ic2: str}, low_memory=False)
        try:
            for ch in reader:
                total += len(ch)
                seen.update(ch[ic1].unique()[:50].tolist())
                m = ch[ic1].isin(want) & ch[ic2].isin(want)
                p1, p2 = ch.loc[m, ip1].to_numpy(np.int64), ch.loc[m, ip2].to_numpy(np.int64)
                ok = (p1 >= start) & (p1 < end) & (p2 >= start) & (p2 < end)
                I.append(p1[ok] // resolution - lo)
                J.append(p2[ok] // resolution - lo)
        except (ValueError, pd.errors.ParserError) as exc:
            raise ContactFileError(f"Could not parse the .pairs file: {exc}") from exc
    if not any(_norm(s) == _norm(chrom) for s in seen):
        raise ContactFileError(f"{chrom} is not in this .pairs file. It has: {', '.join(sorted(map(str, seen))[:12])}.")
    i = np.concatenate(I) if I else np.empty(0, np.int64)
    j = np.concatenate(J) if J else np.empty(0, np.int64)
    return i, j, np.ones(len(i)), f".pairs streamed: {total:,} pairs read, {len(i):,} in the region"
