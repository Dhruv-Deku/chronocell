"""
Coordinate liftover between assemblies with UCSC chain files (Phase B4): hg19 <-> hg38, mm10 <-> mm39.

Chain files are downloaded on demand from UCSC (goldenPath/<from>/liftOver/<from>To<To>.over.chain.gz), checked
against the MD5 in the md5sum.txt UCSC publishes in the same directory, and kept in .chronocell_cache/liftover/.
A position maps through the ungapped block that contains it; positions in gaps or unaligned sequence do not
map (None), and an interval maps only if both ends map to the same chromosome and strand with its length
changed by at most `max_size_change` (as UCSC liftOver's minMatch idea, simplified).

    from chronocell import liftover as LO
    ch = LO.Chain.load("hg19", "hg38")
    ch.position("chr21", 28_000_000)          # -> ("chr21", pos, "+") or None
"""

from __future__ import annotations

import bisect
import gzip
import io
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

PAIRS = {("hg19", "hg38"), ("hg38", "hg19"), ("mm10", "mm39"), ("mm39", "mm10")}
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/{0}/liftOver/"
CACHE = Path(__file__).resolve().parent.parent / ".chronocell_cache" / "liftover"
UA = {"User-Agent": "ChronoCell-5D (research; liftover)"}


def _file(a: str, b: str) -> str:
    return f"{a}To{b[0].upper()}{b[1:]}.over.chain.gz"


def chain_path(a: str, b: str) -> Path:
    """The chain file a -> b, downloaded once and MD5-checked against UCSC md5sum.txt."""
    if (a, b) not in PAIRS:
        raise ValueError(f"No liftover {a} -> {b} (available: {', '.join(f'{x}->{y}' for x, y in sorted(PAIRS))}).")
    name = _file(a, b)
    out = CACHE / name
    if out.exists():
        return out
    from .predict import download_verified
    sums = urllib.request.urlopen(urllib.request.Request(UCSC.format(a) + "md5sum.txt", headers=UA), timeout=60).read()
    md5 = {ln.split()[1].lstrip("*"): ln.split()[0] for ln in sums.decode().splitlines() if len(ln.split()) == 2}
    if name not in md5:
        raise OSError(f"UCSC md5sum.txt does not list {name}; not downloaded.")
    return download_verified(UCSC.format(a) + name, out, md5[name])


@dataclass
class _Blocks:
    t_start: list[int] = field(default_factory=list)      # source (target in UCSC terms) block starts, sorted
    t_end: list[int] = field(default_factory=list)
    q_chrom: list[str] = field(default_factory=list)
    q_start: list[int] = field(default_factory=list)
    q_strand: list[str] = field(default_factory=list)
    q_size: list[int] = field(default_factory=list)


class Chain:
    def __init__(self, blocks: dict[str, _Blocks], source: str = "", target: str = ""):
        self.blocks = blocks
        self.source, self.target = source, target

    @classmethod
    def load(cls, a: str, b: str) -> "Chain":
        return cls.parse(chain_path(a, b).read_bytes(), a, b)

    @classmethod
    def parse(cls, data: bytes, source: str = "", target: str = "") -> "Chain":
        """UCSC chain format: 'chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id'
        followed by 'size dt dq' lines (the last line 'size')."""
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        raw: dict[str, list[tuple]] = {}
        t_name = q_name = q_strand = None
        t_pos = q_pos = q_size = 0
        for line in io.StringIO(data.decode("ascii", "replace")):
            f = line.split()
            if not f:
                continue
            if f[0] == "chain":
                t_name, t_pos = f[2], int(f[5])
                q_name, q_size, q_strand, q_pos = f[7], int(f[8]), f[9], int(f[10])
                continue
            size = int(f[0])
            raw.setdefault(t_name, []).append((t_pos, t_pos + size, q_name, q_pos, q_strand, q_size))
            if len(f) == 3:
                t_pos += size + int(f[1])
                q_pos += size + int(f[2])
        blocks = {}
        for c, rows in raw.items():
            rows.sort()
            b = _Blocks()
            for ts, te, qn, qs, st, qz in rows:
                b.t_start.append(ts)
                b.t_end.append(te)
                b.q_chrom.append(qn)
                b.q_start.append(qs)
                b.q_strand.append(st)
                b.q_size.append(qz)
            blocks[c] = b
        return cls(blocks, source, target)

    def position(self, chrom: str, pos: int) -> tuple[str, int, str] | None:
        """0-based position -> (chrom, 0-based position, strand) in the target assembly, or None."""
        b = self.blocks.get(chrom)
        if b is None:
            return None
        k = bisect.bisect_right(b.t_start, pos) - 1
        if k < 0 or pos >= b.t_end[k]:
            return None
        off = pos - b.t_start[k]
        if b.q_strand[k] == "+":
            return b.q_chrom[k], b.q_start[k] + off, "+"
        return b.q_chrom[k], b.q_size[k] - (b.q_start[k] + off) - 1, "-"

    def interval(self, chrom: str, start: int, end: int, max_size_change: float = 0.5) -> tuple[str, int, int] | None:
        a = self.position(chrom, start)
        z = self.position(chrom, end - 1)
        if a is None or z is None or a[0] != z[0] or a[2] != z[2]:
            return None
        lo, hi = min(a[1], z[1]), max(a[1], z[1]) + 1
        if abs((hi - lo) - (end - start)) > max_size_change * (end - start):
            return None
        return a[0], lo, hi

    def bed(self, text: str) -> tuple[str, int]:
        """Lift a BED / bedGraph text (first three columns); returns (lifted text, lines that did not map)."""
        out, lost = [], 0
        for line in text.splitlines():
            if not line.strip() or line.startswith(("#", "track", "browser")):
                out.append(line)
                continue
            f = line.split("\t")
            r = self.interval(f[0], int(f[1]), int(f[2]))
            if r is None:
                lost += 1
                continue
            out.append("\t".join([r[0], str(r[1]), str(r[2])] + f[3:]))
        return "\n".join(out) + "\n", lost

    def bedpe(self, text: str) -> tuple[str, int]:
        out, lost = [], 0
        for line in text.splitlines():
            if not line.strip() or line.startswith(("#", "track", "browser")):
                out.append(line)
                continue
            f = line.split("\t")
            r1, r2 = self.interval(f[0], int(f[1]), int(f[2])), self.interval(f[3], int(f[4]), int(f[5]))
            if r1 is None or r2 is None:
                lost += 1
                continue
            out.append("\t".join([r1[0], str(r1[1]), str(r1[2]), r2[0], str(r2[1]), str(r2[2])] + f[6:]))
        return "\n".join(out) + "\n", lost
