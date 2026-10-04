"""
Reference genome model for any assembly configured under chronocell/data/genomes/.

An assembly is data, not code: a folder chronocell/data/genomes/<assembly>/ with a manifest.json
(display name, species, main chromosomes, file names, sources, licence) plus the files it names:
chromosome sizes / cytobands / gaps / centromere spans, RefSeq Select genes and a chromosome-alias
table. `python -m chronocell.genome_fetch <assembly>` builds such a folder from the UCSC Genome
Browser. Shipped: hg38 (GRCh38, human; the default) and mm39 (GRCm39, mouse).

All coordinates are 0-based, half-open [start, end) in base pairs, as served by the UCSC Genome
Browser REST API.

`Chrom` bundles one chromosome of one assembly with a bin resolution. Resolution defaults to the
finest of the standard Hi-C resolutions that keeps the chromosome at <= 6,000 beads (hg38 chr22 ->
10 kb, chr1 -> 50 kb); it can also be inferred from the bead count of an uploaded structure.
The module-level names (CHROM, N_BINS, bin_start, ...) describe hg38 chr22 at 10 kb and are kept
for backwards compatibility.
"""

from __future__ import annotations

import gzip
import json
import math
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np

DATA = Path(__file__).with_name("data")
GENOMES = DATA / "genomes"
DEFAULT_ASSEMBLY = "hg38"
RESOLUTIONS = (5_000, 10_000, 20_000, 25_000, 40_000, 50_000, 100_000, 250_000, 500_000, 1_000_000)
MAX_DEFAULT_BEADS = 6_000


@dataclass(frozen=True)
class Band:
    start: int
    end: int
    name: str
    stain: str


@dataclass(frozen=True)
class Gene:
    name: str
    chrom: str
    start: int
    end: int
    strand: str


# ======================================================================================
# Assemblies (data folders with a manifest)
# ======================================================================================
@dataclass(frozen=True)
class Assembly:
    name: str                         # "hg38"
    display: str                      # "GRCh38/hg38"
    species: str                      # "Homo sapiens"
    common: str                       # "human"
    main_chromosomes: tuple[str, ...]
    folder: Path
    files: dict = field(default_factory=dict)
    sources: dict = field(default_factory=dict)
    licence: str = ""
    citation: str = ""

    @property
    def short(self) -> str:
        """'GRCh38' from 'GRCh38/hg38'."""
        return self.display.split("/")[0]

    def path(self, key: str) -> Path | None:
        f = self.files.get(key)
        return (self.folder / f).resolve() if f else None


@lru_cache(maxsize=1)
def assemblies() -> dict[str, Assembly]:
    """Every assembly with a readable manifest under chronocell/data/genomes/."""
    out = {}
    for man in sorted(GENOMES.glob("*/manifest.json")):
        try:
            m = json.loads(man.read_text(encoding="utf-8"))
            out[m["assembly"]] = Assembly(m["assembly"], m.get("display", m["assembly"]), m.get("species", ""),
                                          m.get("common", ""), tuple(m["main_chromosomes"]), man.parent,
                                          m.get("files", {}), m.get("sources", {}), m.get("licence", ""),
                                          m.get("citation", ""))
        except (OSError, ValueError, KeyError):
            continue                      # a broken folder never takes the others down
    return out


def assembly(name: str | None = None) -> Assembly:
    a = assemblies()
    key = name or DEFAULT_ASSEMBLY
    if key not in a:
        for asm in a.values():            # accept the display name too ('GRCh38/hg38', 'GRCm39')
            if key in (asm.display, asm.short):
                return asm
        raise KeyError(f"Unknown genome assembly '{key}' (available: {', '.join(a)}).")
    return a[key]


def main_chromosomes(name: str | None = None) -> tuple[str, ...]:
    return assembly(name).main_chromosomes


@lru_cache(maxsize=8)
def _data(name: str | None = None) -> dict:
    with open(assembly(name).path("chromosomes"), encoding="utf-8") as fh:
        return json.load(fh)


@lru_cache(maxsize=8)
def _aliases(name: str | None = None) -> dict[str, str]:
    """Every alias (UCSC, assembly, Ensembl, GenBank, RefSeq; with and without version) -> UCSC name."""
    p = assembly(name).path("aliases")
    out: dict[str, str] = {}
    if p is None or not p.exists():
        return out
    for line in p.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        ucsc = cols[0]
        for c in cols:
            c = c.strip()
            if c:
                out.setdefault(c, ucsc)
                out.setdefault(c.lower(), ucsc)
                out.setdefault(c.split(".")[0].lower(), ucsc)
    return out


# Backwards-compatible constants (hg38)
ASSEMBLY = "GRCh38/hg38"
MAIN_CHROMOSOMES = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])


def gene(name: str, assembly_name: str | None = None) -> Gene:
    g = _data(assembly_name)["genes"][name]
    return Gene(name, g["chrom"], g["start"], g["end"], g["strand"])


def genes_in(chrom_name: str, start: int, end: int, assembly_name: str | None = None) -> list[Gene]:
    """Annotated anchor genes overlapping [start, end) on chrom_name (sorted by position)."""
    out = [gene(n, assembly_name) for n, g in _data(assembly_name)["genes"].items()
           if g["chrom"] == chrom_name and g["start"] < end and g["end"] > start]
    return sorted(out, key=lambda g: g.start)


def chromosome_size(name: str, assembly_name: str | None = None) -> int:
    return int(_data(assembly_name)["chromosomes"][name]["size"])


def default_resolution(name: str, assembly_name: str | None = None) -> int:
    size = chromosome_size(name, assembly_name)
    for r in RESOLUTIONS:
        if math.ceil(size / r) <= MAX_DEFAULT_BEADS:
            return r
    return RESOLUTIONS[-1]


def resolution_for_beads(name: str, n_beads: int, assembly_name: str | None = None) -> int:
    """Resolution whose bin count equals `n_beads` (standard values first, else size / n)."""
    size = chromosome_size(name, assembly_name)
    for r in RESOLUTIONS:
        if math.ceil(size / r) == n_beads:
            return r
    return max(1, math.ceil(size / max(n_beads, 1)))


@dataclass(frozen=True)
class Chrom:
    name: str
    size: int
    resolution: int
    bands: tuple[Band, ...]
    gaps: tuple[tuple[int, int], ...]
    centromere_model: tuple[int, int] | None
    assembly: str = DEFAULT_ASSEMBLY

    # ---- bins ---------------------------------------------------------------------------
    @property
    def n_bins(self) -> int:
        return math.ceil(self.size / self.resolution)

    @property
    def short(self) -> str:
        return self.name.removeprefix("chr")

    def bin_start(self, i) -> np.ndarray:
        return np.asarray(i, dtype=np.int64) * self.resolution

    def bin_end(self, i) -> np.ndarray:
        return np.minimum((np.asarray(i, dtype=np.int64) + 1) * self.resolution, self.size)

    def bin_lengths(self) -> np.ndarray:
        idx = np.arange(self.n_bins)
        return (self.bin_end(idx) - self.bin_start(idx)).astype(np.int64)

    def interval_to_bins(self, start: int, end: int) -> tuple[int, int]:
        """Half-open bp interval -> half-open bin interval covering it (clipped to the chromosome)."""
        return max(0, start // self.resolution), min(self.n_bins, math.ceil(end / self.resolution))

    def locus(self, i: int) -> str:
        return f"{self.name}:{int(self.bin_start(i)) + 1:,}-{int(self.bin_end(i)):,}"

    def mb(self, i) -> np.ndarray:
        return self.bin_start(i) / 1e6

    # ---- annotation ---------------------------------------------------------------------
    def gap_fraction(self) -> np.ndarray:
        """Fraction of each bin covered by assembly gaps (N), by exact interval overlap."""
        n = self.n_bins
        idx = np.arange(n)
        starts, ends = self.bin_start(idx), self.bin_end(idx)
        covered = np.zeros(n, dtype=np.int64)
        for g0, g1 in self.gaps:
            b0, b1 = self.interval_to_bins(g0, g1)
            if b0 >= b1:
                continue
            sl = slice(b0, b1)
            covered[sl] += np.clip(np.minimum(ends[sl], g1) - np.maximum(starts[sl], g0), 0, None)
        return covered / np.maximum(ends - starts, 1)

    def assembled_mask(self, max_gap: float = 0.5) -> np.ndarray:
        return self.gap_fraction() <= max_gap

    def band_for_bins(self) -> np.ndarray:
        idx = np.arange(self.n_bins)
        mids = (self.bin_start(idx) + self.bin_end(idx)) / 2
        edges = np.array([b.end for b in self.bands])
        return np.searchsorted(edges, mids, side="right").clip(0, len(self.bands) - 1)

    @property
    def acen(self) -> tuple[int, int] | None:
        spans = [(b.start, b.end) for b in self.bands if b.stain == "acen"]
        return (min(s for s, _ in spans), max(e for _, e in spans)) if spans else None

    @property
    def genome(self) -> Assembly:
        return assembly(self.assembly)

    def with_resolution(self, resolution: int) -> "Chrom":
        return Chrom(self.name, self.size, int(resolution), self.bands, self.gaps, self.centromere_model, self.assembly)


@lru_cache(maxsize=128)
def chrom(name: str = "chr22", resolution: int | None = None, assembly_name: str | None = None) -> Chrom:
    asm = assembly(assembly_name).name
    key = name if name in _data(asm)["chromosomes"] else normalize_chrom(name, asm)
    rec = _data(asm)["chromosomes"][key]
    bands = tuple(Band(int(s), int(e), n, st) for s, e, n, st in rec["bands"])
    gaps = tuple((int(s), int(e)) for s, e in rec["gaps"])
    cen = tuple(rec["centromere_model"]) if rec["centromere_model"] else None
    return Chrom(key, int(rec["size"]), int(resolution or default_resolution(key, asm)), bands, gaps, cen, asm)


def chrom_for_beads(name: str, n_beads: int, assembly_name: str | None = None) -> Chrom:
    return chrom(name, resolution_for_beads(name, n_beads, assembly_name), assembly_name)


_REFSEQ_HUMAN = {f"NC_{i:06d}": f"chr{i}" for i in range(1, 23)} | {"NC_000023": "chrX", "NC_000024": "chrY",
                                                                    "NC_012920": "chrM"}


def normalize_chrom(name: str, assembly_name: str | None = None) -> str:
    """Any common spelling of a chromosome name -> the UCSC style used here ('chr9', 'chrX', 'chrM').

    Accepts 'chr9', 'Chr9', 'CHR9', '9', 'X', 'chrx', 'M', 'MT', 'chrMT', and every alias in the
    assembly's alias table (Ensembl, GenBank and RefSeq accessions, with or without version, e.g.
    NC_000009.12 -> chr9 on hg38, CM001012.3 -> chr19 on mm39). Unknown names raise KeyError.
    """
    raw = str(name).strip()
    if not raw:
        raise KeyError("empty chromosome name")
    aliases = _aliases(assembly_name)
    for cand in (raw, raw.lower(), raw.split(".")[0].lower()):
        if cand in aliases:
            return aliases[cand]
    acc = raw.split(".")[0].upper()
    if (assembly_name or DEFAULT_ASSEMBLY) == "hg38" and acc in _REFSEQ_HUMAN:
        return _REFSEQ_HUMAN[acc]
    core = raw[3:] if raw.lower().startswith("chr") else raw
    core = core.strip().upper()
    if core in ("M", "MT"):
        return "chrM"
    if core in ("X", "Y"):
        return f"chr{core}"
    if core.isdigit() and 1 <= int(core) <= 99:
        return f"chr{int(core)}"
    if "_" in core and re.fullmatch(r"([0-9]+|X|Y|UN)_[A-Z0-9]+(V\d+)?(_ALT|_RANDOM|_FIX)?", core):
        rest = raw[3:] if raw.lower().startswith("chr") else raw      # alt / random / unplaced contigs keep their case
        return "chr" + rest
    raise KeyError(f"Unrecognised chromosome name '{name}'.")


def genes_path(assembly_name: str | None = None) -> Path | None:
    return assembly(assembly_name).path("genes")


# ---- hg38 chr22 at 10 kb: backwards-compatible module API ---------------------------------
DEFAULT = chrom("chr22", 10_000)
CHROM = DEFAULT.name
CHROM_SIZE = DEFAULT.size
RESOLUTION = DEFAULT.resolution
N_BINS = DEFAULT.n_bins
CYTOBANDS = DEFAULT.bands
GAPS = DEFAULT.gaps
CENTROMERE_MODEL = DEFAULT.centromere_model
ACEN = DEFAULT.acen


def bin_start(i) -> np.ndarray:
    return DEFAULT.bin_start(i)


def bin_end(i) -> np.ndarray:
    return DEFAULT.bin_end(i)


def bin_lengths(n: int = N_BINS) -> np.ndarray:
    return DEFAULT.bin_lengths()[:n]


def interval_to_bins(start: int, end: int) -> tuple[int, int]:
    return DEFAULT.interval_to_bins(start, end)


def gap_fraction(n: int = N_BINS) -> np.ndarray:
    return DEFAULT.gap_fraction()[:n]


def assembled_mask(n: int = N_BINS, max_gap: float = 0.5) -> np.ndarray:
    return DEFAULT.assembled_mask(max_gap)[:n]


def band_for_bins(n: int = N_BINS) -> np.ndarray:
    return DEFAULT.band_for_bins()[:n]


def locus(i: int) -> str:
    return DEFAULT.locus(i)
