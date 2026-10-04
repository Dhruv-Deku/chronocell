"""
Registry, download and loading of every real dataset used to validate ChronoCell-5D.

Every entry records where the data come from, who produced it, under which licence, and whether
it is a PRACTICE dataset (tuning allowed) or a TEST dataset (run once, after every setting is
frozen). The split is fixed here, in code, before any v4 experiment ran (validation/TUNING.md):

  PRACTICE  Bintu K562 28-30 Mb; HCT116 28-30 Mb untreated and +auxin; HCT116 34-37 Mb untreated
            (all as in v3.3); Su et al. chromosome 2 (main and p-arm replicate), IMR-90.
  TEST      Bintu IMR90 28-30 Mb, IMR90 18-20 Mb, A549 28-30 Mb (as in v3.3); Su et al. chromosome
            21 (both replicates) and the genome-scale DNA-MERFISH set, IMR-90; Bintu HCT116 34-37 Mb
            +auxin (the held-out perturbation test).
  EXCLUDED  Bintu IMR90 28-30 Mb cell-cycle set and IMR90 STORM set (same cell line and region as a
            test set; not used for anything).

Nothing is redistributed: files are downloaded on demand into validation/data/ (git-ignored), and
the SHA-256 of every file actually used is written to validation/data_manifest.json so a result
can be matched to the exact bytes it came from.

    python validation/datasets.py --list
    python validation/datasets.py --fetch bintu_imr90_28_30 su_chr21 ...   (or --fetch all)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
MANIFEST = ROOT / "data_manifest.json"
UA = {"User-Agent": "ChronoCell-5D validation (research; github.com/Sh1voham/ChronoCell-5D)"}

BINTU_BASE = "https://raw.githubusercontent.com/BogdanBintu/ChromatinImaging/master/Data/"
SU_BASE = "https://zenodo.org/records/3928890/files/"

BINTU_CITE = ("Bintu B, Mateo LJ, Su J-H, et al. Super-resolution chromatin tracing reveals domains and cooperative "
              "interactions in single cells. Science 362, eaau1783 (2018). doi:10.1126/science.aau1783")
BINTU_LICENSE = ("No licence file in github.com/BogdanBintu/ChromatinImaging (all rights reserved by default). "
                 "Downloaded on demand for research validation; not redistributed.")
SU_CITE = ("Su J-H, Zheng P, Kinrot SS, Bintu B, Zhuang X. Genome-scale imaging of the 3D organization and "
           "transcriptional activity of chromatin. Cell 182, 1641-1659 (2020). doi:10.1016/j.cell.2020.07.032. "
           "Data: Zenodo record 3928890, doi:10.5281/zenodo.3928890")
SU_LICENSE = "CC-BY-4.0 (Zenodo record 3928890). Downloaded on demand; not redistributed."
RAO2014_CITE = ("Rao SSP, Huntley MH, Durand NC, et al. A 3D map of the human genome at kilobase resolution reveals "
                "principles of chromatin looping. Cell 159, 1665-1680 (2014). doi:10.1016/j.cell.2014.11.021. GEO GSE63525")


@dataclass(frozen=True)
class Entry:
    key: str
    kind: str                    # "bintu_csv" | "su_trace" | "su_genome" | "su_hic"
    role: str                    # "practice" | "test" | "excluded" | "input"
    cell_line: str
    region: str                  # human-readable region (hg38)
    files: tuple[str, ...]       # file name(s) at the source
    base: str
    citation: str
    license: str
    assembly: str = "hg38"
    condition: str = "untreated"
    step_bp: int = 30_000        # spacing between consecutive imaged loci
    locus_bp: int = 30_000       # size of each imaged locus
    region_start_hg38: int | None = None
    region_start_hg19: int | None = None
    paired_hic: str | None = None   # key of the matching Hi-C entry
    notes: str = ""

    @property
    def subdir(self) -> Path:
        return DATA if self.kind == "bintu_csv" else DATA / "su2020"

    def paths(self) -> list[Path]:
        return [self.subdir / f for f in self.files]

    def urls(self) -> list[str]:
        if self.base == SU_BASE:
            return [self.base + urllib.parse.quote(f) + "?download=1" for f in self.files]
        return [self.base + urllib.parse.quote(f) for f in self.files]


def _bintu(key, role, cell, region, fname, start38, start19=None, condition="untreated", notes=""):
    return Entry(key, "bintu_csv", role, cell, region, (fname,), BINTU_BASE, BINTU_CITE, BINTU_LICENSE,
                 condition=condition, region_start_hg38=start38, region_start_hg19=start19, notes=notes)


REGISTRY: dict[str, Entry] = {e.key: e for e in (
    # ---- v3.3 test set (unchanged) ----
    _bintu("bintu_imr90_28_30", "test", "IMR90", "chr21:28.0-29.9 Mb", "IMR90_chr21-28-30Mb.csv", 28_000_071, 29_372_390),
    _bintu("bintu_imr90_18_20", "test", "IMR90", "chr21:18.6-20.6 Mb", "IMR90_chr21-18-20Mb.csv", 18_627_714, 20_000_032),
    _bintu("bintu_a549_28_30", "test", "A549", "chr21:28.0-29.9 Mb", "A549_chr21-28-30Mb.csv", 28_000_071),
    # ---- v3.3 practice set (unchanged) ----
    _bintu("bintu_k562_28_30", "practice", "K562", "chr21:28.0-29.9 Mb", "K562_chr21-28-30Mb.csv", 28_000_071, 29_372_390),
    _bintu("bintu_hct116_28_30", "practice", "HCT116", "chr21:28.0-29.9 Mb", "HCT116_chr21-28-30Mb_untreated.csv",
           28_000_071, 29_372_390),
    _bintu("bintu_hct116_28_30_auxin", "practice", "HCT116", "chr21:28.0-29.9 Mb", "HCT116_chr21-28-30Mb_6h auxin.csv",
           28_000_071, 29_372_390, condition="6 h auxin (RAD21 degraded, cohesin depleted)"),
    _bintu("bintu_hct116_34_37", "practice", "HCT116", "chr21:34.6-37.1 Mb", "HCT116_chr21-34-37Mb_untreated.csv",
           34_628_096, 36_000_395),
    # ---- v4: held-out perturbation test ----
    _bintu("bintu_hct116_34_37_auxin", "test", "HCT116", "chr21:34.6-37.1 Mb", "HCT116_chr21-34-37Mb_6h auxin.csv",
           34_628_096, 36_000_395, condition="6 h auxin (RAD21 degraded, cohesin depleted)",
           notes="Held-out test of the cohesin-depletion perturbation (Pillar 4); never used for tuning."),
    # ---- excluded (same cell line and region as a test set) ----
    _bintu("bintu_imr90_28_30_cellcycle", "excluded", "IMR90", "chr21:28.0-29.9 Mb", "IMR90_chr21-28-30Mb_cell cycle.csv",
           28_000_071, notes="Shares cell line and region with a test set; not used."),
    # ---- Su et al. 2020 (Zhuang lab), IMR-90 ----
    Entry("su_chr2", "su_trace", "practice", "IMR90", "chr2, whole chromosome (935 loci, 250 kb steps)",
          ("chromosome2.tsv",), SU_BASE, SU_CITE, SU_LICENSE, step_bp=250_000, locus_bp=50_000, paired_hic="su_hic_chr2"),
    Entry("su_chr2_parm_rep", "su_trace", "practice", "IMR90", "chr2 p-arm replicate (250 kb steps)",
          ("chromosome2_p-arm_replicate.tsv",), SU_BASE, SU_CITE, SU_LICENSE, step_bp=250_000, locus_bp=50_000,
          paired_hic="su_hic_chr2", notes="Replicate not analysed in the paper."),
    Entry("su_chr21", "su_trace", "test", "IMR90", "chr21, whole chromosome (651 loci, 50 kb steps)",
          ("chromosome21.tsv",), SU_BASE, SU_CITE, SU_LICENSE, step_bp=50_000, locus_bp=50_000, paired_hic="su_hic_chr21",
          notes="Overlaps the Bintu IMR90 test regions, so it is test-only."),
    Entry("su_chr21_rep", "su_trace", "test", "IMR90", "chr21 replicate with cell-cycle markers (651 loci)",
          ("chromosome21-cell_cycle.tsv",), SU_BASE, SU_CITE, SU_LICENSE, step_bp=50_000, locus_bp=50_000,
          paired_hic="su_hic_chr21"),
    Entry("su_genome", "su_genome", "test", "IMR90", "1,041 loci on every chromosome (DNA-MERFISH, 3 replicates)",
          ("genomic-scale.tsv",), SU_BASE, SU_CITE, SU_LICENSE, step_bp=3_000_000, locus_bp=100_000,
          paired_hic="su_hic_genome", notes="chr2 loci overlap the practice set and are reported separately."),
    # ---- Hi-C binned onto the imaged loci by Su et al. from Rao et al. 2014 (IMR-90 in situ Hi-C) ----
    Entry("su_hic_chr21", "su_hic", "input", "IMR90", "Hi-C, chr21 loci (50 kb bins)", ("Hi-C_contacts_chromosome21.tsv",),
          SU_BASE, SU_CITE + " Hi-C reads: " + RAO2014_CITE, SU_LICENSE),
    Entry("su_hic_chr2", "su_hic", "input", "IMR90", "Hi-C, chr2 loci (50 kb bins)", ("Hi-C_contacts_chromosome2.tsv",),
          SU_BASE, SU_CITE + " Hi-C reads: " + RAO2014_CITE, SU_LICENSE),
    Entry("su_hic_genome", "su_hic", "input", "IMR90", "Hi-C, genome-scale loci (500 kb bins)",
          ("Hi-C_contacts_genome-scale.tsv",), SU_BASE, SU_CITE + " Hi-C reads: " + RAO2014_CITE, SU_LICENSE),
)}


def by_role(role: str) -> list[Entry]:
    return [e for e in REGISTRY.values() if e.role == role]


# ======================================================================================
# Download with integrity record
# ======================================================================================
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_manifest() -> dict:
    try:
        return json.loads(MANIFEST.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def record(path: Path, url: str) -> str:
    """SHA-256 of a data file, written to validation/data_manifest.json (committed)."""
    man = _load_manifest()
    rel = path.relative_to(ROOT).as_posix()
    digest = _sha256(path)
    man[rel] = {"sha256": digest, "bytes": path.stat().st_size, "url": url}
    MANIFEST.write_text(json.dumps(dict(sorted(man.items())), indent=1) + "\n", encoding="utf-8")
    return digest


def fetch(entry: Entry, retries: int = 3, verbose: bool = True) -> list[Path]:
    """Download the entry's files if missing (atomic write, resumable by retry) and record their hashes."""
    out = []
    for path, url in zip(entry.paths(), entry.urls()):
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(path.suffix + ".part")
            for attempt in range(1, retries + 1):
                try:
                    t0 = time.time()
                    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r, \
                            tmp.open("wb") as fh:
                        while chunk := r.read(1 << 20):
                            fh.write(chunk)
                    tmp.replace(path)
                    if verbose:
                        print(f"downloaded {path.name} ({path.stat().st_size / 1e6:.1f} MB, {time.time() - t0:.0f} s)",
                              flush=True)
                    break
                except OSError as exc:
                    if attempt == retries:
                        raise
                    print(f"retry {attempt} for {path.name}: {exc}", flush=True)
                    time.sleep(3 * attempt)
        man = _load_manifest().get(path.relative_to(ROOT).as_posix())
        if not man or man.get("bytes") != path.stat().st_size:
            record(path, url)
        out.append(path)
    return out


# ======================================================================================
# Loaders
# ======================================================================================
@dataclass
class Traces:
    """Single-chromosome-copy 3D traces of consecutive loci.

    xyz      (copies, loci, 3) nm, NaN where a locus was not detected
    starts   (loci,) hg38 start of each locus (bp, 0-based)
    chrom    chromosome name (single-chromosome traces) or per-locus names (genome-scale)
    groups   (copies,) replicate / experiment label per copy (genome-scale), else zeros
    """
    key: str
    xyz: np.ndarray
    starts: np.ndarray
    chrom: np.ndarray
    groups: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    meta: dict = field(default_factory=dict)

    @property
    def n_loci(self) -> int:
        return self.xyz.shape[1]

    @property
    def n_copies(self) -> int:
        return self.xyz.shape[0]


def load_bintu(entry: Entry, min_detected: float = 0.5) -> Traces:
    """Bintu et al. 2018 CSV: chromosome index, segment index, Z, X, Y (nm). Copies with fewer than
    `min_detected` of their segments detected are dropped (the v3.3 rule)."""
    path = fetch(entry, verbose=False)[0]
    rows = np.genfromtxt(path, delimiter=",", skip_header=2)
    cell = rows[:, 0].astype(int)
    seg = rows[:, 1].astype(int) - 1
    n_cells, n_seg = cell.max(), seg.max() + 1
    xyz = np.full((n_cells, n_seg, 3), np.nan)
    xyz[cell - 1, seg] = rows[:, 2:5]
    ok = np.isfinite(xyz[..., 0]).mean(axis=1) >= min_detected
    starts = entry.region_start_hg38 + np.arange(n_seg) * 30_000
    return Traces(entry.key, xyz[ok], starts.astype(np.int64), np.full(n_seg, "chr21"), np.zeros(int(ok.sum()), int),
                  {"copies_total": int(n_cells), "min_detected": min_detected})


_LOCUS = re.compile(r"(chr[0-9XYM]+):(\d+)-(\d+)")


def _parse_loci(col: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    names, starts = [], []
    for v in col:
        m = _LOCUS.match(str(v))
        names.append(m.group(1))
        starts.append(int(m.group(2)) - 1)
    return np.asarray(names), np.asarray(starts, np.int64)


def load_su(entry: Entry, min_detected: float = 0.5) -> Traces:
    """Su et al. 2020 sequential-hybridisation traces (one row per locus per chromosome copy).

    The copy boundary is where the locus order restarts; 'Chromosome copy number' alone is not a
    unique copy id. Copies with fewer than `min_detected` loci detected are dropped (as for Bintu).
    """
    import pandas as pd
    path = fetch(entry, verbose=False)[0]
    df = pd.read_csv(path, sep="\t", usecols=[0, 1, 2, 3], low_memory=False)
    df.columns = ["z", "x", "y", "locus"]
    loci = df["locus"].to_numpy()
    uniq, first = np.unique(loci, return_index=True)
    order = uniq[np.argsort(first)]
    n = len(order)
    if len(df) % n:
        raise ValueError(f"{path.name}: {len(df):,} rows are not a whole number of {n}-locus traces.")
    copies = len(df) // n
    if not (loci.reshape(copies, n) == order[None]).all():
        raise ValueError(f"{path.name}: loci are not in the same order in every trace.")
    xyz = df[["x", "y", "z"]].to_numpy(np.float64).reshape(copies, n, 3)
    chrom, starts = _parse_loci(order)
    ok = np.isfinite(xyz[..., 0]).mean(axis=1) >= min_detected
    return Traces(entry.key, xyz[ok], starts, chrom, np.zeros(int(ok.sum()), int),
                  {"copies_total": int(copies), "min_detected": min_detected, "loci": n})


def load_su_genome(entry: Entry, min_detected: float = 0.5) -> Traces:
    """Su et al. 2020 genome-scale DNA-MERFISH: one row per locus per homolog per cell per experiment.
    Each (experiment, cell, homolog) is one copy of one chromosome; returned as copies of the full
    locus list (loci of other chromosomes are NaN for that copy)."""
    import pandas as pd
    path = fetch(entry, verbose=False)[0]
    df = pd.read_csv(path, sep="\t", usecols=[0, 1, 2, 3, 4, 5, 6], low_memory=False)
    df.columns = ["z", "x", "y", "locus", "homolog", "cell", "exp"]
    loci = df["locus"].to_numpy()
    uniq, first = np.unique(loci, return_index=True)
    names, starts = _parse_loci(uniq)
    key = pd.Series(np.arange(len(uniq)), index=uniq)
    li = key.loc[loci].to_numpy()
    chrom_of = names[li]
    copy_id = pd.factorize(pd.Series(list(zip(df["exp"].to_numpy(), df["cell"].to_numpy(), df["homolog"].to_numpy(),
                                              chrom_of))))[0]
    order_loci = np.lexsort((starts, _chrom_rank(names)))
    rank = np.empty(len(uniq), int)
    rank[order_loci] = np.arange(len(uniq))
    xyz = np.full((copy_id.max() + 1, len(uniq), 3), np.nan, dtype=np.float32)
    xyz[copy_id, rank[li]] = df[["x", "y", "z"]].to_numpy(np.float32)
    exp = np.zeros(copy_id.max() + 1, int)
    exp[copy_id] = df["exp"].to_numpy()
    names, starts = names[order_loci], starts[order_loci]
    # a copy only covers its own chromosome: require min_detected of that chromosome's loci
    keep = np.zeros(len(xyz), bool)
    for c in np.unique(names):
        cols = names == c
        has = np.isfinite(xyz[:, cols, 0])
        own = has.any(axis=1)
        keep |= own & (has.mean(axis=1) >= min_detected)
    return Traces(entry.key, xyz[keep], starts, names, exp[keep],
                  {"copies_total": int(len(xyz)), "min_detected": min_detected, "loci": len(uniq)})


def _chrom_rank(names: np.ndarray) -> np.ndarray:
    def r(c: str) -> int:
        s = c.removeprefix("chr")
        return int(s) if s.isdigit() else {"X": 23, "Y": 24, "M": 25}.get(s, 26)
    return np.array([r(c) for c in names])


def load_su_hic(entry: Entry) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Square Hi-C matrix binned on the imaged loci (Rao et al. 2014 reads summed by Su et al.).
    Returns (counts (n, n), chrom (n,), starts (n,))."""
    import pandas as pd
    path = fetch(entry, verbose=False)[0]
    df = pd.read_csv(path, sep="\t", index_col=0)
    chrom, starts = _parse_loci(df.columns.to_numpy())
    mat = df.to_numpy(np.float64)
    if mat.shape[0] != mat.shape[1]:
        raise ValueError(f"{path.name}: Hi-C matrix is not square {mat.shape}.")
    return mat, chrom, starts


def load(key: str, **kw) -> Traces:
    e = REGISTRY[key]
    if e.kind == "bintu_csv":
        return load_bintu(e, **kw)
    if e.kind == "su_trace":
        return load_su(e, **kw)
    if e.kind == "su_genome":
        return load_su_genome(e, **kw)
    raise ValueError(f"{key} is not an imaging dataset.")


def describe() -> list[dict]:
    return [asdict(e) | {"paths": [p.relative_to(ROOT).as_posix() for p in e.paths()]} for e in REGISTRY.values()]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--fetch", nargs="*", default=None, help="dataset keys, or 'all'")
    a = ap.parse_args()
    if a.list or a.fetch is None:
        for e in REGISTRY.values():
            print(f"{e.key:30s} {e.role:9s} {e.cell_line:7s} {e.region}")
    if a.fetch is not None:
        keys = list(REGISTRY) if a.fetch == ["all"] else a.fetch
        for k in keys:
            if k not in REGISTRY:
                sys.exit(f"unknown dataset {k}")
            fetch(REGISTRY[k])


if __name__ == "__main__":
    main()
