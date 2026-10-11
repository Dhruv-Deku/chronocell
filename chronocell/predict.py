"""
No-Hi-C prediction (Pillar 5): the pattern of 3D distances between loci from sequence and CTCF binding
alone, for regions with no contact data. Tested held out in validation/predictor.py (Gate 5): it passed
the pre-registered rule, modestly, so the app and `python -m chronocell.predict` offer it, always
labelled "predicted, not measured".

Per locus (a genome interval):
  fwd, rev   CTCF ChIP-seq peaks in the locus whose best CTCF motif match (a position weight matrix,
             e.g. JASPAR MA0139.1) is on the + / - strand of the reference; unk = peaks with no match
             above the relative-score threshold
  gc         GC fraction of the locus sequence
Per pair i < j (s = genomic separation, bp), with L(x) = log(1 + x):
  conv    L(fwd_i) L(rev_j)                      motifs pointing at each other
  div     L(rev_i) L(fwd_j)                      motifs pointing away from each other
  tandem  L(fwd_i) L(fwd_j) + L(rev_i) L(rev_j)
  anchor  L(n_i) + L(n_j)                        any CTCF peak at either end
  between L(peaks between the loci)              CTCF sites separating them
  gc      (gc_i - m)(gc_j - m)                   compartment-like similarity (m: the region's mean)
  and each of these times log10(s)
Target: the trend-removed log median distance, r_ij = log d_ij - mean of log d over pairs at the same
separation. Ridge regression on standardised features, pooled over training datasets with equal weight
per dataset. Prediction for a new region: log d_ij = a + b log10(s_ij) + x_ij . beta, where a and b are
the training data's distance-separation trend (the no-data genomic baseline uses that trend alone).
Orientation convention: '+' is the strand on which the matrix, as published, matches the reference.
Which orientation pattern matters is learned from the training data, not assumed.
"""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

BASES = b"ACGT"
FEATURES = ("conv", "div", "tandem", "anchor", "between", "gc")
DATA = Path(__file__).with_name("data")
MODEL_PATH = DATA / "predictor.json"                 # frozen copy of validation/predictor_model.json (Gate 5)
JASPAR_PATH = DATA / "jaspar_MA0139.1.jaspar"         # JASPAR 2024, MA0139.1 (CTCF), CC BY 4.0
JASPAR_CITATION = ("Rauluseviciute I et al. JASPAR 2024. Nucleic Acids Res 52, D174-D182 (2024); matrix MA0139.1 "
                   "(CTCF), CC BY 4.0.")


# ======================================================================================
# Sequence and motif
# ======================================================================================
def read_fasta(path: str | Path) -> bytes:
    """One-sequence FASTA (optionally gzipped) as upper-case bytes."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rb") as fh:
        lines = fh.read().split(b"\n")
    return b"".join(ln.strip() for ln in lines if ln and not ln.startswith(b">")).upper()


class GCIndex:
    """GC fraction of any interval in O(1) from cumulative counts (N and other symbols are ignored).
    `offset` = the chromosome coordinate of seq[0], for an index built on part of a chromosome."""

    def __init__(self, seq: bytes, offset: int = 0):
        a = np.frombuffer(seq, dtype=np.uint8)
        gc = (a == ord("G")) | (a == ord("C"))
        acgt = gc | (a == ord("A")) | (a == ord("T"))
        self.offset = int(offset)
        self.gc = np.concatenate([[0], np.cumsum(gc, dtype=np.int64)])
        self.acgt = np.concatenate([[0], np.cumsum(acgt, dtype=np.int64)])

    def fraction(self, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        last = len(self.gc) - 1
        starts = np.clip(np.asarray(starts, np.int64) - self.offset, 0, last)
        ends = np.clip(np.asarray(ends, np.int64) - self.offset, 0, last)
        n = self.acgt[ends] - self.acgt[starts]
        return np.where(n > 0, (self.gc[ends] - self.gc[starts]) / np.maximum(n, 1), np.nan)


def read_jaspar(text: str) -> np.ndarray:
    """JASPAR-format count matrix -> (4, L) counts in A, C, G, T order."""
    rows = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        base, _, rest = line.partition("[")
        rows[base.strip()] = [float(v) for v in rest.rstrip("]").split()]
    return np.array([rows[b] for b in "ACGT"])


def log_odds(counts: np.ndarray, pseudocount: float = 0.25, background=(0.25, 0.25, 0.25, 0.25)) -> np.ndarray:
    p = (counts + pseudocount) / (counts.sum(axis=0, keepdims=True) + 4 * pseudocount)
    return np.log2(p / np.asarray(background)[:, None])


def _onehot_index(seq: bytes) -> np.ndarray:
    lut = np.full(256, 4, np.int64)
    for k, b in enumerate(BASES):
        lut[b] = k
    return lut[np.frombuffer(seq, dtype=np.uint8)]


def best_match(seq: bytes, pwm: np.ndarray) -> tuple[float, int]:
    """Best relative score (0..1) of the matrix on either strand of `seq`, and its strand (+1 / -1);
    windows containing N are skipped. (nan, 0) if the sequence is shorter than the motif."""
    L = pwm.shape[1]
    idx = _onehot_index(seq)
    if len(idx) < L:
        return float("nan"), 0
    ext = np.vstack([pwm, np.full((1, L), -np.inf)])            # row 4 = N
    rc = ext[[3, 2, 1, 0, 4]][:, ::-1]                            # reverse complement
    win = np.lib.stride_tricks.sliding_window_view(idx, L)
    cols = np.arange(L)
    fwd = ext[win, cols].sum(axis=1)
    rev = rc[win, cols].sum(axis=1)
    lo, hi = pwm.min(axis=0).sum(), pwm.max(axis=0).sum()
    f, r = fwd.max(), rev.max()
    best, strand = (f, 1) if f >= r else (r, -1)
    if not np.isfinite(best):
        return float("nan"), 0
    return float((best - lo) / (hi - lo)), strand


def orient_peaks(chrom_seq: bytes, starts: np.ndarray, ends: np.ndarray, summits: np.ndarray, pwm: np.ndarray,
                 min_relative: float = 0.8, half_width: int = 100) -> tuple[np.ndarray, np.ndarray]:
    """Motif strand per peak (+1 / -1, 0 if no match >= min_relative) scanned within half_width of the
    summit (or of the peak centre when no summit is given), and the relative score."""
    strand = np.zeros(len(starts), np.int8)
    score = np.full(len(starts), np.nan)
    for k, (s, e, sm) in enumerate(zip(starts, ends, summits)):
        c = int(s + sm) if sm >= 0 else int((s + e) // 2)
        lo, hi = max(int(s), c - half_width), min(int(e), c + half_width)
        if hi - lo < pwm.shape[1]:
            lo, hi = int(s), int(e)
        rel, st = best_match(chrom_seq[lo:hi], pwm)
        score[k] = rel
        if np.isfinite(rel) and rel >= min_relative:
            strand[k] = st
    return strand, score


# ======================================================================================
# Features
# ======================================================================================
@dataclass
class Loci:
    """Loci of one region with their CTCF and sequence annotations."""
    starts: np.ndarray        # bp, 0-based
    ends: np.ndarray
    fwd: np.ndarray
    rev: np.ndarray
    unk: np.ndarray
    gc: np.ndarray

    @property
    def n(self) -> int:
        return len(self.starts)


def annotate(starts: np.ndarray, ends: np.ndarray, peak_mid: np.ndarray, peak_strand: np.ndarray,
             gc_index: GCIndex) -> Loci:
    starts, ends = np.asarray(starts, np.int64), np.asarray(ends, np.int64)
    order = np.argsort(peak_mid)
    mid, st = np.asarray(peak_mid)[order], np.asarray(peak_strand)[order]
    lo = np.searchsorted(mid, starts, side="left")
    hi = np.searchsorted(mid, ends, side="left")
    cnt = lambda m: np.array([int(m[a:b].sum()) for a, b in zip(lo, hi)])   # noqa: E731
    return Loci(starts, ends, cnt(st == 1), cnt(st == -1), cnt(st == 0), gc_index.fraction(starts, ends))


def pair_features(loci: Loci, peak_mid: np.ndarray, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Feature matrix (pairs x 12) for pairs i < j: the six FEATURES, then each times log10(s)."""
    i, j = np.asarray(i), np.asarray(j)
    L = np.log1p
    n_all = loci.fwd + loci.rev + loci.unk
    conv = L(loci.fwd[i]) * L(loci.rev[j])
    div = L(loci.rev[i]) * L(loci.fwd[j])
    tandem = L(loci.fwd[i]) * L(loci.fwd[j]) + L(loci.rev[i]) * L(loci.rev[j])
    anchor = L(n_all[i]) + L(n_all[j])
    mid = np.sort(np.asarray(peak_mid))
    between = L(np.searchsorted(mid, loci.starts[j], side="left") - np.searchsorted(mid, loci.ends[i], side="left"))
    g = np.nan_to_num(loci.gc - np.nanmean(loci.gc), nan=0.0)
    gc = g[i] * g[j]
    base = np.stack([conv, div, tandem, anchor, between, gc], axis=1)
    s = np.abs((loci.starts[j] + loci.ends[j]) / 2 - (loci.starts[i] + loci.ends[i]) / 2)
    ls = np.log10(np.maximum(s, 1.0))[:, None]
    return np.hstack([base, base * ls])


def feature_names() -> list[str]:
    return list(FEATURES) + [f"{f} x log10 s" for f in FEATURES]


def trend_residual(log_d: np.ndarray, sep: np.ndarray) -> np.ndarray:
    """log d minus the mean of log d over pairs with the same separation (NaN-aware)."""
    out = np.full_like(log_d, np.nan, dtype=np.float64)
    ok = np.isfinite(log_d)
    keys, inv = np.unique(sep, return_inverse=True)
    tot = np.bincount(inv[ok], weights=log_d[ok], minlength=len(keys))
    cnt = np.bincount(inv[ok], minlength=len(keys))
    mean = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
    out[ok] = log_d[ok] - mean[inv[ok]]
    return out


# ======================================================================================
# Model
# ======================================================================================
@dataclass
class Predictor:
    beta: list[float]
    mu: list[float]
    sd: list[float]
    trend: tuple[float, float]               # log d = a + b log10(s / 1 bp), d in nm
    ridge_lambda: float
    names: list[str] = field(default_factory=feature_names)
    trained_on: list[str] = field(default_factory=list)
    motif: str = ""
    notes: str = ""

    def residual(self, X: np.ndarray) -> np.ndarray:
        Z = (X - np.asarray(self.mu)) / np.asarray(self.sd)
        return Z @ np.asarray(self.beta)

    def log_distance(self, X: np.ndarray, sep_bp: np.ndarray) -> np.ndarray:
        a, b = self.trend
        return a + b * np.log10(np.maximum(sep_bp, 1.0)) + self.residual(X)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=1)

    @classmethod
    def from_json(cls, text: str) -> "Predictor":
        d = json.loads(text)
        d["trend"] = tuple(d["trend"])
        return cls(**d)


def fit_trend(log_d: list[np.ndarray], sep_bp: list[np.ndarray]) -> tuple[float, float]:
    """log d = a + b log10 s by least squares, every dataset weighted equally."""
    X, y, w = [], [], []
    for ld, s in zip(log_d, sep_bp):
        ok = np.isfinite(ld) & (s > 0)
        X.append(np.log10(s[ok]))
        y.append(ld[ok])
        w.append(np.full(ok.sum(), 1.0 / max(ok.sum(), 1)))
    X, y, w = np.concatenate(X), np.concatenate(y), np.concatenate(w)
    A = np.stack([np.ones_like(X), X], axis=1)
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
    return float(coef[0]), float(coef[1])


def fit_ridge(Xs: list[np.ndarray], ys: list[np.ndarray], lam: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Weighted ridge on standardised features (each dataset weighs the same); returns beta, mu, sd."""
    X = np.vstack(Xs)
    y = np.concatenate(ys)
    w = np.concatenate([np.full(len(v), 1.0 / max(len(v), 1)) for v in ys])
    w = w / w.sum()
    mu = (w[:, None] * X).sum(axis=0)
    sd = np.sqrt((w[:, None] * (X - mu) ** 2).sum(axis=0))
    sd = np.where(sd > 0, sd, 1.0)
    Z = (X - mu) / sd
    # minimise sum_k w_k (y_k - z_k . beta)^2 + lam |beta|^2 with sum w = 1 (standardised z: lam ~ shrinkage)
    A = Z.T @ (w[:, None] * Z) + lam * np.eye(Z.shape[1])
    beta = np.linalg.solve(A, Z.T @ (w * y))
    return beta, mu, sd


# ======================================================================================
# In the app: the frozen model on a window with no contact data
# ======================================================================================
def load_model() -> tuple["Predictor", dict] | None:
    """The frozen predictor and its settings / validation record (chronocell/data/predictor.json)."""
    try:
        d = json.loads(MODEL_PATH.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    meta = {k: d.pop(k) for k in ("settings", "assembly", "validation", "source_note") if k in d}
    d["trend"] = tuple(d["trend"])
    return Predictor(**d), meta


def load_pwm(pseudocount: float = 0.25) -> np.ndarray:
    return log_odds(read_jaspar(JASPAR_PATH.read_text(encoding="utf-8")), pseudocount)


def read_peaks(data: bytes | str, chrom_name: str, norm=None) -> dict[str, np.ndarray]:
    """ChIP-seq peaks of one chromosome from narrowPeak (summit in column 10) or BED (gzip accepted);
    header, track and browser lines are skipped. norm(name) normalises chromosome names."""
    if isinstance(data, bytes):
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        data = data.decode("utf-8", "replace")
    starts, ends, summits = [], [], []
    for line in data.splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t") if "\t" in line else line.split()
        if len(f) < 3:
            continue
        try:
            c = norm(f[0]) if norm else f[0]
            a, b = int(f[1]), int(f[2])
        except (ValueError, KeyError):
            continue
        if c != chrom_name or b <= a:
            continue
        sm = -1
        if len(f) >= 10:
            try:
                sm = int(f[9])
            except ValueError:
                sm = -1
        starts.append(a)
        ends.append(b)
        summits.append(sm)
    return {"starts": np.asarray(starts, np.int64), "ends": np.asarray(ends, np.int64),
            "summits": np.asarray(summits, np.int64)}


def predict_window(loci_starts: np.ndarray, loci_ends: np.ndarray, seq: bytes, peaks: dict[str, np.ndarray],
                   model: "Predictor", settings: dict, pwm: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    """Predicted median distance (nm, N x N) between loci from sequence + CTCF peaks alone, with the
    frozen model and settings of Gate 5. Returns the matrix and a summary of the inputs used."""
    pwm = load_pwm(settings.get("pseudocount", 0.25)) if pwm is None else pwm
    s, e, sm = peaks["starts"], peaks["ends"], peaks["summits"]
    strand, score = orient_peaks(seq, s, e, sm, pwm, settings.get("motif_min_relative_score", 0.8),
                                 settings.get("summit_half_width_bp", 100))
    mid = np.where(sm >= 0, s + sm, (s + e) // 2)
    lo, hi = int(np.min(loci_starts)), int(np.max(loci_ends))
    loci = annotate(loci_starts, loci_ends, mid, strand, GCIndex(seq[lo:hi], offset=lo))
    i, j = np.triu_indices(loci.n, 1)
    X = pair_features(loci, mid, i, j)
    sep = np.abs((loci.starts[j] + loci.ends[j]) - (loci.starts[i] + loci.ends[i])) / 2
    d = np.zeros((loci.n, loci.n))
    v = np.exp(model.log_distance(X, sep))
    d[i, j] = v
    d[j, i] = v
    info = {"peaks": int(len(s)), "peaks_with_motif": int((strand != 0).sum()),
            "loci_with_peak": int(((loci.fwd + loci.rev + loci.unk) > 0).sum()), "loci": int(loci.n),
            "gc_mean": float(np.nanmean(loci.gc)) if np.isfinite(loci.gc).any() else float("nan")}
    return d, info


def compare_maps(a: np.ndarray, b: np.ndarray, max_pairs: int = 200_000, seed: int = 0) -> dict:
    """Agreement between two median-distance maps of the same loci (e.g. built from contacts vs predicted):
    Spearman rho over the pairs i < j, raw and beyond the separation trend (each distance divided by its
    map's mean at the same |i - j|, as in the validation protocol). Pairs are a seeded random sample when
    there are more than max_pairs. Agreement between two models, not accuracy."""
    from .accuracy import spearman
    a, b = np.asarray(a, np.float64), np.asarray(b, np.float64)
    n = len(a)
    if a.shape != (n, n) or b.shape != (n, n) or n < 3:
        raise ValueError("Need two square maps of the same size (at least 3 loci).")
    i, j = np.triu_indices(n, 1)
    sampled = len(i) > max_pairs
    if sampled:
        pick = np.random.default_rng(seed).choice(len(i), max_pairs, replace=False)
        i, j = i[pick], j[pick]
    x, y = a[i, j], b[i, j]
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x, y, s = x[ok], y[ok], (j - i)[ok]
    _, inv = np.unique(s, return_inverse=True)
    xe = (np.bincount(inv, weights=x) / np.bincount(inv))[inv]
    ye = (np.bincount(inv, weights=y) / np.bincount(inv))[inv]
    return {"spearman": spearman(x, y), "spearman_trend_removed": spearman(x / xe, y / ye), "pairs": int(ok.sum()),
            "sampled": bool(sampled), "median_ratio": float(np.median(y / x)) if ok.any() else float("nan")}


UCSC_CHROMOSOMES = "https://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/chromosomes/"
ENCODE_DOWNLOAD = "https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
SOURCES_PATH = DATA / "validation_sources.json"
USER_AGENT = {"User-Agent": "ChronoCell-5D (research; github.com/Sh1voham/ChronoCell-5D)"}


def _trust_system_certificates() -> None:
    try:
        import truststore                  # the system certificate store, where TLS is inspected
        truststore.inject_into_ssl()
    except ImportError:
        pass


def download_verified(url: str, out: str | Path, md5: str, progress=None, retries: int = 6) -> Path:
    """Download url to out through out.part: interrupted or short transfers resume with an HTTP Range
    request; the file is moved into place only if its MD5 equals the one its source publishes (on a
    mismatch the partial file is deleted and OSError raised). progress(bytes_done, bytes_total)."""
    import hashlib
    import time
    import urllib.request
    _trust_system_certificates()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".part")
    for attempt in range(1, retries + 1):
        have = tmp.stat().st_size if tmp.exists() else 0
        req = urllib.request.Request(url, headers=dict(USER_AGENT, **({"Range": f"bytes={have}-"} if have else {})))
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                resumed = bool(have) and r.status == 206
                total = r.headers.get("Content-Range", "").rpartition("/")[2] if resumed else r.headers.get("Content-Length")
                total = int(total) if total and total.isdigit() else None
                done = have if resumed else 0
                with tmp.open("ab" if resumed else "wb") as fh:
                    while chunk := r.read(1 << 20):
                        fh.write(chunk)
                        done += len(chunk)
                        if progress is not None:
                            progress(done, total)
            if total is not None and tmp.stat().st_size < total:
                raise OSError(f"short read ({tmp.stat().st_size:,} of {total:,} bytes)")
            break
        except OSError:
            if attempt == retries:
                raise
            time.sleep(min(30, 2 * 2 ** attempt))
    h = hashlib.md5()
    with tmp.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != md5:
        tmp.unlink()
        raise OSError(f"{out.name} does not match the MD5 its source publishes; deleted, please retry.")
    tmp.replace(out)
    return out


def fetch_chromosome_fasta(assembly: str, chrom: str, dest_dir: str | Path, progress=None, retries: int = 6) -> Path:
    """One chromosome's sequence from UCSC ({assembly}/chromosomes/{chrom}.fa.gz), kept gzipped in dest_dir.
    Interrupted transfers resume; the file is checked against UCSC's md5sum.txt before it is used (no
    published MD5, no download). progress(bytes_done, bytes_total) is called while downloading."""
    import time
    import urllib.request
    _trust_system_certificates()
    out = Path(dest_dir) / f"{assembly}_{chrom}.fa.gz"
    if out.exists():
        return out
    base = UCSC_CHROMOSOMES.format(assembly=assembly)
    md5 = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(urllib.request.Request(base + "md5sum.txt", headers=USER_AGENT), timeout=60) as r:
                for line in r.read().decode().splitlines():
                    parts = line.split()
                    if len(parts) == 2 and parts[1] == f"{chrom}.fa.gz":
                        md5 = parts[0]
            break
        except OSError:
            if attempt == retries:
                raise
            time.sleep(min(30, 2 * 2 ** attempt))
    if md5 is None:
        raise OSError(f"UCSC's md5sum.txt lists no {chrom}.fa.gz for {assembly}; not downloading an unverifiable file.")
    return download_verified(base + f"{chrom}.fa.gz", out, md5, progress, retries)


def encode_ctcf_sources(path: Path = SOURCES_PATH) -> list[dict]:
    """The ENCODE CTCF ChIP-seq peak files listed in data/validation_sources.json (the Gate 5 inputs):
    cell line, file accession, experiment, the MD5 the ENCODE portal publishes, assembly and citation."""
    try:
        rows = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    out = []
    for r in rows:
        if not str(r.get("key", "")).startswith("encode_ctcf_") or not r.get("md5"):
            continue
        acc = r["url"].rstrip("/").rpartition("/")[2]
        out.append({"cell_line": r["cell_line"], "accession": acc, "experiment": r.get("experiment", ""), "md5": r["md5"],
                    "assembly": r.get("assembly", ""), "citation": r["citation"], "license": r.get("license", ""),
                    "page": r["url"]})
    return out


def fetch_encode_peaks(accession: str, md5: str, dest_dir: str | Path, progress=None, retries: int = 6) -> Path:
    """An ENCODE peak file (bed.gz), downloaded once into dest_dir and checked against the portal's MD5."""
    out = Path(dest_dir) / f"{accession}.bed.gz"
    if out.exists():
        return out
    return download_verified(ENCODE_DOWNLOAD.format(acc=accession), out, md5, progress, retries)


# ======================================================================================
# Command line
# ======================================================================================
CACHE_FASTA = Path(__file__).resolve().parent.parent / ".chronocell_cache" / "fasta"
VALIDATED_LOCUS_BP = (30_000, 50_000)        # locus sizes of the Gate 5 test sets (Bintu 2018; Su 2020 chr21)


def read_fasta_record(path: str | Path, name: str, norm=None) -> bytes:
    """The sequence of one record of a FASTA file (optionally gzipped), upper case: the record whose name
    (first word of the header, normalised by norm) is `name`, or the only record of a one-sequence file."""
    op = gzip.open if str(path).endswith(".gz") else open
    seqs: dict[str, list[bytes]] = {}
    cur = None
    with op(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                head = line[1:].strip().split()
                cur = head[0].decode("utf-8", "replace") if head else ""
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line.strip())
    if not seqs:
        raise ValueError(f"{path}: no FASTA record")
    for k, v in seqs.items():
        try:
            if (norm(k) if norm else k) == name:
                return b"".join(v).upper()
        except (KeyError, ValueError):
            continue
    if len(seqs) == 1:
        return b"".join(next(iter(seqs.values()))).upper()
    raise ValueError(f"{path}: no record named {name} (records: {', '.join(list(seqs)[:5])} ...)")


def _sha256(path: str | Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


MOUSE_ASSEMBLIES = ("mm10", "mm39")


def assembly_supported(assembly_name: str, meta: dict) -> bool:
    """The training assembly (hg38), or a mouse assembly once the pre-registered mouse test (Gate 5m) has passed."""
    if assembly_name == meta.get("assembly"):
        return True
    if assembly_name in MOUSE_ASSEMBLIES:
        from .accuracy import mouse_predictor_evidence
        ev = mouse_predictor_evidence()
        return bool(ev and ev["verdict"] == "pass")
    return False


def run_cli(chrom: str, start: int, end: int, peaks_path: str | Path, out: str | Path, fasta: str | Path | None = None,
            bin_bp: int = 30_000, assembly_name: str = "hg38") -> dict:
    """Predict the median-distance map of [start, end) in bins of bin_bp from CTCF peaks + sequence; write
    out (.npy, nm) and out with .json (the record: inputs with SHA-256, model, validation). Returns the record."""
    import datetime as dt
    from . import accuracy as ACC, genome
    loaded = load_model()
    if loaded is None:
        raise SystemExit(f"The predictor file is missing: {MODEL_PATH}")
    model, meta = loaded
    if not assembly_supported(assembly_name, meta):
        raise SystemExit(f"The predictor was trained and tested on human {meta.get('assembly')} only; "
                         f"{assembly_name} is not supported (validation/RESULTS.md, Gate 5).")
    name = genome.normalize_chrom(chrom, assembly_name)
    size = genome.chromosome_size(name, assembly_name)
    if not (0 <= start < end <= size):
        raise SystemExit(f"Region {name}:{start:,}-{end:,} is outside the chromosome (size {size:,}).")
    if bin_bp <= 0 or (end - start) // bin_bp < 3:
        raise SystemExit("Need at least 3 bins: lower --bin or widen the region.")
    norm = lambda c: genome.normalize_chrom(c, assembly_name)       # noqa: E731
    if fasta is None:
        fasta = fetch_chromosome_fasta(assembly_name, name, CACHE_FASTA,
                                       lambda d, t: print(f"\rdownloading {name}: {d / 1e6:,.1f} MB", end="", flush=True))
        print()
        seq_source = f"UCSC {assembly_name} {name}.fa.gz (MD5 checked)"
    else:
        seq_source = "user file"
    seq = read_fasta_record(fasta, name, norm)
    if len(seq) < end:
        raise SystemExit(f"The sequence of {name} in {fasta} is {len(seq):,} bp, shorter than the region end {end:,}.")
    pk = read_peaks(Path(peaks_path).read_bytes(), name, norm)
    if len(pk["starts"]) == 0:
        raise SystemExit(f"No peaks on {name} in {peaks_path}.")
    edges = np.arange(start, end - bin_bp + 1, bin_bp, dtype=np.int64)
    d, info = predict_window(edges, edges + bin_bp, seq, pk, model, meta["settings"])
    out = Path(out)
    np.save(out, d)
    bench = (ACC.load_benchmark() or {}).get("models", {}).get("predicted_sequence_ctcf") or {}
    record = {
        "status": "predicted, not measured",
        "what": "median 3D distance (nm) across cells between every pair of loci, predicted from CTCF ChIP-seq peaks "
                "(oriented by the JASPAR CTCF motif) and GC content, with no contact data",
        "software": "ChronoCell-5D chronocell.predict",
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "output": {"file": out.name, "sha256": _sha256(out), "shape": list(d.shape), "units": "nm"},
        "region": {"assembly": assembly_name, "chrom": name, "start": int(start), "end": int(edges[-1] + bin_bp),
                   "bin_bp": int(bin_bp), "loci": int(len(edges))},
        "inputs": {"peaks": {"file": str(peaks_path), "sha256": _sha256(peaks_path), **info},
                   "sequence": {"file": str(fasta), "sha256": _sha256(fasta), "source": seq_source}},
        "model": {"file": "chronocell/data/predictor.json", "sha256": _sha256(MODEL_PATH), "trained_on": model.trained_on,
                  "ridge_lambda": model.ridge_lambda, "motif": meta["settings"]["motif"], "motif_citation": JASPAR_CITATION},
        "validation": {**(meta.get("validation") or {}),
                       "overall_percent_of_ceiling": bench.get("overall_percent_of_ceiling"),
                       "per_dataset_percent_of_ceiling": bench.get("per_dataset_percent_of_ceiling"),
                       "cohesin_depleted_control_percent_of_ceiling": bench.get("control_percent_of_ceiling"),
                       "validated_locus_bp": list(VALIDATED_LOCUS_BP),
                       "reading": "a prior, not a measurement; its signal is compartment / insulation level, not CTCF "
                                  "loops (the cohesin-depleted control scored as high); see validation/RESULTS.md, Gate 5"},
    }
    if bin_bp not in VALIDATED_LOCUS_BP:
        record["validation"]["warning"] = (f"bin size {bin_bp:,} bp was not tested; Gate 5 tested "
                                           + " and ".join(f"{b // 1000} kb" for b in VALIDATED_LOCUS_BP) + " loci")
    out.with_suffix(".json").write_text(json.dumps(record, indent=1, default=float) + "\n", encoding="utf-8")
    return record


def main(argv: list[str] | None = None) -> None:
    """python -m chronocell.predict --chrom chr21 --start 28000000 --end 30000000 --peaks ctcf.narrowPeak --out map.npy"""
    import argparse
    ap = argparse.ArgumentParser(prog="python -m chronocell.predict",
                                 description="Predict a median 3D distance map from CTCF peaks + sequence, with no contact "
                                             "data (the frozen Gate 5 model; hg38, and mouse once Gate 5m has passed). Writes OUT (.npy, nm) and a "
                                             "JSON record next to it. The output is predicted, not measured.")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True, help="region start (bp, 0-based)")
    ap.add_argument("--end", type=int, required=True, help="region end (bp)")
    ap.add_argument("--peaks", required=True, help="CTCF ChIP-seq peaks of the cell type: narrowPeak or BED (.gz accepted)")
    ap.add_argument("--fasta", help="the chromosome's sequence (FASTA, .gz accepted); default: download from UCSC once")
    ap.add_argument("--bin", type=int, default=30_000, help="locus size in bp (tested: 30000, 50000; default 30000)")
    ap.add_argument("--assembly", default="hg38")
    ap.add_argument("--out", required=True, help="output .npy (a .json record is written beside it)")
    a = ap.parse_args(argv)
    rec = run_cli(a.chrom, a.start, a.end, a.peaks, a.out, a.fasta, a.bin, a.assembly)
    r = rec["region"]
    print(f"wrote {a.out} ({r['loci']} x {r['loci']} loci, {r['chrom']}:{r['start']:,}-{r['end']:,}) and "
          f"{Path(a.out).with_suffix('.json').name}: PREDICTED, not measured; "
          f"{rec['inputs']['peaks']['peaks']} peaks, {rec['inputs']['peaks']['peaks_with_motif']} with a motif match.")
    if "warning" in rec["validation"]:
        print("warning:", rec["validation"]["warning"])


if __name__ == "__main__":
    main()
