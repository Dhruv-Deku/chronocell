"""
No-Hi-C prediction (Pillar 5): the pattern of 3D distances between loci from sequence and CTCF binding
alone, for regions with no contact data. A hypothesis under test (validation/predictor.py, Gate 5),
not part of the app unless that test passes.

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
    """GC fraction of any interval in O(1) from cumulative counts (N and other symbols are ignored)."""

    def __init__(self, seq: bytes):
        a = np.frombuffer(seq, dtype=np.uint8)
        gc = (a == ord("G")) | (a == ord("C"))
        acgt = gc | (a == ord("A")) | (a == ord("T"))
        self.gc = np.concatenate([[0], np.cumsum(gc, dtype=np.int64)])
        self.acgt = np.concatenate([[0], np.cumsum(acgt, dtype=np.int64)])

    def fraction(self, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        starts, ends = np.asarray(starts, np.int64), np.asarray(ends, np.int64)
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
