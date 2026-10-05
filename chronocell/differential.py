"""
Two-condition differential analysis of contact maps (Phase B2), replicate-aware.

Input: two conditions, each a list of replicate maps of the same region (dense raw counts, same bins).

Pixels (differential contacts)
    1. Library size: each replicate's total within the tested band (pairs min_sep..max_sep bins apart).
    2. log2 counts per million, log2((c + 0.5) / (L + 1) x 1e6).
    3. Distance normalisation (on by default): each replicate's values are centred per diagonal on the mean of
       all replicates on that diagonal, removing sample-wide P(s) differences (as the loess-by-distance step of
       multiHiCcompare does). Turn it off to test P(s) changes as well.
    4. limma-style moderated t (Smyth, Stat Appl Genet Mol Biol 2004): per pixel the condition means and the
       pooled within-condition variance s^2 (d = n_A + n_B - 2), shrunk towards a prior s0^2 with d0 degrees of
       freedom estimated from all pixels (empirical Bayes, method of moments on log s^2, in 10 abundance strata);
       t = (mean_B - mean_A) / sqrt(s~^2 (1/n_A + 1/n_B)) with d + d0 degrees of freedom.
    5. Benjamini-Hochberg false-discovery rate over all tested pixels.
    Pixels with a mean count below `min_count` are not tested. Without at least two replicates in each
    condition no statistic is possible: fold changes only, and the result says so.

Per-bin tracks (insulation, compartment eigenvector) use the same moderated t across replicates of the
per-bin values. Loops, boundaries and compartment switches compare calls made on each condition's pooled map
(chronocell.analysis / chronocell.domains) and attach the pixel or bin statistics where available.
Several resolutions: `aggregate(mat, k)` sums k x k blocks.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import domains as DOM

STRATA = 10


def aggregate(mat: np.ndarray, k: int) -> np.ndarray:
    """Sum k x k blocks (a coarser resolution); the last partial block is dropped."""
    if k == 1:
        return np.asarray(mat, float)
    m = (len(mat) // k) * k
    a = np.asarray(mat, float)[:m, :m]
    return a.reshape(m // k, k, m // k, k).sum(axis=(1, 3))


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    q = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


def squeeze_var(s2: np.ndarray, d: float) -> tuple[float, float]:
    """Empirical Bayes prior (d0, s0^2) for sample variances s2 with d degrees of freedom (Smyth 2004, method
    of moments on log s^2). d0 = inf means no extra variability (all variances alike)."""
    from scipy.special import digamma, polygamma
    s2 = np.asarray(s2, float)
    s2 = s2[np.isfinite(s2) & (s2 > 0)]
    if len(s2) < 3:
        return 0.0, float(np.mean(s2)) if len(s2) else 1.0
    z = np.log(s2)
    e = z - digamma(d / 2) + np.log(d / 2)
    emean = e.mean()
    evar = np.mean((e - emean) ** 2) * len(e) / (len(e) - 1) - polygamma(1, d / 2)
    if evar <= 0:
        return float("inf"), float(np.exp(emean))
    # solve trigamma(d0 / 2) = evar
    lo, hi = 1e-6, 1e6
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if polygamma(1, mid / 2) > evar:
            lo = mid
        else:
            hi = mid
    d0 = float(np.sqrt(lo * hi))
    s0 = float(np.exp(emean + digamma(d0 / 2) - np.log(d0 / 2)))
    return d0, s0


def moderated_t(ya: np.ndarray, yb: np.ndarray, abundance: np.ndarray | None = None) -> dict:
    """ya (n_A, m), yb (n_B, m): per-feature moderated t, p and BH q (abundance-stratified prior)."""
    from scipy.stats import t as tdist
    na, nb = ya.shape[0], yb.shape[0]
    ma, mb = ya.mean(0), yb.mean(0)
    d = na + nb - 2
    lfc = mb - ma
    if d < 1:
        return {"log2_fc": lfc, "t": np.full_like(lfc, np.nan), "p": np.full_like(lfc, np.nan),
                "q": np.full_like(lfc, np.nan), "df_prior": np.nan}
    ss = ((ya - ma) ** 2).sum(0) + ((yb - mb) ** 2).sum(0)
    s2 = ss / d
    post = np.empty_like(s2)
    dfs = np.empty_like(s2)
    ab = abundance if abundance is not None else (ma + mb) / 2
    edges = np.unique(np.quantile(ab, np.linspace(0, 1, STRATA + 1)))
    strat = np.clip(np.searchsorted(edges, ab, side="right") - 1, 0, max(len(edges) - 2, 0))
    d0s = []
    for k in np.unique(strat):
        sel = strat == k
        d0, s0 = squeeze_var(s2[sel], d)
        d0s.append(d0)
        if np.isinf(d0):
            post[sel], dfs[sel] = s0, 1e6
        else:
            post[sel] = (d0 * s0 + d * s2[sel]) / (d0 + d)
            dfs[sel] = d + d0
    se = np.sqrt(np.maximum(post, 1e-12) * (1 / na + 1 / nb))
    t = lfc / se
    p = 2 * tdist.sf(np.abs(t), dfs)
    return {"log2_fc": lfc, "t": t, "p": p, "q": _bh(p), "df_prior": float(np.median(d0s))}


@dataclass
class DiffResult:
    pixels: pd.DataFrame                  # i, j, separation, mean counts, log2_fc, t, p, q
    statistics: bool
    fdr: float
    resolution: int
    notes: list[str] = field(default_factory=list)
    tracks: dict = field(default_factory=dict)       # insulation / compartment per-bin tests
    loops: pd.DataFrame | None = None
    boundaries: pd.DataFrame | None = None
    compartments: pd.DataFrame | None = None

    @property
    def significant(self) -> pd.DataFrame:
        if not self.statistics:
            return self.pixels.iloc[0:0]
        return self.pixels[self.pixels["q"] <= self.fdr]

    def summary(self) -> dict:
        s = self.significant
        return {"resolution": self.resolution, "tested_pixels": int(len(self.pixels)), "statistics": self.statistics,
                "fdr": self.fdr, "significant": int(len(s)), "gained": int((s["log2_fc"] > 0).sum()) if len(s) else 0,
                "lost": int((s["log2_fc"] < 0).sum()) if len(s) else 0}


def _log_cpm(mats: list[np.ndarray], mask: np.ndarray) -> np.ndarray:
    out = []
    for m in mats:
        lib = float(np.asarray(m, float)[mask].sum())
        out.append(np.log2((np.asarray(m, float)[mask] + 0.5) / (lib + 1.0) * 1e6))
    return np.array(out)


def differential_pixels(cond_a: list[np.ndarray], cond_b: list[np.ndarray], resolution: int, fdr: float = 0.05,
                        min_sep: int = 2, max_sep: int | None = None, min_count: float = 5.0,
                        distance_normalise: bool = True) -> DiffResult:
    n = len(cond_a[0])
    for m in list(cond_a) + list(cond_b):
        if np.asarray(m).shape != (n, n):
            raise ValueError("Every replicate must be a square map of the same bins.")
    max_sep = max_sep or n - 1
    sep = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    band = np.triu(np.ones((n, n), bool), min_sep) & (sep <= max_sep)
    allm = list(cond_a) + list(cond_b)
    mean_counts = np.mean([np.asarray(m, float)[band] for m in allm], axis=0)
    keep = mean_counts >= min_count
    ii, jj = np.nonzero(band)
    ii, jj, mc = ii[keep], jj[keep], mean_counts[keep]
    mask = np.zeros((n, n), bool)
    mask[ii, jj] = True
    ya, yb = _log_cpm(cond_a, mask), _log_cpm(cond_b, mask)     # boolean indexing: the row-major order of ii, jj
    notes = []
    if distance_normalise:
        s = jj - ii
        both = np.vstack([ya, yb])
        for d in np.unique(s):
            sel = s == d
            centre = both[:, sel].mean()
            for arr in (ya, yb):
                for r in range(arr.shape[0]):
                    arr[r, sel] += centre - arr[r, sel].mean()
        notes.append("Distance-normalised: each replicate centred per diagonal (P(s) differences removed).")
    stats = len(cond_a) >= 2 and len(cond_b) >= 2
    res = moderated_t(ya, yb) if stats else {"log2_fc": yb.mean(0) - ya.mean(0)}
    if not stats:
        notes.append("Each condition needs at least two replicates for statistics: fold changes only, no p-values.")
    df = pd.DataFrame({"i": ii, "j": jj, "separation_bp": (jj - ii) * resolution, "mean_count": mc,
                       "log2_fc": res["log2_fc"], "t": res.get("t", np.nan), "p": res.get("p", np.nan),
                       "q": res.get("q", np.nan)})
    return DiffResult(df, stats, fdr, resolution, notes)


def _sparse(mat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    m = np.triu(np.asarray(mat, float), 1)
    i, j = np.nonzero(m)
    return i, j, m[i, j]


def track_test(values_a: list[np.ndarray], values_b: list[np.ndarray], name: str) -> pd.DataFrame | None:
    """Per-bin moderated t of a track (insulation, eigenvector) across replicates; None without replicates."""
    A, B = np.array(values_a, float), np.array(values_b, float)
    ok = np.isfinite(A).all(0) & np.isfinite(B).all(0)
    if A.shape[0] < 2 or B.shape[0] < 2 or ok.sum() < 5:
        return None
    r = moderated_t(A[:, ok], B[:, ok], abundance=np.zeros(ok.sum()))
    out = pd.DataFrame({"bin": np.flatnonzero(ok), f"{name}_a": A[:, ok].mean(0), f"{name}_b": B[:, ok].mean(0),
                        "difference": r["log2_fc"], "t": r["t"], "p": r["p"], "q": r["q"]})
    return out


def compare(cond_a: list[np.ndarray], cond_b: list[np.ndarray], resolution: int, fdr: float = 0.05,
            max_sep_bp: int = 2_000_000, min_count: float = 5.0, orient: np.ndarray | None = None,
            loops: bool = True) -> DiffResult:
    """Everything: differential pixels, insulation and compartment tracks, boundary and loop changes."""
    from . import analysis as AN
    n = len(cond_a[0])
    max_sep = min(n - 1, int(max_sep_bp // resolution))
    res = differential_pixels(cond_a, cond_b, resolution, fdr, max_sep=max_sep, min_count=min_count)
    w = DOM.window_for(resolution, n)
    ins_a = [DOM.insulation(*_sparse(m), n, w) for m in cond_a]
    ins_b = [DOM.insulation(*_sparse(m), n, w) for m in cond_b]
    ev_a = [DOM.compartments(*_sparse(m), n, orient)[0] for m in cond_a]
    ev_b = [DOM.compartments(*_sparse(m), n, orient)[0] for m in cond_b]
    if orient is None:      # eigenvector signs are arbitrary per replicate: align each to the first replicate of A
        ref = ev_a[0]
        def align(v):
            ok = np.isfinite(v) & np.isfinite(ref)
            return -v if ok.sum() > 3 and np.corrcoef(v[ok], ref[ok])[0, 1] < 0 else v
        ev_a, ev_b = [align(v) for v in ev_a], [align(v) for v in ev_b]
        res.notes.append("Compartment signs aligned to the first replicate of condition A (no GC / gene track).")
    res.tracks["insulation"] = track_test(ins_a, ins_b, "insulation")
    res.tracks["compartment"] = track_test(ev_a, ev_b, "eigenvector")
    pa, pb = np.sum(cond_a, axis=0), np.sum(cond_b, axis=0)
    ins_pa, ins_pb = DOM.insulation(*_sparse(pa), n, w), DOM.insulation(*_sparse(pb), n, w)
    ba, bb = DOM.boundaries(ins_pa, w), DOM.boundaries(ins_pb, w)
    rows = [{"bin": b, "change": "lost in B"} for b in ba if all(abs(b - x) > 2 for x in bb)]
    rows += [{"bin": b, "change": "gained in B"} for b in bb if all(abs(b - x) > 2 for x in ba)]
    rows += [{"bin": b, "change": "kept"} for b in ba if any(abs(b - x) <= 2 for x in bb)]
    bdf = pd.DataFrame(rows, columns=["bin", "change"])
    if res.tracks["insulation"] is not None and len(bdf):
        t = res.tracks["insulation"].set_index("bin")
        bdf["insulation_q"] = [float(t["q"].get(b, np.nan)) for b in bdf["bin"]]
    res.boundaries = bdf
    eva, evb = np.nanmean(ev_a, axis=0), np.nanmean(ev_b, axis=0)
    sw = np.isfinite(eva) & np.isfinite(evb) & (np.sign(eva) != np.sign(evb)) & (np.abs(evb - eva) > 0.01)
    cdf = pd.DataFrame({"bin": np.flatnonzero(sw), "eigenvector_a": eva[sw], "eigenvector_b": evb[sw],
                        "switch": np.where(evb[sw] > 0, "B to A" if orient is not None else "sign + in B",
                                           "A to B" if orient is not None else "sign - in B")})
    if res.tracks["compartment"] is not None and len(cdf):
        t = res.tracks["compartment"].set_index("bin")
        cdf["q"] = [float(t["q"].get(b, np.nan)) for b in cdf["bin"]]
    res.compartments = cdf
    if loops:
        la = AN.call_loops(*_sparse(pa), n, resolution, max_sep_bp=max_sep_bp)
        lb = AN.call_loops(*_sparse(pb), n, resolution, max_sep_bp=max_sep_bp)
        rows = []
        def near(c, others):
            return any(abs(c.i - o.i) <= 2 and abs(c.j - o.j) <= 2 for o in others)
        for c in la:
            rows.append({"i": c.i, "j": c.j, "in_a": True, "in_b": near(c, lb)})
        for c in lb:
            if not near(c, la):
                rows.append({"i": c.i, "j": c.j, "in_a": False, "in_b": True})
        ldf = pd.DataFrame(rows, columns=["i", "j", "in_a", "in_b"])
        if len(ldf):
            ldf["change"] = np.where(ldf["in_a"] & ldf["in_b"], "kept", np.where(ldf["in_a"], "lost in B", "gained in B"))
            px = res.pixels.set_index(["i", "j"])
            ldf["log2_fc"] = [float(px["log2_fc"].get((i, j), np.nan)) for i, j in zip(ldf["i"], ldf["j"])]
            ldf["q"] = [float(px["q"].get((i, j), np.nan)) for i, j in zip(ldf["i"], ldf["j"])]
        res.loops = ldf
    return res


def to_bedpe(df: pd.DataFrame, chrom: str, start_bp: int, resolution: int) -> str:
    """Differential pixels or loops as BEDPE (0-based, half-open), with log2_fc and q when present."""
    lines = []
    for r in df.itertuples(index=False):
        a, b = start_bp + int(r.i) * resolution, start_bp + int(r.j) * resolution
        extra = "\t".join(str(round(float(getattr(r, k)), 6)) if hasattr(r, k) else "." for k in ("log2_fc", "q"))
        lines.append(f"{chrom}\t{a}\t{a + resolution}\t{chrom}\t{b}\t{b + resolution}\t.\t.\t{extra}")
    return "\n".join(lines) + ("\n" if lines else "")


def plant_changes(mat: np.ndarray, pixels: list[tuple[int, int]], fold: float, rng: np.random.Generator) -> np.ndarray:
    """Spike-in: multiply the expected count of the given pixels by `fold` (binomial thinning for fold < 1,
    added Poisson counts for fold > 1). Used by Gate 7."""
    out = np.asarray(mat, float).copy()
    for i, j in pixels:
        c = out[i, j]
        if fold >= 1:
            out[i, j] = c + rng.poisson(c * (fold - 1))
        else:
            out[i, j] = rng.binomial(int(c), fold)
        out[j, i] = out[i, j]
    return out
