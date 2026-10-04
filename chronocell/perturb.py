"""
Perturbations of a fitted population (ensemble), as predictions that can be checked against data.

1. Cohesin depletion (e.g. RAD21 degradation by auxin)
   From the untreated ensemble's median distance map d (any v3.3 or v4 result):
       log d'_ij = log d_ij + shift(s_ij) - lam * resid_ij
   resid_ij  = log d_ij - mean_{|k-l| = |i-j|} log d_kl     (the domain pattern beyond the trend)
   shift(s)  = c0 + c1 L + c2 L^2,  L = log10(s / 1 bp), clamped to the separations it was fitted on
   lam is the fraction of the domain pattern that cohesin loss removes. lam and c0..c2 are FITTED on
   one practice pair (HCT116 chr21:28-30 Mb, untreated vs 6 h auxin; Bintu et al. 2018) by
   validation/perturbation.py and stored, with their provenance, in
   chronocell/data/perturbation_params.json. They are not first-principles constants. The held-out
   test (HCT116 chr21:34-37 Mb) is reported in validation/RESULTS.md (Gate 4).

2. Structural variants as rearranged pieces of the fitted ensemble (exact, covariance space)
   A rearranged chromosome is a sequence of pieces of the reference chain (kept, reversed, copied, or
   from a partner). Each piece keeps its internal shape law, and the shapes of all native pieces keep
   their joint law (cross-piece correlations included); only the displacement across each junction is
   replaced by an independent bond with the ensemble's median adjacent variance s_b. For beads i in
   piece k and j in piece l > k, with x the reference positions:
       y_j - y_i = (x_exit_k - x_i) + b_k + sum_{k<q<l} [(x_exit_q - x_entry_q) + b_q] + (x_j - x_entry_l)
   a linear map of the reference ensemble plus independent bonds, so the derived ensemble is an exact,
   valid Gaussian computed from the reference pair-variance matrix alone (covariance -1/2 J S J), with
   no matrix inversion. For an ideal chain it reduces to concatenation (s'_ij = s_i,a-1 + s_b + s_b,j
   across a deletion). A tandem copy has the segment's internal law and is independent of everything
   else; a partner segment without data is a homogeneous chain with this ensemble's mean variance per
   separation. Hi-C-like output in reference coordinates sums contacts over copies.
   (A first version rewired the fitted spring network, couplings = -precision; the precision of a
   fitted max-ent ensemble is dominated by inversion noise, |k| beyond 5 beads ~ the backbone k on a
   planted ideal chain, so that approach was dropped: validation/TUNING.md section 7.)
   This is a mechanism simulator. It is tested in validation/perturbation.py where data exist and
   labelled "not validated" everywhere else.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import ensemble as ENS

PARAMS_PATH = Path(__file__).with_name("data") / "perturbation_params.json"


# ======================================================================================
# 1. Cohesin depletion
# ======================================================================================
def trend_residual(log_d: np.ndarray, sep: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(trend, residual) of a log-distance map: per-separation mean over i < j, and what is left."""
    n = len(log_d)
    iu = np.triu_indices(n, 1)
    s = sep[iu]
    _, inv = np.unique(s, return_inverse=True)
    mean = np.bincount(inv, weights=log_d[iu]) / np.bincount(inv)
    trend = np.zeros((n, n))
    trend[iu] = mean[inv]
    trend = trend + trend.T
    resid = log_d - trend
    np.fill_diagonal(resid, 0.0)
    return trend, resid


@dataclass(frozen=True)
class CohesinParams:
    lam: float
    c: tuple[float, float, float]
    log_sep_range: tuple[float, float]
    fitted_on: str = ""
    note: str = ""

    def shift(self, sep_bp: np.ndarray) -> np.ndarray:
        L = np.log10(np.clip(np.asarray(sep_bp, dtype=np.float64), 10 ** self.log_sep_range[0], 10 ** self.log_sep_range[1]))
        return self.c[0] + self.c[1] * L + self.c[2] * L ** 2


def load_params() -> dict:
    try:
        return json.loads(PARAMS_PATH.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def cohesin_params() -> CohesinParams | None:
    p = load_params().get("cohesin_loss")
    if not p:
        return None
    return CohesinParams(float(p["lam"]), tuple(float(v) for v in p["c"]), tuple(float(v) for v in p["log_sep_range"]),
                         p.get("fitted_on", ""), p.get("note", ""))


def cohesin_loss_map(median_nm: np.ndarray, sep_bp: np.ndarray, params: CohesinParams) -> np.ndarray:
    """Predicted median distance map after cohesin depletion (nm), from the untreated map."""
    d = np.asarray(median_nm, dtype=np.float64)
    log_d = np.log(np.where(d > 0, d, np.nan))
    np.fill_diagonal(log_d, 0.0)
    _, resid = trend_residual(np.nan_to_num(log_d), sep_bp)
    out = np.exp(log_d + params.shift(sep_bp) - params.lam * resid)
    np.fill_diagonal(out, 0.0)
    return out


def project_to_ensemble(median_nm: np.ndarray, r_c_nm: float, iterations: int = 600, seed: int = 0):
    """Refit a valid Gaussian ensemble to a target median-distance map (maximum entropy), so a
    perturbed map becomes a population again (members, trajectories, uncertainty)."""
    d = np.asarray(median_nm, dtype=np.float64)
    sigma = d / (ENS.MAXWELL_MEDIAN * r_c_nm)
    f = ENS.contact_probability_from_sigma(np.where(sigma > 0, sigma, 1e-6))
    np.fill_diagonal(f, np.nan)
    return ENS.fit_ensemble(f, 10_000.0, r_c_nm=r_c_nm,
                            cfg=ENS.EnsembleConfig(iterations=iterations, seed=seed, replicas=50, frames=10))


# ======================================================================================
# 2. Structural variants: rearranged pieces of the fitted ensemble
# ======================================================================================
NATIVE, COPY, PARTNER, INVERTED = 0, 1, 2, 3


def pair_variance_from_result(res: ENS.EnsembleResult) -> tuple[np.ndarray, float]:
    """(per-axis pair variance in units of r_c^2, r_c in nm) of any v3.3 / v4 result, exactly."""
    r_c = float(res.config.get("r_c_nm", 150.0))
    if getattr(res, "model", None) is not None:
        return res.model.dense_variance(dtype=np.float64), r_c
    s = (np.asarray(res.median_distance_nm, dtype=np.float64) / (ENS.MAXWELL_MEDIAN * r_c)) ** 2
    np.fill_diagonal(s, 0.0)
    return s, r_c


@dataclass
class Piece:
    beads: np.ndarray        # bead indices in this piece's source, in derived order
    source: str              # "native" | "copy" | "partner"
    kind: int                # NATIVE | COPY | PARTNER | INVERTED


def deletion_pieces(n: int, segments) -> tuple[list[Piece], str]:
    """Kept runs between one or more deleted half-open bead intervals [a, b), in order."""
    drop = np.zeros(n, bool)
    for a, b in segments:
        a, b = max(0, int(a)), min(n, int(b))
        if b > a:
            drop[a:b] = True
    if not drop.any():
        raise ValueError("Empty deletion.")
    if drop.all():
        raise ValueError("The deletion removes every bead.")
    keep = np.flatnonzero(~drop)
    runs = np.split(keep, np.flatnonzero(np.diff(keep) > 1) + 1)
    edges = np.flatnonzero(np.diff(np.r_[0, drop.astype(np.int8), 0]))
    spans = ", ".join(f"{a}-{b - 1}" for a, b in zip(edges[::2], edges[1::2]))
    return [Piece(r, "native", NATIVE) for r in runs], f"deletion of beads {spans}"


def pieces_for(op: str, n: int, params: dict) -> tuple[list[Piece], str]:
    if op == "deletion" and "segments" in params:
        return deletion_pieces(n, params["segments"])
    if op == "deletion":
        a, b = max(1, int(params["a"])), min(n - 1, int(params["b"]))
        if b <= a:
            raise ValueError("Empty deletion.")
        return [Piece(np.arange(a), "native", NATIVE), Piece(np.arange(b, n), "native", NATIVE)], \
            f"deletion of beads {a}-{b - 1}"
    if op == "inversion":
        a, b = max(1, int(params["a"])), min(n - 1, int(params["b"]))
        if b - a < 2:
            raise ValueError("An inversion needs at least two beads.")
        return [Piece(np.arange(a), "native", NATIVE), Piece(np.arange(b - 1, a - 1, -1), "native", INVERTED),
                Piece(np.arange(b, n), "native", NATIVE)], f"inversion of beads {a}-{b - 1}"
    if op == "duplication":
        a, b = max(0, int(params["a"])), min(n, int(params["b"]))
        if b <= a:
            raise ValueError("Empty duplication.")
        ps = [Piece(np.arange(b), "native", NATIVE), Piece(np.arange(a, b), "copy", COPY)]
        if b < n:
            ps.append(Piece(np.arange(b, n), "native", NATIVE))
        return ps, f"tandem duplication of beads {a}-{b - 1}"
    if op == "translocation":
        bp = int(np.clip(params["breakpoint"], 2, n - 1))
        m = int(params.get("partner_beads", 50))
        return [Piece(np.arange(bp), "native", NATIVE), Piece(np.arange(m), "partner", PARTNER)], \
            f"translocation at bead {bp} to a {m}-bead partner segment"
    raise ValueError(f"Unknown operation '{op}'.")


def _homogeneous(S: np.ndarray, m: int) -> np.ndarray:
    """Pair variances of a homogeneous chain with S's mean variance at each separation."""
    n = len(S)
    mean = np.array([np.mean(np.diagonal(S, k)) if k < n else np.nan for k in range(m)])
    if m > n:                                              # beyond the window: extend linearly in separation
        slope = (mean[n - 1] - mean[max(1, n // 2)]) / max(1, n - 1 - max(1, n // 2))
        mean[n:] = mean[n - 1] + slope * (np.arange(n, m) - (n - 1))
    sep = np.abs(np.subtract.outer(np.arange(m), np.arange(m)))
    out = mean[sep]
    np.fill_diagonal(out, 0.0)
    return out


def derive(S: np.ndarray, pieces: list[Piece], partner_S: np.ndarray | None = None,
           s_bond: float | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pair-variance matrix of the rearranged chain (same units as S), with each derived bead's
    reference bead (-1 for partner beads) and kind."""
    n = len(S)
    s_b = float(np.median(np.diag(S, 1))) if s_bond is None else float(s_bond)
    blocks = [-0.5 * _centre(S)]                         # reference covariance, exact from S
    offsets = {"native": 0}
    ext = n
    for p in pieces:
        if p.source == "copy":
            seg = np.unique(p.beads)
            blocks.append(-0.5 * _centre(S[np.ix_(seg, seg)]))
            offsets[id(p)] = ext - int(seg.min())          # copy bead k -> ext + (k - seg.min())
            ext += len(seg)
        elif p.source == "partner":
            Sp = partner_S if partner_S is not None else _homogeneous(S, len(p.beads))
            blocks.append(-0.5 * _centre(Sp[: len(p.beads), : len(p.beads)]))
            offsets[id(p)] = ext
            ext += len(p.beads)
    sigma = np.zeros((ext, ext))
    pos = 0
    for blk in blocks:
        k = len(blk)
        sigma[pos:pos + k, pos:pos + k] = blk
        pos += k

    def col(p: Piece, b: np.ndarray) -> np.ndarray:
        return (b + (0 if p.source == "native" else offsets[id(p)])).astype(np.int64)

    total = sum(len(p.beads) for p in pieces)
    V = np.zeros((total, ext))
    junction = np.zeros(total)
    origin = np.full(total, -1, np.int64)
    kind = np.zeros(total, np.int8)
    base = np.zeros(ext)
    r = 0
    for k, p in enumerate(pieces):
        c = col(p, p.beads)
        m = len(c)
        V[r:r + m] = base[None, :]
        V[np.arange(r, r + m), c] += 1.0
        V[r:r + m, c[0]] -= 1.0
        junction[r:r + m] = k
        if p.source != "partner":
            origin[r:r + m] = p.beads
        kind[r:r + m] = p.kind
        base = base.copy()
        base[c[-1]] += 1.0
        base[c[0]] -= 1.0
        r += m
    G = V @ sigma @ V.T
    g = np.diag(G)
    Sd = np.clip(g[:, None] + g[None, :] - 2.0 * G, 0.0, None) + s_b * np.abs(junction[:, None] - junction[None, :])
    np.fill_diagonal(Sd, 0.0)
    return Sd, origin, kind


def _centre(S: np.ndarray) -> np.ndarray:
    S = np.asarray(S, dtype=np.float64)
    return S - S.mean(0, keepdims=True) - S.mean(1, keepdims=True) + S.mean()


def to_reference(Sd: np.ndarray, origin: np.ndarray, n_ref: int, r_c_nm: float) -> tuple[np.ndarray, np.ndarray]:
    """Contact probability (summed over copies, as Hi-C reads mapped to the reference would be) and
    median distance (averaged over copies) in reference bead coordinates; NaN where a bead is absent."""
    sigma = np.sqrt(Sd)
    p = ENS.contact_probability_from_sigma(np.where(sigma > 0, sigma, 1e-9))
    d = ENS.MAXWELL_MEDIAN * sigma * r_c_nm
    ok = origin >= 0
    P = np.zeros((n_ref, n_ref))
    D = np.zeros((n_ref, n_ref))
    C = np.zeros((n_ref, n_ref))
    ii, jj = np.nonzero(ok[:, None] & ok[None, :])
    keep = ii != jj
    ii, jj = ii[keep], jj[keep]
    np.add.at(P, (origin[ii], origin[jj]), p[ii, jj])
    np.add.at(D, (origin[ii], origin[jj]), d[ii, jj])
    np.add.at(C, (origin[ii], origin[jj]), 1.0)
    present = C > 0
    P = np.where(present, P, np.nan)
    D = np.where(present, D / np.maximum(C, 1), np.nan)
    np.fill_diagonal(P, np.nan)
    return P, D


@dataclass
class VariantImpact:
    description: str
    log2_fc: np.ndarray                     # (n, n) log2 contact-probability change, reference coordinates
    p_before: np.ndarray
    p_after: np.ndarray
    d_before_nm: np.ndarray
    d_after_nm: np.ndarray
    interval: dict | None = None            # bootstrap percentiles of log2_fc, if computed
    notes: list[str] = field(default_factory=list)

    def top_changes(self, k: int = 20, min_sep: int = 2) -> list[dict]:
        n = len(self.log2_fc)
        iu = np.triu_indices(n, min_sep)
        v = self.log2_fc[iu]
        ok = np.isfinite(v)
        order = np.argsort(-np.abs(np.where(ok, v, 0.0)))[:k]
        out = []
        for t in order:
            if not ok[t]:
                continue
            i, j = int(iu[0][t]), int(iu[1][t])
            row = {"i": i, "j": j, "log2_fc": float(self.log2_fc[i, j]), "p_before": float(self.p_before[i, j]),
                   "p_after": float(self.p_after[i, j])}
            if self.interval is not None:
                row["log2_fc_90ci"] = (float(self.interval["lo"][i, j]), float(self.interval["hi"][i, j]))
            out.append(row)
        return out


def variant_impact_from_variance(S: np.ndarray, r_c_nm: float, op: str, params: dict,
                                 partner_S: np.ndarray | None = None) -> VariantImpact:
    n = len(S)
    pieces, desc = pieces_for(op, n, params)
    Sd, origin, _ = derive(S, pieces, partner_S)
    p_after, d_after = to_reference(Sd, origin, n, r_c_nm)
    sig = np.sqrt(np.clip(S, 0.0, None))
    p_before = ENS.contact_probability_from_sigma(np.where(sig > 0, sig, 1e-9))
    np.fill_diagonal(p_before, np.nan)
    d_before = ENS.MAXWELL_MEDIAN * sig * r_c_nm
    with np.errstate(divide="ignore", invalid="ignore"):
        fc = np.log2(p_after / p_before)
    return VariantImpact(desc, fc, p_before, p_after, d_before, d_after)


def variant_impact(res: ENS.EnsembleResult, op: str, params: dict) -> VariantImpact:
    """Apply one rearrangement to a fitted population and compare contacts before and after."""
    S, r_c = pair_variance_from_result(res)
    return variant_impact_from_variance(S, r_c, op, params)


def bootstrap_impact(fit_fn, counts: tuple[np.ndarray, np.ndarray, np.ndarray], op: str, params: dict,
                     reps: int = 8, seed: int = 0, level: float = 0.9) -> dict:
    """Uncertainty of a predicted change from the input data: refit on Poisson-resampled counts
    `reps` times, apply the same variant, and take percentiles of log2 fold change per pair.
    fit_fn(ci, cj, cm) -> EnsembleResult."""
    rng = np.random.default_rng(seed)
    ci, cj, cm = counts
    fcs = []
    for _ in range(reps):
        cm_b = rng.poisson(np.asarray(cm, dtype=np.float64)).astype(np.float64)
        keep = cm_b > 0
        res = fit_fn(ci[keep], cj[keep], cm_b[keep])
        fcs.append(variant_impact(res, op, params).log2_fc)
    arr = np.stack(fcs)
    q = (1 - level) / 2
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return {"lo": np.nanquantile(arr, q, axis=0), "hi": np.nanquantile(arr, 1 - q, axis=0),
                "median": np.nanmedian(arr, axis=0), "reps": reps, "level": level}


# ======================================================================================
# Genes and enhancer-promoter pairs touched by a variant
# ======================================================================================
def enhancer_promoter_pairs(impact: VariantImpact, promoters: list[tuple[str, int]], enhancers: list[int],
                            min_p: float = 0.01, min_abs_log2: float = 0.5) -> list[dict]:
    """Promoter (gene, bead) x enhancer (bead) pairs whose predicted contact changes by at least
    min_abs_log2 (log2) with contact probability >= min_p before or after; pairs that lose a partner
    to a deletion are listed as 'removed'."""
    out = []
    for gene, pb in promoters:
        for eb in enhancers:
            if abs(eb - pb) < 2:
                continue
            fc = impact.log2_fc[pb, eb]
            if not np.isfinite(fc):
                if np.isfinite(impact.p_before[pb, eb]) and np.isnan(impact.p_after[pb, eb]):
                    out.append({"gene": gene, "promoter_bead": pb, "enhancer_bead": eb, "log2_fc": float("-inf"),
                                "p_before": float(impact.p_before[pb, eb]), "p_after": 0.0, "change": "removed"})
                continue
            pmax = max(impact.p_before[pb, eb], impact.p_after[pb, eb])
            if abs(fc) >= min_abs_log2 and pmax >= min_p:
                row = {"gene": gene, "promoter_bead": pb, "enhancer_bead": eb, "log2_fc": float(fc),
                       "p_before": float(impact.p_before[pb, eb]), "p_after": float(impact.p_after[pb, eb]),
                       "change": "gained" if fc > 0 else "lost"}
                if impact.interval is not None:
                    row["log2_fc_90ci"] = (float(impact.interval["lo"][pb, eb]), float(impact.interval["hi"][pb, eb]))
                out.append(row)
    return sorted(out, key=lambda r: -abs(r["log2_fc"]) if math.isfinite(r["log2_fc"]) else -1e9)
