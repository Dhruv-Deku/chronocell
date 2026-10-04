"""
Self-Math PDB State Evaluator: pure mathematics on 3D coordinates.

Input: any set of ordered points (beads of a chromatin model, or the atoms of an uploaded PDB in file
order). Nothing else is used: no contact data, no labels, no stored answer. Every number below is
derived from the coordinates by the formula next to it.

Metrics
  Radius of gyration     R_g = sqrt( (1/N) sum_i |r_i - r_cm|^2 )
  Gyration tensor        S = (1/N) sum_i (r_i - r_cm)(r_i - r_cm)^T,  eigenvalues l1 >= l2 >= l3
  Asphericity            Delta = (3/2) sum_k (l_k - l_mean)^2 / (Tr S)^2      in [0, 1]: 0 sphere, 1 rod
                         (algebraically equal to the relative shape anisotropy kappa^2)
  Acylindricity          c = l2 - l3   (nm^2)
  Packing density        rho = N / V_hull  (beads per um^3; V_hull = convex-hull volume, Qhull)
  Distance-decay slope   gamma_d = OLS slope of log d_ij on log |i - j| over every pair with
                         4 <= |i - j| <= N/10 (a uniform sample of 200,000 when there are more): the
                         scaling window of physics.distance_scaling; below it bond geometry, above it
                         finite-size saturation dominate. For a polymer this is the spatial scaling
                         exponent nu (1/3 compact globule, 1/2 ideal chain).
  Contact-scale exponent gamma_c = alpha * gamma_d with alpha = 3: the contact law d ~ I^(-1/alpha)
                         used everywhere in ChronoCell turns d ~ s^nu into I ~ s^(-alpha nu). The
                         state thresholds below are on this contact scale (fractal globule:
                         gamma_c = 1; ideal chain: 1.5).
  Short-range exponent   gamma_c on pairs with 2 <= |i - j| <= 10
  Scaling break          |gamma_d(lower half) - gamma_d(upper half)| of the scaling window, split at
                         its geometric midpoint: a single power law has no break
  Local density          non-bonded beads within 1.5 b of each bead (b = median bond); robust
                         z-scores (median / 1.4826 MAD); "spikes" are beads with z > 3.5
  Local anomalies        runs of >= 3 consecutive beads whose mean distance to their 10 nearest chain
                         neighbours deviates from the structure's own scaling law by > 3 robust z

Classification (rule set from the ChronoCell v4 specification; every threshold is shown with its
source and is NOT fitted or validated on labelled disease or senescence structures, because no such
labelled 3D data were available; treat the label as a description of the geometry, not a diagnosis)
  Senescent  at least 2 of: short-range gamma_c > 1.3 [specification]; R_g < 0.85 x R_g,ref
             [assumption]; packing density > 1.2 x rho_ref [assumption]
  Diseased   at least 2 of: scaling break > 0.15 [assumption]; density-spike fraction > 1 %
             [assumption]; >= 1 local anomaly run [assumption]
  Normal     gamma_c within [0.75, 1.10] [specification] and fewer than 2 flags of either kind
  Indeterminate  none of the above
  R_g,ref and rho_ref are the values of a compact globule of the same N and bond length at the volume
  fraction 0.19 measured on ChronoCell's reference fractal globule (agent.REFERENCE_PACKING):
  R_g,ref = sqrt(3/5) b (N / (8 x 0.19))^(1/3); rho_ref = N / ((4 pi / 3) (R_g,ref / sqrt(3/5))^3).
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

import numpy as np

ALPHA = 3.0
REF_PACKING = 0.19
WINDOW_MIN = 4          # scaling window [4, N/10] beads, as chronocell.physics.distance_scaling
SHORT_MAX = 10          # "short-range" = 2 ... 10 beads


@dataclass(frozen=True)
class Thresholds:
    gamma_normal: tuple[float, float] = (0.75, 1.10)     # specification
    gamma_senescent: float = 1.30                        # specification
    rg_depressed: float = 0.85                           # assumption
    density_elevated: float = 1.20                       # assumption
    scaling_break: float = 0.15                          # assumption
    spike_fraction: float = 0.01                         # assumption
    spike_z: float = 3.5                                 # assumption
    anomaly_z: float = 3.0                               # assumption
    anomaly_run: int = 3                                 # assumption


SOURCES = {"gamma_normal": "specification", "gamma_senescent": "specification", "rg_depressed": "assumption",
           "density_elevated": "assumption", "scaling_break": "assumption", "spike_fraction": "assumption",
           "spike_z": "assumption", "anomaly_z": "assumption", "anomaly_run": "assumption"}


@dataclass
class Criterion:
    state: str
    name: str
    value: float
    threshold: str
    met: bool
    source: str


@dataclass
class Evaluation:
    n: int
    unit_note: str
    rg_nm: float
    eigenvalues_nm2: tuple[float, float, float]
    asphericity: float
    acylindricity_nm2: float
    hull_volume_um3: float
    packing_density_per_um3: float
    bond_nm: float
    gamma_d: float
    gamma_d_r2: float
    gamma_c: float
    gamma_c_short: float
    gamma_d_short: float
    gamma_d_long: float
    scaling_break: float
    rg_ref_nm: float
    density_ref_per_um3: float
    density_z: np.ndarray
    spike_fraction: float
    anomaly_runs: list[tuple[int, int]]
    decay_s: np.ndarray
    decay_d: np.ndarray
    state: str
    criteria: list[Criterion] = field(default_factory=list)

    def summary(self) -> dict:
        out = {k: v for k, v in asdict(self).items() if k not in ("density_z", "decay_s", "decay_d", "criteria")}
        out["criteria"] = [asdict(c) for c in self.criteria]
        return out


def gyration(x: np.ndarray) -> tuple[float, np.ndarray]:
    d = x - x.mean(axis=0)
    S = d.T @ d / len(x)
    lam = np.sort(np.linalg.eigvalsh(S))[::-1]
    return float(np.sqrt(max(lam.sum(), 0.0))), lam


def asphericity(lam: np.ndarray) -> float:
    tr = float(lam.sum())
    if tr <= 0:
        return float("nan")
    return float(1.5 * np.sum((lam - lam.mean()) ** 2) / tr ** 2)


def hull_volume(x: np.ndarray) -> float:
    """Convex-hull volume (same units cubed); NaN for fewer than 4 points or flat sets."""
    if len(x) < 4:
        return float("nan")
    try:
        from scipy.spatial import ConvexHull
        return float(ConvexHull(x).volume)
    except Exception:          # coplanar / degenerate input, or scipy missing
        return float("nan")


def _pairs_in(n: int, s_min: int, s_max: int, max_pairs: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """All pairs with s_min <= j - i <= s_max, or a uniform sample of max_pairs of them."""
    s_max = min(s_max, n - 1)
    counts = n - np.arange(s_min, s_max + 1)
    total = int(counts.sum()) if counts.size else 0
    if total == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    if total <= max_pairs:
        s = np.repeat(np.arange(s_min, s_max + 1), counts)
        i = np.concatenate([np.arange(c) for c in counts])
        return i, i + s
    pick = rng.choice(total, max_pairs, replace=False)
    edges = np.concatenate([[0], np.cumsum(counts)])
    k = np.searchsorted(edges, pick, side="right") - 1
    i = pick - edges[k]
    return i, i + s_min + k


def decay_fit(x: np.ndarray, max_pairs: int = 200_000, seed: int = 0, s_min: int = 1, s_max: int | None = None):
    """OLS of log d_ij on log |i - j| over pairs with s_min <= |i - j| <= s_max (all of them, or a
    uniform sample), with binned means for plotting."""
    n = len(x)
    rng = np.random.default_rng(seed)
    i, j = _pairs_in(n, int(s_min), int(s_max if s_max is not None else n - 1), max_pairs, rng)
    s = (j - i).astype(np.float64)
    d = np.linalg.norm(x[i] - x[j], axis=1)
    ok = d > 0
    ls, ld = np.log(s[ok]), np.log(d[ok])
    if ls.size < 3 or np.ptp(ls) == 0:
        return float("nan"), float("nan"), np.empty(0), np.empty(0)
    A = np.vstack([ls, np.ones_like(ls)]).T
    coef, *_ = np.linalg.lstsq(A, ld, rcond=None)
    pred = A @ coef
    r2 = 1.0 - float(np.sum((ld - pred) ** 2)) / float(np.sum((ld - ld.mean()) ** 2))
    edges = np.unique(np.round(np.geomspace(1, max(2, s.max()), 30)))
    sb = np.digitize(s[ok], edges)
    xs = np.array([np.exp(ls[sb == k].mean()) for k in np.unique(sb)])
    ys = np.array([np.exp(ld[sb == k].mean()) for k in np.unique(sb)])
    return float(coef[0]), r2, xs, ys


def local_density(x: np.ndarray, r: float) -> np.ndarray:
    from ..physics import local_density as ld
    return ld(x, r)


def _robust_z(v: np.ndarray) -> np.ndarray:
    med = np.median(v)
    mad = np.median(np.abs(v - med)) * 1.4826
    if mad <= 0:
        sd = float(np.std(v))
        return (v - med) / sd if sd > 0 else np.zeros_like(v)
    return (v - med) / mad


def _runs(mask: np.ndarray, min_len: int) -> list[tuple[int, int]]:
    out = []
    i, n = 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            if j - i >= min_len:
                out.append((i, j))
            i = j
        else:
            i += 1
    return out


def evaluate(coords_nm: np.ndarray, thresholds: Thresholds | None = None, unit_note: str = "nm") -> Evaluation:
    """Evaluate one structure. Coordinates must be in nm (see `coords_from_pdb_text` for unit handling)."""
    t = thresholds or Thresholds()
    x = np.asarray(coords_nm, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 3 or len(x) < 8:
        raise ValueError("Need an (N, 3) array with at least 8 points.")
    if not np.isfinite(x).all():
        raise ValueError("Coordinates contain NaN or infinite values.")
    n = len(x)
    rg, lam = gyration(x)
    vol_nm3 = hull_volume(x)
    vol_um3 = vol_nm3 / 1e9 if np.isfinite(vol_nm3) else float("nan")
    density = n / vol_um3 if vol_um3 and np.isfinite(vol_um3) and vol_um3 > 0 else float("nan")
    bonds = np.linalg.norm(np.diff(x, axis=0), axis=1)
    b = float(np.median(bonds)) if bonds.size else float("nan")
    lo_w, hi_w = WINDOW_MIN, max(2 * WINDOW_MIN, n // 10)      # scaling window, as physics.distance_scaling
    g_d, r2, _, _ = decay_fit(x, s_min=lo_w, s_max=hi_w)
    _, _, xs, ys = decay_fit(x)                                   # all separations, for the plot
    g_short, _, _, _ = decay_fit(x, s_min=2, s_max=SHORT_MAX)
    mid = int(round(math.sqrt(lo_w * hi_w)))                     # halves of the window in log separation
    g_lo, _, _, _ = decay_fit(x, s_min=lo_w, s_max=mid)
    g_hi, _, _, _ = decay_fit(x, s_min=mid, s_max=hi_w)
    brk = abs(g_lo - g_hi) if np.isfinite(g_lo) and np.isfinite(g_hi) else float("nan")
    g_c, g_c_short = ALPHA * g_d, ALPHA * g_short
    # compact-globule references at the reference globule's volume fraction
    rg_ref = math.sqrt(3.0 / 5.0) * b * (n / (8.0 * REF_PACKING)) ** (1.0 / 3.0) if np.isfinite(b) else float("nan")
    r_sphere_um = rg_ref / math.sqrt(3.0 / 5.0) / 1000.0
    dens_ref = n / ((4.0 * math.pi / 3.0) * r_sphere_um ** 3) if r_sphere_um > 0 else float("nan")
    # local density spikes
    dens = local_density(x, 1.5 * b) if np.isfinite(b) and b > 0 else np.zeros(n)
    z = _robust_z(dens)
    spike_frac = float(np.mean(z > t.spike_z))
    # local anomalies against the structure's own scaling law
    k = 10
    lo = np.clip(np.arange(n) - k, 0, n - 1)
    hi = np.clip(np.arange(n) + k, 0, n - 1)
    mean_nb = np.array([np.mean(np.linalg.norm(x[lo[q]:hi[q] + 1] - x[q], axis=1)) for q in range(n)])
    if np.isfinite(g_d):
        expect = np.array([np.mean(np.exp(np.log(b) + g_d * np.log(np.maximum(np.abs(np.arange(lo[q], hi[q] + 1) - q), 1))))
                           for q in range(n)])
        dev = np.log(mean_nb / np.maximum(expect, 1e-9))
        za = _robust_z(dev)
        runs = _runs(np.abs(za) > t.anomaly_z, t.anomaly_run)
    else:
        runs = []
    crit = [
        Criterion("Senescent", "short-range contact-scale exponent gamma_c(short)", g_c_short,
                  f"> {t.gamma_senescent}", bool(g_c_short > t.gamma_senescent), SOURCES["gamma_senescent"]),
        Criterion("Senescent", "R_g / R_g,ref (compact globule, same N and b)", rg / rg_ref if rg_ref else float("nan"),
                  f"< {t.rg_depressed}", bool(rg_ref and rg / rg_ref < t.rg_depressed), SOURCES["rg_depressed"]),
        Criterion("Senescent", "packing density / reference", density / dens_ref if dens_ref else float("nan"),
                  f"> {t.density_elevated}", bool(np.isfinite(density) and density / dens_ref > t.density_elevated),
                  SOURCES["density_elevated"]),
        Criterion("Diseased", "scaling break |gamma_d(short) - gamma_d(long)|", brk, f"> {t.scaling_break}",
                  bool(np.isfinite(brk) and brk > t.scaling_break), SOURCES["scaling_break"]),
        Criterion("Diseased", f"density-spike fraction (z > {t.spike_z})", spike_frac, f"> {t.spike_fraction}",
                  bool(spike_frac > t.spike_fraction), SOURCES["spike_fraction"]),
        Criterion("Diseased", f"local anomaly runs (|z| > {t.anomaly_z}, >= {t.anomaly_run} beads)", float(len(runs)),
                  ">= 1", bool(len(runs) >= 1), SOURCES["anomaly_run"]),
        Criterion("Normal", "contact-scale exponent gamma_c", g_c,
                  f"in [{t.gamma_normal[0]}, {t.gamma_normal[1]}]",
                  bool(t.gamma_normal[0] <= g_c <= t.gamma_normal[1]), SOURCES["gamma_normal"]),
    ]
    sen = sum(c.met for c in crit if c.state == "Senescent")
    dis = sum(c.met for c in crit if c.state == "Diseased")
    if sen >= 2:
        state = "Senescent"
    elif dis >= 2:
        state = "Diseased"
    elif crit[-1].met:
        state = "Normal"
    else:
        state = "Indeterminate"
    return Evaluation(n, unit_note, rg, tuple(float(v) for v in lam), asphericity(lam), float(lam[1] - lam[2]),
                      vol_um3, density, b, g_d, r2, g_c, g_c_short, g_short, g_hi, brk, rg_ref, dens_ref, z, spike_frac,
                      runs, xs, ys, state, crit)


def coords_from_pdb_text(text: str) -> tuple[list[np.ndarray], str]:
    """All models of a PDB file as nm coordinates, and a note on the unit used.

    ChronoCell PDBs declare nm (REMARK 250) and are read exactly; any other PDB is in Angstrom by the
    wwPDB convention and is divided by 10. Every ATOM/HETATM record is a point, in file order."""
    models: list[list[tuple[float, float, float]]] = [[]]
    unit_nm = offset = None
    import re
    for line in text.splitlines():
        if line.startswith("REMARK 250  COORDINATE UNIT"):
            m = re.search(r":\s*([0-9.eE+-]+)\s*NM", line)
            unit_nm = float(m.group(1)) if m else None
        elif line.startswith("REMARK 250  ORIGIN OFFSET"):
            try:
                offset = np.array([float(v) for v in line.split(":", 1)[1].split()])
            except ValueError:
                offset = None
        elif line.startswith("ENDMDL"):
            models.append([])
        elif line.startswith(("ATOM  ", "HETATM")):
            try:
                models[-1].append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    models = [np.asarray(m, dtype=np.float64) for m in models if len(m)]
    if not models:
        raise ValueError("No ATOM/HETATM records found.")
    if unit_nm is not None and offset is not None and offset.size == 3:
        return [m * unit_nm - offset for m in models], "nm (declared in REMARK 250)"
    return [m / 10.0 for m in models], "Angstrom (wwPDB convention), converted to nm"


def evaluate_population(members_nm: np.ndarray, thresholds: Thresholds | None = None) -> dict:
    """Evaluate every member of an ensemble: the distribution of each metric and of the state labels."""
    evs = [evaluate(m, thresholds) for m in np.asarray(members_nm, dtype=np.float64)]
    keys = ("rg_nm", "asphericity", "acylindricity_nm2", "packing_density_per_um3", "gamma_d", "gamma_c")
    out = {k: {"median": float(np.nanmedian([getattr(e, k) for e in evs])),
               "p10": float(np.nanpercentile([getattr(e, k) for e in evs], 10)),
               "p90": float(np.nanpercentile([getattr(e, k) for e in evs], 90))} for k in keys}
    states = [e.state for e in evs]
    out["states"] = {s: states.count(s) for s in sorted(set(states))}
    out["members"] = len(evs)
    return out
