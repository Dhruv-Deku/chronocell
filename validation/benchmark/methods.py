"""
Method adapters for the benchmark. Every method gets the same input for a split and returns a
predicted median-distance map (nm); population models also return the per-pair spread (sigma, nm).

Input types
  imaging  contact frequencies of half A at the dataset's contact radius (+ copies observed per pair)
  hic      sequencing counts on the same loci (Rao et al. 2014), the same truth (half B medians)

Methods
  genomic_distance_only  power law in genomic separation fitted on half A's measured medians (the v3.3
                         baseline; note it sees half A's distances, which the other methods never do)
  no_3d                  each pair's frequency inverted on its own (Maxwell law), no 3D model
  v3_2_single            ChronoCell v3.2 single structure (chronocell.egnn.fit_structure)
  v3_3_windowed          ChronoCell v3.3 population model on the fewest equal windows of <= 400 loci
  v4_whole               ChronoCell v4 population model on all loci at once (frozen settings)
  pastis_mds, pastis_pm2 PASTIS 0.4.0 (Varoquaux et al., Bioinformatics 2014), its own code run unmodified
                         from the official source distribution (validation/data/tools; the package
                         __init__, which imports compiled extensions, is skipped). Its output has no
                         length unit: it is scaled so its median adjacent distance equals the adjacent
                         distance implied by the input frequencies (Maxwell inversion), a rule that
                         uses the input only.
Not run (reported with the reason): ShRec3D (MATLAB, no installable implementation), Chrom3D (C++
binary needing Boost and lamina-association data), 3DMax / LorDG (Java, not installed here).
"""

from __future__ import annotations

import sys
import tarfile
import time
import types
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

from chronocell import ensemble as E, population as P        # noqa: E402

TOOLS = ROOT / "data" / "tools"
PASTIS_SDIST = "https://files.pythonhosted.org/packages/64/a0/396cb8b6bbc987a5fdc9d95d2c4bc06972938fea265d18ab744532ea01dd/pastis-0.4.0.tar.gz"
PASTIS_SHA256 = "9c921e3d17f92a44147f31e70ade25dc5cdd0175d85b1410ff36b020202c1fdf"

NOT_RUN = {
    "ShRec3D": "Published as MATLAB code (Lesne et al., Nat Methods 2014); no installable implementation here. "
               "ChronoCell's v3.2 shortest-path MDS start is ShRec3D-style but is not the published tool.",
    "Chrom3D": "C++ program (Paulsen et al., Genome Biol 2017) needing Boost to build and lamina-association data; "
               "no compiler toolchain or LAD data here.",
    "3DMax / LorDG": "Java programs (Oluwadare et al.); no Java runtime on this machine.",
}


@dataclass
class Input:
    kind: str                    # "imaging" | "hic"
    n: int
    freq: np.ndarray             # contact probabilities (imaging: measured; hic: counts -> probability)
    seen: np.ndarray | float     # observations per pair (imaging) or N_eff (hic)
    counts: np.ndarray           # integer-like counts (imaging: copies in contact; hic: reads)
    r_c_nm: float
    median_a: np.ndarray         # half A measured medians (only the genomic baseline uses them)
    sep: np.ndarray
    tiles: list
    seed: int


@dataclass
class Prediction:
    median_nm: np.ndarray
    sigma_nm: np.ndarray | None = None
    seconds: float = 0.0
    contact_fit: float = float("nan")
    notes: list = field(default_factory=list)


def _adjacent_from_input(inp: Input) -> float:
    f = np.diag(inp.freq, 1)
    f = f[np.isfinite(f)]
    return float(np.median(E.gaussian_median_distance(np.clip(f, 1e-4, 0.999), inp.r_c_nm)))


def _dist(x: np.ndarray) -> np.ndarray:
    return np.linalg.norm(x[:, None] - x[None], axis=-1)


def genomic_distance_only(inp: Input) -> Prediction:
    import protocol as PR
    t = time.time()
    m = PR.genomic_baseline(inp.median_a, inp.sep, np.triu(np.ones((inp.n, inp.n), bool), 1))
    return Prediction(m, seconds=time.time() - t, notes=["fitted on half A's measured medians"])


def no_3d(inp: Input) -> Prediction:
    t = time.time()
    seen = np.maximum(np.asarray(inp.seen, dtype=float), 1.0)
    f = np.clip(np.nan_to_num(inp.freq, nan=0.0), 0.5 / seen, 1 - 0.5 / seen)
    m = E.gaussian_median_distance(f, inp.r_c_nm)
    np.fill_diagonal(m, 0.0)
    return Prediction(m, sigma_nm=m / E.MAXWELL_MEDIAN, seconds=time.time() - t)


def v3_2_single(inp: Input) -> Prediction:
    from chronocell import accuracy as ACC, egnn, physics
    t = time.time()
    i, j = np.triu_indices(inp.n, 1)
    c = inp.counts[i, j]
    keep = np.isfinite(c) & (c > 0)
    b0 = physics.bond_length_for(30_000)                   # the v3.2 pipeline's anchor, as in v3.3 validation
    feats = egnn.node_features(np.full(inp.n, 0.42), np.zeros(inp.n), np.ones(inp.n, bool))
    res = egnn.fit_structure(inp.n, feats, i[keep], j[keep], c[keep], egnn.FitConfig(seed=inp.seed), b0=b0)
    d = _dist(res.coords_nm)
    return Prediction(d, seconds=time.time() - t,
                      contact_fit=ACC.contact_fit_structure(res.coords_nm, i[keep], j[keep], c[keep]))


def v3_3_windowed(inp: Input) -> Prediction:
    t = time.time()
    m = np.full((inp.n, inp.n), np.nan)
    fits = []
    for a, b in inp.tiles:
        seen = inp.seen[a:b, a:b] if np.ndim(inp.seen) == 2 else inp.seen
        res = E.fit_ensemble(inp.freq[a:b, a:b], seen, r_c_nm=inp.r_c_nm, cfg=E.EnsembleConfig(seed=inp.seed))
        m[a:b, a:b] = res.median_distance_nm
        fits.append(res.contact_fit)
    return Prediction(m, sigma_nm=m / E.MAXWELL_MEDIAN, seconds=time.time() - t, contact_fit=float(np.mean(fits)))


def v4_whole(inp: Input) -> Prediction:
    from frozen import WHOLE_CHROMOSOME
    t = time.time()
    res = P.fit_population(inp.freq, inp.seen, r_c_nm=inp.r_c_nm,
                           cfg=P.PopulationConfig(seed=inp.seed, replicas=20, frames=2, **WHOLE_CHROMOSOME))
    iu = np.triu_indices(inp.n, 1)
    sig = np.zeros((inp.n, inp.n))
    sig[iu] = P.pair_sigma_nm(res, iu[0], iu[1])
    return Prediction(res.median_distance_nm.astype(np.float64), sigma_nm=sig + sig.T, seconds=time.time() - t,
                      contact_fit=res.contact_fit)


_PASTIS = None


def _pastis():
    """PASTIS 0.4.0's own optimisation modules, imported from the official sdist (verified by SHA-256)."""
    global _PASTIS
    if _PASTIS is not None:
        return _PASTIS
    import hashlib
    import urllib.request
    TOOLS.mkdir(parents=True, exist_ok=True)
    sdist = TOOLS / "pastis-0.4.0.tar.gz"
    if not sdist.exists():
        urllib.request.urlretrieve(PASTIS_SDIST, sdist)
    if hashlib.sha256(sdist.read_bytes()).hexdigest() != PASTIS_SHA256:
        raise RuntimeError("pastis-0.4.0.tar.gz does not match its published SHA-256")
    src = TOOLS / "pastis-0.4.0" / "pastis"
    if not src.exists():
        with tarfile.open(sdist) as tf:
            tf.extractall(TOOLS)
    if "pastis" not in sys.modules:
        pkg = types.ModuleType("pastis")
        pkg.__path__ = [str(src)]
        sys.modules["pastis"] = pkg
    from pastis.optimization import mds, poisson_structure
    _PASTIS = (mds, poisson_structure)
    return _PASTIS


def _pastis_counts(inp: Input) -> np.ndarray:
    c = np.nan_to_num(np.asarray(inp.counts, dtype=float), nan=0.0)
    c = np.round(c)
    np.fill_diagonal(c, 0.0)
    return c


def _scale_to_input(x: np.ndarray, inp: Input) -> np.ndarray:
    adj = float(np.median(np.linalg.norm(np.diff(x, axis=0), axis=1)))
    return x * (_adjacent_from_input(inp) / adj) if adj > 0 else x


def pastis_mds(inp: Input) -> Prediction:
    mds, _ = _pastis()
    t = time.time()
    x = mds.estimate_X(_pastis_counts(inp), alpha=-3.0, beta=1.0, random_state=inp.seed)
    x = _scale_to_input(x, inp)
    return Prediction(_dist(x), seconds=time.time() - t, notes=["PASTIS MDS, alpha = -3 (its default)"])


def pastis_pm2(inp: Input) -> Prediction:
    _, ps = _pastis()
    import contextlib
    import io
    t = time.time()
    model = ps.PM2(alpha=-3.0, beta=1.0, random_state=inp.seed, max_iter=5000)
    with contextlib.redirect_stdout(io.StringIO()):          # PM2 prints alpha / beta every outer iteration
        x = model.fit(_pastis_counts(inp))
    x = _scale_to_input(x, inp)
    alpha = float(np.ravel(model.alpha_)[0])
    return Prediction(_dist(x), seconds=time.time() - t, notes=[f"PASTIS PM2, fitted alpha = {alpha:.3f}"])


def learned_correction(inp: Input) -> Prediction:
    """Phase A4 (Gate 3b): the population model (v3.3 up to 400 loci, the frozen whole-chromosome v4 above,
    as in validation/phase_a.py) followed by the correction frozen in chronocell/data/learned_correction.json."""
    from chronocell import learned_correction as LC
    from frozen import WHOLE_CHROMOSOME
    t = time.time()
    if inp.n <= P.V33_MAX_BEADS:
        res = E.fit_ensemble(inp.freq, inp.seen, r_c_nm=inp.r_c_nm, cfg=E.EnsembleConfig(seed=inp.seed))
    else:
        res = P.fit_population(inp.freq, inp.seen, r_c_nm=inp.r_c_nm,
                               cfg=P.PopulationConfig(seed=inp.seed, **WHOLE_CHROMOSOME))
    base = np.asarray(res.median_distance_nm, dtype=np.float64)
    iu = np.triu_indices(inp.n, 1)
    sig = np.zeros((inp.n, inp.n))
    sig[iu] = P.pair_sigma_nm(res, iu[0], iu[1])
    sig = sig + sig.T
    params = LC.load()
    if not params or "kind" not in params:
        return Prediction(base, sigma_nm=sig, seconds=time.time() - t, contact_fit=res.contact_fit,
                          notes=["no correction frozen (practice chose none)"])
    seen = np.asarray(inp.seen, dtype=float)
    n_eff = float(np.nanmedian(seen[iu])) if seen.ndim == 2 else float(seen)
    step = float(np.median(np.diag(inp.sep, 1)))
    feats = LC.features(base, inp.freq, inp.counts, inp.sep, inp.r_c_nm, inp.kind == "hic", n_eff, step)
    corr = LC.correct(base, feats, params)
    factor = np.where(base > 0, corr / np.where(base > 0, base, 1.0), 1.0)
    return Prediction(corr, sigma_nm=sig * factor, seconds=time.time() - t, contact_fit=res.contact_fit,
                      notes=[f"learned correction ({params['kind']}), frozen on practice data"])


METHODS = {"genomic_distance_only": genomic_distance_only, "no_3d": no_3d, "v3_2_single": v3_2_single,
           "v3_3_windowed": v3_3_windowed, "v4_whole": v4_whole, "pastis_mds": pastis_mds, "pastis_pm2": pastis_pm2,
           "learned_correction": learned_correction}
OURS = ("v3_3_windowed", "v4_whole")
