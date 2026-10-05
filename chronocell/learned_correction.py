"""
Learned correction of a population model's median distances (Phase A4; validation/learned_correction.py
trains it, Gate 3b tests it). Inference only, NumPy: no PyTorch needed to apply a frozen correction.

For every pair i < j the correction predicts r_ij = log(true median / model median) from features
computed from the input and the fitted model only, and returns model * exp(r). Features:

  sep        log10 genomic separation (bp)
  model_res  log model median minus its mean at the same separation (the model's pattern)
  input_res  the same for the direct inversion of the input frequency (no 3D model)
  input_sm   input_res averaged over the 3 x 3 neighbouring pairs (less noisy)
  model_sm   model_res averaged likewise
  evidence   log(1 + contact count of the pair)
  hic        1 for sequencing Hi-C input, 0 for imaging-derived contacts
  depth      log10 of the input's effective depth (n_eff)
  step       log10 locus spacing (bp)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

FEATURES = ("sep", "model_res", "input_res", "input_sm", "model_sm", "evidence", "hic", "depth", "step")
PATH = Path(__file__).with_name("data") / "learned_correction.json"


def _trend_residual(v: np.ndarray, sep: np.ndarray, mask: np.ndarray) -> np.ndarray:
    out = np.full(v.shape, np.nan)
    s = sep[mask]
    x = v[mask]
    ok = np.isfinite(x)
    _, inv = np.unique(s, return_inverse=True)
    mean = np.bincount(inv[ok], weights=x[ok], minlength=inv.max() + 1) / np.maximum(
        np.bincount(inv[ok], minlength=inv.max() + 1), 1)
    out[mask] = x - mean[inv]
    return out


def _smooth(m: np.ndarray) -> np.ndarray:
    n = len(m)
    pad = np.pad(np.where(np.isfinite(m), m, 0.0), 1)
    cnt = np.pad(np.isfinite(m).astype(float), 1)
    s = sum(pad[1 + a:1 + a + n, 1 + b:1 + b + n] for a in (-1, 0, 1) for b in (-1, 0, 1))
    c = sum(cnt[1 + a:1 + a + n, 1 + b:1 + b + n] for a in (-1, 0, 1) for b in (-1, 0, 1))
    return np.where(c > 0, s / np.maximum(c, 1), np.nan)


def features(model_median: np.ndarray, freq: np.ndarray, counts: np.ndarray, sep_bp: np.ndarray, r_c_nm: float,
             hic: bool, n_eff: float, step_bp: float) -> np.ndarray:
    """(n, n, len(FEATURES)) feature array; entries off the upper triangle are NaN."""
    from . import ensemble as E
    n = len(model_median)
    upper = np.triu(np.ones((n, n), bool), 1)
    sep = np.maximum(np.asarray(sep_bp, dtype=np.float64), 1.0)
    logm = np.log(np.where(model_median > 0, model_median, np.nan))
    f = np.asarray(freq, dtype=np.float64)
    inv = np.log(E.gaussian_median_distance(np.clip(np.nan_to_num(f, nan=1e-4), 1e-4, 0.999), r_c_nm))
    inv = np.where(np.isfinite(f), inv, np.nan)
    model_res = _trend_residual(logm, sep, upper)
    input_res = _trend_residual(inv, sep, upper)
    out = np.full((n, n, len(FEATURES)), np.nan)
    out[..., 0] = np.log10(sep)
    out[..., 1] = model_res
    out[..., 2] = input_res
    full_in = np.where(upper, input_res, 0) + np.where(upper, input_res, 0).T
    full_mo = np.where(upper, model_res, 0) + np.where(upper, model_res, 0).T
    out[..., 3] = _smooth(np.where(np.isfinite(full_in) & ~np.eye(n, dtype=bool), full_in, np.nan))
    out[..., 4] = _smooth(np.where(np.isfinite(full_mo) & ~np.eye(n, dtype=bool), full_mo, np.nan))
    out[..., 5] = np.log1p(np.nan_to_num(np.asarray(counts, dtype=np.float64)))
    out[..., 6] = float(hic)
    out[..., 7] = np.log10(max(float(n_eff), 1.0))
    out[..., 8] = np.log10(max(float(step_bp), 1.0))
    out[~upper] = np.nan
    return out


def _mlp(params: dict, x: np.ndarray) -> np.ndarray:
    h = (x - np.asarray(params["mu"])) / np.asarray(params["sd"])
    h = np.where(np.isfinite(h), h, 0.0)
    layers = params["layers"]
    for k, layer in enumerate(layers):
        h = h @ np.asarray(layer["W"]).T + np.asarray(layer["b"])
        if k < len(layers) - 1:
            h = h / (1.0 + np.exp(-h))                     # SiLU
    return h[:, 0]


def load(path: Path = PATH) -> dict | None:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def correct(model_median: np.ndarray, feats: np.ndarray, params: dict) -> np.ndarray:
    """Corrected median distances (symmetric, zero diagonal)."""
    n = len(model_median)
    iu = np.triu_indices(n, 1)
    x = feats[iu]
    if params["kind"] == "linear":
        z = (x - np.asarray(params["mu"])) / np.asarray(params["sd"])
        r = np.where(np.isfinite(z), z, 0.0) @ np.asarray(params["beta"]) + params["intercept"]
    else:
        r = _mlp(params, x)
    out = np.zeros((n, n))
    out[iu] = np.asarray(model_median, dtype=np.float64)[iu] * np.exp(r)
    return out + out.T
