"""Synthetic check (planted ideal Gaussian chains): how large are the fitted ensemble's spring couplings
(precision matrix, K = pinv of the centred covariance) beyond the backbone? For a true ideal chain the
precision is exactly tridiagonal (zero beyond |i-j| = 1).

    python validation/network_check.py        # -> validation/network_surgery_check.json (SYNTHETIC input)
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from chronocell import ensemble as E, perturb as PT
from tests.test_v33 import _population

out = []
for n, seed in ((40, 11), (60, 3)):
    d = _population(n=n, cells=6000, seed=seed)
    res = E.fit_ensemble((d < 1.0).mean(0), 6000, r_c_nm=150.0, cfg=E.EnsembleConfig(seed=0))
    S, _ = PT.pair_variance_from_result(res)
    J = np.eye(n) - 1.0 / n
    C = -0.5 * J @ S @ J
    K = np.linalg.pinv(C)
    sep = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    backbone = np.median(np.abs(K[sep == 1]))
    far = np.abs(K[sep > 5])
    row = {"n": n, "seed": seed, "median_abs_backbone": float(backbone), "median_abs_beyond_5": float(np.median(far)),
           "p95_abs_beyond_5": float(np.percentile(far, 95)), "max_abs_beyond_5": float(far.max()),
           "ratio_p95_to_backbone": float(np.percentile(far, 95) / backbone)}
    out.append(row)
    print(row)
json.dump({"note": "synthetic planted ideal chains; justifies dropping spring-network surgery (TUNING.md section 7)",
           "rows": out}, open(ROOT / "network_surgery_check.json", "w"), indent=1)
