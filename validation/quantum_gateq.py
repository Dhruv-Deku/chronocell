"""
Gate Q (quantum lab): on real data, does the quantum route (QAOA, VQE, a quantum-kernel classifier, all on the
statevector SIMULATOR of chronocell/quantum) do what it claims, and how does it compare with classical methods
on the same input?

Q1  QAOA solves the TAD-boundary QUBO: on held-out windows the best of the shots is a minimum-energy bitstring
    (checked against exact enumeration of all 2^n states).
Q2  The QUBO's boundaries agree with the reference domain calls at least as well as the app's classical callers
    (insulation, TopDom-like), within a margin.
Q3  VQE (UCCSD) reaches the exact (FCI) energy of H2 and HeH+ within chemical accuracy (1.6 mHa) at every
    pre-registered bond length.
Q4  The quantum-kernel SVM (ZZ feature map) predicts whether a gene is expressed from Hi-C features about as
    well as a classical RBF-kernel SVM, trained on one cell line and tested on another.

Data. Hi-C: the ENCODE GRCh38 maps Gate 6 read by region (10 Mb windows at 10 kb, cached): practice GM12878
(ENCFF256UOW) chr1/2/3:100-110 Mb; test K562 (ENCFF616PUW) and IMR-90 (ENCFF188SSH) chr4:100-110, chr7:100-110,
chr11:60-70 Mb. Reference domains: ENCODE's Arrowhead "contact domains" called on those same maps (preferred
default; portal MD5): GM12878 ENCFF531LSJ, K562 ENCFF271SAF, IMR-90 ENCFF166QGX. So Q2 measures agreement with
Arrowhead, the field's standard caller, not biological truth. Gene labels (Q4): GTEx v10 median TPM >= 1 in the
matching cell type (GM12878: Cells - EBV-transformed lymphocytes; IMR-90: Cells - Cultured fibroblasts); K562
has no matching GTEx tissue and is not used for Q4. Genes: protein-coding TSS more than 500 kb inside a window.

    python validation/quantum_gateq.py --practice tad      # settings for Q1/Q2 on GM12878
    python validation/quantum_gateq.py --practice qaoa     # QAOA depth / objective on GM12878 (resumable)
    python validation/quantum_gateq.py --practice qsvm     # kernel settings (cross-validation on GM12878)
    python validation/quantum_gateq.py --practice chem     # literature anchors of the chemistry
    python validation/quantum_gateq.py --test all          # run once (rules in frozen.QUANTUM_GATEQ); resumable
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                                                          # noqa: E402
import loops_gate6 as L6                                                      # noqa: E402,F401  (registers the maps)
from chronocell import analysis as AN, domains as DOM                         # noqa: E402
from chronocell.quantum import chem as CH, kernels as KQ, problems as PR, qubo as QB, sim  # noqa: E402

RES10 = 10_000
WINDOW = 10_000_000
CORE = 20                                     # bins per QUBO window -> 19 qubits
DOMAINS = {"gm12878": ("ENCFF531LSJ", "73a79b71bf4a266d29a27f56dc7282ee"),
           "k562": ("ENCFF271SAF", "80e31b2bf8af3c8bff0e4c71a1a9c8dd"),
           "imr90": ("ENCFF166QGX", "9fe299edc11042270272fbc555678c34")}
PRACTICE = {"gm12878": [("chr1", 100_000_000), ("chr2", 100_000_000), ("chr3", 100_000_000)]}
TEST = {"k562": [("chr4", 100_000_000), ("chr7", 100_000_000), ("chr11", 60_000_000)],
        "imr90": [("chr4", 100_000_000), ("chr7", 100_000_000), ("chr11", 60_000_000)]}
GTEX_TISSUE = {"gm12878": "Cells_EBV-transformed_lymphocytes", "imr90": "Cells_Cultured_fibroblasts"}
EDGE = 500_000
TPM_ON = 1.0
CACHE = D.DATA / "cache" / "gateq"

# practice grids
TAD_GRID = {"res": [40_000, 50_000], "min_size": [2, 3, 4], "boundary_cost": [0.0, 0.5, 1.0, 1.5, 2.0, 3.0],
            "gamma": {"difference": [1.0, 1.5, 2.0, 2.5, 3.0, 4.0], "log": [1.0, 1.5, 2.0, 3.0]}}
CLASSICAL_GRID = {"insulation": [(w, d) for w in (2, 3, 4, 5, 8) for d in (0.05, 0.1, 0.15, 0.25, 0.4)],
                  "topdom": [(w,) for w in (2, 3, 4, 5, 8)]}
QAOA_GRID = [{"p": p, "objective": o} for p in (3, 6) for o in ("cvar", "expectation")]
QAOA_FIXED = {"shots": 4096, "maxiter": 80, "alpha": 0.1}
QSVM_GRID = [{"bandwidth": s, "reps": r, "C": c} for s in (0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8) for r in (1, 2)
             for c in (0.1, 1.0, 10.0, 100.0)]
RBF_GRID = [{"gamma": g, "C": c} for g in (0.005, 0.01, 0.02, 0.05, 0.2, 0.5, 2.0) for c in (0.1, 1.0, 10.0, 100.0)]
FEATURES = ["compartment", "insulation", "coverage", "log_dist_boundary", "gene_density"]


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)),
                   encoding="utf-8")
    tmp.replace(path)


# ======================================================================================
# Data
# ======================================================================================
def region_counts(cell: str, chrom: str, start: int) -> np.ndarray:
    m = D.hic_region(f"encode_loops_{cell}", chrom, start, start + WINDOW, RES10)
    m = np.asarray(m, float)
    return np.maximum(m, m.T)


def reference_boundaries(cell: str, chrom: str, start: int) -> list[int]:
    acc, md5 = DOMAINS[cell]
    p = D.fetch_url(f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bedpe.gz",
                    D.DATA / "encode" / f"{acc}.bedpe.gz", md5)
    out = set()
    with gzip.open(p, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or not f[1].isdigit():
                continue
            c = f[0] if f[0].startswith("chr") else "chr" + f[0]
            if c != chrom:
                continue
            for x in (int(f[1]), int(f[5])):
                if start <= x <= start + WINDOW:
                    out.add(x)
    return sorted(out)


def cores(N: int, span: int) -> list[tuple[int, int]]:
    out = []
    c0 = span
    while c0 + CORE <= N - span:
        out.append((c0, c0 + CORE))
        c0 += CORE
    return out


def interval(start: int, res: int, core: tuple[int, int]) -> tuple[float, float]:
    return start + (core[0] + 0.5) * res, start + (core[1] - 0.5) * res


def in_iv(bps, iv) -> list[float]:
    return [b for b in bps if iv[0] <= b < iv[1]]


def classical_boundaries(C: np.ndarray, start: int, res: int, method: str, param: tuple) -> list[float]:
    N = len(C)
    ci, cj = np.triu_indices(N, 1)
    v = C[ci, cj]
    ok = v > 0
    if method == "insulation":
        w, depth = param
        ins = DOM.insulation(ci[ok], cj[ok], v[ok], N, int(w))
        return [start + (i + 0.5) * res for i in DOM.boundaries(ins, int(w), float(depth))]
    b, _ = AN.topdom(ci[ok], cj[ok], v[ok], N, window=int(param[0]))
    return [start + (i + 1) * res for i in b]


# ======================================================================================
# Q1 / Q2: TAD boundaries
# ======================================================================================
def tad_qubo_for(C: np.ndarray, core, setting: dict, mask=None, exp=None):
    span = 2 * setting["min_size"]
    mask = PR.coverage_mask(C) if mask is None else mask
    exp = PR.expected_by_distance(C, span, mask) if exp is None else exp
    return PR.tad_qubo(C, core, setting["gamma"], setting["min_size"], span, exp, mask,
                       boundary_cost=setting.get("boundary_cost", 0.0), weight=setting.get("weight", "difference"))


def tad_exact(C: np.ndarray, start: int, setting: dict) -> list[dict]:
    """Exact optimum of every core window's QUBO (dynamic programming over the band; equal to enumeration)."""
    res, m = setting["res"], setting["min_size"]
    span = 2 * m
    mask = PR.coverage_mask(C)
    exp = PR.expected_by_distance(C, span, mask)
    out = []
    for core in cores(len(C), span):
        q, pos = tad_qubo_for(C, core, setting, mask, exp)
        ex = QB.banded_exact(q)
        out.append({"core": core, "boundaries": [start + k * res for k in PR.boundaries_from_bits(ex.bits, pos)],
                    "energy": ex.energy, "feasible": PR.feasible(ex.bits, m), "qubits": q.n})
    return out


def score(cell: str, regions: list, res: int, span: int, pred_of) -> dict:
    """Pool TP / FP / FN over every core window of the regions; pred_of(chrom, start, core_index, core) -> bp list."""
    tp = fp = fn = 0
    for chrom, start in regions:
        ref = reference_boundaries(cell, chrom, start)
        N = WINDOW // res
        for k, core in enumerate(cores(N, span)):
            iv = interval(start, res, core)
            a, b, c = PR.match(in_iv(pred_of(chrom, start, k, core), iv), in_iv(ref, iv), res)
            tp, fp, fn = tp + a, fp + b, fn + c
    return PR.f1(tp, fp, fn)


def practice_tad() -> dict:
    cell = "gm12878"
    regs = PRACTICE[cell]
    grid = []
    for res in TAD_GRID["res"]:
        regC = {(c, s): PR.coarsen(region_counts(cell, c, s), res // RES10) for c, s in regs}
        for weight, gammas in TAD_GRID["gamma"].items():
            for m in TAD_GRID["min_size"]:
                for g in gammas:
                    for mu in TAD_GRID["boundary_cost"]:
                        setting = {"res": res, "min_size": m, "gamma": g, "boundary_cost": mu, "weight": weight}
                        sols = {k: tad_exact(C, k[1], setting) for k, C in regC.items()}
                        sc = score(cell, regs, res, 2 * m, lambda c, s, k, core: sols[(c, s)][k]["boundaries"])
                        feas = float(np.mean([w["feasible"] for v in sols.values() for w in v]))
                        grid.append({**setting, **sc, "optimum_feasible": feas,
                                     "windows": sum(len(v) for v in sols.values())})
                        print(f"{setting}: F1 {sc['f1']:.3f} (P {sc['precision']:.2f} R {sc['recall']:.2f}) "
                              f"feasible {feas:.2f}", flush=True)
    best = max(grid, key=lambda r: (r["f1"], -r["boundary_cost"], -r["gamma"]))
    classical = {}
    for method, params in CLASSICAL_GRID.items():
        rows = []
        for prm in params:
            regC = {(c, s): PR.coarsen(region_counts(cell, c, s), best["res"] // RES10) for c, s in regs}
            pred = {k: classical_boundaries(v, k[1], best["res"], method, prm) for k, v in regC.items()}
            sc = score(cell, regs, best["res"], 2 * best["min_size"], lambda c, s, k, core: pred[(c, s)])
            rows.append({"param": list(prm), **sc})
            print(f"{method} {prm}: F1 {sc['f1']:.3f}", flush=True)
        classical[method] = {"grid": rows, "best": max(rows, key=lambda r: r["f1"])}
    out = {"made": _now(), "cell": cell, "regions": regs, "grid": grid,
           "choice": {k: best[k] for k in ("res", "min_size", "gamma", "boundary_cost", "weight")}, "choice_scores": best,
           "classical": classical}
    _json(ROOT / "results_gateq_practice_tad.json", out)
    return out


def qaoa_window(tag: str, cell: str, chrom: str, start: int, k: int, setting: dict, cfg: dict, noise: dict) -> dict:
    """QAOA (+ simulated annealing, SQA, noise) on one core window; cached per window and setting."""
    key = (f"{tag}_{cell}_{chrom}_{start}_{k}_p{cfg['p']}_{cfg['objective']}_r{setting['res']}_m{setting['min_size']}"
           f"_g{setting['gamma']}_c{setting.get('boundary_cost', 0.0)}_{setting.get('weight', 'difference')}")
    path = CACHE / f"{key}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    res, m = setting["res"], setting["min_size"]
    C = PR.coarsen(region_counts(cell, chrom, start), res // RES10)
    span = 2 * m
    core = cores(len(C), span)[k]
    q, pos = tad_qubo_for(C, core, setting)
    e = q.energies()
    ex = QB.exact(q, e)
    dp = QB.banded_exact(q)
    r = sim.qaoa(e, p=cfg["p"], shots=QAOA_FIXED["shots"], objective=cfg["objective"], alpha=QAOA_FIXED["alpha"],
                 maxiter=QAOA_FIXED["maxiter"], seed=k)
    h, J, _ = q.ising()
    counts = sim.qaoa_circuit(h, J, r.gammas, r.betas).counts()
    F = sim.fidelity_estimate(counts, noise)
    rng = np.random.default_rng(1000 + k)
    noisy = sim.sample(r.probs, QAOA_FIXED["shots"], rng, F, noise["readout"], q.n)
    emin = ex.energy
    tol = 1e-9 * max(1.0, abs(emin))
    sa = QB.anneal(q, reads=64, sweeps=600, seed=k)
    sq = QB.sqa(q, reads=16, sweeps=300, trotter=16, seed=k)
    rand = np.random.default_rng(2000 + k).integers(0, 2 ** q.n, QAOA_FIXED["shots"])
    best_bits = QB.to_bits(r.best, q.n)
    out = {"cell": cell, "chrom": chrom, "start": start, "core": core, "qubits": q.n, "setting": setting, "qaoa": cfg,
           "exact_energy": emin, "n_optima": ex.detail["n_optima"],
           "exact_boundaries": [start + j * res for j in PR.boundaries_from_bits(ex.bits, pos)],
           "qaoa_hit": bool(r.hit), "qaoa_p_optimal": r.p_optimal, "qaoa_best_energy": r.best_energy,
           "qaoa_boundaries": [start + j * res for j in PR.boundaries_from_bits(best_bits, pos)],
           "qaoa_approx_ratio": r.approximation_ratio, "qaoa_evaluations": r.evaluations, "qaoa_seconds": r.seconds,
           "device": r.device, "uniform_p_optimal": ex.detail["n_optima"] / 2 ** q.n,
           "random_hit": bool(np.any(np.abs(e[rand] - emin) <= tol)),
           "noisy_fidelity": F, "noisy_hit": bool(np.any(np.abs(e[noisy] - emin) <= tol)),
           "circuit": {k2: counts[k2] for k2 in ("gates", "depth", "cx_equivalent", "one_qubit")},
           "sa_hit": bool(abs(sa.energy - emin) <= tol), "sa_seconds": sa.seconds,
           "sqa_hit": bool(abs(sq.energy - emin) <= tol), "sqa_seconds": sq.seconds, "exact_seconds": ex.seconds,
           "dp_seconds": dp.seconds, "dp_agrees": bool(abs(dp.energy - emin) <= tol),
           "gammas": r.gammas, "betas": r.betas}
    _json(path, out)
    return out


def practice_qaoa(limit_s: float = 6000) -> dict:
    t0 = time.perf_counter()
    tad = json.loads((ROOT / "results_gateq_practice_tad.json").read_text(encoding="utf-8"))
    setting = tad["choice"]
    cell = "gm12878"
    rows = []
    for cfg in QAOA_GRID:
        res = []
        for chrom, start in PRACTICE[cell]:
            n_cores = len(cores(WINDOW // setting["res"], 2 * setting["min_size"]))
            for k in range(n_cores):
                if time.perf_counter() - t0 > limit_s:
                    print("time budget reached; rerun to continue (cached per window)", flush=True)
                    return {}
                w = qaoa_window("practice", cell, chrom, start, k, setting, cfg, sim.NOISE_DEFAULT)
                res.append(w)
                print(f"{cfg} {chrom} core {k}: hit {w['qaoa_hit']} p_opt {w['qaoa_p_optimal']:.4f} "
                      f"({w['qaoa_seconds']:.1f} s, {w['device']})", flush=True)
        rows.append(summarise_qaoa(cfg, res))
    best = max(rows, key=lambda r: (r["hit_rate"], -r["p"], r["mean_p_optimal"]))
    out = {"made": _now(), "setting": setting, "grid": rows, "choice": {"p": best["p"], "objective": best["objective"]},
           "fixed": QAOA_FIXED, "noise": sim.NOISE_DEFAULT}
    _json(ROOT / "results_gateq_practice_qaoa.json", out)
    return out


def summarise_qaoa(cfg: dict, res: list[dict]) -> dict:
    return {**cfg, "windows": len(res), "hit_rate": float(np.mean([w["qaoa_hit"] for w in res])),
            "mean_p_optimal": float(np.mean([w["qaoa_p_optimal"] for w in res])),
            "median_p_optimal": float(np.median([w["qaoa_p_optimal"] for w in res])),
            "mean_uniform_p_optimal": float(np.mean([w["uniform_p_optimal"] for w in res])),
            "random_hit_rate": float(np.mean([w["random_hit"] for w in res])),
            "noisy_hit_rate": float(np.mean([w["noisy_hit"] for w in res])),
            "mean_noisy_fidelity": float(np.mean([w["noisy_fidelity"] for w in res])),
            "sa_hit_rate": float(np.mean([w["sa_hit"] for w in res])),
            "sqa_hit_rate": float(np.mean([w["sqa_hit"] for w in res])),
            "mean_qaoa_seconds": float(np.mean([w["qaoa_seconds"] for w in res])),
            "mean_sa_seconds": float(np.mean([w["sa_seconds"] for w in res])),
            "mean_sqa_seconds": float(np.mean([w["sqa_seconds"] for w in res])),
            "mean_exact_seconds": float(np.mean([w["exact_seconds"] for w in res])),
            "mean_dp_seconds": float(np.mean([w["dp_seconds"] for w in res])),
            "mean_cx": float(np.mean([w["circuit"]["cx_equivalent"] for w in res])),
            "qubits": int(res[0]["qubits"]) if res else 0}


# ======================================================================================
# Q4: genes
# ======================================================================================
def gene_features(cell: str, chrom: str, start: int) -> list[dict]:
    from chronocell import genes as G, annotations as AO
    M = region_counts(cell, chrom, start)
    n = len(M)
    ci, cj = np.triu_indices(n, 1)
    v = M[ci, cj]
    ok = v > 0
    ci, cj, v = ci[ok], cj[ok], v[ok]
    gt = G.table("hg38")
    gt = gt[(gt["chrom"] == chrom) & (gt["biotype"] == "coding")]
    tss_all = gt["tss"].to_numpy()
    dens = np.zeros(n)
    for t in tss_all:
        b = (t - start) // RES10
        if 0 <= b < n:
            dens[b] += 1
    comp, _ = DOM.compartments(ci, cj, v, n, dens, max_bins=100)
    w = 20
    ins = DOM.insulation(ci, cj, v, n, w)
    bnd = np.array(DOM.boundaries(ins, w))
    cov = M.sum(axis=1)
    lc = np.where(cov > 0, np.log1p(cov), np.nan)
    lc = (lc - np.nanmean(lc)) / (np.nanstd(lc) or 1.0)
    tab = AO.gtex_median_tpm(False)
    tissue = GTEX_TISSUE.get(cell)
    tpm = tab[tissue].groupby(level=0).max() if tissue else None
    rows = []
    for r in gt.itertuples(index=False):
        if not (start + EDGE <= r.tss < start + WINDOW - EDGE):
            continue
        b = (r.tss - start) // RES10
        if cov[b] <= 0:
            continue
        near = np.sum(np.abs(tss_all - r.tss) <= 250_000) - 1
        d = float(np.min(np.abs(bnd - b))) * RES10 if len(bnd) else float(WINDOW)
        f = {"compartment": comp[b], "insulation": ins[b], "coverage": lc[b], "log_dist_boundary": math.log10(d + RES10),
             "gene_density": math.log1p(near)}
        if not all(np.isfinite(list(f.values()))):
            continue
        t = float(tpm.get(r.name, np.nan)) if tpm is not None else np.nan
        if not np.isfinite(t):
            continue
        rows.append({"gene": r.name, "chrom": chrom, "tss": int(r.tss), **f, "tpm": t, "label": int(t >= TPM_ON)})
    return rows


def genes_for(cell: str, regions: list) -> list[dict]:
    path = CACHE / f"genes_{cell}_{'_'.join(c + str(s) for c, s in regions)}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    rows = [g for c, s in regions for g in gene_features(cell, c, s)]
    _json(path, rows)
    return rows


def _xy(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    return np.array([[r[f] for f in FEATURES] for r in rows], float), np.array([r["label"] for r in rows], int)


def _fit_predict(kind: str, prm: dict, Xtr, ytr, Xte) -> np.ndarray:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    if kind == "qsvm":
        a, b = KQ.scale_features(Xtr, Xte)
        return KQ.QSVM(prm["C"], prm["bandwidth"], prm["reps"]).fit(a, ytr).decision_function(b)
    sc = StandardScaler().fit(Xtr)
    if kind == "rbf":
        g = prm["gamma"]
        m = SVC(C=prm["C"], kernel="rbf", gamma=g, class_weight="balanced").fit(sc.transform(Xtr), ytr)
        return m.decision_function(sc.transform(Xte))
    m = LogisticRegression(class_weight="balanced", max_iter=2000).fit(sc.transform(Xtr), ytr)
    return m.decision_function(sc.transform(Xte))


def _cv_auc(kind: str, prm: dict, X, y, folds: int = 5, seed: int = 0) -> float:
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    aucs = []
    for tr, te in StratifiedKFold(folds, shuffle=True, random_state=seed).split(X, y):
        aucs.append(roc_auc_score(y[te], _fit_predict(kind, prm, X[tr], y[tr], X[te])))
    return float(np.mean(aucs))


def practice_qsvm() -> dict:
    rows = genes_for("gm12878", PRACTICE["gm12878"])
    X, y = _xy(rows)
    print(f"{len(y)} genes, {y.mean():.2f} expressed", flush=True)
    qs = [{**p, "cv_auc": _cv_auc("qsvm", p, X, y)} for p in QSVM_GRID]
    rb = [{**p, "cv_auc": _cv_auc("rbf", p, X, y)} for p in RBF_GRID]
    lr = _cv_auc("logistic", {}, X, y)
    out = {"made": _now(), "cell": "gm12878", "genes": len(y), "expressed_fraction": float(y.mean()), "features": FEATURES,
           "qsvm_grid": qs, "rbf_grid": rb, "logistic_cv_auc": lr,
           "choice": {"qsvm": {k: v for k, v in max(qs, key=lambda r: r["cv_auc"]).items() if k != "cv_auc"},
                      "rbf": {k: v for k, v in max(rb, key=lambda r: r["cv_auc"]).items() if k != "cv_auc"}},
           "best_cv_auc": {"qsvm": max(r["cv_auc"] for r in qs), "rbf": max(r["cv_auc"] for r in rb), "logistic": lr}}
    print(json.dumps(out["best_cv_auc"]), out["choice"], flush=True)
    _json(ROOT / "results_gateq_practice_qsvm.json", out)
    return out


def _boot_auc(y, scores: dict, reps: int = 1000, seed: int = 0) -> dict:
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    n = len(y)
    out = {k: [] for k in scores}
    diff = []
    for _ in range(reps):
        i = rng.integers(0, n, n)
        if len(set(y[i])) < 2:
            continue
        a = {k: roc_auc_score(y[i], s[i]) for k, s in scores.items()}
        for k in a:
            out[k].append(a[k])
        diff.append(a["qsvm"] - a["rbf"])
    ci = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] for k, v in out.items()}
    ci["qsvm_minus_rbf"] = [float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))]
    return ci


# ======================================================================================
# Q3: chemistry
# ======================================================================================
def chem_rows(geometries: dict) -> list[dict]:
    rows = []
    for name, Rs in geometries.items():
        rows += CH.curve(name, Rs, "uccsd", noise=sim.NOISE_DEFAULT)
        for r in CH.curve(name, Rs, "hea", noise=sim.NOISE_DEFAULT, layers=2):
            rows.append({**r, "ansatz": "hardware-efficient (2 layers)"})
    for r in rows:
        r.setdefault("ansatz", "UCCSD")
    return rows


def practice_chem() -> dict:
    B = CH.BOHR_ANGSTROM
    h2 = CH.molecule("H2", 1.4 * B)
    he = CH.molecule("HeH+", 1.4632 * B)
    out = {"made": _now(), "anchors": {
        "H2 R=1.4 bohr E_HF (Szabo & Ostlund: -1.1167)": h2.e_hf, "H2 R=1.4 bohr E_FCI (Szabo & Ostlund: -1.1373)": h2.e_fci,
        "HeH+ R=1.4632 bohr E_HF (Szabo & Ostlund: -2.86066)": he.e_hf},
        "pauli_terms_h2": h2.paulis, "rows": chem_rows({"H2": [0.735], "HeH+": [0.775]})}
    _json(ROOT / "results_gateq_practice_chem.json", out)
    print(json.dumps(out["anchors"], indent=1), flush=True)
    return out


# ======================================================================================
# Test
# ======================================================================================
def run_test() -> dict:
    import frozen as F
    from sklearn.metrics import roc_auc_score
    G = F.QUANTUM_GATEQ
    setting, cfg = G["tad"], G["qaoa"]
    noise = G["noise"]
    out = {"made": _now(), "rule": G, "q1": {}, "q2": {}, "per_window": {}}
    # ---- Q1 / Q2
    span = 2 * setting["min_size"]
    res = setting["res"]
    q2_rows = {}
    for cell, regs in TEST.items():
        wins = []
        for chrom, start in regs:
            for k in range(len(cores(WINDOW // res, span))):
                w = qaoa_window("test", cell, chrom, start, k, setting, cfg, noise)
                wins.append(w)
                print(f"{cell} {chrom} core {k}: hit {w['qaoa_hit']} p_opt {w['qaoa_p_optimal']:.4f} "
                      f"({w['qaoa_seconds']:.1f} s)", flush=True)
        out["per_window"][cell] = [{k2: w[k2] for k2 in ("chrom", "start", "core", "qaoa_hit", "qaoa_p_optimal",
                                                         "uniform_p_optimal", "noisy_hit", "sa_hit", "sqa_hit",
                                                         "random_hit", "qaoa_seconds", "noisy_fidelity")} for w in wins]
        out["q1"][cell] = summarise_qaoa(cfg, wins)
        by = {(w["chrom"], w["start"], w["core"][0]): w for w in wins}
        cores_list = cores(WINDOW // res, span)

        def from_win(field):
            return lambda c, s, k, core: by[(c, s, cores_list[k][0])][field]
        rows = {"QAOA (simulated quantum)": score(cell, regs, res, span, from_win("qaoa_boundaries")),
                "exact optimum of the QUBO": score(cell, regs, res, span, from_win("exact_boundaries"))}
        for method in ("insulation", "topdom"):
            prm = tuple(G["classical"][method])
            pred = {(c, s): classical_boundaries(PR.coarsen(region_counts(cell, c, s), res // RES10), s, res, method, prm)
                    for c, s in regs}
            rows[f"{method} (classical)"] = score(cell, regs, res, span, lambda c, s, k, core: pred[(c, s)])
        q2_rows[cell] = rows
    hits = [w["qaoa_hit"] for cell in out["per_window"] for w in out["per_window"][cell]]
    out["q1"]["pooled_hit_rate"] = float(np.mean(hits))
    out["q1"]["windows"] = len(hits)
    out["q1"]["pass"] = bool(np.mean(hits) >= G["q1_min_hit_rate"])
    q2_pass = {}
    for cell, rows in q2_rows.items():
        best_cl = max(rows["insulation (classical)"]["f1"], rows["topdom (classical)"]["f1"])
        q2_pass[cell] = bool(rows["QAOA (simulated quantum)"]["f1"] >= best_cl - G["q2_margin"])
    out["q2"] = {"rows": q2_rows, "pass_by_cell": q2_pass, "pass": all(q2_pass.values())}
    _json(ROOT / "results_gateq.json", out)
    # ---- Q3
    rows = chem_rows(G["geometries"])
    ucc = [r for r in rows if r["ansatz"] == "UCCSD"]
    out["q3"] = {"rows": rows, "max_abs_error_mha": float(max(abs(r["error_mha"]) for r in ucc)),
                 "pass": all(r["within_chemical_accuracy"] for r in ucc)}
    _json(ROOT / "results_gateq.json", out)
    # ---- Q4
    tr = genes_for("gm12878", PRACTICE["gm12878"])
    te = genes_for("imr90", TEST["imr90"])
    Xtr, ytr = _xy(tr)
    Xte, yte = _xy(te)
    sc = {"qsvm": _fit_predict("qsvm", G["qsvm"], Xtr, ytr, Xte), "rbf": _fit_predict("rbf", G["rbf"], Xtr, ytr, Xte),
          "logistic": _fit_predict("logistic", {}, Xtr, ytr, Xte)}
    auc = {k: float(roc_auc_score(yte, v)) for k, v in sc.items()}
    ci = _boot_auc(yte, sc)
    # shot noise: the same QSVM with kernel entries estimated from shots
    a, b = KQ.scale_features(Xtr, Xte)
    qs = KQ.QSVM(G["qsvm"]["C"], G["qsvm"]["bandwidth"], G["qsvm"]["reps"], shots=G["qsvm_shots"]).fit(a, ytr)
    auc_shots = float(roc_auc_score(yte, qs.decision_function(b)))
    q4_pass = bool(auc["qsvm"] >= auc["rbf"] - G["q4_margin"] and ci["qsvm"][0] > 0.5)
    out["q4"] = {"train_genes": len(ytr), "test_genes": len(yte), "test_expressed_fraction": float(yte.mean()),
                 "auc": auc, "ci95": ci, "qsvm_auc_with_shots": auc_shots, "pass": q4_pass}
    out["overall"] = {k: out[k]["pass"] for k in ("q1", "q2", "q3", "q4")}
    _json(ROOT / "results_gateq.json", out)
    print(json.dumps(out["overall"]), flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", choices=["tad", "qaoa", "qsvm", "chem"])
    g.add_argument("--test", choices=["all"])
    ap.add_argument("--budget", type=float, default=6000, help="seconds per call before stopping (resumable)")
    a = ap.parse_args()
    if a.practice == "tad":
        practice_tad()
    elif a.practice == "qaoa":
        practice_qaoa(a.budget)
    elif a.practice == "qsvm":
        practice_qsvm()
    elif a.practice == "chem":
        practice_chem()
    else:
        run_test()


if __name__ == "__main__":
    main()
