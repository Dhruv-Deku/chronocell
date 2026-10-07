"""
Round 2 of the quantum gates (October 2026): new methods for the gates that failed (Q2, Q4, Q6), each developed on
data already seen (the failed tests' data become practice), pre-registered in frozen.QUANTUM_ROUND2, and tested ONCE
on data never used before. The original results (results_gateq.json, results_qdrug.json) stay as they are.

Q2b TAD boundaries: the domain QUBO's settings were chosen in Gate Q on three GM12878 windows holding 41 reference
    boundaries, and called too few boundaries on K562 and IMR-90 (recall 0.14-0.17). Q2b: the same QUBO and QAOA,
    settings re-chosen on all nine windows already seen (GM12878, K562, IMR-90; about 400 reference boundaries), the
    classical callers re-tuned on the same windows (same chance); test on two cell lines never used: HMEC and HAP-1
    (ENCODE in situ Hi-C and the Arrowhead domains called on each map).
Q4b Genes: the quantum-kernel SVM predicting whether a gene is active from Hi-C features. Q4 trained on 106 genes of
    one cell line with GTEx labels and lost to the RBF-SVM by 0.034 AUC (margin 0.03). Q4b: labels from ENCODE total
    RNA-seq of the cell line itself (one lab, one protocol: CSHL long total RNA, GENCODE V29 quantifications, TPM >= 1),
    two more Hi-C features, training on every gene of the windows already seen (GM12878, IMR-90, K562); settings by
    leave-one-cell-line-out cross-validation; test on a cell line never used: HMEC.
Q6b Molecules: ADAPT-VQE (chronocell/quantum/molecules.adapt_vqe; Grimsley et al., Nat Commun 2019) instead of the
    fixed-order UCCSD circuit that missed chemical accuracy on stretched N2 and HCN in Q6. Practice: the six Q6 cases
    and nine more (equilibrium and N2 at 2x). Test: eleven molecule / bond-length / active-space cases never run before.

    python validation/quantum_round2.py --practice q2b
    python validation/quantum_round2.py --practice q4b
    python validation/quantum_round2.py --practice q6b
    python validation/quantum_round2.py --test q6b         # run once (rules in frozen.QUANTUM_ROUND2)
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                                                          # noqa: E402
import quantum_gateq as GQ                                                    # noqa: E402  (registers the maps)
from chronocell.quantum import kernels as KQ, molecules as M, problems as PR, sim  # noqa: E402

RESULTS = ROOT / "results_round2.json"
Q6B_PRACTICE = [["H2O", 1.5, [4, 4]], ["NH3", 1.5, [6, 5]], ["N2", 1.5, [6, 6]], ["HF", 2.0, [2, 2]],
                ["CH2O", 1.3, [4, 4]], ["HCN", 1.3, [4, 4]], ["N2", 1.0, [6, 6]], ["CO", 1.0, [6, 6]],
                ["H2O", 1.0, [8, 6]], ["NH3", 1.0, [6, 6]], ["HCN", 1.0, [6, 6]], ["CH2O", 1.0, [6, 6]],
                ["N2", 2.0, [6, 6]], ["LiH", 2.0, [2, 4]], ["CH4", 1.0, [4, 4]]]
Q6B_GRID = [{"pool": p, "grad_tol": 1e-3, "max_operators": 100} for p in ("sd", "gsd")]


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _json(path: Path, obj) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)),
                   encoding="utf-8")
    tmp.replace(path)


# ======================================================================================
# New cell lines, shared by Q2b and Q4b (read only by the tests)
# ======================================================================================
ENC_HIC = "https://www.encodeproject.org/files/{0}/@@download/{0}.hic"
NEW_MAPS = {   # in situ Hi-C map (ENCODE4 pipeline) and the Arrowhead contact domains called on that map (portal MD5)
    "hmec": ("ENCFF943JRY", "ENCFF578NLR", "6fdf3e26d2b051f8be4b4b9f190838b6", "ENCSR711AVS"),
    "hap1": ("ENCFF898HRO", "ENCFF459IDB", "9dccc20b714c5d9c3412475d6fd06a4b", "ENCSR390HMC")}
for _k, (_h, _d, _m, _) in NEW_MAPS.items():
    D.HIC_SOURCES.setdefault(f"encode_loops_{_k}", (ENC_HIC.format(_h), "GRCh38", _k))
    GQ.DOMAINS.setdefault(_k, (_d, _m))
SEEN = {"gm12878": GQ.PRACTICE["gm12878"], "k562": GQ.TEST["k562"], "imr90": GQ.TEST["imr90"]}
NEW_REGIONS = GQ.PRACTICE["gm12878"] + GQ.TEST["k562"]      # six 10 Mb regions, read in the new cell lines only


# ======================================================================================
# Q2b: TAD boundaries (settings re-chosen on every seen window; new cell lines)
# ======================================================================================
Q2B_GRID = {"res": [40_000, 50_000, 60_000, 80_000], "min_size": [2, 3, 4], "gamma": [1.5, 2.0, 3.0, 4.0, 5.0],
            "boundary_cost": [0.0, 0.25, 0.5, 1.0], "weight": ["difference", "log"]}


def _classical_f1(cell: str, regs: list, res: int, span: int, method: str, prm: tuple) -> dict:
    pred = {(c, s): GQ.classical_boundaries(PR.coarsen(GQ.region_counts(cell, c, s), res // GQ.RES10), s, res, method, prm)
            for c, s in regs}
    return GQ.score(cell, regs, res, span, lambda c, s, k, core: pred[(c, s)])


def practice_q2b() -> dict:
    """Every setting's exact QUBO optimum on the nine seen windows; the classical callers' grids on the same windows.
    Choice: the setting whose worst cell line has the largest F1 margin over that cell's better classical caller
    (each classical caller tuned jointly over the three cell lines at the same resolution)."""
    rows, classical = [], {}
    for res in Q2B_GRID["res"]:
        for m in Q2B_GRID["min_size"]:
            span = 2 * m
            cl = {}
            for method, params in GQ.CLASSICAL_GRID.items():
                grid = [{"param": list(p), **{c: _classical_f1(c, regs, res, span, method, p)["f1"] for c, regs in SEEN.items()}}
                        for p in params]
                best = max(grid, key=lambda r: np.mean([r[c] for c in SEEN]))
                cl[method] = {"best": best, "grid": grid}
            classical[f"{res}_{m}"] = cl
            ref = {c: max(cl["insulation"]["best"][c], cl["topdom"]["best"][c]) for c in SEEN}
            for weight in Q2B_GRID["weight"]:
                for g in Q2B_GRID["gamma"]:
                    for mu in Q2B_GRID["boundary_cost"]:
                        st = {"res": res, "min_size": m, "gamma": g, "boundary_cost": mu, "weight": weight}
                        f = {}
                        for c, regs in SEEN.items():
                            sols = {(ch, s): GQ.tad_exact(PR.coarsen(GQ.region_counts(c, ch, s), res // GQ.RES10), s, st)
                                    for ch, s in regs}
                            f[c] = GQ.score(c, regs, res, span, lambda ch, s, k, core: sols[(ch, s)][k]["boundaries"])
                        margin = {c: f[c]["f1"] - ref[c] for c in SEEN}
                        rows.append({**st, "f1": {c: f[c]["f1"] for c in SEEN}, "scores": f, "classical_best": ref,
                                     "margin": margin, "min_margin": min(margin.values())})
                        print(st, {c: round(f[c]["f1"], 3) for c in SEEN}, "min margin", round(rows[-1]["min_margin"], 3),
                              flush=True)
    best = max(rows, key=lambda r: (r["min_margin"], -r["boundary_cost"], -r["gamma"]))
    ch = {k: best[k] for k in ("res", "min_size", "gamma", "boundary_cost", "weight")}
    cl = classical[f"{best['res']}_{best['min_size']}"]
    out = {"made": _now(), "windows": SEEN, "grid_spec": Q2B_GRID, "grid": rows, "choice": ch, "choice_scores": best,
           "classical_at_choice": {k: v["best"] for k, v in cl.items()},
           "classical_choice": {k: v["best"]["param"] for k, v in cl.items()}, "classical": classical}
    _json(ROOT / "results_round2_practice_q2b.json", out)
    print("choice", ch, "classical", out["classical_choice"], flush=True)
    return out


def test_q2b(R: dict) -> dict:
    setting, cfg, noise = R["tad"], R["qaoa"], R["noise"]
    res, span = setting["res"], 2 * setting["min_size"]
    regs = [tuple(x) for x in R["test_regions"]]
    out = {"per_cell": {}, "q1_like": {}}
    for cell in R["test_cells"]:
        wins = []
        cl_ = GQ.cores(GQ.WINDOW // res, span)
        for chrom, start in regs:
            for k in range(len(cl_)):
                w = GQ.qaoa_window("round2", cell, chrom, start, k, setting, cfg, noise)
                wins.append(w)
                print(f"{cell} {chrom} core {k}: hit {w['qaoa_hit']} ({w['qaoa_seconds']:.1f} s)", flush=True)
        by = {(w["chrom"], w["start"], w["core"][0]): w for w in wins}

        def from_win(field):
            return lambda c, s, k, core: by[(c, s, cl_[k][0])][field]
        rows = {"QAOA (simulated quantum)": GQ.score(cell, regs, res, span, from_win("qaoa_boundaries")),
                "exact optimum of the QUBO": GQ.score(cell, regs, res, span, from_win("exact_boundaries"))}
        for method in ("insulation", "topdom"):
            rows[f"{method} (classical)"] = _classical_f1(cell, regs, res, span, method, tuple(R["classical"][method]))
        best_cl = max(rows["insulation (classical)"]["f1"], rows["topdom (classical)"]["f1"])
        out["per_cell"][cell] = {"rows": rows, "pass": bool(rows["QAOA (simulated quantum)"]["f1"] >= best_cl - R["margin"])}
        out["q1_like"][cell] = GQ.summarise_qaoa(cfg, wins)
    out["pass"] = all(v["pass"] for v in out["per_cell"].values())
    return out


# ======================================================================================
# Q4b: genes (ENCODE RNA-seq labels, more training cell lines, a new test cell line)
# ======================================================================================
RNA = {   # CSHL whole-cell long total RNA-seq, GENCODE V29 gene quantifications (accession, portal MD5)
    "gm12878": [("ENCFF910XWA", "9cc4061cdfdebc9071d57cc39c7d9459"), ("ENCFF413MYB", "fce4d85d80294e7af0e973db832d59e5")],
    "imr90": [("ENCFF019KLP", "6ecee6b8b2131e8a39986e812f9bce98"), ("ENCFF268IJX", "92bcdf261a365d2fafb8630b847f95d9")],
    "k562": [("ENCFF136UBW", "9d001bd2413f175a689b06d2fca726e7"), ("ENCFF171FQU", "56e682ec8fefa1f87a1a0848353b86d7")],
    "hmec": [("ENCFF798WGM", "48aa5afd31a58cd92d62626ee913e058")]}
FEATURES_Q4 = list(GQ.FEATURES)
FEATURES_Q4B = FEATURES_Q4 + ["local_oe", "insulation_slope"]
# grids widened once on practice (both kernels alike) when the first run's best settings sat at their edges
Q4B_QSVM_GRID = [{"bandwidth": s, "reps": r, "C": c} for s in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8)
                 for r in (1, 2) for c in (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)]
Q4B_RBF_GRID = [{"gamma": g, "C": c} for g in (0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 2.0)
                for c in (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)]


def encode_tpm(cell: str) -> dict:
    """Gene symbol -> mean TPM over the replicates (Ensembl ids mapped to symbols with the GTEx v10 table)."""
    import gzip
    from chronocell import annotations as AO
    sym = {}
    with gzip.open(AO._gtex_path(False), "rt", encoding="utf-8") as fh:
        fh.readline()
        fh.readline()
        head = fh.readline().rstrip("\n").split("\t")
        i_n, i_d = head.index("Name"), head.index("Description")
        for line in fh:
            f = line.split("\t", max(i_n, i_d) + 1)
            sym[f[i_n].split(".")[0]] = f[i_d]
    reps = []
    for acc, md5 in RNA[cell]:
        p = D.fetch_url(f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.tsv", D.DATA / "encode" / f"{acc}.tsv", md5)
        t: dict = {}
        with p.open(encoding="utf-8") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            i_g, i_t = head.index("gene_id"), head.index("TPM")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                s = sym.get(f[i_g].split(".")[0])
                if s:
                    t[s] = t.get(s, 0.0) + float(f[i_t])
        reps.append(t)
    keys = set.intersection(*(set(r) for r in reps))
    return {k: float(np.mean([r[k] for r in reps])) for k in keys}


def gene_rows(cell: str, chrom: str, start: int, tpm: dict) -> list[dict]:
    """Q4's five features plus local_oe (mean log2 observed / expected contacts of the TSS bin within 100 kb) and
    insulation_slope (insulation 50 kb downstream minus 50 kb upstream of the TSS, strand-aware)."""
    import math
    from chronocell import domains as DOM, genes as G
    Mx = GQ.region_counts(cell, chrom, start)
    n = len(Mx)
    ci, cj = np.triu_indices(n, 1)
    v = Mx[ci, cj]
    ok = v > 0
    ci, cj, v = ci[ok], cj[ok], v[ok]
    gt = G.table("hg38")
    gt = gt[(gt["chrom"] == chrom) & (gt["biotype"] == "coding")]
    tss_all = gt["tss"].to_numpy()
    dens = np.zeros(n)
    for t in tss_all:
        b = (t - start) // GQ.RES10
        if 0 <= b < n:
            dens[b] += 1
    comp, _ = DOM.compartments(ci, cj, v, n, dens, max_bins=100)
    w = 20
    ins = DOM.insulation(ci, cj, v, n, w)
    bnd = np.array(DOM.boundaries(ins, w))
    cov = Mx.sum(axis=1)
    lc = np.where(cov > 0, np.log1p(cov), np.nan)
    lc = (lc - np.nanmean(lc)) / (np.nanstd(lc) or 1.0)
    exp = np.array([np.mean(np.diagonal(Mx, d)) for d in range(11)])
    rows = []
    for r in gt.itertuples(index=False):
        if not (start + GQ.EDGE <= r.tss < start + GQ.WINDOW - GQ.EDGE):
            continue
        b = int((r.tss - start) // GQ.RES10)
        if cov[b] <= 0 or r.name not in tpm:
            continue
        near = np.sum(np.abs(tss_all - r.tss) <= 250_000) - 1
        d = float(np.min(np.abs(bnd - b))) * GQ.RES10 if len(bnd) else float(GQ.WINDOW)
        oe = [math.log2((Mx[b, b + k] + 1) / (exp[abs(k)] + 1)) for k in range(-10, 11) if k and 0 <= b + k < n]
        sgn = 1 if r.strand == "+" else -1
        slope = float(ins[b + 5 * sgn] - ins[b - 5 * sgn]) if 5 <= b < n - 5 else np.nan
        f = {"compartment": comp[b], "insulation": ins[b], "coverage": lc[b], "log_dist_boundary": math.log10(d + GQ.RES10),
             "gene_density": math.log1p(near), "local_oe": float(np.mean(oe)), "insulation_slope": slope}
        if not all(np.isfinite(list(f.values()))):
            continue
        t = tpm[r.name]
        rows.append({"cell": cell, "gene": r.name, "chrom": chrom, "tss": int(r.tss), **f, "tpm": t,
                     "label": int(t >= GQ.TPM_ON)})
    return rows


def genes_q4b(cell: str, regions: list) -> list[dict]:
    path = GQ.CACHE / f"q4b_genes_{cell}_{'_'.join(c + str(s) for c, s in regions)}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    tpm = encode_tpm(cell)
    rows = [g for c, s in regions for g in gene_rows(cell, c, s, tpm)]
    GQ._json(path, rows)
    return rows


def _xy(rows: list[dict], feats: list[str]) -> tuple[np.ndarray, np.ndarray]:
    return np.array([[r[f] for f in feats] for r in rows], float), np.array([r["label"] for r in rows], int)


def _loco(kind: str, prm: dict, by_cell: dict, feats: list[str]) -> dict:
    """Leave-one-cell-line-out AUC: train on the other seen cell lines, score the left-out one."""
    from sklearn.metrics import roc_auc_score
    out = {}
    for held in by_cell:
        tr = [r for c, rows in by_cell.items() if c != held for r in rows]
        Xtr, ytr = _xy(tr, feats)
        Xte, yte = _xy(by_cell[held], feats)
        out[held] = float(roc_auc_score(yte, GQ._fit_predict(kind, prm, Xtr, ytr, Xte)))
    out["mean"] = float(np.mean(list(out.values())))
    return out


def practice_q4b() -> dict:
    by_cell = {c: genes_q4b(c, regs) for c, regs in SEEN.items()}
    out = {"made": _now(), "genes": {c: len(v) for c, v in by_cell.items()},
           "expressed_fraction": {c: float(np.mean([r["label"] for r in v])) for c, v in by_cell.items()},
           "agreement_with_gtex_labels": {}, "feature_sets": {}}
    for c in ("gm12878", "imr90"):                     # the new label source against Q4's GTEx labels
        old = {r["gene"]: r["label"] for r in GQ.genes_for(c, SEEN[c])}
        both = [(old[r["gene"]], r["label"]) for r in by_cell[c] if r["gene"] in old]
        out["agreement_with_gtex_labels"][c] = {"genes": len(both),
                                                "agree": float(np.mean([a == b for a, b in both])) if both else None}
    print(out["genes"], out["expressed_fraction"], out["agreement_with_gtex_labels"], flush=True)
    for name, feats in (("q4", FEATURES_Q4), ("q4b", FEATURES_Q4B)):
        qs = [{**p, **_loco("qsvm", p, by_cell, feats)} for p in Q4B_QSVM_GRID]
        rb = [{**p, **_loco("rbf", p, by_cell, feats)} for p in Q4B_RBF_GRID]
        lr = _loco("logistic", {}, by_cell, feats)
        fs = {"features": feats, "qsvm_best": max(qs, key=lambda r: r["mean"]), "rbf_best": max(rb, key=lambda r: r["mean"]),
              "logistic": lr, "qsvm_grid": qs, "rbf_grid": rb}
        out["feature_sets"][name] = fs
        print(name, "qsvm", fs["qsvm_best"], "rbf", fs["rbf_best"], "logistic", lr, flush=True)
    best = max(out["feature_sets"], key=lambda k: out["feature_sets"][k]["qsvm_best"]["mean"])
    fs = out["feature_sets"][best]
    keep = ("bandwidth", "reps", "C", "gamma")
    out["choice"] = {"features": best, "qsvm": {k: v for k, v in fs["qsvm_best"].items() if k in keep},
                     "rbf": {k: v for k, v in fs["rbf_best"].items() if k in keep}}
    _json(ROOT / "results_round2_practice_q4b.json", out)
    print("choice", out["choice"], flush=True)
    return out


def test_q4b(R: dict) -> dict:
    from sklearn.metrics import roc_auc_score
    feats = FEATURES_Q4B if R["features"] == "q4b" else FEATURES_Q4
    tr = [r for c, regs in SEEN.items() for r in genes_q4b(c, regs)]
    te = genes_q4b(R["test_cell"], [tuple(x) for x in R["test_regions"]])
    Xtr, ytr = _xy(tr, feats)
    Xte, yte = _xy(te, feats)
    sc = {"qsvm": GQ._fit_predict("qsvm", R["qsvm"], Xtr, ytr, Xte), "rbf": GQ._fit_predict("rbf", R["rbf"], Xtr, ytr, Xte),
          "logistic": GQ._fit_predict("logistic", {}, Xtr, ytr, Xte)}
    auc = {k: float(roc_auc_score(yte, v)) for k, v in sc.items()}
    ci = GQ._boot_auc(yte, sc)
    a, b = KQ.scale_features(Xtr, Xte)
    qs = KQ.QSVM(R["qsvm"]["C"], R["qsvm"]["bandwidth"], R["qsvm"]["reps"], shots=R["qsvm_shots"]).fit(a, ytr)
    auc_shots = float(roc_auc_score(yte, qs.decision_function(b)))
    ok = bool(auc["qsvm"] >= auc["rbf"] - R["margin"] and ci["qsvm"][0] > 0.5)
    return {"train_genes": int(len(ytr)), "test_genes": int(len(yte)), "test_expressed_fraction": float(yte.mean()),
            "auc": auc, "ci95": ci, "qsvm_auc_with_shots": auc_shots, "pass": ok}


# ======================================================================================
# Q7b: docking with an empirical (Vina-like) score and rigid-body refinement
# ======================================================================================
Q7B_GRID = [{"score": "vina", "refine_top": k, "maxfev": 300} for k in (0, 10, 30)]


def _complex_pb(cid: str):
    import quantum_drug_gates as QD
    from chronocell.quantum import docking as DK
    zl, zp = QD.pb_zip("ligands"), QD.pb_zip("proteins")
    lig = DK.read_sdf(zl.read(f"posebusters_ligands_256/{cid}_ligand.sdf").decode(), cid)
    prot = DK.read_pdb(zp.read(f"posebusters_spruce_structures_256/{cid}_spruced.pdb").decode(errors="replace"))
    return lig, prot


def dock_q7b(cid: str, loader, settings: dict, qcfg: dict, seed: int, tag: str) -> dict:
    """One complex, three routes on the same graph and score: QAOA's sampled cliques, the exact maximal cliques, and a
    random search; each route's poses are scored with the Vina-like function and the best `refine_top` refined."""
    import hashlib
    import quantum_drug_gates as QD
    from scipy.spatial import cKDTree
    from chronocell.quantum import docking as DK, qubo as QB
    key = hashlib.sha1(json.dumps([settings, qcfg]).encode()).hexdigest()[:10]
    path = QD.CACHE / f"q7b_{tag}_{cid}_{key}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    lig, prot = loader(cid)
    heavy = lig.heavy
    L = lig.xyz[heavy]
    tree = cKDTree(prot.xyz[prot.el != "H"])
    rec = {"cid": cid, "usable": bool(tree.query(L)[0].min() < 4.0)}
    if rec["usable"]:
        cfg = QD.DOCK
        centre = L.mean(0)
        lt = DK.ligand_types(lig)
        pt = DK.protein_types(prot, centre, float(np.linalg.norm(L - centre, axis=1).max()) + 11.0)
        lf = DK.ligand_features(lig, cfg["max_polar"], cfg["max_hydrophobic"])
        hs = DK.hotspots_directional(prot, centre, per_type=cfg["per_type"])
        rec.update({"heavy_atoms": int(heavy.sum()), "crystal_score": DK.vina_score(L, lt, pt)})
        k = settings["refine_top"]

        def route(cliques, g):
            poses = []
            for _, c in cliques:
                bits = np.zeros(len(g.vertices))
                bits[c] = 1
                p = DK.pose_from_clique(lig.xyz, lf, hs, g, bits)
                if p is not None:
                    poses.append((DK.vina_score(p[heavy], lt, pt), p[heavy]))
            if not poses:
                return {"rmsd": float("nan"), "score": float("nan"), "rmsd_unrefined": float("nan"), "poses": 0}
            poses.sort(key=lambda z: z[0])
            out = {"rmsd_unrefined": DK.rmsd(poses[0][1], L), "poses": len(poses)}
            if k:
                best = min((DK.refine(p, lt, pt, settings["maxfev"]) for _, p in poses[:k]), key=lambda z: z[1])
            else:
                best = (poses[0][1], poses[0][0])
            out.update({"rmsd": DK.rmsd(best[0], L), "score": best[1]})
            return out

        if len(hs.kind) >= 3 and len(lf.kind) >= 3:
            full = DK.interaction_graph(lf, hs, tau=cfg["tau"], max_vertices=80)
            g, _ = DK.core_subgraph(full, cfg["qubits"])
            gpath = QD.CACHE / f"q7b_{tag}_{cid}_graph_{hashlib.sha1(json.dumps(qcfg).encode()).hexdigest()[:10]}.json"
            if gpath.exists():                         # QAOA and the exact cliques do not depend on the pose settings
                gr = json.loads(gpath.read_text(encoding="utf-8"))
            else:
                q = DK.clique_qubo(g)
                e = q.energies()
                ex = QB.exact(q, e)
                r = sim.qaoa(e, p=qcfg["p"], shots=4096, objective=qcfg["objective"], maxiter=80, seed=seed)
                gr = {"vertices": q.n, "qaoa_hit": bool(r.hit), "qaoa_p_optimal": r.p_optimal,
                      "uniform_p_optimal": ex.detail["n_optima"] / 2 ** q.n, "qaoa_seconds": r.seconds,
                      "qaoa_cliques": QD.sampled_cliques(g, r, cfg["cliques"]),
                      "exact_cliques": QD.maximal_cliques(g, cfg["cliques"])}
                gpath.parent.mkdir(parents=True, exist_ok=True)
                gpath.write_text(json.dumps(gr), encoding="utf-8")
            rec.update({k2: gr[k2] for k2 in ("vertices", "qaoa_hit", "qaoa_p_optimal", "uniform_p_optimal", "qaoa_seconds")})
            rec["qaoa"] = route(gr["qaoa_cliques"], g)
            rec["classical"] = route(gr["exact_cliques"], g)
        rb = DK.random_search_vina(L, centre, lt, pt, cfg["random_poses"], k, np.random.default_rng(seed),
                                   maxfev=settings["maxfev"])
        rec["random"] = {"rmsd": DK.rmsd(rb[0], L), "score": rb[1]}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec), encoding="utf-8")
    return rec


def summarise_q7b(recs: list[dict]) -> dict:
    use = [r for r in recs if r.get("usable")]
    ok = lambda route: float(np.mean([np.isfinite(r.get(route, {}).get("rmsd", np.nan)) and r[route]["rmsd"] <= 2.0
                                      for r in use])) if use else float("nan")
    qq = [r for r in use if "qaoa_hit" in r]
    return {"complexes": len(recs), "usable": len(use), "success_qaoa": ok("qaoa"), "success_classical": ok("classical"),
            "success_random": ok("random"),
            "success_qaoa_unrefined": float(np.mean([r.get("qaoa", {}).get("rmsd_unrefined", np.nan) <= 2.0 for r in use])),
            "qaoa_hit_rate": float(np.mean([r["qaoa_hit"] for r in qq])) if qq else float("nan"), "graphs": len(qq),
            "crystal_scores_lower_than_qaoa_pose": float(np.mean([r["crystal_score"] <= r["qaoa"]["score"] + 1e-9
                                                                  for r in qq if np.isfinite(r["qaoa"].get("score", np.nan))]))}


def practice_q7b() -> dict:
    """All 256 PoseBusters complexes of the Q7 data (practice and test there; practice now)."""
    import quantum_drug_gates as QD
    import frozen as F
    qcfg = F.QUANTUM_DRUG_GATES["q7"]["qaoa"]
    ids = QD.complexes("practice") + QD.complexes("test")
    out = {"made": _now(), "complexes": len(ids), "qaoa": qcfg, "grid": {}}
    for st in Q7B_GRID:
        recs = []
        t0 = time.perf_counter()
        for i, cid in enumerate(ids):
            recs.append(dock_q7b(cid, _complex_pb, st, qcfg, i, "pb"))
            if i % 20 == 0:
                print(st, i, cid, round(time.perf_counter() - t0), "s", flush=True)
        out["grid"][f"refine{st['refine_top']}"] = {"settings": st, **summarise_q7b(recs)}
        print(st, out["grid"][f"refine{st['refine_top']}"], flush=True)
        _json(ROOT / "results_round2_practice_q7b.json", out)
    best = max(out["grid"].values(), key=lambda g: (g["success_qaoa"], -g["settings"]["refine_top"]))
    out["choice"] = best["settings"]
    _json(ROOT / "results_round2_practice_q7b.json", out)
    return out


# ======================================================================================
# Q6b: ADAPT-VQE
# ======================================================================================
def _model(name: str, scale: float, act) -> tuple[M.QubitModel, M.HF]:
    _, fn, charge = M.LIBRARY[name]
    ints = M.integrals(fn(scale), charge)
    hf = M.rhf(ints)
    asp = M.active_space(ints, hf, act[1], act[0])
    return M.qubit_hamiltonian(asp), hf


def q6b_case(name: str, scale: float, act, settings: dict) -> dict:
    qm, hf = _model(name, scale, act)
    t0 = time.perf_counter()
    u = M.vqe_uccsd(qm)
    t_u = time.perf_counter() - t0
    t0 = time.perf_counter()
    a = M.adapt_vqe(qm, **settings)
    t_a = time.perf_counter() - t0
    row = {"molecule": name, "bond_scale": scale, "active": list(act), "qubits": qm.n_qubits, "hf_converged": hf.converged,
           "e_hf": a.hf, "e_cas_fci": a.fci, "correlation_mEh": 1e3 * (a.fci - a.hf),
           "uccsd_error_mEh": 1e3 * u.error, "uccsd_parameters": u.parameters, "uccsd_seconds": t_u,
           "adapt_energy": a.energy, "adapt_error_mEh": 1e3 * a.error, "adapt_operators": a.parameters,
           "adapt_converged": a.converged, "adapt_final_gradient": a.gradient_norms[-1], "adapt_electrons": a.electrons,
           "adapt_seconds": t_a, "adapt_sequence": a.operators}
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if k != "adapt_sequence"}, flush=True)
    return row


def practice_q6b() -> dict:
    out = {"made": _now(), "cases": Q6B_PRACTICE, "grid": {}}
    for st in Q6B_GRID:
        rows = [q6b_case(n, s, a, st) for n, s, a in Q6B_PRACTICE]
        out["grid"][st["pool"]] = {"settings": st, "rows": rows,
                                   "max_abs_error_mEh": max(abs(r["adapt_error_mEh"]) for r in rows),
                                   "within_chemical_accuracy": int(sum(abs(r["adapt_error_mEh"]) <= 1.6 for r in rows))}
    best = min(out["grid"].values(), key=lambda g: g["max_abs_error_mEh"])
    out["choice"] = best["settings"]
    _json(ROOT / "results_round2_practice_q6b.json", out)
    print("choice", out["choice"], {k: v["max_abs_error_mEh"] for k, v in out["grid"].items()}, flush=True)
    return out


def test_q6b(R: dict) -> dict:
    rows = [q6b_case(n, s, a, R["settings"]) for n, s, a in R["cases"]]
    ok = [abs(r["adapt_error_mEh"]) <= R["tolerance_mEh"] for r in rows]
    return {"rows": rows, "within": int(sum(ok)), "cases": len(rows),
            "max_abs_error_mEh": float(max(abs(r["adapt_error_mEh"]) for r in rows)),
            "uccsd_within": int(sum(abs(r["uccsd_error_mEh"]) <= R["tolerance_mEh"] for r in rows)),
            "pass": bool(all(ok))}


# ======================================================================================
# Test driver (each part once)
# ======================================================================================
TESTS = {"q2b": test_q2b, "q4b": test_q4b, "q6b": test_q6b}
PRACTICE = {"q2b": practice_q2b, "q4b": practice_q4b, "q6b": practice_q6b, "q7b": practice_q7b}


def run_test(part: str) -> dict:
    import frozen as F
    R = F.QUANTUM_ROUND2
    out = json.loads(RESULTS.read_text(encoding="utf-8")) if RESULTS.exists() else {"made": _now()}
    if part in out:
        print(f"{part} has already been run once ({RESULTS.name}); not run again.", flush=True)
        return out
    t0 = time.perf_counter()
    res = TESTS[part](R[part])
    res["rule"] = R[part]
    res["run_at"] = _now()
    res["seconds"] = time.perf_counter() - t0
    out[part] = res
    out["overall"] = {k: out[k]["pass"] for k in TESTS if k in out}
    _json(RESULTS, out)
    print(part, "pass" if res["pass"] else "fail", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", choices=sorted(PRACTICE))
    g.add_argument("--test", choices=sorted(TESTS))
    a = ap.parse_args()
    if a.practice:
        PRACTICE[a.practice]()
    else:
        run_test(a.test)


if __name__ == "__main__":
    main()
