"""
Gate 8: does the Drug lab's mechanism simulator (chronocell/therapy.py) predict how a real drug changes 3D chromatin?

Data: MINA chromatin tracing of an 840 kb region of chrX (77.66-78.50 Mb, GRCh38; 30 kb loci) in IMR-90 cells
(Cheng et al. 2021, Wang lab; 4DN Data Portal, FOF-CT core + trace tables, MD5 from the portal, files from 4DN
Open Data on AWS): untreated (9 pooled replicates) and, one experiment each, cells treated with GSK126 (EZH2
inhibitor), Trichostatin A + sodium butyrate (HDAC inhibitors), 5-aza-2'-deoxycytidine (DNMT inhibitor), DMOG
(JmjC demethylase / 2-OG enzyme inhibitor) and alpha-amanitin (RNA polymerase II inhibitor). Each trace is labelled
active (Xa) or inactive (Xi) X.

Prediction (no fitting to treated data): the untreated median distance matrix of an allele is embedded in 3D
(classical MDS), the IMR-90 H3K27ac peaks (ENCODE ENCFF730BVO) give each locus its activity signal, and the drug
lab simulates the drug's class at full dose (efficacy 0.8, mechanism-only mode: no healthy baseline). Predicted
change of each pair: log(treated / untreated distance) in the simulated fold. Measured change: log of the ratio of
median distances (treated traces / untreated traces). Both are centred (global size changes removed: batch effects
between imaging experiments make the overall scale unreliable), and compared by Spearman rho over all locus pairs.
Uncertainty: 200 bootstrap resamples of the traces of both conditions. Permutation check: the same simulation with
the signal track shuffled across loci (200 permutations).

    python validation/drug_gate8.py --practice    # untreated split halves (noise) and alpha-amanitin (pipeline)
    python validation/drug_gate8.py --test        # run once (rule in frozen.DRUG_GATE8)
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                                     # noqa: E402
from chronocell import therapy as TH                     # noqa: E402

S3 = "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/"
SETS = {   # condition -> (4DN set, (core acc, uuid, md5), (trace acc, uuid, md5), drug-lab class)
    "untreated": ("4DNESDQV4VNJ", ("4DNFI7PBQK6G", "64986b2b-2765-4302-8ba1-a7085ae93c08", "3da41c14e199f5c50886115959a20eb0"),
                  ("4DNFIJ62XTT2", "a7eb90c4-a60f-46f9-920e-1ceb71190f27", "44c53f9c1f9b5969b01c84b06b22ccb8"), None),
    "GSK126 (EZH2 inhibitor)": ("4DNES73TYNAK", ("4DNFIZ9MVTME", "c31ea110-e76a-47bc-a272-564cb64da228", "cf444aa3455a2a852fd95033b936d6b8"),
                                ("4DNFIR3U55FW", "bb1e6e92-4382-431b-8cd5-d530e34edadc", "3628708bf9764e3b6272d7ce3e9d7cf7"), "ezh2"),
    "TSA + NaBu (HDAC inhibitors)": ("4DNESUDN4NP6", ("4DNFIEWN8S38", "9e496dd3-04eb-4a82-9119-0daefc5664be", "889be95da0049e304d4420a3635b6e64"),
                                     ("4DNFIZQ3ZIJQ", "071a9cfb-cdf6-4a9a-9f2e-6011754f9f6a", "2ce966afe7dacef4a3d044e2f842d130"), "hdac"),
    "5-aza-dC (DNMT inhibitor)": ("4DNES5XKLWKW", ("4DNFISVSJ9I2", "e1b5bb99-1e6b-4b25-a3e6-2cddc85f4d54", "44edd31178296726398727a05dd3d61b"),
                                  ("4DNFIVU9TG94", "5649415d-027e-48b9-8a81-3adb3f6c79ca", "7f37765c0df6d64e0148fe7f551e3a66"), "dnmt"),
    "DMOG (demethylase blocker)": ("4DNES4Y519FW", ("4DNFIYSCNG5N", "43d3e656-7a74-4b96-8672-e9a3c588bb99", "f50f06a1068e1571ba585b2188165584"),
                                   ("4DNFINVSXJOR", "ce63d67f-7139-462c-9d6b-9d78d5121b0f", "82cbcbf7bb85933712f8e301256ba83d"), "kdm"),
    "alpha-amanitin (transcription inhibitor)": ("4DNES82IH3OG", ("4DNFIN8AVPBA", "03125b53-6dae-4470-9830-322c3df7207d", "f49870e456de5f3f32d88ba0efb372a2"),
                                                 ("4DNFIKXYNRQ1", "c691a581-2307-4d1a-bfbf-74e7c5654895", "4fb171d8c02c9197b1343a0217b51749"), "txn"),
}
PRACTICE_DRUG = "alpha-amanitin (transcription inhibitor)"
TEST_DRUGS = ["GSK126 (EZH2 inhibitor)", "TSA + NaBu (HDAC inhibitors)", "5-aza-dC (DNMT inhibitor)", "DMOG (demethylase blocker)"]
MIN_DETECTED = 0.7
EFFICACY = 0.8
BOOT = 200
PERM = 200


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def fetch(acc: str, uuid: str, md5: str) -> Path:
    return D.fetch_url(f"{S3}{uuid}/{acc}.csv", D.DATA / "4dn_mina" / f"{acc}.csv", md5)


def load(condition: str) -> dict:
    """{allele: (traces (n, L, 3) nm, loci starts)} with traces detected at >= MIN_DETECTED of loci."""
    import pandas as pd
    _, core, trace, _ = SETS[condition]
    cp, tp = fetch(*core), fetch(*trace)
    rows = [ln.split(",")[:8] for ln in cp.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    df = pd.DataFrame(rows, columns=["spot", "trace", "x", "y", "z", "chrom", "start", "end"])
    for c in ("x", "y", "z"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["trace"] = df["trace"].astype(int)
    df["start"] = df["start"].astype(int)
    status = {}
    for ln in tp.read_text(encoding="utf-8").splitlines():
        if ln and not ln.startswith("#"):
            f = ln.split(",")
            if len(f) >= 2 and f[0].strip().isdigit():
                status[int(f[0])] = f[1].strip()
    loci = np.sort(df["start"].unique())
    out = {}
    for allele in ("Xa", "Xi"):
        ids = [t for t in df["trace"].unique() if status.get(int(t)) == allele]
        sub = df[df["trace"].isin(ids)]
        tid = pd.Index(sorted(ids))
        xyz = np.full((len(tid), len(loci), 3), np.nan)
        xyz[tid.get_indexer(sub["trace"]), np.searchsorted(loci, sub["start"].to_numpy())] = sub[["x", "y", "z"]].to_numpy(float)
        step = np.nanmedian(np.linalg.norm(np.diff(xyz, axis=1), axis=-1))
        xyz *= 1.0 if step > 20 else 1000.0                       # micron -> nm
        ok = np.isfinite(xyz[..., 0]).mean(axis=1) >= MIN_DETECTED
        out[allele] = (xyz[ok], loci)
    return out


def median_distances(xyz: np.ndarray) -> np.ndarray:
    d = np.linalg.norm(xyz[:, :, None, :] - xyz[:, None, :, :], axis=-1)
    with np.errstate(all="ignore"):
        return np.nanmedian(d, axis=0)


def mds(D: np.ndarray) -> np.ndarray:
    n = len(D)
    Dm = np.where(np.isfinite(D), D, np.nanmedian(D))
    J = np.eye(n) - 1.0 / n
    Bm = -0.5 * J @ (Dm ** 2) @ J
    w, V = np.linalg.eigh(Bm)
    idx = np.argsort(w)[::-1][:3]
    return V[:, idx] * np.sqrt(np.maximum(w[idx], 0))


def h3k27ac(loci: np.ndarray, width: int = 30_000) -> np.ndarray:
    p = D.mark_peaks_path("h3k27ac", "IMR90")
    sig = np.zeros(len(loci))
    with gzip.open(p, "rt") as fh:
        for ln in fh:
            f = ln.split("\t")
            if f[0] != "chrX":
                continue
            a, b, v = int(f[1]), int(f[2]), float(f[6])
            for k, s in enumerate(loci):
                ov = min(b, s + width) - max(a, s)
                if ov > 0:
                    sig[k] += v * ov / (b - a)
    return sig


def predicted(x0: np.ndarray, sig: np.ndarray, drug: str) -> np.ndarray:
    b0 = float(np.median(np.linalg.norm(np.diff(x0, axis=0), axis=1)))
    r = TH.simulate_treatment(x0, sig, np.ones(len(x0), bool), b0, drug, None, EFFICACY, np.array([0.0, 1.0]))
    x1 = r.frames[-1]
    iu = np.triu_indices(len(x0), 1)
    d0 = np.linalg.norm(x0[:, None] - x0[None], axis=-1)[iu]
    d1 = np.linalg.norm(x1[:, None] - x1[None], axis=-1)[iu]
    p = np.log(d1 / d0)
    return p - p.mean()


def measured(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    iu = np.triu_indices(a.shape[1], 1)
    m = np.log(median_distances(b)[iu] / median_distances(a)[iu])
    return m - np.nanmean(m)


def _rho(p: np.ndarray, m: np.ndarray) -> float:
    from scipy.stats import spearmanr
    ok = np.isfinite(p) & np.isfinite(m)
    return float(spearmanr(p[ok], m[ok]).statistic) if ok.sum() > 10 else float("nan")


def score(untreated: np.ndarray, treated: np.ndarray, loci: np.ndarray, drug: str, seed: int = 0) -> dict:
    x0 = mds(median_distances(untreated))
    sig = h3k27ac(loci)
    p = predicted(x0, sig, drug)
    m = measured(untreated, treated)
    rho = _rho(p, m)
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(BOOT):
        a = untreated[rng.integers(0, len(untreated), len(untreated))]
        b = treated[rng.integers(0, len(treated), len(treated))]
        boot.append(_rho(p, measured(a, b)))
    perm = []
    for _ in range(PERM):
        perm.append(_rho(predicted(x0, rng.permutation(sig), drug), m))
    boot = np.array(boot)
    perm = np.array(perm)
    return {"drug_class": drug, "rho": rho, "ci95": [float(np.nanpercentile(boot, 2.5)), float(np.nanpercentile(boot, 97.5))],
            "perm_mean": float(np.nanmean(perm)), "perm_p": float((np.sum(perm >= rho) + 1) / (len(perm) + 1)),
            "n_untreated": int(len(untreated)), "n_treated": int(len(treated)), "loci": int(len(loci)),
            "measured_global_log_change": float(np.nanmean(np.log(median_distances(treated)[np.triu_indices(len(loci), 1)] /
                                                                    median_distances(untreated)[np.triu_indices(len(loci), 1)]))),
            "predicted_direction": int(TH.ALL_DRUGS[drug].direction)}


def run_practice() -> dict:
    un = load("untreated")
    out = {"made": _now(), "split_halves": {}, "practice_drug": {}}
    rng = np.random.default_rng(1)
    for allele, (xyz, loci) in un.items():
        idx = rng.permutation(len(xyz))
        a, b = xyz[idx[: len(xyz) // 2]], xyz[idx[len(xyz) // 2:]]
        out["split_halves"][allele] = {k: score(a, b, loci, k) for k in ("ezh2", "hdac", "dnmt", "kdm", "txn")}
        print(allele, {k: round(v["rho"], 3) for k, v in out["split_halves"][allele].items()}, flush=True)
    tr = load(PRACTICE_DRUG)
    for allele in ("Xa", "Xi"):
        out["practice_drug"][allele] = score(un[allele][0], tr[allele][0], un[allele][1], SETS[PRACTICE_DRUG][3])
        print(PRACTICE_DRUG, allele, out["practice_drug"][allele], flush=True)
    (ROOT / "results_gate8_practice.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


def run_test() -> dict:
    import frozen as F
    R = F.DRUG_GATE8
    un = load("untreated")
    out = {"made": _now(), "rule": R, "drugs": {}}
    for cond in TEST_DRUGS:
        tr = load(cond)
        out["drugs"][cond] = {allele: score(un[allele][0], tr[allele][0], un[allele][1], SETS[cond][3])
                              for allele in ("Xa", "Xi")}
        x = out["drugs"][cond][R["allele"]]
        x["passes"] = bool(x["rho"] >= R["min_rho"] and x["ci95"][0] > 0 and x["perm_p"] <= R["max_perm_p"])
        print(cond, {a: round(v["rho"], 3) for a, v in out["drugs"][cond].items()}, "pass" if x["passes"] else "fail", flush=True)
    n = sum(v[R["allele"]]["passes"] for v in out["drugs"].values())
    out["passing_drugs"] = int(n)
    out["verdict"] = "pass" if n >= R["min_drugs"] else "fail"
    (ROOT / "results_gate8.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(out["verdict"], flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    run_practice() if a.practice else run_test()


if __name__ == "__main__":
    main()
