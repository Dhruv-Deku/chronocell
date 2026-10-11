"""
Gate 5m: the frozen human Gate 5 predictor (sequence + CTCF, validation/predictor_model.json) on mouse
ES-cell chromatin tracing. Pre-registered in validation/frozen.py (PREDICTOR_MOUSE) before any mouse data
were read; nothing here is fitted on mouse data.

    python validation/predictor_mouse.py --test      # run once -> validation/results_predictor_mouse.json

Truth: ORCA traces of Hafner et al., Mol Cell 83:1377 (2023), 4DN FOF-CT "DNA-spot/trace core" tables. The
4DN portal's own download links need an account key; the same files are public on the 4DN AWS Open Data
bucket (each file's `open_data_url` in the portal's public metadata), and are checked against the MD5 the
portal publishes. Sets with several experiments (replicates) are pooled. Traces with fewer than half of their
loci detected are dropped (the rule used for every tracing set here).

Inputs: ENCODE CTCF IDR thresholded peaks of mouse ES-Bruce4 (mm10) and the UCSC mm10 chromosome sequence
(MD5 from UCSC md5sum.txt). Scoring: Gate 5's, unchanged (half-B medians as truth, all pairs, 3 splits,
95 % interval from split 0 with 100 resamples, as for the other 30 kb tracing sets).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import urllib.request
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D                              # noqa: E402
import protocol as PR                             # noqa: E402
from frozen import PREDICTOR as CFG, PREDICTOR_MOUSE as RULE   # noqa: E402
from chronocell import predict as PD              # noqa: E402

OUT = ROOT / "results_predictor_mouse.json"
DATA = D.DATA / "4dn"
S3 = "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/"
# 4DN experiment set -> (description, [(core table accession, open-data path, MD5 published by the 4DN portal)])
SETS = {
    "4DNESD28H8O7": ("chr6 locus, untreated", [("4DNFILYJ6RJ7", "64be8671-c1b1-430d-8675-d78ab53c02a8", "ff7537b67e63631c36414c174fb504da")]),
    "4DNESWDXDZSE": ("chr3 locus, untreated", [("4DNFIA63SK83", "f96751cd-86a7-4248-8e46-abfacc5d1a2a", "9a83829564859f148bfe1b2e696586bb")]),
    "4DNESJ3TXVIR": ("chr6 locus, CTCF-AID untreated", [
        ("4DNFIZWUQTXF", "644f7b9a-f4e0-4bbe-b8c3-bff0c819ad9d", "0858436e2533e530721d23eca5aeca6e"),
        ("4DNFIKVZTN9K", "d6372243-7901-4bc5-88e8-d250a1b202e9", "d2ef9f30e65949a11d0413abe4af5927"),
        ("4DNFID7HZ8LX", "9a169672-dcdf-4f8f-829b-fade994ce38f", "c58d1e324be448ec7f6c96a97182dde8")]),
    "4DNESQ49IXDU": ("chr6 locus, RAD21-AID untreated", [
        ("4DNFIRKSEW93", "ccb0c740-d19f-4bf3-b73b-1bca319ca598", "9a60058a12b0bb10b1dbbee99fd398c1"),
        ("4DNFIJIZ3RQH", "cbc1a19c-576f-48ae-aff1-7e0f31dab043", "6e66e8e564ff4b2ce55eed0b96d0178a"),
        ("4DNFIFW72AJF", "a7f1bba4-03a6-4021-ba35-4e9cc4e672bc", "b1e8e28f44718fdd40831a6b8bb00369")]),
    "4DNESTNG39BO": ("chr3 locus, CTCF-AID untreated", [
        ("4DNFI2OG3IPB", "30a0322b-b661-42bb-b908-073d14bf8b5a", "c2abd0b8cce1c4e2c47b1ff94465a9a9"),
        ("4DNFIU5X9MXT", "3db16c1a-91c8-40ca-80b2-556a065ed01f", "1236ad81e534c3b2c464582318810eaf"),
        ("4DNFIIQ8YIG4", "87667821-d014-46c9-9f6e-be12e64fe0c7", "15ec7ef14f6f00748cb014ebb7720611")]),
    "4DNESLRTOSQT": ("chr3 locus, RAD21-AID untreated", [
        ("4DNFI4SYERQD", "14cfb096-f08d-4c3b-8d48-71218da6af9b", "13e2c9b7072959ce9d112f6e6edd509b"),
        ("4DNFILHN6ZHZ", "15df1440-77c4-4964-a011-82dc296d1f89", "5a70533e9e8815a1ca0849921e40dac0"),
        ("4DNFIPUTK6GO", "27bcd217-ed49-441d-9364-dd85616512b6", "d2ba5f83f65d8ccc00f14638aef67b6f")]),
    "4DNES2KX6HQ5": ("chr6 locus, CTCF-AID + auxin", [
        ("4DNFIJYFR8FT", "fe3c5350-6d7a-4566-914c-e7b9bb9c305b", "99941665ab6ab50f653b4536bb23948b"),
        ("4DNFIIOVHG1I", "469582c7-8d3d-4a42-a179-10f4a7803b6d", "288d404eb132f8e77baa260db06619fd"),
        ("4DNFI4VLY3WU", "35e54b5f-8945-4a2b-a554-316e84bad5ae", "9064ac673d857cec433b36f47be7834c")]),
    "4DNESMN7RCSB": ("chr6 locus, RAD21-AID + auxin", [
        ("4DNFICN9NE8U", "bf0df72e-52b2-49cb-b5ff-c2b81d924f18", "566054f18f489f0caf64cc116580826d"),
        ("4DNFIBI4PIOM", "d520c455-dc51-4f07-a58a-ca52ee2b1c90", "32f3e88bdf8ff02d23182942d371b984"),
        ("4DNFI4AKVCVH", "8ad2bb7c-9e59-4fe9-b74d-8738f7b630e8", "ccb9fa556fbc6aca325118e3b1b34323")]),
    "4DNESG62SAVA": ("chr3 locus, CTCF-AID + auxin", [
        ("4DNFILD6TDR6", "fdaaf048-629f-4f3f-97af-16a47a9deb03", "26c6dc23119e7e6470b41d6440d387aa"),
        ("4DNFIKXR227H", "15d17583-a7ca-41ab-a6e8-f22e857c0c8e", "132da6ac223a6342079a73b6e50516e0"),
        ("4DNFIW3V9CRO", "0f9367da-7264-424e-9864-6688d68dede4", "575fcfc7a0a9bb4b5840276da68eb292")]),
    "4DNESBH54BG2": ("chr3 locus, RAD21-AID + auxin", [
        ("4DNFI1MBPSTC", "aab3e7d5-6122-46d5-af06-7d23e9e4f87f", "8d18b40b566af53168ddaae2bfb7ecc0"),
        ("4DNFIHA5PA9U", "b02e02c0-6bce-495a-ba3f-084fafe0cc9c", "16b684f66e319a61d656309e2a55f86c"),
        ("4DNFIDNY8W7L", "78348eef-7fb5-4486-9b88-2b79e890fbab", "71410803e444f2c4e7a88a048ac4bb22")]),
}
CITE = ("Hafner A, Park M, Berger SE, Murphy SE, Nora EP, Boettiger AN. Loop stacking organizes genome folding from TADs "
        "to chromosomes. Mol Cell 83, 1377-1392 (2023). 4DN Data Portal (4DN Open Data on AWS).")
UCSC_MM10 = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/"
BOOT = 100
R_C = 150.0
MIN_DETECTED = 0.5


def core_paths(set_acc: str) -> list[Path]:
    return [D.fetch_url(f"{S3}{uuid}/{acc}.csv", DATA / f"{acc}.csv", md5) for acc, uuid, md5 in SETS[set_acc][1]]


def read_fofct(paths: list[Path], min_detected: float = MIN_DETECTED) -> D.Traces:
    """FOF-CT core tables (Spot_ID, Trace_ID, X, Y, Z, Chrom, Chrom_Start, Chrom_End, Cell_ID) -> traces in nm.
    Trace ids are made unique per file. The coordinates are in nm whatever the XYZ_unit header line says when
    the median step between consecutive loci exceeds 20 units (ORCA files carry 'micron' with nm values)."""
    import pandas as pd
    frames, unit = [], ""
    for k, p in enumerate(paths):
        lines = p.read_text(encoding="utf-8").splitlines()
        unit = next((ln.split("=", 1)[1].split(",")[0].strip() for ln in lines if ln.startswith("##XYZ_unit")), unit)
        rows = [ln for ln in lines if ln.strip() and not ln.lstrip('"').startswith("#")]
        df = pd.DataFrame([r.split(",")[:8] for r in rows],      # Cell_ID (9th column) is absent in some files
                          columns=["spot", "trace", "x", "y", "z", "chrom", "start", "end"])
        df["trace"] = df["trace"].astype(np.int64) + k * 10_000_000
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    for c in ("x", "y", "z"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["start"] = df["start"].astype(np.int64)
    df["end"] = df["end"].astype(np.int64)
    chroms = df["chrom"].unique()
    if len(chroms) != 1:
        raise ValueError(f"expected one chromosome, found {list(chroms)}")
    # replicate files give the same locus starts 1-2 bp apart: loci are matched to the nearest kb
    df["locus"] = np.rint(df["start"] / 1000.0).astype(np.int64)
    keys = np.sort(df["locus"].unique())
    loci = df.groupby("locus")["start"].min().reindex(keys).to_numpy()
    ends = df.groupby("locus")["end"].max().reindex(keys).to_numpy()
    traces = df["trace"].unique()
    li = np.searchsorted(keys, df["locus"].to_numpy())
    ti = pd.Index(traces).get_indexer(df["trace"])
    xyz = np.full((len(traces), len(loci), 3), np.nan)
    xyz[ti, li] = df[["x", "y", "z"]].to_numpy(float)
    step = np.nanmedian(np.linalg.norm(np.diff(xyz, axis=1), axis=-1))
    scale = 1.0 if step > 20 else 1000.0
    xyz *= scale
    ok = np.isfinite(xyz[..., 0]).mean(axis=1) >= min_detected
    return D.Traces("+".join(p.stem for p in paths), xyz[ok], loci.astype(np.int64), np.full(len(loci), str(chroms[0])),
                    np.zeros(int(ok.sum()), int),
                    {"copies_total": int(len(traces)), "min_detected": min_detected, "xyz_unit_header": unit,
                     "scale_to_nm": scale, "locus_bp": int(np.median(ends - loci)), "files": [p.name for p in paths]})


def mm10_fasta(chrom: str) -> Path:
    sums = urllib.request.urlopen(urllib.request.Request(UCSC_MM10 + "md5sum.txt", headers=D.UA), timeout=60).read()
    md5 = {ln.split()[1]: ln.split()[0] for ln in sums.decode().splitlines() if len(ln.split()) == 2}[f"{chrom}.fa.gz"]
    return D.fetch_url(UCSC_MM10 + f"{chrom}.fa.gz", D.DATA / "mm10" / f"{chrom}.fa.gz", md5)


def mouse_peaks(chrom: str, seq: bytes) -> dict:
    import pandas as pd
    acc, _, md5 = RULE["peaks"]
    path = D.fetch_url(f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz", D.DATA / "encode" / f"{acc}.bed.gz", md5)
    df = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    df = df[df[0] == chrom]
    s, e, summit = df[1].to_numpy(np.int64), df[2].to_numpy(np.int64), df[9].to_numpy(np.int64)
    pwm = PD.log_odds(PD.read_jaspar(D.jaspar_ctcf_path().read_text()), CFG["pseudocount"])
    strand, score = PD.orient_peaks(seq, s, e, summit, pwm, CFG["motif_min_relative_score"], CFG["summit_half_width_bp"])
    mid = np.where(summit >= 0, s + summit, (s + e) // 2)
    return {"mid": mid, "strand": strand, "score": score}


def design(tr: D.Traces) -> dict:
    chrom = str(tr.chrom[0])
    seq = PD.read_fasta(mm10_fasta(chrom))
    pk = mouse_peaks(chrom, seq)
    starts = np.asarray(tr.starts, np.int64)
    ends = starts + int(tr.meta["locus_bp"])
    loci = PD.annotate(starts, ends, pk["mid"], pk["strand"], PD.GCIndex(seq))
    i, j = np.triu_indices(loci.n, 1)
    X = PD.pair_features(loci, pk["mid"], i, j)
    sep = np.rint(np.abs((starts[j] + ends[j]) - (starts[i] + ends[i])) / 2).astype(np.int64)
    return {"loci": loci, "i": i, "j": j, "X": X, "sep": sep, "chrom": chrom,
            "motif_rate": float(np.mean(pk["strand"] != 0)) if len(pk["strand"]) else float("nan"),
            "peaks_on_chrom": int(len(pk["mid"]))}


def _matrix(n: int, i: np.ndarray, j: np.ndarray, v: np.ndarray) -> np.ndarray:
    m = np.zeros((n, n))
    m[i, j] = v
    return m + m.T


def test() -> None:
    from validation.benchmark import methods as M, run as BR
    if OUT.exists():
        sys.exit(f"{OUT.name} exists: Gate 5m runs once.")
    model = PD.Predictor.from_json((ROOT / "predictor_model.json").read_text())
    roles = {k: "test" for k in RULE["test"]} | {k: "secondary" for k in RULE["secondary"]} | \
            {k: "control" for k in RULE["control"]}
    rows, info = [], {}
    for k, role in roles.items():
        tr = read_fofct(core_paths(k))
        d = design(tr)
        n = d["loci"].n
        sep_m = PR.separation(n, tr.starts)
        masks = BR._masks(n, PR.tiles(n, 400))
        pred = np.exp(_matrix(n, d["i"], d["j"], model.log_distance(d["X"], d["sep"].astype(float))))
        a, b = model.trend
        base = np.exp(_matrix(n, d["i"], d["j"], a + b * np.log10(np.maximum(d["sep"], 1.0))))
        np.fill_diagonal(pred, 0.0)
        np.fill_diagonal(base, 0.0)
        info[k] = {"role": role, "description": SETS[k][0], "loci": n, "copies": int(tr.n_copies),
                   "chrom": d["chrom"], "region": f"{d['chrom']}:{int(tr.starts[0]):,}-{int(tr.starts[-1]) + tr.meta['locus_bp']:,}",
                   "ctcf_peaks_on_chrom": d["peaks_on_chrom"], "motif_rate": d["motif_rate"], **tr.meta}
        for split in range(RULE["splits"]):
            st = BR._stats(f"mouse_{k}", tr.xyz, split, R_C)
            xyz_b = tr.xyz[st["b_idx"]]
            for name, mat in (("predictor", pred), ("genomic_trend_no_data", base)):
                sc = BR._score_pred(M.Prediction(mat), st, sep_m, masks, xyz_b, BOOT, None, split == 0, split)
                rows.append({"dataset": k, "role": role, "split": split, "method": name, "loci": n, "status": "ok",
                             "input": "sequence + CTCF", **sc})
                ap = sc.get("all_pairs", {})
                print(f"{k} ({SETS[k][0]:32s}) s{split} {name:22s} {ap.get('percent_of_ceiling', np.nan):6.1f}% raw "
                      f"{ap.get('raw_spearman', np.nan):.3f}", flush=True)
    summ = BR.summarise(rows)
    passes = {}
    for k in RULE["test"]:
        s = summ.get(k, {}).get("sequence + CTCF", {}).get("all_pairs", {})
        p, g = s.get("predictor"), s.get("genomic_trend_no_data")
        ci = p.get("ci95", {}).get("percent_of_ceiling", [np.nan, np.nan])
        passes[k] = {"i": bool(p["percent_of_ceiling"] > 0 and ci[0] > 0), "ii": bool(p["raw_spearman"] > g["raw_spearman"]),
                     "percent_of_ceiling": p["percent_of_ceiling"], "ci95": ci, "raw_spearman": p["raw_spearman"],
                     "baseline_raw_spearman": g["raw_spearman"]}
    verdict = "pass" if all(v["i"] and v["ii"] for v in passes.values()) else "fail"
    OUT.write_text(json.dumps({"mode": "test (run once; human model unchanged)", "rule": RULE, "sets": info,
                               "gate5m": passes, "verdict": verdict, "summary": summ, "rows": rows, "citation": CITE,
                               "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")},
                              indent=1, default=float), encoding="utf-8")
    print(json.dumps(passes, indent=1, default=float))
    print(f"Gate 5m: {verdict.upper()} -> {OUT}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", action="store_true", required=True)
    ap.parse_args()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        test()


if __name__ == "__main__":
    main()
