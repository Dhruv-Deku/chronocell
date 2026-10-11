"""
Exports for other viewers (Phase B6). Coordinates are 0-based, half-open (BED convention).

IGV          BED (boundaries, domains), BEDPE (loops, differential pixels), bedGraph (insulation, compartment
             eigenvector). bigWig needs a compiled writer (pyBigWig or UCSC bedGraphToBigWig), neither of which
             installs on this Windows machine: convert the bedGraph with bedGraphToBigWig where available.
Juicebox     2D annotation files (loops / domains, Juicebox's own text format, loaded with "Load 2D annotations")
             and contacts in Juicer "short with score" text, which Juicer Tools `pre` turns into a .hic.
HiGlass      a multi-resolution .mcool (cooler schema v3, written with h5py; the coarser resolutions are sums),
             which HiGlass serves directly; beddb needs clodius and is not written here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def write_text(path: Path, text: str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def bed(df: pd.DataFrame, name_col: str | None = None) -> str:
    if df is None or not len(df):
        return ""
    lines = []
    for r in df.itertuples(index=False):
        row = [str(r.chrom), str(int(r.start)), str(int(r.end))]
        if name_col:
            row.append(str(getattr(r, name_col)))
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def bedgraph(df: pd.DataFrame, value: str, name: str | None = None) -> str:
    ok = df[np.isfinite(df[value].astype(float))]
    head = f'track type=bedGraph name="{name or value}"\n'
    return head + "".join(f"{r.chrom}\t{int(r.start)}\t{int(r.end)}\t{float(getattr(r, value)):.6g}\n"
                          for r in ok.itertuples(index=False))


def loops_bedpe(df: pd.DataFrame) -> str:
    if df is None or not len(df):
        return ""
    cols = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    extra = [c for c in df.columns if c not in cols]
    return df[cols + extra].to_csv(sep="\t", index=False, header=False)


def juicebox_2d(loops: pd.DataFrame | None = None, domains: pd.DataFrame | None = None) -> str:
    """Juicebox 2D annotation text (chr1 x1 x2 chr2 y1 y2 color comment)."""
    lines = ["chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\tcomment"]
    if loops is not None:
        for r in loops.itertuples(index=False):
            lines.append(f"{str(r.chrom1).removeprefix('chr')}\t{int(r.start1)}\t{int(r.end1)}\t"
                         f"{str(r.chrom2).removeprefix('chr')}\t{int(r.start2)}\t{int(r.end2)}\t0,255,255\tloop")
    if domains is not None:
        for r in domains.itertuples(index=False):
            c = str(r.chrom).removeprefix("chr")
            lines.append(f"{c}\t{int(r.start)}\t{int(r.end)}\t{c}\t{int(r.start)}\t{int(r.end)}\t255,255,0\tdomain")
    return "\n".join(lines) + "\n"


def juicer_short(ci: np.ndarray, cj: np.ndarray, cm: np.ndarray, chrom: str, start: int, resolution: int) -> str:
    """Juicer 'short with score' format: str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 score (one line per pixel;
    positions are bin starts). `java -jar juicer_tools.jar pre file.txt out.hic hg38` converts it."""
    c = chrom.removeprefix("chr")
    return "".join(f"0 {c} {start + int(i) * resolution} 0 0 {c} {start + int(j) * resolution} 1 {float(v):g}\n"
                   for i, j, v in zip(ci, cj, cm))


def write_mcool(path: Path, chroms: list[tuple[str, int]], pixels: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
                base_resolution: int, factors: tuple[int, ...] = (1, 2, 5, 10, 25), weights: dict | None = None) -> Path:
    """Multi-resolution cooler. pixels: chrom -> (bin_i, bin_j, count), chromosome-local bins at base_resolution,
    upper triangle (cis only). Coarser resolutions are block sums. weights: chrom -> per-bin balancing weights
    at the base resolution (stored as the 'weight' column, cooler's convention)."""
    import h5py
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.attrs["format"] = "HDF5::MCOOL"
        f.attrs["format-version"] = 2
        res_grp = f.create_group("resolutions")
        for k in factors:
            res = base_resolution * k
            g = res_grp.create_group(str(res))
            g.attrs.update({"format": "HDF5::Cooler", "format-version": 3, "bin-type": "fixed", "bin-size": res,
                            "storage-mode": "symmetric-upper", "generated-by": "ChronoCell-5D"})
            names = np.array([c for c, _ in chroms], dtype="S")
            lengths = np.array([L for _, L in chroms], dtype=np.int64)
            cg = g.create_group("chroms")
            cg.create_dataset("name", data=names)
            cg.create_dataset("length", data=lengths)
            nb = [int(-(-L // res)) for _, L in chroms]
            offs = np.concatenate([[0], np.cumsum(nb)])
            bchrom = np.concatenate([np.full(m, i, np.int32) for i, m in enumerate(nb)])
            bstart = np.concatenate([np.arange(m, dtype=np.int64) * res for m in nb])
            bend = np.concatenate([np.minimum((np.arange(m, dtype=np.int64) + 1) * res, L) for m, (_, L) in zip(nb, chroms)])
            bg = g.create_group("bins")
            bg.create_dataset("chrom", data=bchrom)
            bg.create_dataset("start", data=bstart)
            bg.create_dataset("end", data=bend)
            if weights is not None and k == 1:
                w = np.concatenate([np.asarray(weights.get(c, np.full(m, np.nan)), float)[:m] for (c, _), m in zip(chroms, nb)])
                bg.create_dataset("weight", data=w)
            B1, B2, C = [], [], []
            for ci_, (c, _) in enumerate(chroms):
                if c not in pixels:
                    continue
                i, j, v = (np.asarray(x) for x in pixels[c])
                a, b = np.minimum(i, j) // k, np.maximum(i, j) // k
                key = a * nb[ci_] + b
                u, inv = np.unique(key, return_inverse=True)
                s = np.bincount(inv, weights=np.asarray(v, float))
                B1.append(offs[ci_] + u // nb[ci_])
                B2.append(offs[ci_] + u % nb[ci_])
                C.append(s)
            b1 = np.concatenate(B1) if B1 else np.empty(0, np.int64)
            b2 = np.concatenate(B2) if B2 else np.empty(0, np.int64)
            cnt = np.concatenate(C) if C else np.empty(0)
            order = np.lexsort((b2, b1))
            b1, b2, cnt = b1[order], b2[order], cnt[order]
            pg = g.create_group("pixels")
            pg.create_dataset("bin1_id", data=b1.astype(np.int64))
            pg.create_dataset("bin2_id", data=b2.astype(np.int64))
            pg.create_dataset("count", data=np.rint(cnt).astype(np.int32) if np.allclose(cnt, np.rint(cnt)) else cnt)
            ig = g.create_group("indexes")
            ig.create_dataset("chrom_offset", data=offs.astype(np.int64))
            ig.create_dataset("bin1_offset", data=np.searchsorted(b1, np.arange(offs[-1] + 1)).astype(np.int64))
            g.attrs["nbins"], g.attrs["nchroms"], g.attrs["nnz"] = int(offs[-1]), len(chroms), int(len(b1))
    return path
