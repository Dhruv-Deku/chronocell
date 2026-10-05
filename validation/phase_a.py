"""
Phase A (accuracy core): the units every Phase A gate scores, built once and cached.

A unit is one (dataset, input, Hi-C source, depth, split): the population model fitted to the input
(imaging-derived contacts of half A, or sequencing Hi-C on the imaged loci) and the held-out truth
(half B of the imaged copies). Gates 1c (sizes), 2d (intervals), 2e (reliability) and 3b (learned
correction) all read these units, so every gate scores the same fits.

Hi-C sources (all raw counts, MAPQ >= 30):
  rao2014                Rao et al. 2014 in situ Hi-C, as the earlier gates use it (validation/datasets.py)
  encode_<cell>_<prot>   ENCODE "mapping quality thresholded contact matrix" .hic files, GRCh38, read
                         remotely by HTTP range (only the blocks of the region are fetched; the MD5 the
                         portal publishes covers the whole file, so it cannot be checked on a range read;
                         the bytes fetched are recorded with each cached region, as for rao2014)
Depth: `thin` < 1 keeps each read with that probability (binomial thinning of the counts, seeded), a
lower-depth library of the same cells.

    python validation/phase_a.py --build practice      # cache every practice unit (resumable)
    python validation/phase_a.py --build test          # cache every test unit (only gates read them)
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import datasets as D            # noqa: E402
import protocol as PR           # noqa: E402
from chronocell import ensemble as E, physics, population as P   # noqa: E402

CACHE = ROOT / "data" / "cache" / "phaseA"
ENCODE_HIC = {   # key: (file accession, cell line, protocol, MD5 published by the ENCODE portal); GRCh38, MAPQ >= 30
    "encode_k562_intact": ("ENCFF621AIY", "K562", "intact", "71e5dc072beddbdcda11c9016696e142"),
    "encode_hct116_insitu": ("ENCFF750AOC", "HCT116", "in situ", "7301909c0600d813ebfb73bc8d37b4bb"),
    "encode_hct116_intact": ("ENCFF573OPJ", "HCT116", "intact", "b74eec26bd44e3867812185fe3ca3cfb"),
    "encode_imr90_intact": ("ENCFF281ILS", "IMR90", "intact", "588e6252f5ddf47fa78e3c22ca01aca2"),
    "encode_imr90_dilution": ("ENCFF636CEM", "IMR90", "dilution", "3a5574993799b4e6ea2d5175df9c610d"),
    "encode_a549_insitu": ("ENCFF689CUX", "A549", "in situ", "8bee09203e1d4dbc482f13354a8dec45"),
}
for _k, (_acc, _cell, _prot, _md5) in ENCODE_HIC.items():
    D.HIC_SOURCES.setdefault(_k, (f"https://www.encodeproject.org/files/{_acc}/@@download/{_acc}.hic", "GRCh38", _cell))
PROTOCOL = {"rao2014": "in situ", **{k: v[2] for k, v in ENCODE_HIC.items()}}


@dataclass(frozen=True)
class UnitSpec:
    dataset: str                 # validation/datasets.py key, or "su_genome*:chrN" for one chromosome
    input: str                   # "imaging" | "hic"
    source: str = "rao2014"      # Hi-C source (hic input only)
    thin: float = 1.0            # fraction of reads kept (hic input only)
    split: int = 0

    @property
    def name(self) -> str:
        tail = "" if self.input == "imaging" else f"__{self.source}__t{self.thin:g}"
        return f"{self.dataset.replace(':', '-')}__{self.input}{tail}__s{self.split}"


def _hic_counts(spec: UnitSpec, tr, chrom_cols) -> np.ndarray:
    base = spec.dataset.split(":")[0]
    e = D.REGISTRY[base]
    if spec.source == "rao2014":
        if e.kind == "bintu_csv":
            return D.bintu_hic(base)[0]
        H, hch, hst = D.load_su_hic(D.REGISTRY[e.paired_hic])
        starts = tr.starts if chrom_cols is None else tr.starts[chrom_cols]
        chrom = None if chrom_cols is None else str(np.asarray(tr.chrom)[chrom_cols][0])
        sel = np.ones(len(hst), bool) if chrom is None else (hch == chrom)
        idx = np.flatnonzero(sel)[np.searchsorted(hst[sel], starts)]
        if not np.array_equal(hst[idx], starts):
            raise KeyError(f"{spec.dataset}: the Hi-C loci do not cover the imaged loci")
        return H[np.ix_(idx, idx)]
    if e.kind != "bintu_csv":
        raise KeyError(f"{spec.source} is only paired with the Bintu regions")
    if D.HIC_SOURCES[spec.source][2] != e.cell_line:
        raise KeyError(f"{spec.source} is {D.HIC_SOURCES[spec.source][2]} Hi-C, {base} is {e.cell_line}")
    n = len(tr.starts)
    starts38 = e.region_start_hg38 + np.arange(n) * 30_000          # GRCh38 coordinates of the 30 kb segments
    return D.hic_for_segments(spec.source, "chr21", starts38, 30_000)


def thin_counts(counts: np.ndarray, frac: float, seed: int = 0) -> np.ndarray:
    """Keep each read with probability frac (binomial thinning), symmetric, diagonal kept as is."""
    if frac >= 1.0:
        return counts
    c = np.rint(np.asarray(counts, dtype=np.float64)).astype(np.int64)
    rng = np.random.default_rng(seed)
    iu = np.triu_indices(len(c), 1)
    out = np.zeros_like(c)
    out[iu] = rng.binomial(c[iu], frac)
    out = out + out.T
    np.fill_diagonal(out, np.diag(c))
    return out.astype(np.float64)


@lru_cache(maxsize=3)
def _load_traces(base: str):
    return D.load(base)          # genome-scale files are hundreds of MB: parse once per process


def _traces(spec: UnitSpec):
    base, _, chrom = spec.dataset.partition(":")
    tr = _load_traces(base)
    if not chrom:
        return tr, None, tr.xyz.astype(np.float64)
    cols = np.asarray(tr.chrom) == chrom
    own = np.isfinite(tr.xyz[:, cols, 0]).mean(axis=1) >= 0.5
    return tr, cols, tr.xyz[own][:, cols].astype(np.float64)


def _fit(freq, seen, r_c: float, n: int, split: int):
    if n <= P.V33_MAX_BEADS:
        return E.fit_ensemble(freq, seen, r_c_nm=r_c, cfg=E.EnsembleConfig(seed=split))
    from frozen import WHOLE_CHROMOSOME
    return P.fit_population(freq, seen, r_c_nm=r_c, cfg=P.PopulationConfig(seed=split, **WHOLE_CHROMOSOME))


def build(spec: UnitSpec, force: bool = False) -> Path:
    """Fit the unit's model and store it with the held-out truth (cached; returns the cache path)."""
    path = CACHE / f"{spec.name}.npz"
    if path.exists() and not force:
        return path
    from benchmark.run import _hic_input
    from frozen import GATE1
    t0 = time.time()
    tr, cols, xyz = _traces(spec)
    base = spec.dataset.split(":")[0]
    e = D.REGISTRY[base]
    n = xyz.shape[1]
    starts = tr.starts if cols is None else tr.starts[cols]
    a, b = PR.split(len(xyz), spec.split)
    r_img = 150.0 if e.kind == "bintu_csv" else None
    hb = PR.half_stats(xyz[b], r_img)
    meta = {"spec": spec.__dict__, "dataset_role": e.role, "cell_line": e.cell_line, "loci": n,
            "copies": int(len(xyz)), "step_bp": int(e.step_bp if e.kind != "bintu_csv" else 30_000)}
    if spec.input == "imaging":
        ha = PR.half_stats(xyz[a], r_img)
        r_c = ha.adjacent_median if r_img is None else r_img
        freq, seen, counts = ha.freq, ha.seen, ha.freq * ha.seen
        n_eff = float(np.nanmedian(ha.seen[np.triu_indices(n, 1)]))
        meta.update(protocol="imaging")
    else:
        raw = _hic_counts(spec, tr, cols)
        raw = thin_counts(raw, spec.thin, seed=1000 + spec.split)
        freq, n_eff, counts = _hic_input(raw)
        seen = n_eff
        r_c = physics.bond_length_for(meta["step_bp"]) / float(E.gaussian_median_distance(float(GATE1["p_adjacent"]), 1.0))
        meta.update(protocol=PROTOCOL[spec.source], thin=spec.thin, reads_in_window=float(np.triu(raw, 1).sum()))
    res = _fit(freq, seen, r_c, n, spec.split)
    iu = np.triu_indices(n, 1)
    sigma = np.zeros((n, n), np.float32)
    sigma[iu] = P.pair_sigma_nm(res, iu[0], iu[1])
    sigma = sigma + sigma.T
    meta.update(model=res.config.get("model", "ensemble_v3_3"), device=res.config.get("device_used", "cpu"),
                contact_fit=float(res.contact_fit), n_eff=float(n_eff), r_c_nm=float(r_c),
                seconds=round(time.time() - t0, 1), utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, median=np.asarray(res.median_distance_nm, np.float32), sigma=sigma,
                        freq=np.asarray(freq, np.float32), counts=np.asarray(counts, np.float32),
                        seen=np.broadcast_to(np.asarray(seen, np.float32), (n, n)).copy(),
                        truth=np.asarray(hb.median, np.float32), sep=PR.separation(n, starts), starts=np.asarray(starts),
                        b_idx=np.asarray(b), meta=np.array(json.dumps(meta)))
    return path


def load(spec: UnitSpec) -> dict:
    z = np.load(build(spec), allow_pickle=False)
    out = {k: z[k] for k in z.files if k != "meta"}
    out["meta"] = json.loads(str(z["meta"]))
    return out


def single_copies(spec: UnitSpec) -> np.ndarray:
    """Half B's single-copy coordinates of the unit (the truth for interval coverage)."""
    _, _, xyz = _traces(spec)
    return xyz[load(spec)["b_idx"]]


def genome_units(key: str, input_kind: str, source: str = "rao2014", thin: float = 1.0, split: int = 0,
                 skip=("chr2", "chrY")) -> list[UnitSpec]:
    """One unit per chromosome of a genome-scale set (chr2 overlaps the practice chr2 traces; Y is too sparse)."""
    tr = _load_traces(key)
    chroms = [c for c in dict.fromkeys(np.asarray(tr.chrom).tolist()) if c not in skip]
    return [UnitSpec(f"{key}:{c}", input_kind, source, thin, split) for c in chroms]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", choices=("practice", "test"), required=True)
    a = ap.parse_args()
    import phase_a_units as U
    specs = U.PRACTICE if a.build == "practice" else U.TEST
    for s in specs:
        t0 = time.time()
        try:
            p = build(s)
            print(f"{s.name:70s} ok ({time.time() - t0:.0f} s) {p.name}", flush=True)
        except Exception as exc:                          # reported, never dropped silently
            print(f"{s.name:70s} FAILED {type(exc).__name__}: {exc}", flush=True)


if __name__ == "__main__":
    main()
