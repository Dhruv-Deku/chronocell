"""
Gate Q8b data split (October 2026): a scaffold split of eight TDC datasets that no ChronoCell test has used (19
endpoints; Tox21 has twelve assays), written to validation/admet_q8b_split.json before any model is fitted.

Runs in the isolated quantum environment (RDKit):

    .chronocell_cache/quantum-venv/Scripts/python validation/admet_q8b_split.py

Each molecule's Bemis-Murcko scaffold (RDKit MurckoScaffoldSmiles, chirality ignored; acyclic molecules share the
empty scaffold) defines a group; groups are split as TDC's create_scaffold_split does (groups larger than half the
test size first, then the rest, each list shuffled with the seed), filling train to 80 % and test with the rest.
Only the SMILES are read here; labels are not looked at. Tox21: one split over all its molecules, used by every
assay (each assay then keeps its labelled rows). Molecules RDKit cannot read are left out and counted.
"""

from __future__ import annotations

import csv
import hashlib
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "tdc" / "q8b"
OUT = ROOT / "admet_q8b_split.json"
SEED, TRAIN_FRAC = 42, 0.8
# file -> (Dataverse file id, MD5 published by Dataverse, SMILES column)
FILES = {"cyp1a2_veith": (4259573, "e5eeb84ca332cd059c73b816f7964193", "Drug"),
         "cyp2c19_veith": (4259576, "b4d3da8e8365f63fc920621b0be8ef32", "Drug"),
         "pampa_ncats": (6695858, "535f90234e93a5d83c258f221c1532be", "Drug"),
         "hydrationfreeenergy_freesolv": (4259594, "fbbad23184e31f8487c0a59790fa8d16", "Drug"),
         "skin_reaction": (4259609, "8435dbc2a8c509a40d5befa952c50c83", "Drug"),
         "carcinogens_lagunin": (4259570, "d523b34868ce3eb352607c44dbbae114", "Drug"),
         "clintox": (4259572, "4b0aa740b54f6096775c6fd88c625aea", "Drug"),
         "tox21": (4259612, "6f926279d60d413f0524894fdcb9ba5e", "X")}


def smiles_of(name: str) -> list[str]:
    fid, md5, col = FILES[name]
    p = DATA / f"{name}.csv"
    got = hashlib.md5(p.read_bytes()).hexdigest()
    if got != md5:
        sys.exit(f"{p}: MD5 {got} does not match Dataverse's {md5}")
    with p.open(encoding="utf-8", newline="") as f:
        return [r[col] for r in csv.DictReader(f)]


def scaffold_split(smiles: list[str]) -> tuple[list[int], list[int], int]:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    groups: dict[str, list[int]] = defaultdict(list)
    bad = 0
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            bad += 1
            continue
        groups[MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)].append(i)
    n = sum(len(g) for g in groups.values())
    n_train = int(round(TRAIN_FRAC * n))
    half_test = (n - n_train) / 2
    sets = sorted(groups.values(), key=lambda g: g[0])                  # deterministic before shuffling
    big = [g for g in sets if len(g) > half_test]
    small = [g for g in sets if len(g) <= half_test]
    rng = random.Random(SEED)
    rng.shuffle(big)
    rng.shuffle(small)
    train, test = [], []
    for g in big + small:
        (train if len(train) + len(g) <= n_train else test).extend(g)
    return sorted(train), sorted(test), bad


def main() -> None:
    out = {"seed": SEED, "train_fraction": TRAIN_FRAC, "files": {k: {"dataverse_id": v[0], "md5": v[1]} for k, v in FILES.items()},
           "splits": {}}
    for name in FILES:
        tr, te, bad = scaffold_split(smiles_of(name))
        out["splits"][name] = {"train": tr, "test": te, "unreadable": bad}
        print(f"{name:30s} train {len(tr):6d}  test {len(te):6d}  unreadable {bad}")
    OUT.write_text(json.dumps(out), encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
