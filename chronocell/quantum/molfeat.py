"""
Molecular descriptors from a SMILES string, in pure Python (no RDKit), for the quantum-kernel safety classifier.

A small SMILES reader (organic subset B C N O P S F Cl Br I, bracket atoms with isotopes, charges, explicit H,
aromatic lower case, bonds - = # : / \\, branches, ring closures incl. %nn, dot-separated fragments) builds the
heavy-atom graph and fills implicit hydrogens with the default valences. From the graph:

heavy_atoms, mw (average masses), n_N, n_O, n_S, n_halogen, n_aromatic_atoms, rings (cycle rank), aromatic_rings
(SSSR rings made only of aromatic atoms), hbd (N-H and O-H, Lipinski), hba (N + O, Lipinski), rot_bonds (single,
non-ring bonds between two heavy atoms each with another heavy neighbour, not C-N in an amide, not to a terminal
triple-bond atom), tpsa (Ertl et al. 2000 contributions for N and O; S and P ignored, as RDKit's default),
fsp3 (sp3 carbons / carbons), charge (formal), basic_n (an aliphatic amine N: not aromatic, no multiple bond, not
attached to C=O / C=S / S(=O) / aromatic ring), logp_proxy (a crude atom-count estimate; not Crippen).

The values are checked against RDKit (validation/quantum_drug_gates.py, in the isolated environment); the
agreement is reported in RESULTS.md. Small differences (unusual valences, aromaticity of unusual rings) are
expected and do not matter much for a coarse classifier.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

MASS = {"H": 1.008, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Na": 22.990, "Mg": 24.305,
        "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "K": 39.098, "Ca": 40.078, "Fe": 55.845, "Zn": 65.38,
        "Se": 78.971, "Br": 79.904, "I": 126.904, "Li": 6.94, "Pt": 195.08, "As": 74.922, "Sn": 118.71, "Hg": 200.59,
        "Co": 58.933, "Cu": 63.546, "Ga": 69.723, "Ge": 72.630, "Al": 26.982, "Te": 127.60, "Bi": 208.98, "Gd": 157.25,
        "Ti": 47.867, "Ag": 107.87, "Au": 196.97, "Mn": 54.938, "Cr": 51.996, "Ni": 58.693, "V": 50.942, "Ba": 137.33,
        "Sb": 121.76, "Tc": 98.0, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Cd": 112.41, "In": 114.82, "Cs": 132.91,
        "Sr": 87.62, "Rb": 85.468, "Zr": 91.224, "Mo": 95.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22,
        "Tl": 204.38, "Pb": 207.2, "La": 138.91, "Y": 88.906, "Sc": 44.956, "Nb": 92.906, "Ta": 180.95, "Hf": 178.49,
        "U": 238.03, "Th": 232.04, "Eu": 151.96, "Sm": 150.36, "Nd": 144.24, "Ce": 140.12, "Lu": 174.97, "Yb": 173.05,
        "Er": 167.26, "Ho": 164.93, "Dy": 162.50, "Tb": 158.93, "Pr": 140.91, "Tm": 168.93, "Be": 9.0122, "Ra": 226.0,
        "Xe": 131.29, "Kr": 83.798, "Ar": 39.948, "Ne": 20.180, "He": 4.0026, "At": 210.0, "Po": 209.0, "Fr": 223.0,
        "Ac": 227.0, "Pa": 231.04, "Np": 237.0, "Pu": 244.0, "Am": 243.0}
VALENCE = {"B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5), "S": (2, 4, 6), "F": (1,), "Cl": (1,),
           "Br": (1,), "I": (1,)}
ORGANIC = ("Cl", "Br", "B", "C", "N", "O", "P", "S", "F", "I", "b", "c", "n", "o", "p", "s")
BOND = {"-": 1.0, "=": 2.0, "#": 3.0, "$": 4.0, ":": 1.5, "/": 1.0, "\\": 1.0}
TOKEN = re.compile(r"(\[[^\]]+\]|Br|Cl|[BCNOPSFI]|[bcnops]|\(|\)|\.|=|#|\$|-|\+|/|\\|:|%\d\d|\d)")
BRACKET = re.compile(r"\[(\d*)([A-Z][a-z]?|[a-z][a-z]?|\*)(@*)(H\d*)?([+-]+\d*|[+-]\d+)?(:\d+)?\]")
FEATURES = ["heavy_atoms", "mw", "n_N", "n_O", "n_S", "n_halogen", "n_aromatic_atoms", "rings", "aromatic_rings", "hbd",
            "hba", "rot_bonds", "tpsa", "fsp3", "charge", "basic_n", "logp_proxy"]


class SmilesError(ValueError):
    pass


@dataclass
class Atom:
    el: str
    aromatic: bool
    charge: int = 0
    hx: int | None = None          # explicit H count (bracket atoms)
    h: int = 0                     # total H after filling
    nbrs: list = field(default_factory=list)   # (atom index, bond order)


@dataclass
class Mol:
    atoms: list
    smiles: str

    def bonds(self):
        for i, a in enumerate(self.atoms):
            for j, o in a.nbrs:
                if i < j:
                    yield i, j, o


def parse(smiles: str) -> Mol:
    s = smiles.strip().split()[0] if smiles.strip() else ""
    if not s:
        raise SmilesError("empty SMILES")
    toks = TOKEN.findall(s)
    if "".join(toks) != s:
        raise SmilesError(f"cannot read {smiles!r}")
    atoms: list[Atom] = []
    prev = None
    stack = []
    bond = None
    rings: dict = {}
    for t in toks:
        if t == "(":
            stack.append(prev)
        elif t == ")":
            prev = stack.pop()
        elif t == ".":
            prev = None
        elif t in BOND:
            bond = BOND[t]
        elif t[0] == "%" or t.isdigit():
            k = t
            if k in rings:
                j, b0 = rings.pop(k)
                o = bond or b0 or (1.5 if atoms[prev].aromatic and atoms[j].aromatic else 1.0)
                atoms[prev].nbrs.append((j, o))
                atoms[j].nbrs.append((prev, o))
            else:
                rings[k] = (prev, bond)
            bond = None
        else:
            if t.startswith("["):
                m = BRACKET.fullmatch(t)
                if not m:
                    raise SmilesError(f"bracket atom {t}")
                sym = m.group(2)
                arom = sym[0].islower()
                el = sym.capitalize() if arom else sym
                h = m.group(4)
                hx = 0 if h is None else (1 if h == "H" else int(h[1:]))
                ch = m.group(5) or ""
                charge = 0
                if ch:
                    sign = 1 if ch[0] == "+" else -1
                    rest = ch[1:]
                    charge = sign * (int(rest) if rest.isdigit() else len(ch))
                a = Atom(el, arom, charge, hx)
            else:
                arom = t[0].islower()
                a = Atom(t.capitalize() if arom else t, arom)
            atoms.append(a)
            i = len(atoms) - 1
            if prev is not None:
                o = bond or (1.5 if atoms[prev].aromatic and a.aromatic else 1.0)
                atoms[prev].nbrs.append((i, o))
                a.nbrs.append((prev, o))
            prev = i
            bond = None
    if rings:
        raise SmilesError("unclosed ring")
    for a in atoms:
        if a.hx is not None:
            a.h = a.hx
            continue
        if a.el not in VALENCE or (a.aromatic and a.el in ("O", "S")):
            a.h = 0
            continue
        used = sum(o for _, o in a.nbrs)
        if a.aromatic:
            used = sum(1 for _, o in a.nbrs if o == 1.5) + sum(o for _, o in a.nbrs if o != 1.5) + \
                (1 if any(o == 1.5 for _, o in a.nbrs) else 0)
            if a.el in ("N", "P") and sum(1 for _, o in a.nbrs) == 3:
                used = 3
        used = int(round(used))
        target = next((v for v in VALENCE[a.el] if v >= used), used)
        a.h = max(0, target - used)
    return Mol(atoms, s)


def _cycle_rank(m: Mol) -> int:
    n = len(m.atoms)
    e = sum(1 for _ in m.bonds())
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, j, _ in m.bonds():
        parent[find(i)] = find(j)
    comps = len({find(i) for i in range(n)})
    return e - n + comps


def _ring_bonds(m: Mol) -> set:
    """Bonds on some cycle (a bond is in a ring if removing it keeps its ends connected)."""
    adj = {i: [j for j, _ in a.nbrs] for i, a in enumerate(m.atoms)}
    out = set()
    for i, j, _ in m.bonds():
        seen = {i}
        stack = [i]
        found = False
        while stack and not found:
            x = stack.pop()
            for y in adj[x]:
                if (x == i and y == j) or (x == j and y == i):
                    continue
                if y == j:
                    found = True
                    break
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        if found:
            out.add((i, j))
    return out


def _sssr_aromatic_rings(m: Mol, max_size: int = 7) -> int:
    """Count smallest rings (size <= max_size) made only of aromatic atoms (approximate SSSR)."""
    arom = {i for i, a in enumerate(m.atoms) if a.aromatic}
    adj = {i: [j for j, _ in m.atoms[i].nbrs if j in arom] for i in arom}
    rings = set()
    for start in arom:
        # BFS for short cycles through start
        paths = [[start]]
        for _ in range(max_size):
            new = []
            for p in paths:
                for y in adj[p[-1]]:
                    if y == start and len(p) >= 3:
                        rings.add(frozenset(p))
                    elif y not in p and y > start:
                        new.append(p + [y])
            paths = new
    # keep rings that are not unions of smaller ones: greedy by size with independent edges
    chosen = []
    edges_used: list[set] = []
    for r in sorted(rings, key=len):
        es = set()
        rl = list(r)
        for a in rl:
            for b in adj[a]:
                if b in r:
                    es.add(frozenset((a, b)))
        if any(es <= u for u in edges_used):
            continue
        chosen.append(r)
        edges_used.append(es)
    return min(len(chosen), _cycle_rank(m))


def descriptors(smiles: str) -> dict:
    m = parse(smiles)
    A = m.atoms
    heavy = [a for a in A if a.el != "H"]
    mw = sum(MASS.get(a.el, 0.0) + a.h * MASS["H"] for a in A)
    n_N = sum(a.el == "N" for a in A)
    n_O = sum(a.el == "O" for a in A)
    n_S = sum(a.el == "S" for a in A)
    n_hal = sum(a.el in ("F", "Cl", "Br", "I") for a in A)
    n_ar = sum(a.aromatic for a in A)
    rings = _cycle_rank(m)
    hbd = sum(a.h for a in A if a.el in ("N", "O"))
    hba = n_N + n_O
    ring_b = _ring_bonds(m)

    def carbonyl_like(i):
        a = A[i]
        return a.el in ("C", "S") and any(A[j].el in ("O", "S") and o == 2.0 for j, o in a.nbrs)
    rot = 0
    for i, j, o in m.bonds():
        if o != 1.0 or (i, j) in ring_b:
            continue
        ai, aj = A[i], A[j]
        if len(ai.nbrs) < 2 or len(aj.nbrs) < 2:
            continue
        if any(oo == 3.0 for _, oo in ai.nbrs) or any(oo == 3.0 for _, oo in aj.nbrs):
            continue
        if (ai.el == "N" and carbonyl_like(j)) or (aj.el == "N" and carbonyl_like(i)):
            continue
        rot += 1
    tpsa = 0.0
    for a in A:
        nh = a.h
        heavy_n = len(a.nbrs)
        dbl = sum(1 for _, o in a.nbrs if o == 2.0)
        tri = sum(1 for _, o in a.nbrs if o == 3.0)
        if a.el == "N":
            if a.aromatic:
                tpsa += {0: 12.89, 1: 15.79}.get(nh, 12.89) if a.charge == 0 else 4.10
            elif a.charge > 0:
                tpsa += {0: 0.0, 1: 4.44, 2: 8.38, 3: 13.97}.get(nh, 3.01 if dbl else 0.0)
            elif tri:
                tpsa += 23.79
            elif dbl:
                tpsa += {0: 12.36 if heavy_n == 2 else 11.68, 1: 23.85}.get(nh, 12.36)
            else:
                tpsa += {0: 3.24, 1: 12.03, 2: 26.02}.get(nh, 26.02)
        elif a.el == "O":
            if a.aromatic:
                tpsa += 13.14
            elif a.charge < 0:
                tpsa += 23.06
            elif dbl:
                tpsa += 17.07
            else:
                tpsa += {0: 9.23, 1: 20.23}.get(nh, 20.23)
    carbons = [a for a in A if a.el == "C"]
    fsp3 = (sum(1 for a in carbons if not a.aromatic and all(o == 1.0 for _, o in a.nbrs)) / len(carbons)) if carbons else 0.0
    charge = sum(a.charge for a in A)

    def basic(i):
        a = A[i]
        if a.el != "N" or a.aromatic or any(o != 1.0 for _, o in a.nbrs):
            return False
        for j, _ in a.nbrs:
            b = A[j]
            if b.aromatic or carbonyl_like(j) or b.el in ("S", "N", "O") or any(oo >= 2.0 for _, oo in b.nbrs):
                return False
        return True
    basic_n = sum(basic(i) for i in range(len(A)))
    logp = 0.5 * len([a for a in carbons]) + 0.6 * n_hal - 1.0 * hbd - 0.6 * (n_N + n_O - 0) + 0.3 * n_ar / 6.0
    return {"heavy_atoms": len(heavy), "mw": mw, "n_N": n_N, "n_O": n_O, "n_S": n_S, "n_halogen": n_hal,
            "n_aromatic_atoms": n_ar, "rings": rings, "aromatic_rings": _sssr_aromatic_rings(m), "hbd": hbd, "hba": hba,
            "rot_bonds": rot, "tpsa": tpsa, "fsp3": fsp3, "charge": charge, "basic_n": basic_n, "logp_proxy": logp}


def matrix(smiles_list: list[str]) -> tuple[np.ndarray, list[int]]:
    """Descriptor matrix (rows for readable SMILES) and the indices that were read."""
    rows, ok = [], []
    for k, s in enumerate(smiles_list):
        try:
            d = descriptors(s)
        except (SmilesError, KeyError, IndexError):
            continue
        rows.append([d[f] for f in FEATURES])
        ok.append(k)
    return np.array(rows, float), ok


def lipinski(d: dict) -> dict:
    """Rule of five (Lipinski 2001) and Veber's oral-bioavailability rules, from the descriptors."""
    ro5 = {"MW <= 500": d["mw"] <= 500, "H-bond donors <= 5": d["hbd"] <= 5, "H-bond acceptors (N+O) <= 10": d["hba"] <= 10}
    veber = {"rotatable bonds <= 10": d["rot_bonds"] <= 10, "TPSA <= 140": d["tpsa"] <= 140}
    return {"rule_of_five": ro5, "veber": veber, "ro5_violations": sum(not v for v in ro5.values())}
