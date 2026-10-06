"""
Rigid re-docking as a maximum-weight clique (the formulation of Kuhl et al. 1984 and of the quantum docking
study of Banchi et al., Sci Adv 2020), solved by QAOA on the simulator and exactly.

1. Pocket hot-spots, from the protein only: grid points (1 A) within `box` A of the site centre that sit 2.6-4.6 A
   from protein heavy atoms and are buried; each typed by what a ligand atom there would do:
   A (accept an H-bond: 1.6-2.6 A from a protein polar hydrogen), D (donate one: 2.5-3.3 A from a protein O, or an N
   carrying no H), H (hydrophobic: no protein N/O within 3.6 A, >= 4 protein carbons within 4.5 A). The best point
   per neighbourhood is kept (non-maximum suppression), a few per type.
2. Ligand features, from its 2D graph and one 3D conformer: acceptors (O; N without H that is aromatic or triple
   bonded), donors (N or O with H), hydrophobic carbons (no N/O neighbour), thinned by farthest-point sampling.
3. Binding-interaction graph: a vertex per type-compatible (feature, hot-spot) pair, weighted by the hot-spot's
   burial; an edge when two vertices use different features and different hot-spots whose distances agree within
   tau. A clique is a set of matches that one rigid placement can satisfy at once.
4. QUBO for the maximum-weight clique: minimise -sum w_v x_v + P sum_{non-edges} x_u x_v (P > max weight).
5. Pose: the ligand's features are superposed on their matched hot-spots (Kabsch); several good cliques are turned
   into poses and the pose with the best simple score (matches minus clashes) is kept.

The site centre comes from the crystal ligand (site-specific docking, as a docking box is placed in practice), and
the ligand is the crystal conformer (rigid re-docking), so this measures placement, not conformer search.
Success: heavy-atom RMSD to the crystal pose <= 2 A (no symmetry correction).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .qubo import QUBO

POLAR = ("N", "O")


# ======================================================================================
# Reading
# ======================================================================================
@dataclass
class Ligand:
    name: str
    el: list
    xyz: np.ndarray
    bonds: list                 # (i, j, order)
    h: np.ndarray = None        # implicit + explicit H per heavy atom

    @property
    def heavy(self) -> np.ndarray:
        return np.array([e != "H" for e in self.el])


def read_sdf(text: str, name: str = "ligand") -> Ligand:
    lines = text.splitlines()
    na, nb = int(lines[3][:3]), int(lines[3][3:6])
    el, xyz = [], []
    for ln in lines[4:4 + na]:
        xyz.append([float(ln[0:10]), float(ln[10:20]), float(ln[20:30])])
        el.append(ln[31:34].strip())
    bonds = []
    for ln in lines[4 + na:4 + na + nb]:
        bonds.append((int(ln[0:3]) - 1, int(ln[3:6]) - 1, int(ln[6:9])))
    lig = Ligand(name, el, np.array(xyz), bonds)
    lig.h = _implicit_h(lig)
    return lig


VAL = {"C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "B": 3}


def _implicit_h(lig: Ligand) -> np.ndarray:
    n = len(lig.el)
    used = np.zeros(n)
    explicit = np.zeros(n, int)
    for i, j, o in lig.bonds:
        oo = 1.5 if o == 4 else o
        used[i] += oo
        used[j] += oo
        if lig.el[i] == "H":
            explicit[j] += 1
        if lig.el[j] == "H":
            explicit[i] += 1
    h = np.zeros(n, int)
    for k, e in enumerate(lig.el):
        if e in VAL:
            v = VAL[e]
            if e in ("S", "P") and used[k] > v:
                v = 6 if e == "S" else 5
            h[k] = max(0, int(round(v - used[k])))
        h[k] = max(h[k], 0)
    # explicit hydrogens already counted in `used` when present
    return h


@dataclass
class Protein:
    el: np.ndarray
    xyz: np.ndarray
    name: np.ndarray
    res: np.ndarray


def read_pdb(text: str) -> Protein:
    el, xyz, nm, rs = [], [], [], []
    for ln in text.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            try:
                x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
            except ValueError:
                continue
            e = ln[76:78].strip() or ln[12:16].strip()[0]
            if ln[17:20] in ("HOH", "WAT"):
                continue
            el.append(e.capitalize())
            xyz.append((x, y, z))
            nm.append(ln[12:16].strip())
            rs.append(ln[17:20].strip())
    return Protein(np.array(el), np.array(xyz), np.array(nm), np.array(rs))


# ======================================================================================
# Features and hot-spots
# ======================================================================================
@dataclass
class Points:
    xyz: np.ndarray
    kind: list                 # "A" | "D" | "H"
    weight: np.ndarray
    atom: list = field(default_factory=list)   # ligand atom index (ligand features)


def _fps(xyz: np.ndarray, k: int) -> list[int]:
    if len(xyz) <= k:
        return list(range(len(xyz)))
    chosen = [int(np.argmax(np.linalg.norm(xyz - xyz.mean(0), axis=1)))]
    d = np.linalg.norm(xyz - xyz[chosen[0]], axis=1)
    while len(chosen) < k:
        nxt = int(np.argmax(d))
        chosen.append(nxt)
        d = np.minimum(d, np.linalg.norm(xyz - xyz[nxt], axis=1))
    return chosen


def ligand_features(lig: Ligand, max_polar: int = 5, max_hydrophobic: int = 3) -> Points:
    nb = {i: [] for i in range(len(lig.el))}
    for i, j, o in lig.bonds:
        nb[i].append((j, o))
        nb[j].append((i, o))
    acc, don, hyd = [], [], []
    for k, e in enumerate(lig.el):
        if e == "O":
            acc.append(k)
            if lig.h[k] > 0:
                don.append(k)
        elif e == "N":
            if lig.h[k] > 0:
                don.append(k)
            elif any(o in (3, 4) for _, o in nb[k]) or (len(nb[k]) == 2 and any(o == 2 for _, o in nb[k])):
                acc.append(k)
        elif e == "C" and not any(lig.el[j] in POLAR for j, _ in nb[k]):
            hyd.append(k)
    pts, kinds, idx = [], [], []
    polar = [(k, "A") for k in acc] + [(k, "D") for k in don]
    if polar:
        sel = _fps(np.array([lig.xyz[k] for k, _ in polar]), max_polar)
        for s in sel:
            pts.append(lig.xyz[polar[s][0]])
            kinds.append(polar[s][1])
            idx.append(polar[s][0])
    if hyd:
        sel = _fps(lig.xyz[hyd], max_hydrophobic)
        for s in sel:
            pts.append(lig.xyz[hyd[s]])
            kinds.append("H")
            idx.append(hyd[s])
    w = np.array([1.0 if k != "H" else 0.6 for k in kinds])
    return Points(np.array(pts), kinds, w, idx)


def hotspots(prot: Protein, centre: np.ndarray, box: float = 8.0, spacing: float = 1.0, per_type: int = 4) -> Points:
    from scipy.spatial import cKDTree
    heavy = prot.el != "H"
    hx = prot.xyz[heavy]
    he = prot.el[heavy]
    tree = cKDTree(hx)
    g = np.arange(-box, box + 1e-9, spacing)
    grid = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3) + centre
    d, _ = tree.query(grid)
    grid = grid[(d >= 2.6) & (d <= 4.6)]
    if not len(grid):
        return Points(np.zeros((0, 3)), [], np.zeros(0))
    bur = np.array([len(x) for x in tree.query_ball_point(grid, 8.0)], float)
    keep = bur >= np.percentile(bur, 40)
    grid, bur = grid[keep], bur[keep]
    # polar hydrogens: H within 1.15 A of an N or O
    hmask = prot.el == "H"
    polar_heavy = np.isin(prot.el, POLAR)
    ptree = cKDTree(prot.xyz[polar_heavy])
    hxyz = prot.xyz[hmask]
    dpol, _ = ptree.query(hxyz) if len(hxyz) else (np.zeros(0), None)
    polarH = hxyz[dpol < 1.15] if len(hxyz) else np.zeros((0, 3))
    # acceptor atoms: O, or N with no H within 1.15 A
    acc_atoms = []
    htree = cKDTree(hxyz) if len(hxyz) else None
    for k in np.flatnonzero(polar_heavy):
        if prot.el[k] == "O":
            acc_atoms.append(prot.xyz[k])
        elif htree is not None and not htree.query_ball_point(prot.xyz[k], 1.15):
            acc_atoms.append(prot.xyz[k])
    acc_atoms = np.array(acc_atoms) if acc_atoms else np.zeros((0, 3))
    kinds = np.array([""] * len(grid), dtype=object)
    if len(polarH):
        dh, _ = cKDTree(polarH).query(grid)
        kinds[(dh >= 1.6) & (dh <= 2.6)] = "A"
    if len(acc_atoms):
        da, _ = cKDTree(acc_atoms).query(grid)
        kinds[(kinds == "") & (da >= 2.5) & (da <= 3.3)] = "D"
    polar_tree = cKDTree(hx[np.isin(he, POLAR)]) if np.isin(he, POLAR).any() else None
    carb_tree = cKDTree(hx[he == "C"])
    dp = polar_tree.query(grid)[0] if polar_tree is not None else np.full(len(grid), 99.0)
    nc = np.array([len(x) for x in carb_tree.query_ball_point(grid, 4.5)])
    kinds[(kinds == "") & (dp > 3.6) & (nc >= 4)] = "H"
    score = (bur - bur.min()) / (np.ptp(bur) or 1.0)
    pts, kk, ww = [], [], []
    for t, cap in (("A", per_type), ("D", per_type), ("H", per_type + 2)):
        idx = np.flatnonzero(kinds == t)
        order = idx[np.argsort(-score[idx])]
        chosen = []
        for i in order:
            if all(np.linalg.norm(grid[i] - grid[j]) >= 2.0 for j in chosen):
                chosen.append(i)
            if len(chosen) >= cap:
                break
        for i in chosen:
            pts.append(grid[i])
            kk.append(t)
            ww.append(0.5 + score[i])
    return Points(np.array(pts) if pts else np.zeros((0, 3)), kk, np.array(ww))


# ======================================================================================
# Graph, QUBO, pose
# ======================================================================================
COMPATIBLE = {"A": "A", "D": "D", "H": "H"}


@dataclass
class Graph:
    vertices: list             # (feature index, hotspot index)
    weight: np.ndarray
    adj: np.ndarray            # bool


def interaction_graph(lf: Points, hs: Points, tau: float = 1.5, max_vertices: int = 20) -> Graph:
    verts, w = [], []
    for i, ki in enumerate(lf.kind):
        for s, ks in enumerate(hs.kind):
            if COMPATIBLE[ki] == ks:
                verts.append((i, s))
                w.append(lf.weight[i] * hs.weight[s])
    w = np.array(w)
    if len(verts) > max_vertices:                  # keep the strongest, at most 4 hot-spots per feature
        order = np.argsort(-w)
        per = {}
        keep = []
        for v in order:
            f = verts[v][0]
            if per.get(f, 0) < 4:
                keep.append(v)
                per[f] = per.get(f, 0) + 1
            if len(keep) >= max_vertices:
                break
        keep = sorted(keep)
        verts = [verts[v] for v in keep]
        w = w[keep]
    n = len(verts)
    adj = np.zeros((n, n), bool)
    for a in range(n):
        for b in range(a + 1, n):
            (i, s), (j, t) = verts[a], verts[b]
            if i == j or s == t:
                continue
            dl = np.linalg.norm(lf.xyz[i] - lf.xyz[j])
            dp = np.linalg.norm(hs.xyz[s] - hs.xyz[t])
            if abs(dl - dp) <= tau:
                adj[a, b] = adj[b, a] = True
    return Graph(verts, w, adj)


def clique_qubo(g: Graph, penalty: float | None = None) -> QUBO:
    n = len(g.vertices)
    P = penalty if penalty is not None else 2.0 * float(g.weight.max(initial=1.0))
    quad = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            if not g.adj[a, b]:
                quad[a, b] = P
    return QUBO(-g.weight, quad, 0.0, [f"f{i}-s{s}" for i, s in g.vertices], "Max-weight clique")


def is_clique(g: Graph, bits: np.ndarray) -> bool:
    on = np.flatnonzero(np.asarray(bits) > 0)
    return all(g.adj[a, b] for k, a in enumerate(on) for b in on[k + 1:])


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t minimising |R P + t - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def pose_from_clique(lig_xyz: np.ndarray, lf: Points, hs: Points, g: Graph, bits: np.ndarray) -> np.ndarray | None:
    on = np.flatnonzero(np.asarray(bits) > 0)
    if len(on) < 3:
        return None
    P = np.array([lf.xyz[g.vertices[v][0]] for v in on])
    Q = np.array([hs.xyz[g.vertices[v][1]] for v in on])
    if np.linalg.matrix_rank(P - P.mean(0), tol=0.5) < 2:
        return None
    R, t = kabsch(P, Q)
    return lig_xyz @ R.T + t


def pose_score(pose: np.ndarray, prot: Protein, tree=None) -> float:
    """Simple score: protein heavy atoms 3.3-5 A from ligand atoms (contacts) minus 10 x clashes (< 2.6 A)."""
    from scipy.spatial import cKDTree
    if tree is None:
        tree = cKDTree(prot.xyz[prot.el != "H"])
    d, _ = tree.query(pose)
    clashes = int(np.sum(d < 2.6))
    contacts = sum(len(x) for x in tree.query_ball_point(pose, 5.0)) - sum(len(x) for x in tree.query_ball_point(pose, 3.3))
    return float(contacts - 10 * clashes)


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                     [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
                     [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d]])


def random_search(lig_xyz: np.ndarray, centre: np.ndarray, prot: Protein, n: int, rng: np.random.Generator,
                  shift: float = 2.0) -> np.ndarray:
    """Classical baseline with the same score: n random rigid placements around the site centre, best kept."""
    from scipy.spatial import cKDTree
    tree = cKDTree(prot.xyz[prot.el != "H"])
    c0 = lig_xyz.mean(0)
    best, best_s = None, -np.inf
    for _ in range(n):
        R = random_rotation(rng)
        pose = (lig_xyz - c0) @ R.T + centre + rng.normal(0, shift, 3)
        s = pose_score(pose, prot, tree)
        if s > best_s:
            best, best_s = pose, s
    return best


def hotspots_directional(prot: Protein, centre: np.ndarray, box: float = 9.0, per_type: int = 6) -> Points:
    """Pharmacophore site points placed along protein H-bond vectors: a ligand acceptor 1.9 A beyond each polar
    hydrogen (on the X-H axis), a ligand donor 2.8 A beyond each acceptor O (opposite its heavy neighbours);
    hydrophobic sites from buried, carbon-lined grid points (as `hotspots`). Sites must lie within `box` of the
    centre, not clash with protein heavy atoms (>= 2.4 A) and be buried."""
    from scipy.spatial import cKDTree
    heavy = prot.el != "H"
    hx, he = prot.xyz[heavy], prot.el[heavy]
    tree = cKDTree(hx)
    near = np.flatnonzero(np.linalg.norm(prot.xyz - centre, axis=1) <= box + 3.0)
    sites, kinds = [], []
    for k in near:
        e = prot.el[k]
        if e == "H":
            d, j = tree.query(prot.xyz[k])
            if d < 1.15 and he[j] in POLAR:
                u = prot.xyz[k] - hx[j]
                sites.append(prot.xyz[k] + 1.9 * u / np.linalg.norm(u))
                kinds.append("A")
        elif e == "O":
            nb = [j for j in tree.query_ball_point(prot.xyz[k], 1.65) if np.linalg.norm(hx[j] - prot.xyz[k]) > 0.1]
            if not nb:
                continue
            u = sum((prot.xyz[k] - hx[j]) / np.linalg.norm(prot.xyz[k] - hx[j]) for j in nb)
            if np.linalg.norm(u) < 1e-6:
                continue
            sites.append(prot.xyz[k] + 2.8 * u / np.linalg.norm(u))
            kinds.append("D")
    sites = np.array(sites) if sites else np.zeros((0, 3))
    kinds = np.array(kinds, dtype=object)
    if len(sites):
        ok = (np.linalg.norm(sites - centre, axis=1) <= box) & (tree.query(sites)[0] >= 2.4)
        sites, kinds = sites[ok], kinds[ok]
    grid_h = hotspots(prot, centre, box=box, per_type=per_type)
    hmask = np.array([k == "H" for k in grid_h.kind], bool)
    pts = list(sites) + list(grid_h.xyz[hmask]) if len(grid_h.kind) else list(sites)
    kk = list(kinds) + ["H"] * int(hmask.sum())
    if not pts:
        return Points(np.zeros((0, 3)), [], np.zeros(0))
    pts = np.array(pts)
    bur = np.array([len(x) for x in tree.query_ball_point(pts, 8.0)], float)
    score = (bur - bur.min()) / (np.ptp(bur) or 1.0)
    out_p, out_k, out_w = [], [], []
    for t in ("A", "D", "H"):
        idx = [i for i in range(len(pts)) if kk[i] == t]
        idx.sort(key=lambda i: -score[i])
        chosen = []
        for i in idx:
            if all(np.linalg.norm(pts[i] - pts[j]) >= 1.5 for j in chosen):
                chosen.append(i)
            if len(chosen) >= per_type:
                break
        for i in chosen:
            out_p.append(pts[i])
            out_k.append(t)
            out_w.append(0.5 + score[i])
    return Points(np.array(out_p), out_k, np.array(out_w))


def core_subgraph(g: Graph, k: int = 20) -> tuple[Graph, list[int]]:
    """The k vertices most connected to the rest (weighted degree): where large cliques live. Used to fit the
    interaction graph into a qubit budget."""
    if len(g.vertices) <= k:
        return g, list(range(len(g.vertices)))
    deg = (g.adj * g.weight[None, :]).sum(1) + g.weight
    keep = sorted(np.argsort(-deg)[:k].tolist())
    return Graph([g.vertices[i] for i in keep], g.weight[keep], g.adj[np.ix_(keep, keep)]), keep


# ======================================================================================
# App side: PoseBusters examples on demand, one docking run with everything the panel shows
# ======================================================================================
PB_URL = "https://zenodo.org/api/records/13851241/files/{}/content"
PB_FILES = {"ligands": ("posebusters_ligands_256.zip", "fdeac0fa6ff4aeb46ce662f3ff8753dd"),
            "proteins": ("posebusters_spruce_structures_256.zip", "faef0d4a67d05dbfddbeba61fee5df6e")}
DEFAULT = {"tau": 1.0, "max_polar": 8, "max_hydrophobic": 4, "per_type": 8, "qubits": 20, "cliques": 30,
           "random_poses": 1000}


def posebusters(kind: str, progress=None):
    import zipfile
    from pathlib import Path
    from ..predict import download_verified
    name, md5 = PB_FILES[kind]
    p = Path(__file__).resolve().parent.parent.parent / ".chronocell_cache" / "posebusters" / name
    if not p.exists():
        download_verified(PB_URL.format(name), p, md5, progress)
    return zipfile.ZipFile(p)


def list_complexes(zl) -> list[str]:
    return sorted(n.split("/")[-1].replace("_ligand.sdf", "") for n in zl.namelist()
                  if n.startswith("posebusters_ligands_256/") and n.endswith("_ligand.sdf"))


@dataclass
class DockRun:
    cid: str
    lig: Ligand
    prot: Protein
    features: Points
    sites: Points
    graph: Graph
    qubo: QUBO
    exact_bits: np.ndarray
    qaoa: object
    pose_qaoa: np.ndarray | None
    pose_classical: np.ndarray | None
    pose_random: np.ndarray
    rmsd_qaoa: float
    rmsd_classical: float
    rmsd_random: float
    usable: bool = True


def dock(zl, zp, cid: str, cfg: dict | None = None, p: int = 3, objective: str = "expectation", seed: int = 0) -> DockRun:
    from scipy.spatial import cKDTree
    from . import sim, qubo as QB
    cfg = {**DEFAULT, **(cfg or {})}
    lig = read_sdf(zl.read(f"posebusters_ligands_256/{cid}_ligand.sdf").decode(), cid)
    prot = read_pdb(zp.read(f"posebusters_spruce_structures_256/{cid}_spruced.pdb").decode(errors="replace"))
    tree = cKDTree(prot.xyz[prot.el != "H"])
    heavy = lig.heavy
    L = lig.xyz[heavy]
    usable = tree.query(L)[0].min() < 4.0
    centre = L.mean(0)
    lf = ligand_features(lig, cfg["max_polar"], cfg["max_hydrophobic"])
    hs = hotspots_directional(prot, centre, per_type=cfg["per_type"])
    full = interaction_graph(lf, hs, tau=cfg["tau"], max_vertices=80)
    g, _ = core_subgraph(full, cfg["qubits"])
    q = clique_qubo(g)
    e = q.energies()
    ex = QB.exact(q, e)
    r = sim.qaoa(e, p=p, shots=4096, objective=objective, maxiter=80, seed=seed)

    def best(cliques):
        poses = []
        for _, c in cliques:
            bits = np.zeros(len(g.vertices))
            bits[c] = 1
            pz = pose_from_clique(lig.xyz, lf, hs, g, bits)
            if pz is not None:
                poses.append((pose_score(pz[heavy], prot, tree), pz))
        return max(poses, key=lambda z: z[0])[1] if poses else None

    seen, qc = set(), []
    for idx in list(np.argsort(-r.probs)[:20000]) + list(r.samples):
        bits = QB.to_bits(int(idx), q.n)
        on = tuple(np.flatnonzero(bits).tolist())
        if len(on) >= 3 and on not in seen and is_clique(g, bits):
            seen.add(on)
            qc.append((float(g.weight[list(on)].sum()), list(on)))
    qc = sorted(qc, reverse=True)[:cfg["cliques"]]
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(range(q.n))
    for a, b in zip(*np.nonzero(np.triu(g.adj, 1))):
        G.add_edge(int(a), int(b))
    cc = sorted(((float(g.weight[c].sum()), sorted(c)) for c in nx.find_cliques(G) if len(c) >= 3), reverse=True)[:cfg["cliques"]]
    pq, pc = best(qc), best(cc)
    prand = random_search(L, centre, prot, cfg["random_poses"], np.random.default_rng(seed))
    rq = rmsd(pq[heavy], L) if pq is not None else float("nan")
    rc = rmsd(pc[heavy], L) if pc is not None else float("nan")
    return DockRun(cid, lig, prot, lf, hs, g, q, ex.bits, r, pq, pc, prand, rq, rc, rmsd(prand, L), bool(usable))
