"""
Information about the Drug lab's drug classes (general, educational; not medical advice) and on-demand molecule
data from PubChem.

CLASSES holds, per drug class: the protein it acts on, what that protein does to chromatin, regulatory status as of
2025 (check current labels), example compounds with their status, safety themes noted on labels or in trials, and
the compound names used to look up structures. Statements are kept general on purpose; dates are first approvals by
the US FDA unless stated otherwise.

PubChem (https://pubchem.ncbi.nlm.nih.gov, PUG REST) is queried only when the user asks for a molecule: computed
properties, the 2D depiction (PNG) and a 3D conformer (SDF). Responses are cached under .chronocell_cache/pubchem/
with their SHA-256 recorded (PubChem publishes no checksums). Nothing about the user's data is sent: only the
compound name typed or picked.
"""

from __future__ import annotations

import hashlib
import json
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

CACHE = Path(__file__).resolve().parent.parent / ".chronocell_cache" / "pubchem"
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
UA = {"User-Agent": "ChronoCell-5D (research software; contact via repository)"}
PROPS = ("MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,"
         "HeavyAtomCount,Complexity,IsomericSMILES,IUPACName,Charge")


@dataclass(frozen=True)
class ClassInfo:
    key: str
    target_protein: str
    what_it_does: str          # the protein's normal job on chromatin
    status: str
    examples: tuple            # (compound, status)
    safety: str                # recurring safety themes (general)
    lookup: tuple              # PubChem names, first is the default


CLASSES: dict[str, ClassInfo] = {c.key: c for c in (
    ClassInfo("ezh2", "EZH2 (with EED and SUZ12: Polycomb repressive complex 2)",
              "Writes H3K27me3, the Polycomb mark that keeps developmental genes silent and compacted.",
              "Approved: tazemetostat (2020; epithelioid sarcoma, follicular lymphoma). Valemetostat (EZH1/2) approved "
              "in Japan (2022, adult T-cell leukaemia-lymphoma). EED inhibitors in trials.",
              (("tazemetostat", "approved"), ("valemetostat", "approved in Japan"), ("GSK126", "research tool"),
               ("EPZ-6438", "= tazemetostat")),
              "Labels note secondary malignancies; cytopenias.",
              ("tazemetostat", "GSK126", "valemetostat")),
    ClassInfo("hdac", "Histone deacetylases (class I/II HDACs)",
              "Remove acetyl groups from histones, closing chromatin and quieting genes.",
              "Approved: vorinostat (2006) and romidepsin (2009) for cutaneous T-cell lymphoma; belinostat (2014) for "
              "peripheral T-cell lymphoma; panobinostat (2015, myeloma; later withdrawn in the US); tucidinostat "
              "(China).",
              (("vorinostat", "approved"), ("romidepsin", "approved"), ("belinostat", "approved"),
               ("panobinostat", "approved, later withdrawn in the US"), ("trichostatin A", "research tool"),
               ("sodium butyrate", "research tool")),
              "Cytopenias, fatigue, gastrointestinal effects; ECG / QT changes are noted for several members.",
              ("vorinostat", "romidepsin", "panobinostat", "trichostatin A")),
    ClassInfo("bet", "BET bromodomain proteins (BRD2, BRD3, BRD4, BRDT)",
              "Read acetylated histones and recruit transcription machinery to enhancers and super-enhancers.",
              "No approved BET inhibitor; several in trials (e.g. pelabresib in myelofibrosis).",
              (("JQ1", "research tool"), ("pelabresib", "clinical trials"), ("I-BET762 (molibresib)", "clinical trials")),
              "Thrombocytopenia and gastrointestinal effects are the recurring class findings in trials.",
              ("JQ1", "pelabresib", "molibresib")),
    ClassInfo("ctcf", "CTCF and cohesin (loop extrusion)",
              "Cohesin extrudes loops until it meets CTCF, fencing the genome into neighbourhoods.",
              "No drug stabilises loops; the class is hypothetical (it mimics stronger cohesin retention, e.g. WAPL loss "
              "in research).",
              (("none", "hypothetical"),),
              "Not applicable (no compound).",
              ()),
    ClassInfo("dnmt", "DNA methyltransferases (DNMT1, DNMT3A/B)",
              "Copy and write DNA methylation, which keeps genes and repeats silenced.",
              "Approved: azacitidine (2004) and decitabine (2006) for myelodysplastic syndromes; oral decitabine-"
              "cedazuridine and oral azacitidine (2020).",
              (("azacitidine", "approved"), ("decitabine", "approved"), ("5-aza-2'-deoxycytidine", "= decitabine")),
              "Myelosuppression (low blood counts) is the main dose-limiting effect.",
              ("decitabine", "azacitidine")),
    ClassInfo("kdm", "JmjC-domain histone demethylases (and other 2-oxoglutarate-dependent enzymes)",
              "Remove methyl marks from histones, including repressive H3K9 and H3K27 methylation.",
              "No approved chromatin drug; DMOG is a research tool. Related 2-oxoglutarate-site inhibitors of HIF prolyl "
              "hydroxylases (e.g. daprodustat, roxadustat) are used for anaemia of chronic kidney disease.",
              (("dimethyloxalylglycine (DMOG)", "research tool"),),
              "Research compound; broad 2-oxoglutarate-enzyme inhibition also stabilises HIF.",
              ("dimethyloxalylglycine",)),
    ClassInfo("lsd1", "LSD1 (KDM1A)",
              "Removes H3K4me1/2 at enhancers, decommissioning them (part of CoREST complexes).",
              "No approved LSD1 drug for cancer; iadademstat and bomedemstat in trials; tranylcypromine (an approved "
              "antidepressant) is a weak LSD1 inhibitor.",
              (("iadademstat", "clinical trials"), ("bomedemstat", "clinical trials"), ("tranylcypromine", "approved (depression)")),
              "Thrombocytopenia is a recurring finding in trials.",
              ("iadademstat", "tranylcypromine")),
    ClassInfo("dot1l", "DOT1L",
              "Writes H3K79 methylation on transcribed genes; MLL-fusion leukaemias depend on it.",
              "No approved DOT1L drug; pinometostat reached early clinical trials with limited activity.",
              (("pinometostat", "clinical trials"),),
              "Early-phase data only.",
              ("pinometostat",)),
    ClassInfo("menin", "Menin (MEN1) bound to KMT2A (MLL)",
              "Tethers KMT2A complexes to target genes such as HOXA9 and MEIS1 that sustain leukaemic growth.",
              "Approved: revumenib (2024, relapsed or refractory acute leukaemia with a KMT2A translocation); other "
              "menin inhibitors in trials.",
              (("revumenib", "approved"),),
              "Differentiation syndrome (boxed warning on the revumenib label); QT prolongation monitored.",
              ("revumenib",)),
    ClassInfo("p300", "p300 / CBP (histone acetyltransferases)",
              "Write H3K27ac at enhancers and co-activate transcription factors.",
              "No approved drug; A-485 (catalytic inhibitor) is a research tool; inobrodib (bromodomain) in trials.",
              (("A-485", "research tool"), ("inobrodib", "clinical trials"), ("SGC-CBP30", "research tool")),
              "Early-phase data only.",
              ("A-485", "inobrodib")),
    ClassInfo("bet_degrader", "BRD4 (degraded through an E3 ligase: PROTAC)",
              "Same reader as BET inhibitors, but the drug marks BRD4 for destruction rather than blocking it.",
              "Research tools; degrader drugs for other targets are in trials.",
              (("dBET6", "research tool"), ("ARV-771", "research tool"), ("MZ1", "research tool")),
              "Research compounds.",
              ("dBET6", "ARV-771", "MZ1")),
    ClassInfo("txn", "RNA polymerase II and its kinase CDK9",
              "Transcribes genes; ongoing transcription helps keep active chromatin open and organised.",
              "Approved: dactinomycin (actinomycin D, 1964; paediatric and other cancers). Flavopiridol (alvocidib) and "
              "triptolide derivatives in trials; alpha-amanitin is a research toxin (amanitin antibody-drug conjugates "
              "in early trials).",
              (("dactinomycin", "approved"), ("alpha-amanitin", "research tool"), ("flavopiridol", "clinical trials"),
               ("triptolide", "research tool")),
              "Dactinomycin: myelosuppression, mucositis, tissue damage on extravasation.",
              ("dactinomycin", "flavopiridol", "alpha-amanitin")),
)}


# ======================================================================================
# PubChem (on demand)
# ======================================================================================
class PubChemError(RuntimeError):
    pass


def _get(url: str, binary: bool = False, timeout: int = 30):
    CACHE.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(url.encode()).hexdigest()[:24]
    path = CACHE / (key + (".bin" if binary else ".txt"))
    if path.exists():
        data = path.read_bytes()
    else:
        try:
            req = urllib.request.Request(url, headers=UA)
            data = urllib.request.urlopen(req, timeout=timeout).read()
        except Exception as exc:              # network or 404
            raise PubChemError(f"PubChem request failed: {exc}") from exc
        path.write_bytes(data)
        man = CACHE / "manifest.json"
        rec = json.loads(man.read_text(encoding="utf-8")) if man.exists() else {}
        rec[path.name] = {"url": url, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
        man.write_text(json.dumps(rec, indent=1), encoding="utf-8")
    return data if binary else data.decode("utf-8", errors="replace")


def properties(name: str) -> dict:
    q = urllib.parse.quote(name)
    d = json.loads(_get(f"{PUG}/compound/name/{q}/property/{PROPS}/JSON"))
    try:
        p = d["PropertyTable"]["Properties"][0]
    except (KeyError, IndexError) as exc:
        raise PubChemError(f"no PubChem record for {name!r}") from exc
    p["source"] = f"PubChem CID {p.get('CID')}"
    return p


def image_png(cid: int, size: int = 300) -> bytes:
    return _get(f"{PUG}/compound/cid/{int(cid)}/PNG?image_size={size}x{size}", binary=True)


def conformer_3d(cid: int) -> tuple[list, np.ndarray, list]:
    """(elements, coordinates in angstrom, bonds (i, j, order)) of PubChem's 3D conformer."""
    sdf = _get(f"{PUG}/compound/cid/{int(cid)}/SDF?record_type=3d")
    lines = sdf.splitlines()
    try:
        na, nb = int(lines[3][:3]), int(lines[3][3:6])
    except (IndexError, ValueError) as exc:
        raise PubChemError("no 3D conformer") from exc
    el, xyz, bonds = [], [], []
    for ln in lines[4:4 + na]:
        xyz.append([float(ln[0:10]), float(ln[10:20]), float(ln[20:30])])
        el.append(ln[31:34].strip())
    for ln in lines[4 + na:4 + na + nb]:
        bonds.append((int(ln[0:3]) - 1, int(ln[3:6]) - 1, int(ln[6:9])))
    return el, np.array(xyz), bonds


def drug_likeness(p: dict) -> dict:
    """Lipinski's rule of five and Veber's rules from PubChem's computed properties."""
    mw = float(p.get("MolecularWeight", 0) or 0)
    xlogp = p.get("XLogP")
    rules = {"MW <= 500": mw <= 500, "XLogP <= 5": (xlogp is None) or float(xlogp) <= 5,
             "H-bond donors <= 5": int(p.get("HBondDonorCount", 0)) <= 5,
             "H-bond acceptors <= 10": int(p.get("HBondAcceptorCount", 0)) <= 10,
             "rotatable bonds <= 10 (Veber)": int(p.get("RotatableBondCount", 0)) <= 10,
             "TPSA <= 140 (Veber)": float(p.get("TPSA", 0) or 0) <= 140}
    ro5 = sum(not rules[k] for k in ("MW <= 500", "XLogP <= 5", "H-bond donors <= 5", "H-bond acceptors <= 10"))
    return {"rules": rules, "ro5_violations": ro5}
