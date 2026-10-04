"""
Add a genome assembly to ChronoCell-5D as data, not code.

    python -m chronocell.genome_fetch mm39 --species "Mus musculus" --common mouse --display GRCm39/mm39
    python -m chronocell.genome_fetch hg38 --aliases-only        # add the alias table to an existing genome

writes chronocell/data/genomes/<assembly>/
    manifest.json      assembly, display name, species, main chromosomes, file names, sources, licences
    chromosomes.json   per chromosome: size, cytobands, assembly gaps, centromere span (UCSC REST API)
    genes.json.gz      one RefSeq Select transcript per gene (UCSC track ncbiRefSeqSelect)
    chromAlias.tsv     UCSC / Ensembl / GenBank / RefSeq names of every sequence (UCSC bigZips)

Any folder with a valid manifest.json is picked up by chronocell.genome.assemblies(); nothing else
needs to change. Main chromosomes = sequences without '_' (no alt / random / unplaced contigs) except
chrM. Data: UCSC Genome Browser (Kent et al., Genome Res 2002; Perez et al., NAR 2025); RefSeq
(O'Leary et al., NAR 2016). UCSC data are free for academic, non-profit and personal use; RefSeq is
public NCBI data. HTTPS uses the operating system's trust store when the `truststore` package is
present (Python's own store can lack the intermediate certificate UCSC serves).
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

GENOMES = Path(__file__).with_name("data") / "genomes"
API = "https://api.genome.ucsc.edu"
UA = {"User-Agent": "ChronoCell-5D genome_fetch (research)"}


def _get(url: str, tries: int = 5) -> bytes:
    try:
        import truststore
        truststore.inject_into_ssl()
    except ImportError:
        pass
    last = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
                return r.read()
        except OSError as exc:
            last = exc
            time.sleep(min(30, 2 * 2 ** k))
    raise OSError(f"could not fetch {url}: {last}")


def _track(genome: str, track: str, chrom: str) -> list[dict]:
    q = urllib.parse.urlencode({"genome": genome, "track": track, "chrom": chrom}).replace("&", ";")
    data = json.loads(_get(f"{API}/getData/track?{q}"))
    return data.get(track, []) if isinstance(data.get(track), list) else []


def fetch(assembly: str, species: str, common: str, display: str, genes: bool = True) -> Path:
    out = GENOMES / assembly
    out.mkdir(parents=True, exist_ok=True)
    chroms = json.loads(_get(f"{API}/list/chromosomes?genome={assembly}"))["chromosomes"]
    main = sorted((c for c in chroms if "_" not in c and c != "chrM"), key=_order)
    rec = {}
    gene_rows = []
    for c in main:
        bands = [[int(b["chromStart"]), int(b["chromEnd"]), b.get("name", ""), b.get("gieStain", "")]
                 for b in _track(assembly, "cytoBandIdeo", c)]
        gaps = [[int(g["chromStart"]), int(g["chromEnd"])] for g in _track(assembly, "gap", c)]
        acen = [b for b in bands if b[3] == "acen"]
        rec[c] = {"size": int(chroms[c]), "bands": sorted(bands) or [[0, int(chroms[c]), "", "gneg"]],
                  "gaps": sorted(gaps), "centromere_model": [min(b[0] for b in acen), max(b[1] for b in acen)] if acen else None}
        if genes:
            for g in _track(assembly, "ncbiRefSeqSelect", c):
                acc = str(g.get("name", ""))
                gene_rows.append([g.get("name2") or acc, c, int(g["txStart"]), int(g["txEnd"]), g.get("strand", "+"),
                                  "coding" if acc.startswith("NM_") else "noncoding"])
        print(f"{c}: {len(rec[c]['bands'])} bands, {len(rec[c]['gaps'])} gaps", flush=True)
    stamp = dt.date.today().isoformat()
    (out / "chromosomes.json").write_text(json.dumps(
        {"assembly": display, "source": f"UCSC Genome Browser REST API: chromosomes, cytoBandIdeo, gap (retrieved {stamp})",
         "chromosomes": rec, "genes": {}}, separators=(",", ":")), encoding="utf-8")
    files = {"chromosomes": "chromosomes.json"}
    if genes and gene_rows:
        seen, uniq = set(), []
        for r in sorted(gene_rows, key=lambda r: (_order(r[1]), r[2])):
            if r[0] not in seen:
                seen.add(r[0])
                uniq.append(r)
        with gzip.open(out / "genes.json.gz", "wt", encoding="utf-8") as fh:
            json.dump({"source": f"UCSC Genome Browser REST API, {assembly} track ncbiRefSeqSelect (RefSeq Select, one "
                                 "transcript per gene)", "fields": ["name", "chrom", "start", "end", "strand", "biotype"],
                       "coordinates": "0-based half-open", "retrieved": stamp, "genes": uniq}, fh)
        files["genes"] = "genes.json.gz"
    fetch_aliases(assembly, out)
    files["aliases"] = "chromAlias.tsv"
    manifest = {"assembly": assembly, "display": display, "species": species, "common": common,
                "main_chromosomes": main, "files": files,
                "sources": {"annotation": "UCSC Genome Browser REST API (api.genome.ucsc.edu)",
                            "genes": "NCBI RefSeq Select via UCSC ncbiRefSeqSelect",
                            "aliases": f"hgdownload.soe.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}.chromAlias.txt"},
                "licence": "UCSC Genome Browser data: free for academic, non-profit and personal use; RefSeq: public NCBI data.",
                "citation": "Kent WJ et al., Genome Res 12:996 (2002); O'Leary NA et al., Nucleic Acids Res 44:D733 (2016)",
                "retrieved": stamp}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    return out


def fetch_aliases(assembly: str, out: Path) -> None:
    text = _get(f"https://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}.chromAlias.txt").decode()
    (out / "chromAlias.tsv").write_text(text, encoding="utf-8")


def _order(c: str) -> tuple[int, str]:
    s = c.removeprefix("chr")
    return (int(s), "") if s.isdigit() else ({"X": 100, "Y": 101, "M": 102}.get(s, 200), s)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("assembly")
    ap.add_argument("--species", default="")
    ap.add_argument("--common", default="")
    ap.add_argument("--display", default="")
    ap.add_argument("--no-genes", action="store_true")
    ap.add_argument("--aliases-only", action="store_true")
    a = ap.parse_args()
    if a.aliases_only:
        out = GENOMES / a.assembly
        out.mkdir(parents=True, exist_ok=True)
        fetch_aliases(a.assembly, out)
        print(f"-> {out / 'chromAlias.tsv'}")
        return
    print(f"-> {fetch(a.assembly, a.species, a.common, a.display or a.assembly, not a.no_genes)}")


if __name__ == "__main__":
    main()
