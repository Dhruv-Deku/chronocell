"""
Gene annotations for the variant impact report and the gene tables (Phase B1 / B5). Information only:
none of these is a diagnosis, and the app labels them so.

ClinVar   gene_specific_summary.txt (NCBI, public domain), downloaded on demand and checked against the
          .md5 file NCBI publishes beside it: per gene, the number of alleles reported pathogenic or likely
          pathogenic.
GTEx      median gene TPM per tissue (GTEx Analysis v10, GTEx portal open-access data), downloaded on demand
          from the portal's Google Cloud bucket and checked against the MD5 the bucket publishes for the object
          (x-goog-hash header).
COSMIC    the Cancer Gene Census needs a COSMIC account and licence, so it is never downloaded here; a user can
          load their own CGC export (CSV with a 'Gene Symbol' column) with `cosmic_from_file`.
"""

from __future__ import annotations

import base64
import gzip
import io
import urllib.request
from functools import lru_cache
from pathlib import Path

import pandas as pd

CACHE = Path(__file__).resolve().parent.parent / ".chronocell_cache" / "annotations"
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/gene_specific_summary.txt"
GTEX_URL = ("https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
            "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")
UA = {"User-Agent": "ChronoCell-5D (research; annotations)"}
SOURCES = {
    "clinvar": {"url": CLINVAR_URL, "checksum": "NCBI .md5 file", "licence": "NCBI ClinVar: public domain (US government work); cite "
                "Landrum MJ et al., Nucleic Acids Res 46, D1062 (2018)."},
    "gtex": {"url": GTEX_URL, "checksum": "Google Cloud object MD5 (x-goog-hash)", "licence": "GTEx Portal open-access data; "
             "cite the GTEx Consortium, Science 369, 1318 (2020)."},
    "cosmic": {"url": "https://cancer.sanger.ac.uk/census", "checksum": "user-supplied file", "licence": "COSMIC licence "
               "(registration required); not downloaded by ChronoCell."},
}


def _clinvar_path(download: bool) -> Path | None:
    out = CACHE / "gene_specific_summary.txt"
    if out.exists():
        return out
    if not download:
        return None
    from .predict import download_verified
    md5 = urllib.request.urlopen(urllib.request.Request(CLINVAR_URL + ".md5", headers=UA), timeout=60).read().decode().split()[0]
    return download_verified(CLINVAR_URL, out, md5)


@lru_cache(maxsize=2)
def clinvar_genes(download: bool = False) -> dict[str, int] | None:
    """Gene symbol -> alleles reported pathogenic / likely pathogenic in ClinVar (None if not downloaded)."""
    path = _clinvar_path(download)
    if path is None:
        return None
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.startswith("#Symbol") or ln.startswith("Symbol"))
    df = pd.read_csv(io.StringIO("\n".join(lines[head:]).lstrip("#")), sep="\t", low_memory=False)
    col = next(c for c in df.columns if "Pathogenic" in c and "Alleles" in c)
    return dict(zip(df["Symbol"].astype(str), pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)))


def _gtex_path(download: bool) -> Path | None:
    out = CACHE / Path(GTEX_URL).name
    if out.exists():
        return out
    if not download:
        return None
    from .predict import download_verified
    req = urllib.request.Request(GTEX_URL, method="HEAD", headers=UA)
    h = urllib.request.urlopen(req, timeout=60).headers.get("x-goog-hash", "")
    b64 = next((p.split("=", 1)[1] for p in h.split(",") if p.strip().startswith("md5=")), None)
    if not b64:
        raise OSError("The GTEx bucket did not publish an MD5 for the file; not downloaded.")
    return download_verified(GTEX_URL, out, base64.b64decode(b64).hex())


@lru_cache(maxsize=1)
def gtex_median_tpm(download: bool = False) -> pd.DataFrame | None:
    """Genes x tissues median TPM (index: gene symbol), or None if not downloaded."""
    path = _gtex_path(download)
    if path is None:
        return None
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        fh.readline()
        fh.readline()
        df = pd.read_csv(fh, sep="\t")
    df = df.drop(columns=[c for c in ("Name", "id") if c in df.columns]).groupby("Description").max()
    df.index.name = "gene"
    return df


def gtex_for(genes: list[str], tissue: str | None = None, download: bool = False) -> pd.DataFrame:
    tab = gtex_median_tpm(download)
    if tab is None:
        return pd.DataFrame(columns=["gene"])
    sub = tab.reindex([g for g in genes if g in tab.index])
    if tissue:
        sub = sub[[c for c in sub.columns if c == tissue]]
    return sub.reset_index()


def cosmic_from_file(data: bytes, name: str = "") -> set[str]:
    """Gene symbols of a user-supplied COSMIC Cancer Gene Census export (CSV or TSV)."""
    text = data.decode("utf-8", "replace")
    sep = "\t" if name.endswith((".tsv", ".txt")) or text.count("\t") > text.count(",") else ","
    df = pd.read_csv(io.StringIO(text), sep=sep)
    col = next((c for c in df.columns if c.lower().replace("_", " ") in ("gene symbol", "gene", "symbol")), None)
    if col is None:
        raise ValueError("No 'Gene Symbol' column in the COSMIC file.")
    return set(df[col].astype(str).str.strip())
