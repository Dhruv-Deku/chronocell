"""
Structural-variant files: VCF (SVTYPE / symbolic ALT / breakend ALT) and BEDPE.

Every variant becomes a `Variant` with 0-based, half-open reference coordinates:
    DEL / DUP / INV : chrom, start, end
    BND / TRA       : chrom, start (this side's breakpoint), mate_chrom, mate_pos, orientation
Chromosome names are normalised by the caller-supplied function (chronocell.genome.normalize_chrom
by default), so 'chr9', '9', 'Chr9' and RefSeq/GenBank accessions all work. Lines that cannot be
read are returned as `skipped` with the reason; nothing is silently dropped.
"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field
from typing import Callable

KINDS = ("DEL", "DUP", "INV", "BND")


@dataclass
class Variant:
    kind: str                   # DEL | DUP | INV | BND
    chrom: str
    start: int                  # 0-based
    end: int                    # half-open (= start + 1 for a breakend)
    mate_chrom: str | None = None
    mate_pos: int | None = None
    orientation: str = ""       # BEDPE strands or VCF bracket form, informative only
    vid: str = ""
    source_line: int = 0

    @property
    def length(self) -> int:
        return self.end - self.start if self.kind != "BND" else 0


@dataclass
class VariantFile:
    variants: list[Variant] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    format: str = ""


def _text(data: bytes | str) -> str:
    if isinstance(data, str):
        return data
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8", "replace")


def _info(s: str) -> dict[str, str]:
    out = {}
    for kv in s.split(";"):
        k, _, v = kv.partition("=")
        out[k.strip()] = v.strip()
    return out


_BND = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def read_vcf(data: bytes | str, norm: Callable[[str], str] | None = None) -> VariantFile:
    norm = norm or _default_norm
    vf = VariantFile(format="VCF")
    for ln, line in enumerate(_text(data).splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            vf.skipped.append(f"line {ln}: fewer than 8 VCF columns")
            continue
        chrom, pos, vid, ref, alt, _, _, info = f[:8]
        try:
            chrom = norm(chrom)
            p0 = int(pos) - 1
        except (ValueError, KeyError) as exc:
            vf.skipped.append(f"line {ln}: {exc}")
            continue
        inf = _info(info)
        kind = inf.get("SVTYPE", "").upper()
        m = re.match(r"<([A-Z]+)", alt)
        if not kind and m:
            kind = m.group(1)
        bnd = _BND.search(alt)
        if bnd and kind in ("", "BND", "TRA"):
            kind = "BND"
        kind = {"TRA": "BND", "CNV": "", "DUP:TANDEM": "DUP"}.get(kind, kind)
        if kind.startswith("DUP"):
            kind = "DUP"
        if kind not in KINDS:
            vf.skipped.append(f"line {ln}: unsupported SV type '{kind or alt[:20]}'")
            continue
        try:
            if kind == "BND":
                if bnd:
                    mate_chrom, mate_pos = norm(bnd.group(1)), int(bnd.group(2)) - 1
                else:
                    mate_chrom, mate_pos = norm(inf["CHR2"]), int(inf["END"]) - 1
                vf.variants.append(Variant("BND", chrom, p0, p0 + 1, mate_chrom, mate_pos, alt, vid, ln))
            else:
                if "END" in inf:
                    end = int(inf["END"])
                elif "SVLEN" in inf:
                    end = p0 + 1 + abs(int(inf["SVLEN"].split(",")[0]))
                else:
                    raise ValueError("no END or SVLEN")
                # symbolic alleles: POS is the base before the event (VCF 4.x), so the event starts at POS
                start = p0 + 1 if alt.startswith("<") else p0
                if end <= start:
                    raise ValueError(f"END {end} is not after POS")
                vf.variants.append(Variant(kind, chrom, start, end, vid=vid, source_line=ln))
        except (ValueError, KeyError) as exc:
            vf.skipped.append(f"line {ln}: {exc}")
    return vf


def read_bedpe(data: bytes | str, norm: Callable[[str], str] | None = None) -> VariantFile:
    """chrom1 start1 end1 chrom2 start2 end2 [name score strand1 strand2 [type]]. Type, when absent,
    comes from the strands: same chromosome (+,-) deletion, (-,+) duplication, (+,+)/(-,-) inversion;
    different chromosomes: translocation (BND)."""
    norm = norm or _default_norm
    vf = VariantFile(format="BEDPE")
    for ln, line in enumerate(_text(data).splitlines(), 1):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = re.split(r"\t|\s+", line.strip())
        if len(f) < 6:
            vf.skipped.append(f"line {ln}: fewer than 6 BEDPE columns")
            continue
        try:
            c1, s1, e1, c2, s2, e2 = norm(f[0]), int(f[1]), int(f[2]), norm(f[3]), int(f[4]), int(f[5])
        except (ValueError, KeyError) as exc:
            if ln == 1:
                continue                                    # header row
            vf.skipped.append(f"line {ln}: {exc}")
            continue
        name = f[6] if len(f) > 6 else f"sv{ln}"
        st1, st2 = (f[8], f[9]) if len(f) > 9 else ("", "")
        kind = f[10].upper() if len(f) > 10 else ""
        kind = {"TRA": "BND", "TRANSLOCATION": "BND", "DELETION": "DEL", "DUPLICATION": "DUP",
                "INVERSION": "INV"}.get(kind, kind)
        if c1 != c2:
            vf.variants.append(Variant("BND", c1, (s1 + e1) // 2, (s1 + e1) // 2 + 1, c2, (s2 + e2) // 2,
                                       f"{st1}{st2}", name, ln))
            continue
        if kind not in KINDS:
            kind = {("+", "-"): "DEL", ("-", "+"): "DUP", ("+", "+"): "INV", ("-", "-"): "INV"}.get((st1, st2), "")
        if kind not in ("DEL", "DUP", "INV"):
            vf.skipped.append(f"line {ln}: type unknown (give strands or an 11th type column)")
            continue
        a, b = sorted(((s1 + e1) // 2, (s2 + e2) // 2))
        if b <= a:
            vf.skipped.append(f"line {ln}: breakpoints coincide")
            continue
        vf.variants.append(Variant(kind, c1, a, b, orientation=f"{st1}{st2}", vid=name, source_line=ln))
    return vf


def read_any(data: bytes, name: str, norm: Callable[[str], str] | None = None) -> VariantFile:
    low = name.lower()
    text = _text(data)
    if low.endswith((".vcf", ".vcf.gz")) or text.lstrip().startswith("##fileformat=VCF"):
        return read_vcf(text, norm)
    return read_bedpe(text, norm)


def _default_norm(name: str) -> str:
    try:
        from .genome import normalize_chrom
        return normalize_chrom(name)
    except ImportError:          # pragma: no cover
        n = name.strip()
        return n if n.startswith("chr") else f"chr{n}"
