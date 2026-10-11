"""
The `chronocell` command (Phase B7). After `pip install -e .`:

    chronocell analyze FILE --region chr21:28000000-30000000 --res 10000 [--norm kr] --out DIR
    chronocell diff --a A1.mcool A2.mcool --b B1.mcool B2.mcool --region ... --res 10000 [--fdr 0.05] --out DIR
    chronocell impact FILE --region ... --res 10000 (--variants SV.vcf | --joins "A:120:L-A:180:R" |
                      --segments "A:0-120 + A:180-300" | --cnv CNV.bed | --cnv-from-coverage)
                      [--partner FILE2 --partner-region chr22:...] [--zygosity heterozygous] [--refits 8] --out DIR
    chronocell predict ...            (the existing `python -m chronocell.predict`, unchanged)
    chronocell batch SHEET.csv --out DIR          columns: sample,path,region,resolution[,condition][,normalization]
    chronocell report DIR [--pdf]                 an HTML (and PDF) report of an analyze / diff / impact folder

Every command also runs as `python -m chronocell.cli ...`; the existing `python -m chronocell.<module>` commands
keep working. Outputs: tables (TSV / CSV), BED / BEDPE / bedGraph, Juicebox 2D annotations, summary.json and a
run record (run.json: inputs with SHA-256, versions, parameters).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np


def parse_region(text: str) -> tuple[str, int, int]:
    m = re.match(r"^\s*([^:\s]+):([\d,_]+)-([\d,_]+)\s*$", text or "")
    if not m:
        raise SystemExit(f"Region '{text}' should look like chr21:28000000-30000000.")
    a, b = int(m.group(2).replace(",", "").replace("_", "")), int(m.group(3).replace(",", "").replace("_", ""))
    if b <= a:
        raise SystemExit(f"Region '{text}' is empty.")
    return m.group(1), a, b


def _sha256(path: str) -> str | None:
    p = Path(path)
    if not p.exists():
        return None
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_record(out: Path, command: str, args: dict, inputs: list[str]) -> Path:
    from . import provenance as PV
    rec = {"command": command, "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "parameters": {k: (v if isinstance(v, (int, float, str, bool, type(None), list)) else str(v)) for k, v in args.items()},
           "inputs": [{"path": p, "sha256": _sha256(p)} for p in inputs], "software": PV.software_versions()}
    out.mkdir(parents=True, exist_ok=True)
    (out / "run.json").write_text(json.dumps(rec, indent=1, default=str), encoding="utf-8")
    return out / "run.json"


def _read(path: str, region: str, res: int, norm: str = "none", assembly: str | None = None):
    from . import contacts_io as IO
    chrom, a, b = parse_region(region)
    try:
        return IO.read_region(path, chrom, a, b, res, norm, assembly)
    except IO.ContactFileError as exc:
        raise SystemExit(f"{Path(path).name}: {exc}") from exc


def cmd_analyze(a) -> dict:
    from . import pipelines as PL
    c = _read(a.file, a.region, a.res, a.norm, a.assembly)
    out = Path(a.out)
    r = PL.analyze(c, out, loops=not a.no_loops)
    run_record(out, "analyze", vars(a), [a.file])
    print(json.dumps(r["summary"], indent=1, default=float))
    return r


def cmd_diff(a) -> dict:
    from . import pipelines as PL
    A = [_read(p, a.region, a.res, "none", a.assembly) for p in a.a]
    B = [_read(p, a.region, a.res, "none", a.assembly) for p in a.b]
    out = Path(a.out)
    r = PL.diff(A, B, a.fdr, a.max_sep, a.min_count, out, loops=not a.no_loops)
    run_record(out, "diff", vars(a), list(a.a) + list(a.b))
    print(json.dumps(r["summary"], indent=1, default=float))
    return r


def _karyotype(a, srcs, contacts):
    from . import svio, sv_engine as SV
    sizes = {s.name: s.n for s in srcs}
    windows = {s.name: (s.chrom, s.bin0, s.n, s.resolution) for s in srcs}
    out = []
    if a.joins:
        joins = []
        for part in a.joins.split(","):
            m = re.match(r"^\s*(\w+):(\d+):([LR])\s*-\s*(\w+):(\d+):([LR])\s*$", part)
            if not m:
                raise SystemExit(f"Join '{part}' should look like A:120:L-B:40:R.")
            joins.append(SV.Join(m.group(1), int(m.group(2)), m.group(3), m.group(4), int(m.group(5)), m.group(6)))
        out.append(("joins " + a.joins, SV.karyotype_from_joins(sizes, joins, a.zygosity)))
    if a.segments:
        derivs = [SV.parse_segments(t) for t in a.segments.split(";")]
        out.append(("segments " + a.segments, SV.karyotype_from_segments(sizes, derivs)))
    if a.cnv or a.cnv_from_coverage:
        s0, c0 = srcs[0], contacts[0]
        if a.cnv:
            segs = SV.read_cnv_bed(Path(a.cnv).read_bytes(), s0.chrom, s0.bin0, s0.n, s0.resolution)
            label = f"copy number from {Path(a.cnv).name}"
        else:
            _, segs = SV.estimate_copy_number(c0.ci, c0.cj, c0.cm, c0.n)
            label = "copy number estimated from Hi-C coverage (not validated)"
        if segs:
            out.append((label, SV.karyotype_from_copy_number(s0.n, segs, s0.name)))
    if a.variants:
        data = Path(a.variants).read_bytes()
        vf = svio.read_any(data, a.variants)
        gts = SV.vcf_genotypes(data) if vf.format == "VCF" else {}
        for v in vf.variants:
            js = SV.joins_for_variant(v, windows)
            label = svio.label(v)
            if isinstance(js, str):
                print(f"skipped {label}: {js}", file=sys.stderr)
                continue
            zyg = SV.zygosity(gts.get(v.source_line, {}).get("gt")) if gts else a.zygosity
            try:
                out.append((label, SV.karyotype_from_joins(sizes, js, zyg)))
            except ValueError as exc:
                print(f"skipped {label}: {exc}", file=sys.stderr)
    if not out:
        raise SystemExit("Nothing to apply: give --variants, --joins, --segments, --cnv or --cnv-from-coverage.")
    return out


def cmd_impact(a) -> dict:
    from . import pipelines as PL
    contacts = [_read(a.file, a.region, a.res, "none", a.assembly)]
    names = ["A"]
    if a.partner:
        contacts.append(_read(a.partner, a.partner_region, a.res, "none", a.assembly))
        names.append("B")
    fits = [PL.fit_source(nm, c) for nm, c in zip(names, contacts)]
    srcs = [f[0] for f in fits]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    results = []
    for label, kt in _karyotype(a, srcs, contacts):
        refit = None
        if a.refits:
            rng = np.random.default_rng(0)

            def refit(r, cs=contacts, nm=names):
                from . import contacts_io as IO
                return [PL.fit_source(n, IO.Contacts(c.ci, c.cj, rng.poisson(c.cm).astype(float), c.n, c.chrom, c.start,
                                                     c.resolution))[0] for n, c in zip(nm, cs)]
        r = PL.impact(srcs, kt, a.assembly, refit, a.refits)
        results.append((label, r))
    table = PL.rank_variants(results)
    table.to_csv(out / "variant_ranking.csv", index=False)
    rows = []
    for label, r in results:
        for g in r["genes"]:
            rows.append({"variant": label, **{k: (", ".join(v) if isinstance(v, list) else v) for k, v in g.items()}})
    import pandas as pd
    pd.DataFrame(rows).to_csv(out / "genes.csv", index=False)
    pd.DataFrame([{"variant": label, **b} for label, r in results for b in r["boundaries"]]).to_csv(out / "boundaries.csv", index=False)
    pd.DataFrame([{"variant": label, **e} for label, r in results for e in r["ep"]]).to_csv(out / "ep_pairs.csv", index=False)
    (out / "summary.json").write_text(json.dumps({"variants": [{"variant": l, **r["summary"]} for l, r in results],
                                                  "standing": PL.STANDING["impact"]}, indent=1, default=float),
                                      encoding="utf-8")
    run_record(out, "impact", vars(a), [a.file] + ([a.partner] if a.partner else []) + ([a.variants] if a.variants else []))
    print(table.to_string(index=False))
    return {"ranking": table, "results": results}


def cmd_batch(a) -> list:
    import pandas as pd
    sheet = pd.read_csv(a.sheet)
    need = {"sample", "path", "region", "resolution"}
    if not need <= set(sheet.columns):
        raise SystemExit(f"The sample sheet needs columns {', '.join(sorted(need))} (has {', '.join(sheet.columns)}).")
    out = Path(a.out)
    done = []
    for r in sheet.itertuples(index=False):
        dest = out / str(r.sample)
        if (dest / "summary.json").exists() and not a.force:
            print(f"{r.sample}: done earlier (resume)")
            done.append(dest)
            continue
        ns = argparse.Namespace(file=str(r.path), region=str(r.region), res=int(r.resolution),
                                norm=str(getattr(r, "normalization", "none") or "none"), assembly=a.assembly,
                                out=str(dest), no_loops=False)
        print(f"{r.sample}: analyse")
        cmd_analyze(ns)
        done.append(dest)
    if "condition" in sheet.columns:
        conds = list(dict.fromkeys(sheet["condition"].astype(str)))
        if len(conds) >= 2:
            ca, cb = conds[:2]
            sa, sb = sheet[sheet["condition"].astype(str) == ca], sheet[sheet["condition"].astype(str) == cb]
            regions = set(sa["region"]) | set(sb["region"])
            if len(regions) == 1:
                ns = argparse.Namespace(a=list(sa["path"]), b=list(sb["path"]), region=regions.pop(),
                                        res=int(sheet["resolution"].iloc[0]), fdr=a.fdr, max_sep=2_000_000, min_count=5.0,
                                        assembly=a.assembly, out=str(out / f"diff_{ca}_vs_{cb}"), no_loops=False)
                print(f"differential: {ca} vs {cb}")
                cmd_diff(ns)
    return done


def cmd_report(a) -> Path:
    from . import report_html as RH
    path = RH.report_folder(Path(a.folder), pdf=a.pdf)
    print(path)
    return path


def cmd_predict(rest: list[str]) -> None:
    from . import predict as PD
    sys.argv = ["chronocell predict"] + rest
    PD.main()


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="chronocell", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("analyze", help="loops, domains, compartments, insulation, P(s) of one map")
    p.add_argument("file")
    p.add_argument("--region", required=True)
    p.add_argument("--res", type=int, required=True)
    p.add_argument("--norm", choices=("none", "ice", "kr"), default="none")
    p.add_argument("--assembly")
    p.add_argument("--no-loops", action="store_true")
    p.add_argument("--out", required=True)
    p = sub.add_parser("diff", help="two conditions with replicates: differential contacts, loops, boundaries, compartments")
    p.add_argument("--a", nargs="+", required=True)
    p.add_argument("--b", nargs="+", required=True)
    p.add_argument("--region", required=True)
    p.add_argument("--res", type=int, required=True)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--max-sep", type=int, default=2_000_000)
    p.add_argument("--min-count", type=float, default=5.0)
    p.add_argument("--assembly")
    p.add_argument("--no-loops", action="store_true")
    p.add_argument("--out", required=True)
    p = sub.add_parser("impact", help="variant impact engine v2 (mechanism simulator, not validated)")
    p.add_argument("file")
    p.add_argument("--region", required=True)
    p.add_argument("--res", type=int, required=True)
    p.add_argument("--partner")
    p.add_argument("--partner-region")
    p.add_argument("--variants")
    p.add_argument("--joins")
    p.add_argument("--segments")
    p.add_argument("--cnv")
    p.add_argument("--cnv-from-coverage", action="store_true")
    p.add_argument("--zygosity", choices=("heterozygous", "homozygous"), default="heterozygous")
    p.add_argument("--refits", type=int, default=0)
    p.add_argument("--assembly")
    p.add_argument("--out", required=True)
    p = sub.add_parser("batch", help="the analysis over every sample of a CSV sample sheet (resumable)")
    p.add_argument("sheet")
    p.add_argument("--out", required=True)
    p.add_argument("--assembly")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--force", action="store_true")
    p = sub.add_parser("report", help="HTML (and PDF) report of an output folder")
    p.add_argument("folder")
    p.add_argument("--pdf", action="store_true")
    sub.add_parser("predict", help="prediction without contact data (python -m chronocell.predict)", add_help=False)
    return ap


def main(argv: list[str] | None = None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "predict":
        return cmd_predict(argv[1:])
    a = build_parser().parse_args(argv)
    {"analyze": cmd_analyze, "diff": cmd_diff, "impact": cmd_impact, "batch": cmd_batch, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
