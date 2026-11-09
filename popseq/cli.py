"""The `popseq` command. A thin shell over popseq.query; no logic lives here.

  popseq species --min-isolates 100 --format tsv
  popseq panel --taxid 746128 --n 100 --min-depth 30 --format samples
  popseq describe --taxid 1047171
  popseq index

Exit 0 on success, 2 when the species is not found or a name is ambiguous,
3 when the filters leave nothing. A caller branches on that without parsing.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys

from popseq import query as q

RUN_FORMATS = ["csv", "tsv", "json", "accessions", "samples"]


def _write_rows(rows, fields, fmt):
    w = csv.DictWriter(sys.stdout, fieldnames=fields, extrasaction="ignore",
                       delimiter="\t" if fmt == "tsv" else ",")
    w.writeheader()
    for r in rows:
        w.writerow(r)


def _emit_runs(result, fmt):
    if fmt == "json":
        json.dump(result, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return
    # The metadata goes to stderr so stdout stays a clean table a pipe can read.
    for k, v in result.items():
        if k == "runs":
            continue
        if k == "caveats":
            # One line each, not a repr of the list.
            for c in v:
                print(f"# caveat: {c}", file=sys.stderr)
            continue
        print(f"# {k}: {v}", file=sys.stderr)
    if fmt == "accessions":
        for r in result["runs"]:
            print(r["run_accession"])
    elif fmt == "samples":
        _write_rows(q.to_sample_sheet(result), ["sample", "fastq_1", "fastq_2"], "csv")
        dropped = q.sample_sheet_dropped(result)
        if dropped:
            print(f"# caveat: {dropped} run(s) also have an unpaired FASTQ file; the sheet "
                  "uses only _1 and _2 (the full list is in fastq_ftp, --format json)",
                  file=sys.stderr)
    else:
        _write_rows(result["runs"], q.RUN_FIELDS, fmt)


def _target(a):
    return {"taxid": a.taxid, "name": a.species}


def _filters(a):
    return dict(platform=a.platform, layout=a.layout, strategy=a.strategy,
                min_depth=a.min_depth, min_gbases=a.min_gbases, country=a.country,
                host=a.host,
                since=a.since, include_excluded=a.include_excluded,
                only_exclusion=a.only_exclusion, require_fastq=a.require_fastq)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="popseq", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("species", help="list species in the catalog")
    s.add_argument("--min-isolates", type=int, default=1)
    s.add_argument("--contains")
    s.add_argument("--limit", type=int)
    s.add_argument("--format", choices=["csv", "tsv", "json"], default="csv")

    def add_target(sp):
        g = sp.add_mutually_exclusive_group(required=True)
        g.add_argument("--species")
        g.add_argument("--taxid", type=int)

    def add_filters(sp):
        sp.add_argument("--platform", default="ILLUMINA",
                        help="the catalog holds ILLUMINA runs only (issue #24)")
        sp.add_argument("--layout", default="PAIRED",
                        help="PAIRED, SINGLE, or '' for any. Matches what ENA serves "
                             "(ena_layout: _1 and _2 files = PAIRED), not ENA's label")
        sp.add_argument("--strategy", default="WGS",
                        help="the catalog holds WGS GENOMIC runs only (issue #24)")
        sp.add_argument("--min-depth", type=float,
                        help="estimated x; a run with no estimate is dropped")
        sp.add_argument("--min-gbases", type=float)
        sp.add_argument("--country", help="ISO code (CH), name, or raw string")
        sp.add_argument("--host", help="host taxid (4565) or name (Triticum aestivum)")
        sp.add_argument("--since", help="first_public on or after, YYYY-MM-DD")
        sp.add_argument("--include-excluded", action="store_true",
                        help="keep mutant/evolved/cross/non-WGS rows")
        sp.add_argument("--only-exclusion",
                        choices=["mutant", "evolution", "cross", "not_wgs"],
                        help="keep ONLY that label, for a mutant or evolved panel")
        sp.add_argument("--require-fastq", action="store_true", default=True)
        sp.add_argument("--no-require-fastq", dest="require_fastq",
                        action="store_false")

    r = sub.add_parser("runs", help="all runs of one species, filtered")
    add_target(r)
    add_filters(r)
    r.add_argument("--limit", type=int)
    r.add_argument("--format", choices=RUN_FORMATS, default="csv")

    pa = sub.add_parser("panel", help="select N strains spread across studies")
    add_target(pa)
    add_filters(pa)
    pa.add_argument("--n", type=int, default=100)
    pa.add_argument("--max-per-study", type=int)
    pa.add_argument("--seed", type=int, default=0)
    pa.add_argument("--format", choices=RUN_FORMATS, default="csv")

    d = sub.add_parser("describe", help="everything known about one species, JSON")
    add_target(d)

    sub.add_parser("index", help="the catalog manifest, JSON")

    a = p.parse_args(argv)
    try:
        if a.cmd == "species":
            rows = q.species(min_isolates=a.min_isolates, contains=a.contains,
                             limit=a.limit)
            if not rows:
                raise q.EmptyResult("no species match")
            if a.format == "json":
                json.dump({"popseq_schema": q.SCHEMA_VERSION, "species": rows},
                          sys.stdout, indent=2)
                sys.stdout.write("\n")
            else:
                _write_rows(rows, q.SPECIES_FIELDS, a.format)
        elif a.cmd == "runs":
            _emit_runs(q.runs(**_target(a), limit=a.limit, **_filters(a)), a.format)
        elif a.cmd == "panel":
            _emit_runs(q.panel(**_target(a), n=a.n, max_per_study=a.max_per_study,
                               seed=a.seed, **_filters(a)), a.format)
        elif a.cmd == "describe":
            json.dump(q.describe(**_target(a)), sys.stdout, indent=2, default=str)
            sys.stdout.write("\n")
        elif a.cmd == "index":
            json.dump(q.index(), sys.stdout, indent=2)
            sys.stdout.write("\n")
    except q.SpeciesNotFound as exc:
        print(f"popseq: {exc}", file=sys.stderr)
        return 2
    except q.EmptyResult as exc:
        print(f"popseq: {exc}", file=sys.stderr)
        return 3
    except q.CatalogError as exc:
        print(f"popseq: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
