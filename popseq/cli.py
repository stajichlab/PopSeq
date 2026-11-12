"""The `popseq` command. A thin shell over popseq.query; no logic lives here.

  popseq species --min-isolates 100 --format tsv
  popseq panel --taxid 746128 --n 100 --min-depth 30 --format samples
  popseq describe --taxid 1047171
  popseq index
  popseq fetch SRR1553610 --tier coverage --coverage 30 --out "$SCRATCH/reads"

Exit 0 on success, 2 when the species is not found or a name is ambiguous,
3 when the filters leave nothing. A caller branches on that without parsing.

`popseq fetch` has its own exit codes (spec 002): 0 ok, 2 usage, 10 unavailable_at_ncbi,
11 unpaired, 12 verify_failed, 13 transient, 14 hook_failed. With a manifest the exit code
is the highest of its runs.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

from popseq import query as q

RUN_FORMATS = ["csv", "tsv", "json", "accessions", "samples"]


def _write_rows(rows, fields, fmt):
    w = csv.DictWriter(
        sys.stdout,
        fieldnames=fields,
        extrasaction="ignore",
        delimiter="\t" if fmt == "tsv" else ",",
    )
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
            print(
                f"# caveat: {dropped} run(s) also have an unpaired FASTQ file; the sheet "
                "uses only _1 and _2 (the full list is in fastq_ftp, --format json)",
                file=sys.stderr,
            )
    else:
        _write_rows(result["runs"], q.RUN_FIELDS, fmt)


def _target(a):
    return {"taxid": a.taxid, "name": a.species}


def _filters(a):
    return dict(
        platform=a.platform,
        layout=a.layout,
        strategy=a.strategy,
        min_depth=a.min_depth,
        min_gbases=a.min_gbases,
        country=a.country,
        host=a.host,
        since=a.since,
        include_excluded=a.include_excluded,
        only_exclusion=a.only_exclusion,
        require_fastq=a.require_fastq,
    )


def _cmd_fetch(a) -> int:
    """popseq fetch. Imported here so the query commands do not load the download code."""
    from popseq import fetch as fx
    from popseq import profile as prof
    from popseq.lookup import open_table

    try:
        if os.path.exists(a.run) or a.run.endswith((".csv", ".tsv", ".gz", ".zst")):
            runs = fx.read_manifest(a.run)
        else:
            runs = [a.run]
        plan_rows = {}
        if a.plan:
            with open_table(a.plan) as fh:
                head = fh.readline()
            delim = "\t" if "\t" in head else ","
            with open_table(a.plan) as fh:
                plan_rows = {r["run"]: r for r in csv.DictReader(fh, delimiter=delim)}
        opt = fx.Options(
            tier=a.tier,
            sampling=a.sampling,
            out=a.out,
            coverage=a.coverage,
            genome_mb=a.genome_mb,
            windows=a.windows,
            spots=a.spots,
            require_pairs=a.require_pairs,
            exec_cmd=a.exec_cmd,
            keep=a.keep,
            threads=a.threads,
            concurrency=a.concurrency,
        )
        crows = fx.catalog_rows(runs)
        fetcher = fx.Fetcher()
        worst = 0
        if a.header:
            print("\t".join(fx.FETCH_FIELDS))
        results_fh = None
        if a.results:
            new = not os.path.exists(a.results) or os.path.getsize(a.results) == 0
            results_fh = open(a.results, "a")
            if new:
                results_fh.write("\t".join(fx.FETCH_FIELDS) + "\n")
        for run in runs:
            try:
                res = fetcher.fetch(run, opt, plan_row=plan_rows.get(run), crow=crows.get(run))
            except fx.Usage as exc:
                print(f"popseq fetch: {run}: {exc}", file=sys.stderr)
                worst = max(worst, fx.EXIT["usage"])
                continue
            if a.profile and res.profile:
                prof.append_profile(a.profile, res, res.profile, a.threads)
            line = "\t".join(res.row())
            print(line)
            if results_fh:
                results_fh.write(line + "\n")
                results_fh.flush()
            if res.message:
                print(f"popseq fetch: {run}: {res.message}", file=sys.stderr)
            worst = max(worst, res.exit_code)
        if results_fh:
            results_fh.close()
        return worst
    except fx.Usage as exc:
        print(f"popseq fetch: {exc}", file=sys.stderr)
        return fx.EXIT["usage"]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="popseq", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
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
        sp.add_argument(
            "--platform",
            default="ILLUMINA",
            help="ILLUMINA (default), BGISEQ, DNBSEQ, or '' for any",
        )
        sp.add_argument(
            "--layout",
            default="PAIRED",
            help="PAIRED, SINGLE, or '' for any. Matches what ENA serves "
            "(ena_layout: _1 and _2 files = PAIRED), not ENA's label",
        )
        sp.add_argument(
            "--strategy", default="WGS", help="the catalog holds WGS GENOMIC runs only (issue #24)"
        )
        sp.add_argument(
            "--min-depth", type=float, help="estimated x; a run with no estimate is dropped"
        )
        sp.add_argument("--min-gbases", type=float)
        sp.add_argument("--country", help="ISO code (CH), name, or raw string")
        sp.add_argument("--host", help="host taxid (4565) or name (Triticum aestivum)")
        sp.add_argument("--since", help="first_public on or after, YYYY-MM-DD")
        sp.add_argument(
            "--include-excluded", action="store_true", help="keep mutant/evolved/cross/non-WGS rows"
        )
        sp.add_argument(
            "--only-exclusion",
            choices=["mutant", "evolution", "cross", "not_wgs"],
            help="keep ONLY that label, for a mutant or evolved panel",
        )
        sp.add_argument("--require-fastq", action="store_true", default=True)
        sp.add_argument("--no-require-fastq", dest="require_fastq", action="store_false")

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

    f = sub.add_parser(
        "fetch",
        help="download one run (or a manifest) at a depth tier, verified (spec 002)",
        description="Fetch reads from the archives, verify them against the archive that served "
        "them, and write a <RUN>.json sidecar. Temporary files go in $SCRATCH. One run per call "
        "is the supported mode for a Slurm array or a Nextflow task.",
    )
    f.add_argument("run", help="a run accession, or a CSV/TSV manifest with a run_accession column")
    f.add_argument(
        "--tier",
        choices=["head", "windows", "stride", "full", "coverage"],
        default="windows",
        help="depth tier (default: windows)",
    )
    f.add_argument(
        "--out",
        help="output directory (required unless --tier full, which goes to cache/<taxid>/<RUN>/)",
    )
    f.add_argument("--coverage", type=float, help="tier coverage: target x of the genome size")
    f.add_argument(
        "--sampling",
        choices=["stride", "windows"],
        default="stride",
        help="tier coverage: `stride` (default) keeps evenly spaced reads of the verified full run "
        "and reaches every part of the genome even if the run is sorted by genome position; "
        "`windows` cuts contiguous windows (cheaper, but covers only 25/depth of the genome on a "
        "position-sorted run). With --exec, a hook that exits 75 after a windowed sample makes "
        "fetch resample the same transfer by stride",
    )
    f.add_argument(
        "--genome-mb",
        type=float,
        help="tier coverage: genome size in Mb when catalog/genome_sizes.csv has none",
    )
    f.add_argument("--windows", type=int, default=8, help="windows for tiers windows/coverage")
    f.add_argument(
        "--spots",
        type=int,
        help="spots for head (default 500000), windows (default 2000000) or stride (target)",
    )
    f.add_argument(
        "--require-pairs",
        action="store_true",
        help="fail with exit 11 rather than fetch a run as single-end when no archive serves pairs",
    )
    f.add_argument(
        "--exec",
        dest="exec_cmd",
        metavar="CMD",
        help="run `CMD RUN FASTQ_1 [FASTQ_2]` on the verified reads before they are "
        "deleted (reads are not kept unless --keep, or tier full)",
    )
    f.add_argument("--keep", action="store_true", help="keep the reads even with --exec")
    f.add_argument(
        "--plan",
        help="a plan CSV (benchmark 7.1 format: run, ncbi_spots, "
        "ncbi_spots_with_mates, fastq_md5, ...) instead of live lookups",
    )
    f.add_argument("--results", help="append the result line(s) to this TSV")
    f.add_argument(
        "--profile",
        help="append one line per fetch to this TSV: seconds per stage (plan, transfer, convert, "
        "count, md5, cut, postcut, check, hook, finish, other) and peak memory of this process, "
        "of its largest child and, on Linux, of the whole process tree",
    )
    f.add_argument("--threads", type=int, default=4)
    f.add_argument(
        "--concurrency",
        type=int,
        default=20,
        help="how many fetches the caller runs at once; recorded in every "
        "downloads.csv row and the sidecar, not enforced (default 20; no failures "
        "up to 40 in the 2026-11-11 sweep, issue #27)",
    )
    f.add_argument("--header", action="store_true", help="print the result header line first")

    sub.add_parser("index", help="the catalog manifest, JSON")

    a = p.parse_args(argv)
    if a.cmd == "fetch":
        return _cmd_fetch(a)
    try:
        if a.cmd == "species":
            rows = q.species(min_isolates=a.min_isolates, contains=a.contains, limit=a.limit)
            if not rows:
                raise q.EmptyResult("no species match")
            if a.format == "json":
                json.dump(
                    {"popseq_schema": q.SCHEMA_VERSION, "species": rows}, sys.stdout, indent=2
                )
                sys.stdout.write("\n")
            else:
                _write_rows(rows, q.SPECIES_FIELDS, a.format)
        elif a.cmd == "runs":
            _emit_runs(q.runs(**_target(a), limit=a.limit, **_filters(a)), a.format)
        elif a.cmd == "panel":
            _emit_runs(
                q.panel(
                    **_target(a), n=a.n, max_per_study=a.max_per_study, seed=a.seed, **_filters(a)
                ),
                a.format,
            )
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
