#!/usr/bin/env python3
"""popseq query: ask the catalog what strains are available to profile.

This is the advertised interface. Other tools call this, not the tables, so the
table layout can change without breaking them. Everything it prints is derived
from committed catalog files; it never touches the network.

Four verbs:

  species   which species exist, with counts. The discovery entry point.
  runs      the runs of one species, filtered. One row per run.
  panel     a selected set of runs, one per isolate, spread across studies.
            This is what you want for "give me N strains to profile".
  describe  everything known about one species, as JSON.

Examples:

  # what is worth profiling at all
  popseq_query.py species --min-isolates 100 --format tsv

  # 100 Aspergillus fumigatus strains at 30x or better, no two from the same
  # study if that can be avoided, written as an nf-core sample sheet
  popseq_query.py panel --species "Aspergillus fumigatus" \\
      --n 100 --min-depth 30 --format samples > samples.csv

  # the same selection as JSON, for a program
  popseq_query.py panel --taxid 746128 --n 100 --min-depth 30 --format json

Exit status is 0 on success, 2 when the species is not in the catalog, and 3 when
the filters leave nothing. A tool can branch on that without parsing the output.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(HERE, "catalog", "runs", "fungi_all_runs.csv.gz")
SUMMARY = os.path.join(HERE, "catalog", "kingdom", "fungi_species_summary.csv.gz")
CLONALITY = os.path.join(HERE, "catalog", "summary", "clonality.csv.gz")
COVERAGE = os.path.join(HERE, "catalog", "summary", "coverage.csv.gz")
SCHEMA_VERSION = "1.0.0"

SAMPLE_FIELDS = ["sample", "fastq_1", "fastq_2"]
PANEL_FIELDS = ["run_accession", "sample_accession", "study_accession",
                "species_taxid", "species_name", "strain", "library_layout",
                "read_count", "base_count", "est_depth_x", "country",
                "collection_date", "first_public", "exclusion", "ncbi_status",
                "fastq_ftp", "fastq_aspera"]


def opener(path):
    return gzip.open(path, "rt", newline="") if path.endswith(".gz") else open(path, newline="")


def read_table(path):
    with opener(path) as fh:
        if not fh.readline().startswith("#"):
            fh.seek(0)
        yield from csv.DictReader(fh)


def die(msg, code):
    print(f"popseq_query: {msg}", file=sys.stderr)
    raise SystemExit(code)


def resolve(args):
    """Turn --species or --taxid into a species taxid using the catalog only."""
    rows = list(read_table(SUMMARY))
    if args.taxid:
        hit = [r for r in rows if int(r["species_taxid"]) == args.taxid]
        if not hit:
            die(f"taxid {args.taxid} is not in the catalog", 2)
        return int(hit[0]["species_taxid"]), hit[0]["species_name"], hit[0]
    want = args.species.strip().lower()
    exact = [r for r in rows if r["species_name"].lower() == want]
    if len(exact) == 1:
        return int(exact[0]["species_taxid"]), exact[0]["species_name"], exact[0]
    near = [r for r in rows if want in r["species_name"].lower()]
    if len(near) == 1:
        return int(near[0]["species_taxid"]), near[0]["species_name"], near[0]
    if not near:
        die(f'no species matching "{args.species}". '
            f"Try: popseq_query.py species --contains <word>", 2)
    names = ", ".join(sorted(r["species_name"] for r in near)[:8])
    die(f'"{args.species}" matches {len(near)} species: {names}'
        f'{" ..." if len(near) > 8 else ""}. Use --taxid.', 2)


def load_runs(taxid):
    for r in read_table(RUNS):
        if r["species_taxid"] and int(r["species_taxid"]) == taxid:
            yield r


def apply_filters(rows, a):
    """Filter runs. Every rule is a column test, so it is auditable from the table."""
    out = []
    for r in rows:
        if a.platform and (r["instrument_platform"] or "").upper() != a.platform.upper():
            continue
        if a.layout and (r["library_layout"] or "").upper() != a.layout.upper():
            continue
        if a.strategy and (r["library_strategy"] or "").upper() != a.strategy.upper():
            continue
        if not a.include_excluded and r["exclusion"]:
            continue
        if a.only_exclusion and r["exclusion"] != a.only_exclusion:
            continue
        if a.min_depth:
            try:
                if float(r["est_depth_x"]) < a.min_depth:
                    continue
            except (TypeError, ValueError):
                continue        # no estimate means it cannot be shown to pass
        if a.min_bases:
            try:
                if float(r["base_count"]) < a.min_bases * 1e9:
                    continue
            except (TypeError, ValueError):
                continue
        if a.country and (r["country"] or "").lower() != a.country.lower():
            continue
        if a.since and (r["first_public"] or "") < a.since:
            continue
        if a.require_fastq and not r["fastq_ftp"]:
            continue
        out.append(r)
    return out


def depth_of(r):
    try:
        return float(r["est_depth_x"])
    except (TypeError, ValueError):
        return -1.0


def pick_panel(rows, n, max_per_study, seed):
    """One run per isolate, then spread across studies by round-robin.

    Within an isolate the deepest run wins. Across studies the selection takes one
    run from each study in turn, so a species dominated by one outbreak study does
    not fill the whole panel with that outbreak unless there is nothing else. This
    is a sampling rule, not a claim of independence; see figure 04.
    """
    best = {}
    for r in rows:
        key = r["sample_accession"] or r["run_accession"]
        if key not in best or depth_of(r) > depth_of(best[key]):
            best[key] = r

    by_study = defaultdict(list)
    for r in best.values():
        by_study[r["study_accession"]].append(r)
    rng = random.Random(seed)
    for v in by_study.values():
        v.sort(key=lambda r: (-depth_of(r), r["run_accession"]))

    studies = sorted(by_study, key=lambda s: (-len(by_study[s]), s))
    rng.shuffle(studies)
    picked, i = [], 0
    while len(picked) < n:
        added = False
        for s in studies:
            if i < len(by_study[s]) and (max_per_study is None or i < max_per_study):
                picked.append(by_study[s][i])
                added = True
                if len(picked) >= n:
                    break
        if not added:
            break
        i += 1
    return picked, len(best), len(by_study)


def emit(rows, fmt, meta):
    if fmt == "accessions":
        for r in rows:
            print(r["run_accession"])
    elif fmt == "samples":
        # nf-core convention. The FASTQ columns are ENA URLs, not local paths: a
        # run that has not been fetched has no local path, and inventing one would
        # be a lie. Point them at the fetcher, or at the local panel once fetched.
        w = csv.DictWriter(sys.stdout, fieldnames=SAMPLE_FIELDS)
        w.writeheader()
        for r in rows:
            parts = (r["fastq_ftp"] or "").split(";")
            w.writerow({
                "sample": r["sample_accession"] or r["run_accession"],
                "fastq_1": "https://" + parts[0] if parts and parts[0] else "",
                "fastq_2": "https://" + parts[1] if len(parts) > 1 else "",
            })
    elif fmt == "json":
        json.dump({"popseq_schema": SCHEMA_VERSION, **meta,
                   "runs": [{k: r.get(k, "") for k in PANEL_FIELDS} for r in rows]},
                  sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        w = csv.DictWriter(sys.stdout, fieldnames=PANEL_FIELDS,
                           delimiter="\t" if fmt == "tsv" else ",",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def cmd_species(a):
    rows = [r for r in read_table(SUMMARY) if int(r["isolates"]) >= a.min_isolates]
    if a.contains:
        rows = [r for r in rows if a.contains.lower() in r["species_name"].lower()]
    rows.sort(key=lambda r: -int(r["isolates"]))
    if a.limit:
        rows = rows[:a.limit]
    if not rows:
        die("no species match", 3)
    cols = ["species_taxid", "species_name", "runs", "isolates", "studies",
            "est_download_tb", "first_year", "last_year"]
    if a.format == "json":
        json.dump({"popseq_schema": SCHEMA_VERSION,
                   "species": [{k: r[k] for k in cols} for r in rows]},
                  sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        w = csv.DictWriter(sys.stdout, fieldnames=cols,
                           delimiter="\t" if a.format == "tsv" else ",",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def cmd_runs(a):
    taxid, name, _ = resolve(a)
    rows = apply_filters(list(load_runs(taxid)), a)
    if a.limit:
        rows = rows[:a.limit]
    if not rows:
        die(f"{name}: no runs pass the filters", 3)
    emit(rows, a.format, {"species_taxid": taxid, "species_name": name,
                          "selected": len(rows)})


def cmd_panel(a):
    taxid, name, srow = resolve(a)
    rows = apply_filters(list(load_runs(taxid)), a)
    if not rows:
        die(f"{name}: no runs pass the filters", 3)
    picked, isolates, studies = pick_panel(rows, a.n, a.max_per_study, a.seed)
    if not picked:
        die(f"{name}: filters passed {len(rows)} runs but none survived selection", 3)
    used = len({r["study_accession"] for r in picked})
    meta = {
        "species_taxid": taxid, "species_name": name,
        "requested": a.n, "selected": len(picked),
        "runs_passing_filters": len(rows), "isolates_available": isolates,
        "studies_available": studies, "studies_in_panel": used,
        "note": ("Spread across studies by round-robin, deepest run per isolate. "
                 "Study spread is a sampling rule, not a guarantee of independence."),
    }
    if a.format in ("json",):
        emit(picked, a.format, meta)
    else:
        for k, v in meta.items():
            print(f"# {k}: {v}", file=sys.stderr)
        emit(picked, a.format, meta)


def cmd_describe(a):
    taxid, name, srow = resolve(a)
    clon = next((r for r in read_table(CLONALITY)
                 if int(r["species_taxid"]) == taxid), None)
    cov = next((r for r in read_table(COVERAGE)
                if int(r["species_taxid"]) == taxid), None)
    runs = list(load_runs(taxid))
    excl = defaultdict(int)
    for r in runs:
        excl[r["exclusion"] or "none"] += 1
    out = {
        "popseq_schema": SCHEMA_VERSION,
        "species_taxid": taxid, "species_name": name,
        "summary": srow, "clonality": clon, "coverage": cov,
        "runs_by_exclusion_label": dict(excl),
        "ncbi_status": {k: sum(1 for r in runs if r["ncbi_status"] == k)
                        for k in ("in_ncbi", "ena_only", "not_checked")},
        "caveats": [
            "Counts are a floor. ENA holds runs with an empty tax_id that no "
            "taxonomy query reaches (issue #2).",
            "est_depth_x is base_count over a reference assembly size, not a "
            "measurement on the reads.",
            "A run count is not a strain count. Check clonality.largest_study.",
        ],
    }
    json.dump(out, sys.stdout, indent=2, default=str)
    sys.stdout.write("\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("species", help="list species in the catalog")
    s.add_argument("--min-isolates", type=int, default=1)
    s.add_argument("--contains", help="substring of the species name")
    s.add_argument("--limit", type=int)
    s.add_argument("--format", choices=["csv", "tsv", "json"], default="csv")
    s.set_defaults(func=cmd_species)

    def common(q):
        g = q.add_mutually_exclusive_group(required=True)
        g.add_argument("--species")
        g.add_argument("--taxid", type=int)
        q.add_argument("--platform", default="ILLUMINA")
        q.add_argument("--layout", default="PAIRED", help="PAIRED, SINGLE or '' for any")
        q.add_argument("--strategy", default="WGS")
        q.add_argument("--min-depth", type=float,
                       help="estimated x; runs with no estimate are dropped")
        q.add_argument("--min-bases", type=float, help="Gbase")
        q.add_argument("--country")
        q.add_argument("--since", help="first_public on or after, YYYY-MM-DD")
        q.add_argument("--include-excluded", action="store_true",
                       help="keep mutant/evolved/cross/non-WGS rows")
        q.add_argument("--only-exclusion",
                       choices=["mutant", "evolution", "cross", "not_wgs"],
                       help="keep ONLY this label, for a mutant or evolved panel")
        q.add_argument("--require-fastq", action="store_true", default=True)
        q.add_argument("--no-require-fastq", dest="require_fastq",
                       action="store_false")

    r = sub.add_parser("runs", help="all runs of one species, filtered")
    common(r)
    r.add_argument("--limit", type=int)
    r.add_argument("--format", choices=["csv", "tsv", "json", "accessions", "samples"],
                   default="csv")
    r.set_defaults(func=cmd_runs)

    q = sub.add_parser("panel", help="select N strains spread across studies")
    common(q)
    q.add_argument("--n", type=int, default=100)
    q.add_argument("--max-per-study", type=int)
    q.add_argument("--seed", type=int, default=0)
    q.add_argument("--format", choices=["csv", "tsv", "json", "accessions", "samples"],
                   default="csv")
    q.set_defaults(func=cmd_panel)

    d = sub.add_parser("describe", help="everything known about one species, JSON")
    g = d.add_mutually_exclusive_group(required=True)
    g.add_argument("--species")
    g.add_argument("--taxid", type=int)
    d.set_defaults(func=cmd_describe)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
