#!/usr/bin/env python3
"""Build catalog/runs/fungi_all_runs.csv.gz: the one run-level table.

This promotes the ENA snapshot out of cache/ and into the catalog, because the
catalog is what a consumer gets. A clone of this repository must be able to answer
"what strains of species X can I profile" without re-querying an archive.

Columns added on top of the ENA fields:

    species_taxid, species_name   rolled up through the NCBI taxonomy, never the
                                  deposited name (issue #1: synonyms split a taxon)
    exclusion                     "", mutant, evolution, cross or not_wgs. This
                                  LABELS; nothing is deleted (spec 001)
    est_depth_x                   base_count / reference assembly size, blank when
                                  no reference size resolved
    ncbi_status                   in_ncbi, ena_only, or not_checked

    uv run python scripts/popseq_runs_table.py --ena cache/fungi_runs_ena_2026-11-06.tsv.gz
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))
from popseq_survey import EXCLUDE_PATTERNS, RESCUE_PATTERN  # noqa: E402

ENA_FIELDS = [
    "run_accession", "study_accession", "sample_accession", "scientific_name",
    "tax_id", "instrument_platform", "instrument_model", "library_layout",
    "library_strategy", "library_source", "read_count", "base_count",
    "first_public", "center_name", "country", "collection_date", "strain",
    "isolate", "sample_title", "study_title", "fastq_ftp", "fastq_aspera",
]
OUT_FIELDS = ENA_FIELDS + ["species_taxid", "species_name", "exclusion",
                           "est_depth_x", "ncbi_status"]


def opener(path, mode="rt"):
    return gzip.open(path, mode, newline="") if path.endswith(".gz") else open(path, mode, newline="")


def label(row):
    text = " ".join(filter(None, (row.get("study_title", ""), row.get("sample_title", ""),
                                  row.get("strain", ""), row.get("isolate", ""))))
    if not text or re.search(RESCUE_PATTERN, text, re.I):
        return ""
    for category, pattern in EXCLUDE_PATTERNS.items():
        if re.search(pattern, text, re.I):
            return category
    return ""


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ena", required=True)
    p.add_argument("--taxmap", default=os.path.join(HERE, "catalog", "taxonomy_species_map.csv.gz"))
    p.add_argument("--genome-sizes", default=os.path.join(HERE, "catalog", "genome_sizes.csv"))
    p.add_argument("--ena-only", default=os.path.join(HERE, "catalog", "kingdom",
                                                      "fungi_byspecies.ena_only_runs.csv.gz"))
    p.add_argument("--summary", default=os.path.join(HERE, "catalog", "kingdom",
                                                     "fungi_species_summary.csv.gz"))
    p.add_argument("--out", default=os.path.join(HERE, "catalog", "runs", "fungi_all_runs.csv.gz"))
    args = p.parse_args()

    taxmap, names = {}, {}
    with opener(args.taxmap) as fh:
        for r in csv.DictReader(fh):
            if r.get("species_taxid"):
                taxmap[int(r["taxid"])] = int(r["species_taxid"])
                names[int(r["species_taxid"])] = r.get("species_name") or ""
    print(f"taxonomy: {len(taxmap)} taxids", file=sys.stderr)

    sizes = {}
    with open(args.genome_sizes, newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                sizes[int(r["species_taxid"])] = float(r["genome_mb"]) * 1e6
            except (TypeError, ValueError, KeyError):
                continue
    print(f"genome sizes: {len(sizes)} species", file=sys.stderr)

    # Which species were reconciled against NCBI at all. A run of a species that was
    # never checked must not be reported as "absent from NCBI".
    checked = set()
    with opener(args.summary) as fh:
        if not fh.readline().startswith("#"):
            fh.seek(0)
        for r in csv.DictReader(fh):
            if int(r["isolates"]) >= 100:
                checked.add(int(r["species_taxid"]))
    ena_only = set()
    with opener(args.ena_only) as fh:
        for r in csv.DictReader(fh):
            ena_only.add(r["run_accession"])
    print(f"reconciled species: {len(checked)}; ena-only runs: {len(ena_only)}", file=sys.stderr)

    n = 0
    counts = {"in_ncbi": 0, "ena_only": 0, "not_checked": 0}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    tmp = args.out + ".tmp.gz"
    with opener(args.ena) as src, gzip.open(tmp, "wt", newline="") as dst:
        w = csv.DictWriter(dst, fieldnames=OUT_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in csv.DictReader(src, delimiter="\t"):
            if not (r.get("run_accession") or "").strip():
                continue
            try:
                raw = int(str(r.get("tax_id") or "").strip())
            except ValueError:
                raw = None
            sp = taxmap.get(raw) if raw is not None else None
            r["species_taxid"] = sp or ""
            r["species_name"] = names.get(sp, "") if sp else ""
            r["exclusion"] = label(r)
            gsize = sizes.get(sp)
            try:
                bases = float(r.get("base_count") or 0)
            except ValueError:
                bases = 0
            r["est_depth_x"] = f"{bases / gsize:.1f}" if gsize and bases else ""
            if sp in checked:
                status = "ena_only" if r["run_accession"] in ena_only else "in_ncbi"
            else:
                status = "not_checked"
            r["ncbi_status"] = status
            counts[status] += 1
            w.writerow(r)
            n += 1
    os.replace(tmp, args.out)
    print(f"wrote {args.out}: {n} runs", file=sys.stderr)
    for k, v in counts.items():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
