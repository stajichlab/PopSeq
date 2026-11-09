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
    country_iso, country_name,    `country` normalised (issue #23). `country` itself
    country_status                stays exactly as harvested. See popseq/normalise.py
    host                          raw host from the BioSample table (first of host,
                                  host_plant, specific hosts, host_common_name)
    host_taxid, host_name,        `host` resolved through catalog/vocab/host_taxa.csv.gz,
    host_status                   built by scripts/popseq_vocab.py
    ena_layout                    PAIRED, SINGLE or "" from the files ENA serves
                                  (popseq.query.ena_layout); library_layout is the
                                  submitter's label and is kept as metadata only

    uv run python scripts/popseq_runs_table.py --ena cache/fungi_runs_ena_2026-11-06.tsv.gz \\
        cache/fungi_runs_ena_bgi_2026-11-09.tsv.gz
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
from popseq_vocab import biosample_hosts  # noqa: E402
from popseq_species_summary import summarise  # noqa: E402

sys.path.insert(0, HERE)
from popseq import normalise  # noqa: E402
from popseq.query import ena_layout  # noqa: E402

ENA_FIELDS = [
    "run_accession", "study_accession", "sample_accession", "scientific_name",
    "tax_id", "instrument_platform", "instrument_model", "library_layout",
    "library_strategy", "library_source", "read_count", "base_count",
    "first_public", "center_name", "country", "collection_date", "strain",
    "isolate", "sample_title", "study_title", "fastq_ftp", "fastq_aspera",
]
OUT_FIELDS = ENA_FIELDS + ["species_taxid", "species_name", "exclusion",
                           "est_depth_x", "ncbi_status",
                           "country_iso", "country_name", "country_status",
                           "host", "host_taxid", "host_name", "host_status",
                           "ena_layout"]


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
    p.add_argument("--ena", nargs="+", required=True,
                   help="ENA snapshot(s); a run in more than one is taken from the first")
    p.add_argument("--reconciled-ena",
                   help="the snapshot the NCBI reconciliation compared against "
                        "(default: the first --ena). Runs not in it are not_checked")
    p.add_argument("--taxmap", default=os.path.join(HERE, "catalog", "taxonomy_species_map.csv.gz"))
    p.add_argument("--genome-sizes", default=os.path.join(HERE, "catalog", "genome_sizes.csv"))
    p.add_argument("--ena-only", default=os.path.join(HERE, "catalog", "kingdom",
                                                      "fungi_byspecies.ena_only_runs.csv.gz"))
    p.add_argument("--biosamples", default=os.path.join(HERE, "catalog", "kingdom",
                                                        "fungi_biosamples.csv.gz"))
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

    # Which runs were reconciled against NCBI at all. A run that was never checked must
    # not be reported as in_ncbi or as absent from NCBI. The reconciliation covered the
    # species with 100 or more isolates in the snapshot it compared against, and only
    # the runs in that snapshot. A run from a later snapshot is not_checked.
    reconciled_path = args.reconciled_ena or args.ena[0]
    with opener(reconciled_path) as fh:
        rec_rows = list(csv.DictReader(fh, delimiter="\t"))
    reconciled_runs = {r["run_accession"] for r in rec_rows}
    checked = {int(d["species_taxid"]) for d in
               summarise(rec_rows, {str(k): (str(v), "") for k, v in taxmap.items()})
               if d["isolates"] >= 100}
    del rec_rows
    ena_only = set()
    with opener(args.ena_only) as fh:
        for r in csv.DictReader(fh):
            ena_only.add(r["run_accession"])
    print(f"reconciled species: {len(checked)}; ena-only runs: {len(ena_only)}", file=sys.stderr)

    hosts = biosample_hosts(args.biosamples)
    print(f"biosample hosts: {sum(1 for h, _ in hosts.values() if h)} of {len(hosts)}",
          file=sys.stderr)

    n = 0
    counts = {"in_ncbi": 0, "ena_only": 0, "not_checked": 0}
    cstat, hstat, lstat = {}, {}, {}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    tmp = args.out + ".tmp.gz"
    def snapshot_rows():
        for path in args.ena:
            with opener(path) as src:
                yield from csv.DictReader(src, delimiter="\t")

    seen = set()
    with gzip.open(tmp, "wt", newline="") as dst:
        w = csv.DictWriter(dst, fieldnames=OUT_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in snapshot_rows():
            if not (r.get("run_accession") or "").strip() or r["run_accession"] in seen:
                continue
            seen.add(r["run_accession"])
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
            if sp in checked and r["run_accession"] in reconciled_runs:
                status = "ena_only" if r["run_accession"] in ena_only else "in_ncbi"
            else:
                status = "not_checked"
            r["ncbi_status"] = status
            counts[status] += 1
            r.update(normalise.country(r.get("country") or ""))
            # sample_accession can hold several ;-joined BioSamples (issue #3).
            raw_hosts = [hosts.get(s.strip(), ("", ""))[0]
                         for s in (r.get("sample_accession") or "").split(";")]
            r["host"] = ";".join(raw_hosts) if any(raw_hosts) else ""
            r.update(normalise.host(r["host"], organism_taxids=(raw, sp)))
            r["ena_layout"] = ena_layout(r.get("fastq_ftp") or "")
            k = (r.get("library_layout") or "", r["ena_layout"])
            lstat[k] = lstat.get(k, 0) + 1
            cstat[r["country_status"]] = cstat.get(r["country_status"], 0) + 1
            hstat[r["host_status"]] = hstat.get(r["host_status"], 0) + 1
            w.writerow(r)
            n += 1
    os.replace(tmp, args.out)
    print(f"wrote {args.out}: {n} runs", file=sys.stderr)
    for k, v in counts.items():
        print(f"  {k}: {v}", file=sys.stderr)
    print("library_layout (label) x ena_layout (served): " + ", ".join(
        f"{a or '(empty)'}->{b or '(none)'} {v}" for (a, b), v in sorted(lstat.items())),
        file=sys.stderr)
    for label_, d in (("country_status", cstat), ("host_status", hstat)):
        print(f"{label_}: " + ", ".join(f"{k} {v}" for k, v in
                                        sorted(d.items(), key=lambda x: -x[1])),
              file=sys.stderr)


if __name__ == "__main__":
    main()
