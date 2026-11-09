#!/usr/bin/env python3
"""Take a dated ENA snapshot of fungal WGS runs: cache/fungi_runs_ena_<label>_<date>.tsv.gz.

Everything downstream (scripts/popseq_runs_table.py, popseq_species_summary.py) is built
from snapshots, so a snapshot is the unit of provenance. This script records the exact
query beside the output, in <out>.query.txt.

The query is the catalog scope: tax_tree(4751), library_strategy WGS, library_source
GENOMIC, and one or more instrument platforms. The fields are the 22 columns of the
2026-11-06 snapshot, in the same order.

    # The BGISEQ/DNBSEQ supplement (owner decision 2026-11-08):
    uv run python scripts/popseq_snapshot.py --platform BGISEQ DNBSEQ \\
        --out cache/fungi_runs_ena_bgi_2026-11-09.tsv.gz

Network: one POST to the ENA portal API per platform. Run it where the internet is
reachable (the login node can; compute nodes may not).
"""
from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import io
import os
import sys
import urllib.parse
import urllib.request

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = ["run_accession", "study_accession", "sample_accession", "scientific_name",
          "tax_id", "instrument_platform", "instrument_model", "library_layout",
          "library_strategy", "library_source", "read_count", "base_count", "first_public",
          "center_name", "country", "collection_date", "strain", "isolate", "sample_title",
          "study_title", "fastq_ftp", "fastq_aspera"]
SCOPE = ('tax_tree({taxid}) AND library_strategy="WGS" AND library_source="GENOMIC"'
         ' AND instrument_platform="{platform}"')


def query(taxid: int, platform: str) -> str:
    return SCOPE.format(taxid=taxid, platform=platform)


def fetch(q: str, timeout: int = 1800) -> list[dict]:
    body = urllib.parse.urlencode({"result": "read_run", "query": q,
                                   "fields": ",".join(FIELDS), "format": "tsv",
                                   "limit": 0}).encode()
    req = urllib.request.Request(ENA_SEARCH, data=body,
                                 headers={"User-Agent": "PopSeq-snapshot"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        text = resp.read().decode()
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if rows and list(rows[0]) != FIELDS:
        raise SystemExit(f"ENA returned unexpected columns: {list(rows[0])}")
    return rows


def write(rows: list[dict], out: str) -> None:
    tmp = out + ".tmp"
    with gzip.open(tmp, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, out)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--platform", nargs="+", required=True,
                   help="ENA instrument_platform values, e.g. ILLUMINA BGISEQ DNBSEQ")
    p.add_argument("--taxid", type=int, default=4751)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    rows, queries, seen = [], [], set()
    for plat in a.platform:
        q = query(a.taxid, plat)
        got = fetch(q)
        queries.append(f"{q}\t{len(got)} runs")
        print(f"{q}: {len(got)} runs", file=sys.stderr)
        for r in got:
            if r["run_accession"] not in seen:
                seen.add(r["run_accession"])
                rows.append(r)
    write(rows, a.out)
    with open(a.out + ".query.txt", "w") as fh:
        fh.write(f"# {ENA_SEARCH} result=read_run, taken "
                 f"{datetime.datetime.now(datetime.timezone.utc):%Y-%m-%dT%H:%M:%SZ}\n")
        fh.write("# fields: " + ",".join(FIELDS) + "\n")
        fh.write("\n".join(queries) + "\n")
    print(f"wrote {a.out}: {len(rows)} runs", file=sys.stderr)


if __name__ == "__main__":
    main()
