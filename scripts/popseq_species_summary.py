#!/usr/bin/env python3
"""Build catalog/kingdom/fungi_species_summary.csv.gz from ENA run snapshots.

One row per species taxid (NCBI species-rank ancestor of the run's tax_id, through
catalog/taxonomy_species_map.csv.gz). A run whose tax_id has no species-rank ancestor is
not counted here (issue #22).

Columns:
    species_taxid, species_name
    runs              runs in the snapshot(s)
    isolates          distinct BioSamples (a ';'-joined sample_accession is split)
    studies           distinct study_accession values
    total_gbase       sum of base_count / 1e9
    est_download_tb   total bases / 3.5 / 1e12 (gzip FASTQ estimate; divisor UNMEASURED,
                      spec 002)
    paired_runs, single_runs   by ENA's library_layout label, as harvested
    first_year, last_year      range of first_public

Rows are sorted by isolates, descending (a stable sort on first appearance in the input).
Lines end in CRLF, Python's csv default, as in the first committed build.

    uv run python scripts/popseq_species_summary.py \\
        --ena cache/fungi_runs_ena_2026-11-06.tsv.gz cache/fungi_runs_ena_bgi_2026-11-09.tsv.gz
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIELDS = ["species_taxid", "species_name", "runs", "isolates", "studies", "total_gbase",
          "est_download_tb", "paired_runs", "single_runs", "first_year", "last_year"]


def opener(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    if path.endswith(".zst"):
        import zstandard
        return io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(open(path, "rb")),
                                newline="")
    return open(path, newline="")


def load_taxmap(path):
    out = {}
    with opener(path) as fh:
        for r in csv.DictReader(fh):
            if r.get("species_taxid"):
                out[r["taxid"]] = (r["species_taxid"], r["species_name"])
    return out


def summarise(rows, taxmap):
    acc = defaultdict(lambda: {"runs": 0, "samples": set(), "studies": set(), "bases": 0,
                               "paired": 0, "single": 0, "years": set(), "name": ""})
    seen = set()
    for r in rows:
        run = (r.get("run_accession") or "").strip()
        if not run or run in seen:
            continue
        seen.add(run)
        e = taxmap.get((r.get("tax_id") or "").strip())
        if not e:
            continue
        a = acc[e[0]]
        a["name"] = e[1]
        a["runs"] += 1
        # A multi-sample run lists its BioSamples ';'-joined (issue #3); count each.
        a["samples"].update(x.strip() for x in (r.get("sample_accession") or "").split(";")
                            if x.strip())
        if r.get("study_accession"):
            a["studies"].add(r["study_accession"])
        try:
            a["bases"] += int(r.get("base_count") or 0)
        except ValueError:
            pass
        lay = r.get("library_layout") or ""
        a["paired"] += lay == "PAIRED"
        a["single"] += lay == "SINGLE"
        y = (r.get("first_public") or "")[:4]
        if y.isdigit():
            a["years"].add(y)
    out = []
    for sp, a in acc.items():
        out.append({
            "species_taxid": sp, "species_name": a["name"], "runs": a["runs"],
            "isolates": len(a["samples"]), "studies": len(a["studies"]),
            "total_gbase": round(a["bases"] / 1e9, 1),
            "est_download_tb": round(a["bases"] / 3.5 / 1e12, 4),
            "paired_runs": a["paired"], "single_runs": a["single"],
            "first_year": min(a["years"]) if a["years"] else "",
            "last_year": max(a["years"]) if a["years"] else "",
        })
    out.sort(key=lambda d: -d["isolates"])
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ena", nargs="+", required=True, help="ENA run snapshot(s), TSV")
    p.add_argument("--taxmap", default=os.path.join(HERE, "catalog", "taxonomy_species_map.csv.gz"))
    p.add_argument("--out", default=os.path.join(HERE, "catalog", "kingdom",
                                                  "fungi_species_summary.csv.gz"))
    a = p.parse_args()
    taxmap = load_taxmap(a.taxmap)

    def rows():
        for path in a.ena:
            with opener(path) as fh:
                yield from csv.DictReader(fh, delimiter="\t")

    out = summarise(rows(), taxmap)
    with gzip.open(a.out, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out)
    print(f"wrote {a.out}: {len(out)} species", file=sys.stderr)


if __name__ == "__main__":
    main()
