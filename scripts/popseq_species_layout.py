#!/usr/bin/env python3
"""Build catalog/species_layout.csv.gz: which species can be worked up from paired reads,
and which only from single-end reads.

Derived from the canonical run table, catalog/runs/fungi_all_runs.csv.gz (AGENTS.md §7).
Pairing is ena_layout, what ENA serves, not the submitter's library_layout (spec 001).
Rows with an exclusion label (mutant, evolution, cross, not_wgs) are left out, as in the
default panels. An isolate is a BioSample; a ';'-joined sample_accession is split.

Columns:
    species_taxid, species_name
    isolates                  isolates with at least one run
    isolates_paired           isolates with at least one PAIRED run
    isolates_single_only      isolates whose runs are all SINGLE: usable only as single-end
    isolates_no_fastq         isolates for which ENA serves no FASTQ at all
    single_only_median_depth_x   median est_depth_x over single-only isolates (the deepest
                              run per isolate); blank when no estimate exists
    single_only_ge_10x        single-only isolates whose deepest run is 10x or more
    runs_label_mismatch       runs whose library_layout disagrees with ena_layout

NCBI can hold a run unpaired that ENA serves paired; that is not visible here (spec 002).

    uv run python scripts/popseq_species_layout.py
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import statistics
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIELDS = ["species_taxid", "species_name", "isolates", "isolates_paired",
          "isolates_single_only", "isolates_no_fastq", "single_only_median_depth_x",
          "single_only_ge_10x", "runs_label_mismatch"]


def build(rows):
    lay = defaultdict(lambda: defaultdict(set))   # sp -> isolate -> {PAIRED, SINGLE, ""}
    depth = defaultdict(lambda: defaultdict(float))
    mismatch = defaultdict(int)
    names = {}
    for r in rows:
        sp = r.get("species_taxid") or ""
        if not sp or r.get("exclusion"):
            continue
        names[sp] = r.get("species_name") or ""
        served = r.get("ena_layout") or ""
        if served and (r.get("library_layout") or "") != served:
            mismatch[sp] += 1
        try:
            d = float(r.get("est_depth_x") or "")
        except ValueError:
            d = None
        for s in (x.strip() for x in (r.get("sample_accession") or "").split(";")):
            if not s:
                continue
            lay[sp][s].add(served)
            if d is not None and d > depth[sp][s]:
                depth[sp][s] = d
    out = []
    for sp, isolates in lay.items():
        paired = [s for s, v in isolates.items() if "PAIRED" in v]
        single = [s for s, v in isolates.items() if "SINGLE" in v and "PAIRED" not in v]
        none = [s for s, v in isolates.items() if v == {""}]
        dep = [depth[sp][s] for s in single if s in depth[sp]]
        out.append({
            "species_taxid": sp, "species_name": names[sp], "isolates": len(isolates),
            "isolates_paired": len(paired), "isolates_single_only": len(single),
            "isolates_no_fastq": len(none),
            "single_only_median_depth_x": f"{statistics.median(dep):.1f}" if dep else "",
            "single_only_ge_10x": sum(1 for x in dep if x >= 10),
            "runs_label_mismatch": mismatch[sp],
        })
    out.sort(key=lambda d: (-d["isolates_single_only"], -d["isolates"], int(d["species_taxid"])))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs", default=os.path.join(HERE, "catalog", "runs", "fungi_all_runs.csv.gz"))
    p.add_argument("--out", default=os.path.join(HERE, "catalog", "species_layout.csv.gz"))
    a = p.parse_args()
    with gzip.open(a.runs, "rt", newline="") as fh:
        out = build(csv.DictReader(fh))
    with gzip.open(a.out, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    only = [d for d in out if d["isolates_single_only"] and not d["isolates_paired"]]
    print(f"wrote {a.out}: {len(out)} species; {len(only)} with single-end isolates only; "
          f"{sum(d['isolates_single_only'] for d in out)} single-only isolates", file=sys.stderr)


if __name__ == "__main__":
    main()
