#!/usr/bin/env python3
"""Render catalog/api/: one static JSON file per species, plus an index.

For a tool that cannot clone this repository or run Python. Each file is fetchable
over raw.githubusercontent.com with nothing but an HTTP client:

    https://raw.githubusercontent.com/stajichlab/PopSeq/main/catalog/api/index.json
    https://raw.githubusercontent.com/stajichlab/PopSeq/main/catalog/api/746128.json

A species file carries the full `describe` payload plus a default panel: up to 200
strains, one run per isolate, spread across studies, paired Illumina WGS with the
mutant/evolved/cross/non-WGS rows dropped. That is the common request answered in
advance. Anything else needs the CLI, the package or the MCP server, because the
filter space is too large to pre-render.

Only species with 100 or more isolates get a file. Below that a species is a
handful of genomes, not a population, and 6,373 nearly-empty files would make the
directory useless.

    uv run python scripts/popseq_api.py
    uv run python scripts/popseq_api.py --out /some/dir   # render elsewhere

tests/test_api_fresh.py renders into a temporary directory and fails when the
result differs from the committed catalog/api/ (issue #21). After any change to
catalog/runs/, catalog/summary/, catalog/qc/, catalog/labels/ or the query code,
rerun this script and commit catalog/api/ with the change.
"""
from __future__ import annotations

import argparse
import json
import os
from datetime import date

import popseq

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "catalog", "api")
RAW = "https://raw.githubusercontent.com/stajichlab/PopSeq/main/catalog/api"

PANEL_N = 100

# A compact projection. species_taxid and species_name are on the parent document,
# and fastq_aspera is the same path as fastq_ftp against the fasp endpoint, so
# repeating either on 100 runs x 113 species only makes the files slower to fetch.
PANEL_RUN_FIELDS = ["run_accession", "sample_accession", "study_accession",
                    "strain", "read_count", "base_count", "est_depth_x",
                    "country", "collection_date", "first_public", "ncbi_status",
                    "fastq_ftp"]
PANEL_FILTERS = dict(platform="ILLUMINA", layout="PAIRED", strategy="WGS",
                     include_excluded=False, require_fastq=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--min-isolates", type=int, default=100)
    p.add_argument("--panel-n", type=int, default=PANEL_N)
    p.add_argument("--out", default=OUT, help=f"output directory (default {OUT})")
    a = p.parse_args(argv)
    out = a.out

    os.makedirs(out, exist_ok=True)
    sp = popseq.species(min_isolates=a.min_isolates)
    taxids = [int(r["species_taxid"]) for r in sp]
    held = popseq.preload(taxids)
    print(f"{len(sp)} species, {held:,} runs held in memory")

    listing, written, bytes_total = [], 0, 0
    for r in sp:
        tx = int(r["species_taxid"])
        doc = popseq.describe(taxid=tx)
        try:
            pan = popseq.panel(taxid=tx, n=a.panel_n, **PANEL_FILTERS)
        except popseq.EmptyResult:
            pan = None
        doc["default_panel"] = {
            "filters": {**PANEL_FILTERS, "n": a.panel_n},
            "selection_rule": pan["selection_rule"] if pan else None,
            "selected": pan["selected"] if pan else 0,
            "isolates_available": pan["isolates_available"] if pan else 0,
            "studies_in_panel": pan["studies_in_panel"] if pan else 0,
            "runs": [{k: r[k] for k in PANEL_RUN_FIELDS} for r in pan["runs"]] if pan else [],
            "run_fields": PANEL_RUN_FIELDS,
            "note": ("A pre-rendered answer to the common request. For any other "
                     "filter use the CLI, the popseq package or the MCP server. "
                     "fastq_ftp has no scheme: prepend https://. The Aspera path "
                     "is the same path against fasp.sra.ebi.ac.uk."),
        }
        doc["self"] = f"{RAW}/{tx}.json"
        path = os.path.join(out, f"{tx}.json")
        with open(path, "w") as fh:
            json.dump(doc, fh, indent=1, default=str)
            fh.write("\n")
        size = os.path.getsize(path)
        bytes_total += size
        written += 1
        listing.append({
            "species_taxid": tx, "species_name": r["species_name"],
            "isolates": int(r["isolates"]), "runs": int(r["runs"]),
            "studies": int(r["studies"]),
            "panel_available": doc["default_panel"]["selected"],
            "url": f"{RAW}/{tx}.json", "bytes": size,
        })

    idx = {
        "popseq_schema": popseq.SCHEMA_VERSION,
        "generated": date.today().isoformat(),
        "what": ("One static JSON file per fungal species with 100 or more "
                 "isolates. Each holds the full description of the species and a "
                 "pre-selected panel of up to "
                 f"{a.panel_n} strains to profile."),
        "scope": popseq.index()["scope"],
        "archive_snapshot": popseq.index()["archive_snapshot"],
        "richer_queries": popseq.index()["query_interface"],
        "caveats": popseq.index()["caveats"],
        "species_count": len(listing),
        "species": sorted(listing, key=lambda d: -d["isolates"]),
    }
    with open(os.path.join(out, "index.json"), "w") as fh:
        json.dump(idx, fh, indent=1)
        fh.write("\n")
    print(f"wrote {written} species files + index.json, "
          f"{bytes_total / 1e6:.1f} MB total")


if __name__ == "__main__":
    main()
