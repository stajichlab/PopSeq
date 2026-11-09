#!/usr/bin/env python3
"""Write catalog/index.json: the one file another tool reads first.

A consumer should not have to know this repository's layout. It reads the index,
learns the schema version, where every table is, what each column means, and what
the catalog currently holds. Everything else follows from that.

    uv run python scripts/popseq_index.py
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from datetime import date

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "catalog", "index.json")
SCHEMA_VERSION = "1.0.0"
SNAPSHOT = "2026-11-06"
REPO = "https://github.com/stajichlab/PopSeq"
SHARED_PATH = "/bigdata/stajichlab/shared/projects/Fungi/PopSeq"

TABLES = {
    "runs": dict(
        path="catalog/runs/fungi_all_runs.csv.gz",
        grain="one row per sequencing run",
        key="run_accession",
        built_by="scripts/popseq_runs_table.py"),
    "species_summary": dict(
        path="catalog/kingdom/fungi_species_summary.csv.gz",
        grain="one row per species taxid",
        key="species_taxid",
        built_by="scripts/popseq_summarize.py"),
    "biosamples": dict(
        path="catalog/kingdom/fungi_biosamples.csv.gz",
        grain="one row per BioSample, 783 sparse attribute columns",
        key="sample_accession",
        built_by="scripts/popseq_biosamples.py"),
    "taxonomy_species_map": dict(
        path="catalog/taxonomy_species_map.csv.gz",
        grain="one row per fungal taxid, with its species-rank ancestor",
        key="taxid",
        built_by="scripts/popseq_taxonomy.py"),
    "genome_sizes": dict(
        path="catalog/genome_sizes.csv",
        grain="one row per species with a resolved NCBI reference assembly",
        key="species_taxid",
        built_by="scripts/popseq_genome_sizes.py"),
    "summary_by_year": dict(path="catalog/summary/by_year.csv.gz",
                            grain="species x year", key="species_taxid,year",
                            built_by="scripts/popseq_summarize.py"),
    "summary_clonality": dict(path="catalog/summary/clonality.csv.gz",
                              grain="one row per species", key="species_taxid",
                              built_by="scripts/popseq_summarize.py"),
    "summary_coverage": dict(path="catalog/summary/coverage.csv.gz",
                             grain="one row per species with a genome size",
                             key="species_taxid",
                             built_by="scripts/popseq_summarize.py"),
    "summary_geography": dict(path="catalog/summary/geography.csv.gz",
                              grain="species x country", key="species_taxid,country",
                              built_by="scripts/popseq_summarize.py"),
    "summary_host_source": dict(path="catalog/summary/host_source.csv.gz",
                                grain="species x field x value",
                                key="species_taxid,field,value",
                                built_by="scripts/popseq_summarize.py"),
    "vocab_country_aliases": dict(path="catalog/vocab/country_aliases.csv",
                                  grain="one row per curated country spelling",
                                  key="alias", built_by="curated by hand"),
    "vocab_host_taxa": dict(path="catalog/vocab/host_taxa.csv.gz",
                            grain="one row per distinct normalised host string",
                            key="host_key", built_by="scripts/popseq_vocab.py hosts"),
    "vocab_host_overrides": dict(path="catalog/vocab/host_overrides.csv",
                                 grain="one row per curated host string",
                                 key="host_key", built_by="curated by hand"),
    "ena_only_runs": dict(path="catalog/kingdom/fungi_byspecies.ena_only_runs.csv.gz",
                          grain="runs present at ENA and absent from the filtered "
                                "NCBI query", key="run_accession",
                          built_by="scripts/popseq_reconcile.py"),
    "ncbi_only_runs": dict(path="catalog/kingdom/fungi_byspecies.ncbi_only_runs.csv.gz",
                           grain="runs present at NCBI and absent from ENA",
                           key="run_accession",
                           built_by="scripts/popseq_reconcile.py"),
}

COLUMNS = {
    "species_taxid": "NCBI species-rank taxid. The join key. Never join on a name: "
                     "NCBI labels the same taxon with more than one name.",
    "species_name": "The name of species_taxid in the NCBI taxonomy, not the name "
                    "the submitter deposited.",
    "scientific_name": "The name the submitter deposited. May be a synonym, a "
                       "strain-level name, or a rank above species.",
    "exclusion": "'' or mutant, evolution, cross, not_wgs. A LABEL from the study "
                 "and sample titles. Nothing is deleted. Filter it yourself.",
    "est_depth_x": "base_count divided by the NCBI reference assembly size for the "
                   "species. An estimate from metadata, not a measurement on reads. "
                   "Blank when no reference size resolved.",
    "ncbi_status": "in_ncbi, ena_only, or not_checked. Only the 113 species with "
                   "100 or more isolates were reconciled against NCBI.",
    "fastq_ftp": "';'-joined ENA FTP paths, no scheme. Prepend https://.",
    "fastq_aspera": "';'-joined Aspera paths for the ENA fasp endpoint.",
    "country": "Exactly as the submitter wrote it. Free text. Use country_iso to group.",
    "country_iso": "ISO 3166-1 alpha-2 (XK for Kosovo). Blank unless country_status "
                   "is ok. Parsed from country by popseq/normalise.py (issue #23).",
    "country_status": "ok, absent (field empty), missing (an explicit INSDC missing "
                      "term such as 'not collected'), region (ocean, continent), "
                      "historical (e.g. Yugoslavia), conflicting (multi-sample run "
                      "whose samples disagree), unresolved (add it to "
                      "catalog/vocab/country_aliases.csv).",
    "host": "Raw host from the BioSample: first non-empty of host, host_plant, "
            "specific hosts, host_common_name.",
    "host_taxid": "NCBI taxid of host. Blank unless host_status is ok or "
                  "same_as_organism. Resolved via catalog/vocab/host_taxa.csv.gz.",
    "host_status": "ok, absent, missing, not_a_host (a lab or environment term such "
                   "as 'soil' or 'culture'), same_as_organism (the host field names "
                   "the sequenced organism itself), ambiguous, conflicting, "
                   "unresolved.",
    "sample_accession": "BioSample. In the raw ENA export this field is sometimes "
                        "';'-joined for multi-sample runs (issue #3).",
}

CAVEATS = [
    "Every count is a floor. ENA holds runs with an empty tax_id that no taxonomy "
    "query reaches, so the true per-species totals are higher (issue #2).",
    "A run count is not a strain count and an isolate count is not an "
    "independent-sample count. 52 of the 113 largest species take half or more of "
    "their isolates from one study. Read catalog/summary/clonality.csv.gz before "
    "treating a species as a population.",
    "est_depth_x comes from metadata, not from the reads.",
    "Only 113 species were reconciled against NCBI. For every other species "
    "ncbi_status is not_checked, which is not the same as absent.",
    "Nothing here has been downloaded or verified. These are archive records.",
    "Scope: the run table holds only ENA runs with library_strategy=WGS, "
    "library_source=GENOMIC and instrument_platform=ILLUMINA "
    "(scripts/popseq_survey.py, QUERY_RUNS). A --strategy or --platform filter for "
    "anything else returns nothing because those runs were never harvested, not because "
    "none exist (issue #24).",
]


def sha(path, n=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while (b := fh.read(n)):
            h.update(b)
    return h.hexdigest()


def main():
    rows = []
    p = os.path.join(HERE, "catalog", "kingdom", "fungi_species_summary.csv.gz")
    with gzip.open(p, "rt", newline="") as fh:
        if not fh.readline().startswith("#"):
            fh.seek(0)
        rows = list(csv.DictReader(fh))

    tables = {}
    for name, meta in TABLES.items():
        full = os.path.join(HERE, meta["path"])
        if not os.path.exists(full):
            continue
        with (gzip.open(full, "rt", newline="") if full.endswith(".gz")
              else open(full, newline="")) as fh:
            first = fh.readline()
            if first.startswith("#"):
                first = fh.readline()
            cols = next(csv.reader([first]))
        tables[name] = {**meta, "columns": cols,
                        "bytes": os.path.getsize(full), "sha256": sha(full)}

    index = {
        "popseq_schema": SCHEMA_VERSION,
        "generated": date.today().isoformat(),
        "archive_snapshot": SNAPSHOT,
        "scope": "public whole-genome Illumina sequencing runs of fungal isolates, "
                 "as ENA and NCBI describe them",
        "repository": REPO,
        "shared_filesystem_path": SHARED_PATH,
        "query_interface": {
            "command": "uv run popseq",
            "python": "import popseq; popseq.panel(taxid=746128, n=100)",
            "mcp": "uv run --extra mcp python -m popseq.mcp_server (stdio)",
            "static_json": "catalog/api/index.json and catalog/api/<taxid>.json, "
                           "fetchable over raw.githubusercontent.com with no clone",
            "legacy_command": "uv run python scripts/popseq_query.py (shim)",
            "verbs": {
                "species": "list species with counts; the discovery entry point",
                "runs": "all runs of one species, filtered",
                "panel": "select N strains, one per isolate, spread across studies",
                "describe": "everything known about one species, as JSON",
            },
            "exit_codes": {"0": "ok", "2": "species not found or ambiguous",
                           "3": "filters left nothing"},
            "stable": "Verb names, flag names and the JSON key names are the "
                      "contract. Table paths and columns may change underneath.",
        },
        "totals": {
            "runs": sum(int(r["runs"]) for r in rows),
            "isolates": sum(int(r["isolates"]) for r in rows),
            "species_taxids": len(rows),
            "species_with_100_isolates": sum(1 for r in rows
                                             if int(r["isolates"]) >= 100),
            "estimated_tb_one_run_per_isolate":
                round(sum(float(r["est_download_tb"]) for r in rows), 1),
        },
        "tables": tables,
        "column_notes": COLUMNS,
        "caveats": CAVEATS,
    }
    with open(OUT, "w") as fh:
        json.dump(index, fh, indent=2)
        fh.write("\n")
    print(f"wrote {OUT}: {len(tables)} tables, "
          f"{index['totals']['runs']:,} runs")


if __name__ == "__main__":
    main()
