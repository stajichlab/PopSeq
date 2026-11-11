#!/usr/bin/env python3
"""Build a many-species panel: up to N isolates for every species with enough isolates.

Writes panels/<name>/ per spec 003: panel.yaml, samples.csv, metadata.csv. Issue #36:
GmG builds PCoA and PCA pages from read sketches of this panel (owner decision
2026-11-09: at least 10 paired isolates, at most 300 per species).

Selection, per species, is popseq.panel() with the default filters: Illumina, PAIRED by
ena_layout (what ENA serves), WGS, exclusion labels dropped, FASTQ required. That is one
run per isolate, the deepest run wins, then round-robin across studies. A species is
kept when at least --min-isolates isolates pass those filters. Candidates come from
catalog/species_layout.csv.gz (isolates_paired), which counts every platform, so it is
a superset of what passes.

samples.csv holds ENA URLs, not local paths: nothing has been fetched yet. Rewrite the
FASTQ columns after fetching (popseq fetch, issue #13).

    uv run python scripts/popseq_panel_set.py                # pop_sra_ge10_cap300
    uv run python scripts/popseq_panel_set.py --min-isolates 30 --n 100 --name pop_sra_ge30_cap100
"""
from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

import yaml  # noqa: E402

import popseq  # noqa: E402
from popseq import query as q  # noqa: E402

FILTERS = dict(platform="ILLUMINA", layout="PAIRED", strategy="WGS",
               include_excluded=False, require_fastq=True)
METADATA = ["sample", "biosample", "taxid", "scientific_name", "species_taxid",
            "species_name", "run_accession", "study_accession", "strain", "isolate",
            "country", "country_iso", "country_status", "host", "host_name",
            "collection_date", "first_public", "instrument_model", "read_count",
            "base_count", "est_depth_x", "ncbi_status"]
SAFE_ID = re.compile(r"^[A-Za-z0-9_-][A-Za-z0-9_.-]*$")


def candidates(layout_rows, min_isolates: int) -> list[int]:
    """Species taxids with at least min_isolates paired isolates, largest first."""
    keep = [r for r in layout_rows if int(r["isolates_paired"] or 0) >= min_isolates]
    keep.sort(key=lambda r: (-int(r["isolates_paired"]), int(r["species_taxid"])))
    return [int(r["species_taxid"]) for r in keep]


def sample_id(run: dict) -> str:
    """The BioSample accession, or the run accession when that is not a usable id
    (empty, or several samples joined with ';')."""
    s = (run.get("sample_accession") or "").strip()
    return s if s and SAFE_ID.match(s) else run["run_accession"]


def rows_for(panel: dict, full: dict[str, dict]) -> tuple[list[dict], list[dict]]:
    """samples.csv and metadata.csv rows for one species' panel() result.

    `full` maps run_accession to the complete run-table row, for the columns panel()
    does not return (scientific_name, tax_id, isolate, instrument_model)."""
    sheet = q.to_sample_sheet(panel)
    samples, meta = [], []
    for run, s in zip(panel["runs"], sheet):
        f = full.get(run["run_accession"], {})
        sid = sample_id(run)
        samples.append({"sample": sid, "fastq_1": s["fastq_1"], "fastq_2": s["fastq_2"]})
        m = {k: run.get(k, "") or f.get(k, "") for k in METADATA}
        m.update(sample=sid, biosample=run.get("sample_accession") or "",
                 taxid=f.get("tax_id", ""), scientific_name=f.get("scientific_name", ""),
                 species_taxid=panel["species_taxid"], species_name=panel["species_name"])
        meta.append(m)
    return samples, meta


def catalog_commit() -> str:
    try:
        return subprocess.run(["git", "-C", HERE, "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--min-isolates", type=int, default=10)
    p.add_argument("--n", type=int, default=300, help="isolates per species at most")
    p.add_argument("--seed", type=int, default=0, help="study order shuffle, as panel()")
    p.add_argument("--name", help="default: pop_sra_ge<min>_cap<n>")
    p.add_argument("--layout-table",
                   default=os.path.join(HERE, "catalog", "species_layout.csv.gz"))
    a = p.parse_args()
    name = a.name or f"pop_sra_ge{a.min_isolates}_cap{a.n}"
    out = os.path.join(HERE, "panels", name)

    with gzip.open(a.layout_table, "rt", newline="") as fh:
        taxids = candidates(csv.DictReader(fh), a.min_isolates)
    held = popseq.preload(taxids)
    print(f"{len(taxids)} candidate species, {held:,} runs held", file=sys.stderr)

    samples, meta, per_species, dropped = [], [], [], 0
    for tx in taxids:
        try:
            pan = popseq.panel(taxid=tx, n=a.n, seed=a.seed, **FILTERS)
        except popseq.EmptyResult:
            continue
        if pan["isolates_available"] < a.min_isolates:
            continue
        full = {r["run_accession"]: r for r in q._species_runs(tx)}
        s, m = rows_for(pan, full)
        samples += s
        meta += m
        dropped += q.sample_sheet_dropped(pan)
        per_species.append({"species_taxid": tx, "species_name": pan["species_name"],
                            "selected": pan["selected"],
                            "isolates_available": pan["isolates_available"],
                            "studies_in_panel": pan["studies_in_panel"]})

    ids = [s["sample"] for s in samples]
    if len(ids) != len(set(ids)):
        raise SystemExit("duplicate sample ids across species; refusing to write")
    no_depth = sum(1 for m in meta if not m["est_depth_x"])

    os.makedirs(out, exist_ok=True)
    for fname, rows, cols in (("samples.csv", samples, ["sample", "fastq_1", "fastq_2"]),
                              ("metadata.csv", meta, METADATA)):
        with open(os.path.join(out, fname), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    doc = {
        "name": name,
        "title": f"Up to {a.n} isolates of every fungal species with at least "
                 f"{a.min_isolates} paired Illumina isolates",
        "taxid": "many: see species_taxid in metadata.csv",
        "created": datetime.date.today().isoformat(),
        "created_by": "scripts/popseq_panel_set.py",
        "question": ("Population structure per species from read sketches: GmG static "
                     "PCoA and PCA pages (PopSeq issue #36). Each species is a separate "
                     "population; do not pool species in one ordination."),
        "selection": {
            "source_catalog_commit": catalog_commit(),
            "criteria": (f"popseq.panel(n={a.n}, seed={a.seed}) with filters {FILTERS}; "
                         "one run per isolate, deepest wins, round-robin across studies; "
                         f"species kept when >= {a.min_isolates} isolates pass the filters"),
            "rule_caveat": "a sampling rule, not a guarantee of independence",
            "n_species": len(per_species),
            "n_samples": len(samples),
            "n_runs": len(samples),
            "samples_without_depth_estimate": no_depth,
            "runs_with_an_unpaired_third_file": dropped,
        },
        "fastq_paths": ("ENA HTTPS URLs, not local paths: nothing has been fetched. "
                        "Rewrite after fetching (issue #13)."),
        "depth_tier": "not fetched",
        "species": per_species,
    }
    with open(os.path.join(out, "panel.yaml"), "w") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, allow_unicode=True, width=88)
    print(f"wrote {out}: {len(per_species)} species, {len(samples)} samples, "
          f"{no_depth} without a depth estimate", file=sys.stderr)


if __name__ == "__main__":
    main()
