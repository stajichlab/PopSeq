#!/usr/bin/env python3
"""Build catalog/summary/unplaced_runs.csv.gz: every run with no species taxid (issue #22).

Issue #1 keyed the catalog on the species-rank taxid. A run whose `tax_id` does not
resolve to a species gets an empty `species_taxid` in the run table, so `popseq species`,
`runs`, `panel` and every per-species table miss it. This table lists those runs, one
row each, so none is dropped without a record.

The decision (issue #22): list them in their own table. Do NOT roll them up to a ranked
ancestor. A run filed as "Aspergillus" or "Fungi" is not evidence for any one species,
and putting it in a species row would make exactly the claim the taxid rekeying removed.

Reads only committed files, no network:

    catalog/runs/fungi_all_runs.csv.gz      the run table (species_taxid empty = unplaced)
    catalog/taxonomy_species_map.csv.gz     taxid -> name, rank, species ancestor, resolution

Columns

    run_accession, study_accession,   from the run table, unchanged
    sample_accession, scientific_name,
    tax_id, instrument_platform,
    first_public
    placement            why the run has no species:
                           above_species         tax_id is a ranked taxon above species
                           multiple_taxids       tax_id is a ';'-joined list (a multi-sample run)
                           not_in_taxonomy_map   tax_id is not in the fungal taxonomy map
                           unresolved_in_map     the map holds tax_id but found no ranked ancestor
                           no_taxid              tax_id is empty
                           species_in_map        the map now gives a species: the run table
                                                 is older than the map and must be rebuilt
    resolved_rank        the rank tax_id resolved to (genus, section, kingdom, ...), or
                         `unresolved` when it did not resolve
    ranked_taxid,        the nearest ranked taxon at or above tax_id, from the map. Blank
    ranked_name          when the map does not record it: the map has no parent column, so
                         an unranked tax_id whose walk stopped at a ranked ancestor cannot
                         name that ancestor (none in the 2026-11-18 catalog)
    part_species_taxids  for multiple_taxids only: the distinct species taxids the parts
                         resolve to, ';'-joined. Blank when no part resolves

The output is gzip with mtime 0, so a rebuild from the same inputs is byte-identical.

    uv run python scripts/popseq_unplaced.py
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(HERE, "catalog", "runs", "fungi_all_runs.csv.gz")
TAXMAP = os.path.join(HERE, "catalog", "taxonomy_species_map.csv.gz")
OUT = os.path.join(HERE, "catalog", "summary", "unplaced_runs.csv.gz")

RUN_COLUMNS = ["run_accession", "study_accession", "sample_accession", "scientific_name",
               "tax_id", "instrument_platform", "first_public"]
FIELDS = RUN_COLUMNS + ["placement", "resolved_rank", "ranked_taxid", "ranked_name",
                        "part_species_taxids"]
PLACEMENTS = ["above_species", "multiple_taxids", "not_in_taxonomy_map",
              "unresolved_in_map", "no_taxid", "species_in_map"]


def opener(path, mode="rt"):
    return gzip.open(path, mode, newline="") if path.endswith(".gz") else open(path, mode, newline="")


def load_taxmap(path) -> dict[str, dict]:
    with opener(path) as fh:
        return {r["taxid"]: r for r in csv.DictReader(fh)}


def classify(tax_id: str, taxmap: dict[str, dict]) -> dict:
    """Why a tax_id gives no species. Returns the derived columns of one row."""
    out = dict(placement="", resolved_rank="unresolved", ranked_taxid="", ranked_name="",
               part_species_taxids="")
    t = (tax_id or "").strip()
    if not t:
        out["placement"] = "no_taxid"
        return out
    if ";" in t:
        out["placement"] = "multiple_taxids"
        parts = {p.strip() for p in t.split(";") if p.strip()}
        sp = {taxmap[p]["species_taxid"] for p in parts
              if p in taxmap and taxmap[p].get("species_taxid")}
        out["part_species_taxids"] = ";".join(sorted(sp, key=int))
        return out
    m = taxmap.get(t)
    if m is None:
        out["placement"] = "not_in_taxonomy_map"
        return out
    if m.get("species_taxid"):
        out["placement"] = "species_in_map"
        out["resolved_rank"] = "species"
        out["ranked_taxid"], out["ranked_name"] = m["species_taxid"], m.get("species_name", "")
        return out
    res = m.get("resolution") or ""
    if res.startswith("above_species:"):
        rank = res.split(":", 1)[1]
        out["placement"] = "above_species"
        out["resolved_rank"] = rank
        if m.get("rank") == rank:          # tax_id is itself the ranked taxon
            out["ranked_taxid"], out["ranked_name"] = t, m.get("name", "")
        return out
    out["placement"] = "unresolved_in_map"
    return out


def unplaced(run_rows, taxmap: dict[str, dict]) -> list[dict]:
    """One row per run with an empty species_taxid, sorted for a stable diff."""
    rows = []
    for r in run_rows:
        if (r.get("species_taxid") or "").strip():
            continue
        row = {k: r.get(k, "") for k in RUN_COLUMNS}
        row.update(classify(r.get("tax_id", ""), taxmap))
        rows.append(row)
    order = {p: i for i, p in enumerate(PLACEMENTS)}
    rows.sort(key=lambda x: (order.get(x["placement"], 99), x["resolved_rank"],
                             x["run_accession"]))
    return rows


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while (b := fh.read(1 << 20)):
            h.update(b)
    return h.hexdigest()


def write(rows: list[dict], out: str, provenance: str) -> None:
    os.makedirs(os.path.dirname(out), exist_ok=True)
    buf = io.StringIO()
    buf.write(provenance + "\n")
    w = csv.DictWriter(buf, fieldnames=FIELDS, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    tmp = out + ".tmp"
    with open(tmp, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw,
                                               mtime=0) as gz:
        gz.write(buf.getvalue().encode())
    os.replace(tmp, out)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs", default=RUNS)
    p.add_argument("--taxmap", default=TAXMAP)
    p.add_argument("--out", default=OUT)
    a = p.parse_args(argv)

    taxmap = load_taxmap(a.taxmap)
    total = 0

    def run_rows():
        nonlocal total
        with opener(a.runs) as fh:
            for r in csv.DictReader(fh):
                total += 1
                yield r

    rows = unplaced(run_rows(), taxmap)
    prov = (f"# generated by scripts/popseq_unplaced.py from {os.path.relpath(a.runs, HERE)} "
            f"(sha256 {sha256(a.runs)[:12]}) and {os.path.relpath(a.taxmap, HERE)} "
            f"(sha256 {sha256(a.taxmap)[:12]})")
    write(rows, a.out, prov)
    print(f"wrote {a.out}: {len(rows)} of {total} runs have no species taxid", file=sys.stderr)
    for k, v in Counter(r["placement"] for r in rows).most_common():
        print(f"  {k:<22} {v}", file=sys.stderr)
    print("  by resolved_rank: " + ", ".join(
        f"{k} {v}" for k, v in Counter(r["resolved_rank"] for r in rows).most_common()),
        file=sys.stderr)
    print("  by platform: " + ", ".join(
        f"{k or '(empty)'} {v}" for k, v in
        Counter(r["instrument_platform"] for r in rows).most_common()), file=sys.stderr)
    no_anc = sum(1 for r in rows if r["placement"] == "above_species" and not r["ranked_taxid"])
    print(f"  above_species rows with no ranked taxon in the map: {no_anc}", file=sys.stderr)
    stale = sum(1 for r in rows if r["placement"] == "species_in_map")
    if stale:
        print(f"WARNING: {stale} runs resolve to a species in the map but not in the run "
              "table; rebuild the run table", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
