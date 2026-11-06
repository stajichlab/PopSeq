#!/usr/bin/env python3
"""Fetch reference genome sizes from NCBI for the species in the catalog.

Every coverage number depends on a genome size, so a guessed size silently biases
every depth estimate for that species. This fetches the NCBI reference assembly and
records its accession, so the figure is traceable.

A species with no reference assembly gets an empty size and a reason, never a guess.
Issue #16 covers inheriting a relative's size for those, explicitly flagged.

Resumable: an existing row is kept unless --update.

Usage:
    uv run python scripts/popseq_genome_sizes.py --min-isolates 100
"""

import argparse
import csv
import gzip
import json
import os
import sys
import time
import urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon"
OUT = "catalog/genome_sizes.csv"
COLS = ["species_taxid", "species_name", "genome_mb", "accession", "assembly_level", "note"]


def get(url, tries=3):
    key = os.environ.get("NCBI_API_KEY")
    if key:
        url += ("&" if "?" in url else "?") + "api_key=" + key
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PopSeq/0.1"})
            with urllib.request.urlopen(req, timeout=120) as fh:
                return fh.read().decode("utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            time.sleep(2 * (attempt + 1))
    return None


def reference(taxid):
    for flt in ("&filters.reference_only=true", "&filters.assembly_level=chromosome", ""):
        body = get(f"{API}/{taxid}/dataset_report?page_size=1{flt}")
        if not body:
            continue
        try:
            reports = json.loads(body).get("reports") or []
        except json.JSONDecodeError:
            continue
        if not reports:
            continue
        r = reports[0]
        size = r.get("assembly_stats", {}).get("total_sequence_length")
        if size:
            return (round(int(size) / 1e6, 2), r.get("accession"),
                    r.get("assembly_info", {}).get("assembly_level"), "")
    return (None, None, None, "no reference assembly at NCBI")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", default="catalog/kingdom/fungi_species_summary.csv.gz")
    p.add_argument("--min-isolates", type=int, default=100)
    p.add_argument("--out", default=OUT)
    p.add_argument("--update", action="store_true")
    args = p.parse_args()

    have = {}
    if os.path.exists(args.out):
        with open(args.out, newline="") as fh:
            have = {r["species_taxid"]: r for r in csv.DictReader(fh)}
        print(f"{len(have)} already recorded", file=sys.stderr)

    with gzip.open(args.summary, "rt", newline="") as fh:
        spp = [r for r in csv.DictReader(fh)
               if int(r["isolates"]) >= args.min_isolates]
    todo = [r for r in spp if args.update or r["species_taxid"] not in have]
    print(f"{len(spp)} species at/above {args.min_isolates} isolates, "
          f"{len(todo)} to fetch", file=sys.stderr)

    for n, r in enumerate(todo, 1):
        mb, acc, lvl, note = reference(int(r["species_taxid"]))
        have[r["species_taxid"]] = {
            "species_taxid": r["species_taxid"], "species_name": r["species_name"],
            "genome_mb": mb if mb else "", "accession": acc or "",
            "assembly_level": lvl or "", "note": note,
        }
        if n % 20 == 0 or note:
            print(f"  [{n}/{len(todo)}] {r['species_name']}: "
                  f"{mb or note}", file=sys.stderr)
        time.sleep(0.11)
        if n % 25 == 0:
            with open(args.out, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=COLS)
                w.writeheader()
                w.writerows(have.values())

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(sorted(have.values(), key=lambda x: x["species_name"]))
    ok = sum(1 for v in have.values() if v["genome_mb"])
    print(f"\nWrote {args.out}: {len(have)} species, {ok} with a genome size, "
          f"{len(have) - ok} without", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
