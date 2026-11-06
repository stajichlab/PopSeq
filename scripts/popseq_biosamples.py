#!/usr/bin/env python3
"""Harvest NCBI BioSample attributes into the PopSeq catalog.

A rewrite of ~/projects/biosample_metadata/scripts/biosample2table.py for catalog
scale. The attribute-flattening logic and the incremental-update behaviour of that
script are correct and are kept. Three things are changed:

1. Batched requests. The original calls esearch then efetch once per accession. At
   the unauthenticated limit of 3 requests per second, 10,000 BioSamples take about
   two hours. epost plus one efetch handles BATCH accessions per round trip.
2. API key support. NCBI_API_KEY in the environment raises the limit from 3 to 10
   requests per second. The key is read from the environment and never written out.
3. Resume. Entrez fails intermittently. Results are flushed after every batch, so a
   harvest that dies at 9,000 of 10,000 resumes rather than starting over.

Output is a wide, sparse CSV: one row per BioSample, one column per attribute seen,
using NCBI's harmonized_name where one exists. Empty means unknown; NA is never
written (spec 001).

Usage:
    uv run python scripts/popseq_biosamples.py -e you@ucr.edu \
        --input biosamples.txt --output catalog/biosamples/109871.Bd.biosample.csv.gz

    uv run python scripts/popseq_biosamples.py -e you@ucr.edu \
        --from-runs catalog/runs/109871.Bd.runs.csv.gz \
        --output catalog/biosamples/109871.Bd.biosample.csv.gz
"""

import argparse
import csv
import gzip
import io
import os
import sys
import time
import xml.etree.ElementTree as ET

from Bio import Entrez

BATCH = 100  # an OR query of 100 accessions is about 1.8 kB, well inside URL limits
# Floor between requests. 10/s with a key, 3/s without. Stay just inside both.
SLEEP_WITH_KEY = 0.11
SLEEP_NO_KEY = 0.34
ID_COL = "biosample"


def open_maybe_gzip(path, mode="rt"):
    """Read or write plain or gzipped text. Spec: every script reads compressed input."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, newline="" if "t" in mode else None)
    return open(path, mode, newline="" if "t" in mode else None)


def load_existing(path):
    """Return (rows_by_id, column_order). Missing file is not an error."""
    if not path or not os.path.exists(path):
        return {}, []
    with open_maybe_gzip(path, "rt") as fh:
        reader = csv.DictReader(fh)
        cols = list(reader.fieldnames or [])
        rows = {r[ID_COL]: r for r in reader if r.get(ID_COL)}
    return rows, cols


def write_table(path, rows, columns):
    """Atomic write: build beside the target, then rename."""
    ordered = [ID_COL] + [c for c in columns if c != ID_COL]
    tmp = f"{path}.tmp"
    with open_maybe_gzip(tmp, "wt") as fh:
        writer = csv.DictWriter(fh, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        for key in sorted(rows):
            writer.writerow({c: rows[key].get(c, "") for c in ordered})
    os.replace(tmp, path)


def read_accessions(args):
    """Accessions from --sample, --input, --from-runs, or stdin."""
    acc = []
    if args.sample:
        acc = list(args.sample)
    elif args.from_runs:
        with open_maybe_gzip(args.from_runs, "rt") as fh:
            for row in csv.DictReader(fh):
                # ENA calls it sample_accession; it holds the BioSample.
                value = row.get("sample_accession") or row.get("biosample") or ""
                if value:
                    acc.append(value.strip())
    else:
        handle = open(args.input, "rt") if args.input else sys.stdin
        with handle as fh:
            acc = [line.strip() for line in fh if line.strip()]
    # Preserve order, drop duplicates.
    return list(dict.fromkeys(acc))


def fetch_batch(accessions, debug=False):
    """Resolve a batch of accessions to UIDs, fetch them, return {biosample: {attr: value}}.

    epost on the biosample database needs NCBI's internal numeric UIDs; it rejects
    SAMN accessions ("IDs contain invalid characters"). Measured 2026-11-06. So the
    accession-to-UID step of the original script is still required. The difference
    here is that it runs once per batch, as a single esearch over an OR query held on
    the history server, instead of once per accession.
    """
    term = " OR ".join(f"{a}[accn]" for a in accessions)
    handle = Entrez.esearch(db="biosample", term=term, retmax=len(accessions), usehistory="y")
    found = Entrez.read(handle)
    handle.close()
    if not found.get("IdList") and int(found.get("Count", 0)) == 0:
        return {}

    handle = Entrez.efetch(
        db="biosample",
        webenv=found["WebEnv"],
        query_key=found["QueryKey"],
        retmax=len(accessions),
    )
    raw = handle.read()
    handle.close()
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", errors="replace")

    out = {}
    root = ET.parse(io.StringIO(raw)).getroot()
    for sample in root.iter("BioSample"):
        acc = sample.attrib.get("accession")
        if not acc:
            continue
        record = {ID_COL: acc}

        organism = sample.find(".//Description/Organism")
        if organism is not None:
            if organism.attrib.get("taxonomy_id"):
                record["taxid"] = organism.attrib["taxonomy_id"]
            if organism.attrib.get("taxonomy_name"):
                record["scientific_name"] = organism.attrib["taxonomy_name"]

        # Scope to this sample, not the whole document. The original script iterated
        # from the root inside the per-sample loop, which copies one sample's
        # attributes onto every sample in a multi-sample response.
        for attribute in sample.iterfind(".//Attributes/Attribute"):
            name = attribute.attrib.get("harmonized_name") or attribute.attrib.get("attribute_name")
            if not name:
                continue
            text = (attribute.text or "").strip()
            # Spec 001: empty means unknown. Never store the string NA.
            if text and text.upper() not in {"NA", "N/A", "MISSING", "NOT APPLICABLE", "NOT COLLECTED"}:
                record[name] = text

        if debug:
            print(f"  {acc}: {len(record) - 1} attributes", file=sys.stderr)
        out[acc] = record
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-e", "--email", required=True, help="NCBI Entrez email (required by NCBI)")
    p.add_argument("-o", "--output", required=True, help="output CSV, .csv or .csv.gz")
    p.add_argument("-i", "--input", help="file of accessions, one per line")
    p.add_argument("-s", "--sample", nargs="+", help="accessions on the command line")
    p.add_argument("--from-runs", help="a catalog runs CSV; reads its sample_accession column")
    p.add_argument("-u", "--update", action="store_true", help="re-query accessions already present")
    p.add_argument("--batch", type=int, default=BATCH, help=f"accessions per request (default {BATCH})")
    p.add_argument("--debug", action="store_true")
    args = p.parse_args()

    Entrez.email = args.email
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key
        sleep = SLEEP_WITH_KEY
        print("Using NCBI_API_KEY from the environment (10 requests/s).", file=sys.stderr)
    else:
        sleep = SLEEP_NO_KEY
        print("No NCBI_API_KEY set; limited to 3 requests/s.", file=sys.stderr)

    rows, columns = load_existing(args.output)
    if rows:
        print(f"Loaded {len(rows)} existing records from {args.output}", file=sys.stderr)

    wanted = read_accessions(args)
    todo = wanted if args.update else [a for a in wanted if a not in rows]
    print(f"{len(wanted)} requested, {len(todo)} to fetch", file=sys.stderr)
    if not todo:
        return 0

    seen = set(columns)
    failures = []
    start = time.time()

    for i in range(0, len(todo), args.batch):
        chunk = todo[i : i + args.batch]
        label = f"batch {i // args.batch + 1} ({len(chunk)} accessions)"
        try:
            fetched = fetch_batch(chunk, debug=args.debug)
        except Exception as exc:  # noqa: BLE001 - Entrez raises many unrelated types
            print(f"{label} FAILED: {exc}", file=sys.stderr)
            failures.extend(chunk)
            time.sleep(5)
            continue

        for acc, record in fetched.items():
            rows[acc] = record
            seen.update(record.keys())

        missing = [a for a in chunk if a not in fetched]
        if missing:
            # Report it. A silently absent accession is how a panel loses samples.
            print(f"{label}: {len(missing)} returned no record, first {missing[:3]}", file=sys.stderr)
            failures.extend(missing)

        # Flush after every batch so a crash resumes instead of restarting.
        write_table(args.output, rows, sorted(seen))
        done = min(i + args.batch, len(todo))
        rate = done / max(time.time() - start, 1e-9)
        print(f"{label}: {len(fetched)} fetched, {done}/{len(todo)} done, {rate:.1f}/s", file=sys.stderr)
        time.sleep(sleep)

    write_table(args.output, rows, sorted(seen))
    print(f"Wrote {len(rows)} records, {len(seen)} columns, to {args.output}", file=sys.stderr)

    if failures:
        path = f"{args.output}.failed.txt"
        with open(path, "wt") as fh:
            fh.write("\n".join(dict.fromkeys(failures)) + "\n")
        print(f"{len(set(failures))} accessions failed; listed in {path}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
