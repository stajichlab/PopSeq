#!/usr/bin/env python3
"""Import per-run quality flags from a producer table into catalog/qc/ (specification 008).

    uv run python scripts/popseq_qc_import.py --source-id matpredict_reads_screen_2026-11-08 --format matpredict_screen \\
        --input contamination_flags.tsv --runs cache/fungi_runs_ena_2026-11-06.tsv.gz \\
        --project MATPredict --tool "GmG 2026.10.5 + screen_reads.py" --check contamination \\
        --method "..." --detection-limit "..." --reference "..." --date 2026-11-08

Each run is resolved through the PopSeq kingdom run table (run -> sample, tax_id) and catalog/taxonomy_species_map.csv.gz
(tax_id -> species). The rows of --source-id are replaced; rows of other sources and every curated_* column are kept.
Input may be CSV or TSV, plain, .gz or .zst. Standard library only (zstd through the `zstd` binary).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

QC_COLS = ["run_accession", "sample_accession", "tax_id", "species_taxid", "species_name", "source_id", "check", "status",
           "detected_genus", "second_taxon", "second_hashes", "second_share", "deeper_check", "subsample", "strain",
           "source_note", "curated_status", "curated_note"]
SOURCE_COLS = ["source_id", "project", "tool", "check", "method", "detection_limit", "reference", "date", "rows", "notes"]
CURATED = ("curated_status", "curated_note")

# producer format -> {qc column: producer column}; a value that is a constant is written as given when prefixed with "="
FORMATS = {
    "matpredict_screen": {"run_accession": "run", "status": "screen", "detected_genus": "screen_genus", "second_taxon": "second_genus",
                          "second_hashes": "second_hashes", "second_share": "second_share", "deeper_check": "deeper_screen",
                          "strain": "strain", "subsample": "=2000000 spots, 8 windows (6000000 when deeper_check=yes)",
                          "source_note": "species"},
}


def open_text(path: str):
    """A text handle for a plain, .gz or .zst file."""
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), newline="")
    if path.endswith(".zst"):
        p = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True)
        return io.StringIO(p.stdout.decode(), newline="")
    return open(path, newline="")


def read_table(path: str) -> list[dict]:
    with open_text(path) as fh:
        text = fh.read()
    first = text.split("\n", 1)[0]
    return list(csv.DictReader(io.StringIO(text), delimiter="\t" if first.count("\t") > first.count(",") else ","))


def species_file(species_taxid: str, species_name: str) -> str:
    if not species_taxid:
        return "unresolved.qc.csv.gz"
    return f"{species_taxid}.{re.sub(r'[^A-Za-z0-9]+', '_', species_name).strip('_')}.qc.csv.gz"


def load_qc(path: str) -> list[dict]:
    return read_table(path) if os.path.exists(path) else []


def write_qc(path: str, rows: list[dict]) -> None:
    rows = sorted(rows, key=lambda r: (r["run_accession"], r["source_id"], r["check"]))
    with gzip.open(path, "wt", newline="") as fh:     # gzip header without a name and with mtime 0 would be nicer, but readers do not care
        w = csv.DictWriter(fh, QC_COLS, extrasaction="ignore", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--source-id", required=True); ap.add_argument("--format", required=True, choices=sorted(FORMATS))
    ap.add_argument("--input", required=True, help="the producer table")
    ap.add_argument("--runs", required=True, help="PopSeq kingdom run table (run_accession, sample_accession, tax_id, ...)")
    ap.add_argument("--taxonomy", default="catalog/taxonomy_species_map.csv.gz")
    ap.add_argument("--qc-dir", default="catalog/qc")
    ap.add_argument("--project", required=True); ap.add_argument("--tool", required=True); ap.add_argument("--check", default="contamination")
    ap.add_argument("--method", default=""); ap.add_argument("--detection-limit", default=""); ap.add_argument("--reference", default="")
    ap.add_argument("--date", required=True); ap.add_argument("--notes", default="")
    a = ap.parse_args(argv)

    fmt = FORMATS[a.format]
    want = {r[fmt["run_accession"]] for r in read_table(a.input)}
    runinfo = {}
    with open_text(a.runs) as fh:
        first = fh.readline(); delim = "\t" if first.count("\t") > first.count(",") else ","
        cols = next(csv.reader([first], delimiter=delim))
        for r in csv.DictReader(fh, fieldnames=cols, delimiter=delim):
            if r["run_accession"] in want:
                runinfo[r["run_accession"]] = (r.get("sample_accession", ""), r.get("tax_id", ""))
    sp = {r["taxid"]: (r["species_taxid"], r["species_name"]) for r in read_table(a.taxonomy)}
    new = defaultdict(list); unresolved = 0
    for r in read_table(a.input):
        run = r[fmt["run_accession"]]
        sample, taxid = runinfo.get(run, ("", ""))
        st, sn = sp.get(taxid, ("", ""))
        if not st: unresolved += 1
        row = {"run_accession": run, "sample_accession": sample, "tax_id": taxid, "species_taxid": st, "species_name": sn,
               "source_id": a.source_id, "check": a.check}
        for qc, src in fmt.items():
            if qc == "run_accession": continue
            row[qc] = src[1:] if src.startswith("=") else r.get(src, "")
        new[species_file(st, sn)].append(row)

    os.makedirs(a.qc_dir, exist_ok=True)
    touched = set(new)
    for f in os.listdir(a.qc_dir):          # a species file that only held rows of this source must lose them too
        if f.endswith(".qc.csv.gz"): touched.add(f)
    total = Counter()
    for f in sorted(touched):
        path = os.path.join(a.qc_dir, f)
        old = load_qc(path)
        cur = {(r["run_accession"], r["source_id"], r["check"]): {k: r.get(k, "") for k in CURATED} for r in old}
        keep = [r for r in old if r["source_id"] != a.source_id]
        add = [dict(r, **cur.get((r["run_accession"], r["source_id"], r["check"]), {k: "" for k in CURATED})) for r in new.get(f, [])]
        rows = keep + add
        if rows:
            write_qc(path, rows)
        elif os.path.exists(path):
            os.remove(path)
        for r in rows: total[r["source_id"]] += 1

    sp_path = os.path.join(a.qc_dir, "sources.csv")
    srcs = [r for r in (read_table(sp_path) if os.path.exists(sp_path) else []) if r["source_id"] != a.source_id]
    srcs.append({"source_id": a.source_id, "project": a.project, "tool": a.tool, "check": a.check, "method": a.method,
                 "detection_limit": a.detection_limit, "reference": a.reference, "date": a.date, "rows": "", "notes": a.notes})
    for s in srcs: s["rows"] = str(total.get(s["source_id"], 0))
    with open(sp_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, SOURCE_COLS, extrasaction="ignore", lineterminator="\n"); w.writeheader(); w.writerows(sorted(srcs, key=lambda s: s["source_id"]))
    n = sum(len(v) for v in new.values())
    print(f"{a.source_id}: {n} rows in {len([f for f in new if f != 'unresolved.qc.csv.gz'])} species files; {unresolved} unresolved; "
          f"{len(runinfo)} of {len(want)} runs found in the run table")
    return 0


if __name__ == "__main__":
    sys.exit(main())
