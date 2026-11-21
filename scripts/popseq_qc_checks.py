#!/usr/bin/env python3
"""Run-QC checks that PopSeq evaluates itself from catalog data, imported into catalog/qc/ (specification 008).

Each subcommand writes a findings table (format `popseq_findings`) and imports it with scripts/popseq_qc_import.py, so the
rules of spec 008 hold: the source replaces only its own rows and curated_* columns survive.

    # base_count / genome size below a threshold, for the runs named in an issue (#37, #49). Genome sizes from
    # catalog/genome_sizes.csv. The run list is required; the check does not sweep the catalog unless --catalog-wide is given.
    uv run python scripts/popseq_qc_checks.py low-coverage --run-list catalog/qc/findings/issue_named_runs.csv --date 2026-11-21

    # archive bytes per base above a threshold (issue #49). The archive sizes are measured by a consumer and kept in a table
    # (run_accession, archive_bytes, archive_note); the catalog does not hold archive sizes.
    uv run python scripts/popseq_qc_checks.py archive-size --sizes catalog/qc/findings/archive_sizes.csv \\
        --source-id gmg_archive_size_2026-11-21 --project GmG --date 2026-11-21

Inputs may be plain, .gz or .zst. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("popseq_qc_import", os.path.join(HERE, "popseq_qc_import.py"))
qi = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(qi)

FINDING_COLS = ["run_accession", "status", "source_note"]


def _fmt_x(x: float) -> str:
    return f"{x:g}"


def genome_sizes(path: str) -> dict[str, float]:
    """species_taxid -> genome size in bases, for the species that have one."""
    out = {}
    for r in qi.read_table(path):
        try:
            out[r["species_taxid"]] = float(r["genome_mb"]) * 1e6
        except (KeyError, TypeError, ValueError):
            continue
    return out


def iter_runs(path: str):
    with qi.open_text(path) as fh:
        first = fh.readline(); delim = "\t" if first.count("\t") > first.count(",") else ","
        cols = next(csv.reader([first], delimiter=delim))
        yield from csv.DictReader(fh, fieldnames=cols, delimiter=delim)


def read_run_list(path: str) -> set[str]:
    """Run accessions of an explicit list: a table with a run_accession column, or one accession per line."""
    with qi.open_text(path) as fh:
        lines = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    if lines and "run_accession" in lines[0]:
        first = lines[0]; delim = "\t" if first.count("\t") > first.count(",") else ","
        return {r["run_accession"].strip() for r in csv.DictReader(lines, delimiter=delim) if r.get("run_accession", "").strip()}
    return {l.split(",")[0].split("\t")[0].strip() for l in lines}


def low_coverage(runs, sizes: dict[str, float], below: float, only: set[str] | None = None) -> tuple[list[dict], dict]:
    """Runs whose base_count / genome size is below `below`. With `only` (the runs an issue names), no other run is
    examined or flagged; `only=None` is the catalog-wide pass. A run with no genome size for its species, or no positive
    base_count, is not examined (counted, not flagged)."""
    status = f"below_{_fmt_x(below)}x"
    rows, n = [], {"examined": 0, "no_genome_size": 0, "no_base_count": 0, "flagged": 0}
    for r in runs:
        if only is not None and r["run_accession"] not in only:
            continue
        g = sizes.get(r.get("species_taxid", ""))
        if not g:
            n["no_genome_size"] += 1; continue
        try:
            bases = int(r.get("base_count") or 0)
        except ValueError:
            bases = 0
        if bases <= 0:
            n["no_base_count"] += 1; continue
        n["examined"] += 1
        depth = bases / g
        if depth < below:
            n["flagged"] += 1
            rows.append({"run_accession": r["run_accession"], "status": status,
                         "source_note": f"depth {depth:.4f}x = base_count {bases} / genome {g / 1e6:g} Mb (catalog/genome_sizes.csv)"})
    return rows, n


def archive_size(runs, sizes_table: list[dict], max_bpb: float) -> tuple[list[dict], dict]:
    """Runs whose archive bytes per base (base_count from the run table) is above `max_bpb`. Only the runs in the sizes
    table are examined."""
    want = {r["run_accession"]: r for r in sizes_table}
    bases = {}
    for r in runs:
        if r["run_accession"] in want:
            try:
                bases[r["run_accession"]] = int(r.get("base_count") or 0)
            except ValueError:
                pass
    rows, n = [], {"in_sizes_table": len(want), "examined": 0, "flagged": 0}
    for acc, s in want.items():
        b = bases.get(acc, 0)
        if b <= 0 or not s.get("archive_bytes"):
            continue
        n["examined"] += 1
        bpb = int(s["archive_bytes"]) / b
        if bpb > max_bpb:
            n["flagged"] += 1
            note = f"{bpb:.2f} archive bytes per base = {int(s['archive_bytes'])} bytes / base_count {b}; threshold {_fmt_x(max_bpb)}"
            if s.get("archive_note"):
                note += f"; {s['archive_note']}"
            rows.append({"run_accession": acc, "status": "oversized", "source_note": note})
    return rows, n


def _import(rows: list[dict], a, check: str, method: str, detection_limit: str, exclude: str) -> int:
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "findings.csv")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, FINDING_COLS, lineterminator="\n"); w.writeheader(); w.writerows(rows)
        return qi.main(["--source-id", a.source_id, "--format", "popseq_findings", "--input", path, "--runs", a.runs,
                        "--taxonomy", a.taxonomy, "--qc-dir", a.qc_dir, "--project", a.project, "--tool", "scripts/popseq_qc_checks.py " + a.cmd,
                        "--check", check, "--method", method, "--detection-limit", detection_limit, "--reference", a.reference,
                        "--date", a.date, "--notes", a.notes, "--exclude-statuses", exclude])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("low-coverage", "archive-size"):
        p = sub.add_parser(name)
        p.add_argument("--runs", default="catalog/runs/fungi_all_runs.csv.gz")
        p.add_argument("--taxonomy", default="catalog/taxonomy_species_map.csv.gz")
        p.add_argument("--qc-dir", default="catalog/qc")
        p.add_argument("--date", required=True)
        p.add_argument("--reference", default="specifications/008-run-qc-flags.md")
        p.add_argument("--notes", default="")
        p.add_argument("--no-exclude", action="store_true", help="record the findings without making runs and panel skip them")
    lc = sub.choices["low-coverage"]
    lc.add_argument("--genome-sizes", default="catalog/genome_sizes.csv")
    grp = lc.add_mutually_exclusive_group(required=True)
    grp.add_argument("--run-list", help="examine only the runs in this file (a table with run_accession, or one accession per line): "
                     "the runs named in the linked issues. The default way to run the check")
    grp.add_argument("--catalog-wide", action="store_true", help="opt in: examine every run of the catalog (not done by default; "
                     "the owner limited below_1x to the runs named in issues #37 and #49)")
    lc.add_argument("--below", type=float, default=1.0, help="flag runs below this depth (x); default 1")
    lc.add_argument("--source-id", default="popseq_low_coverage"); lc.add_argument("--project", default="PopSeq")
    asz = sub.choices["archive-size"]
    asz.add_argument("--sizes", required=True, help="table with run_accession, archive_bytes and optionally archive_note")
    asz.add_argument("--max-bytes-per-base", type=float, default=2.0, help="flag runs above this; default 2.0 (a placeholder)")
    asz.add_argument("--source-id", required=True); asz.add_argument("--project", required=True)
    a = ap.parse_args(argv)

    if a.cmd == "low-coverage":
        only = None if a.catalog_wide else read_run_list(a.run_list)
        rows, n = low_coverage(iter_runs(a.runs), genome_sizes(a.genome_sizes), a.below, only)
        if only is not None:
            n["listed"] = len(only)
        status = f"below_{_fmt_x(a.below)}x"
        method = (f"depth = base_count (run table) / genome_mb (catalog/genome_sizes.csv) of the run's species_taxid; "
                  f"{status} when depth < {_fmt_x(a.below)}")
        scope = "every run of the catalog (--catalog-wide)" if only is None else f"only the {len(only)} runs of --run-list (named in the linked issues)"
        method += f"; examined: {scope}"
        limit = ("metadata only, not the reads: base_count as the archive reports it; runs of a species with no genome size, "
                 "or with no positive base_count, are not examined; "
                 + ("" if only is None else "no run outside the list is examined, so other low-depth runs are UNMEASURED"))
    else:
        rows, n = archive_size(iter_runs(a.runs), qi.read_table(a.sizes), a.max_bytes_per_base)
        status = "oversized"
        method = (f"archive bytes per base = archive object size (from --sizes) / base_count (run table); "
                  f"oversized when above {_fmt_x(a.max_bytes_per_base)}")
        limit = "only the runs in the sizes table are examined; the catalog holds no archive sizes"
    check = "low_coverage" if a.cmd == "low-coverage" else "oversized_archive_object"
    print(f"{a.cmd}: " + ", ".join(f"{k} {v}" for k, v in n.items()), file=sys.stderr)
    return _import(rows, a, check, method, limit, "" if a.no_exclude else status)


if __name__ == "__main__":
    sys.exit(main())
