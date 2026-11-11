#!/usr/bin/env python3
"""Download concurrency sweep (issue #27): plan the runs, then summarise the results.

The question: how many `popseq fetch` downloads can run at once before throughput per
stream falls or failures rise? `popseq fetch` does not cap concurrency; the caller does.
The sweep runs the same pool of runs at several levels N, each level a SLURM array with
`%N` (scripts/concurrency_sweep.slurm), one level after another.

`%N` is a ceiling, not a measurement: the queue decides how many tasks actually run
together. So the summary computes the concurrency each fetch really saw, from the
start and end time of every task: the time-weighted mean number of fetches in flight
during that fetch's window, itself included.

    # 1. choose the pool and look each run up once (ENA, NCBI): pool.csv and plan.csv.
    #    The fetches then read plan.csv (`popseq fetch --plan`) and make no lookups, so
    #    only the transfers compete, not 40 NCBI queries at once.
    uv run python scripts/popseq_concurrency.py plan --n 80 --out $W/pool.csv --plan $W/plan.csv
    # 2. run scripts/concurrency_sweep.slurm (see its header)
    # 3. summarise
    uv run python scripts/popseq_concurrency.py summary $W/tasks.tsv \\
        --out benchmarks/2026-11-11_concurrency_sweep.csv \\
        --by-node-out benchmarks/2026-11-11_concurrency_by_node.csv

tasks.tsv, one line per task, written by the SLURM script:
    level  run  taxid  node  start_epoch  end_epoch  exit  bytes  seconds_fetch
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import random
import statistics
import sys
from collections import defaultdict

PLAN_FIELDS = ["run", "taxid", "species", "study", "read_count", "base_count",
               "fastq_bytes", "fastq_md5", "fastq_ftp", "fastq_aspera", "role",
               "ncbi_spots", "ncbi_spots_with_mates"]
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TASK_FIELDS = ["level", "run", "taxid", "node", "start_epoch", "end_epoch", "exit",
               "bytes", "seconds_fetch"]
OK_EXIT = 0
# popseq fetch exit codes that are download outcomes (spec 002). Any other non-zero exit
# (a traceback is 1) failed before or outside the transfer: reported as `crashed`, not
# counted as a download failure, because it says nothing about concurrency.
DOWNLOAD_FAIL_EXITS = {10, 11, 12, 13}


def plan(rows, n: int, lo_gbases: float, hi_gbases: float, seed: int = 0) -> list[dict]:
    """n runs held at NCBI, paired, within a size band, at most one per species first.

    The size band keeps per-run time similar so that throughput differences come from
    concurrency, not from run size. Species are spread so that no single study or
    submitter's storage dominates."""
    ok = [r for r in rows
          if r.get("ncbi_status") == "in_ncbi" and r.get("fastq_ftp")
          and lo_gbases * 1e9 <= float(r.get("base_count") or 0) <= hi_gbases * 1e9]
    rng = random.Random(seed)
    rng.shuffle(ok)
    by_sp: dict[str, list[dict]] = defaultdict(list)
    for r in ok:
        by_sp[r["species_taxid"]].append(r)
    out: list[dict] = []
    i = 0
    while len(out) < n and any(i < len(v) for v in by_sp.values()):
        for sp in sorted(by_sp):
            if i < len(by_sp[sp]) and len(out) < n:
                out.append(by_sp[sp][i])
        i += 1
    return out


def in_flight(tasks: list[dict]) -> list[float]:
    """Time-weighted mean number of tasks running during each task's window.

    Tasks of different levels never overlap by design, but the computation does not
    assume it: it counts every task in the list."""
    iv = [(float(t["start_epoch"]), float(t["end_epoch"])) for t in tasks]
    out = []
    for s, e in iv:
        if e <= s:
            out.append(1.0)
            continue
        overlap = sum(max(0.0, min(e, e2) - max(s, s2)) for s2, e2 in iv)
        out.append(overlap / (e - s))
    return out


def summarise(tasks: list[dict]) -> list[dict]:
    """One row per level: realised concurrency, per-stream MB/s, aggregate, failures."""
    by_level: dict[int, list[dict]] = defaultdict(list)
    for t in tasks:
        by_level[int(t["level"])].append(t)
    rows = []
    for lvl in sorted(by_level):
        ts = by_level[lvl]
        conc = in_flight(ts)
        ok = [t for t in ts if int(t["exit"]) == OK_EXIT
              and float(t["seconds_fetch"] or 0) > 0]
        mbps = sorted(float(t["bytes"]) / float(t["seconds_fetch"]) / 1e6 for t in ok)
        span = (max(float(t["end_epoch"]) for t in ts)
                - min(float(t["start_epoch"]) for t in ts))
        q = statistics.quantiles(mbps, n=4) if len(mbps) >= 2 else [None] * 3
        rows.append({
            "level": lvl, "tasks": len(ts), "ok": len(ok),
            "failed": sum(1 for t in ts if int(t["exit"]) in DOWNLOAD_FAIL_EXITS),
            "crashed": sum(1 for t in ts if int(t["exit"]) not in DOWNLOAD_FAIL_EXITS
                           and int(t["exit"]) != OK_EXIT),
            "realised_concurrency_median": round(statistics.median(conc), 1),
            "realised_concurrency_max": round(max(conc), 1),
            "nodes": len({t["node"] for t in ts}),
            "stream_mb_s_median": round(statistics.median(mbps), 2) if mbps else "",
            "stream_mb_s_q1": round(q[0], 2) if q[0] is not None else "",
            "stream_mb_s_q3": round(q[2], 2) if q[2] is not None else "",
            "aggregate_mb_s": round(sum(float(t["bytes"]) for t in ok) / span / 1e6, 2)
            if span > 0 else "",
            "span_s": round(span),
        })
    return rows


NODE_BINS = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 1000)]


def by_node(tasks: list[dict]) -> list[dict]:
    """Per-stream MB/s binned by how many fetches shared the task's node.

    SLURM packs array tasks onto few nodes, so a level's `%N` mixes two things: fetches
    cluster-wide and fetches per node. Each successful task gets its time-weighted
    in-flight count among all tasks (failures included) of its own level on its own node, then the rates are
    binned on that count. Levels are pooled: the bin, not the level, is the variable."""
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for t in tasks:
        groups[(t["level"], t["node"])].append(t)   # failures still used the node
    pairs = []
    for g in groups.values():
        for t, c in zip(g, in_flight(g)):
            if int(t["exit"]) == OK_EXIT and float(t["seconds_fetch"] or 0) > 0:
                pairs.append((c, float(t["bytes"]) / float(t["seconds_fetch"]) / 1e6))
    rows = []
    for lo, hi in NODE_BINS:
        v = sorted(r for c, r in pairs if lo <= c < hi)
        if not v:
            continue
        q = statistics.quantiles(v, n=4) if len(v) >= 2 else [v[0], v[0], v[0]]
        rows.append({"same_node_in_flight": f"{lo}-{hi}" if hi < 1000 else f">={lo}",
                     "n": len(v), "stream_mb_s_median": round(statistics.median(v), 2),
                     "stream_mb_s_q1": round(q[0], 2), "stream_mb_s_q3": round(q[2], 2)})
    return rows


def _write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def _open(path):
    return gzip.open(path, "rt", newline="") if path.endswith(".gz") else open(path, newline="")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("plan", help="choose the run pool from the catalog")
    a.add_argument("--runs", default=os.path.join(HERE, "catalog", "runs",
                                                  "fungi_all_runs.csv.gz"))
    a.add_argument("--panel", default=os.path.join(HERE, "panels", "pop_sra_ge10_cap300",
                                                   "metadata.csv"),
                   help="draw runs from this panel's run_accession list")
    a.add_argument("--n", type=int, default=80)
    a.add_argument("--lo-gbases", type=float, default=0.8)
    a.add_argument("--hi-gbases", type=float, default=1.6)
    a.add_argument("--seed", type=int, default=0)
    a.add_argument("--out", required=True)
    a.add_argument("--plan", required=True,
                   help="also write the fetch plan (benchmark 7.1 format) here")
    b = sub.add_parser("summary", help="summarise tasks.tsv per level")
    b.add_argument("tasks")
    b.add_argument("--out", required=True)
    b.add_argument("--by-node-out", help="also write the per-node binning here")
    args = p.parse_args()

    if args.cmd == "plan":
        with _open(args.panel) as fh:
            want = {r["run_accession"] for r in csv.DictReader(fh)}
        with _open(args.runs) as fh:
            rows = [r for r in csv.DictReader(fh) if r["run_accession"] in want]
        # Draw a quarter more than needed: runs NCBI holds unpaired, or that NCBI cannot
        # answer for, are dropped, so every fetch can take the AWS route.
        pool = plan(rows, args.n + args.n // 4, args.lo_gbases, args.hi_gbases, args.seed)
        sys.path.insert(0, HERE)
        from popseq import lookup
        ena = lookup.ena_lookup([r["run_accession"] for r in pool])
        prow, kept = [], []
        for r in pool:
            if len(kept) == args.n:
                break
            spots, mates = lookup.ncbi_spots(r["run_accession"])
            e = ena.get(r["run_accession"])
            if not e or not spots or mates != spots:
                continue
            kept.append(r)
            prow.append({"run": r["run_accession"], "taxid": r["species_taxid"],
                         "species": r["species_name"], "study": r.get("study_accession", ""),
                         "read_count": e.get("read_count", ""), "base_count": e.get("base_count", ""),
                         "fastq_bytes": e.get("fastq_bytes", ""), "fastq_md5": e.get("fastq_md5", ""),
                         "fastq_ftp": e.get("fastq_ftp", ""), "fastq_aspera": e.get("fastq_aspera", ""),
                         "role": "concurrency", "ncbi_spots": spots,
                         "ncbi_spots_with_mates": mates})
        pool = kept
        with open(args.plan, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=PLAN_FIELDS, lineterminator="\n")
            w.writeheader()
            w.writerows(prow)
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["run_accession", "species_taxid", "species_name", "base_count"])
            for r in pool:
                w.writerow([r["run_accession"], r["species_taxid"], r["species_name"],
                            r["base_count"]])
        print(f"wrote {args.out}: {len(pool)} runs from "
              f"{len({r['species_taxid'] for r in pool})} species", file=sys.stderr)
    else:
        with _open(args.tasks) as fh:
            tasks = list(csv.DictReader(fh, delimiter="\t"))
        rows = summarise(tasks)
        _write(args.out, rows)
        if args.by_node_out:
            _write(args.by_node_out, by_node(tasks))
        for r in rows:
            print(r, file=sys.stderr)


if __name__ == "__main__":
    main()
