#!/usr/bin/env python3
"""Is a run's FASTQ in genome-position order? A cheap probe from the first reads.

Why: tier `coverage` with `--sampling windows` keeps contiguous blocks of reads. In a run
sorted by genome position those blocks reach only `kept bases / depth` of the genome (GmG
phase 1 pilot: two C. auris runs gave a 2.2 and 2.3 Mb genome against 12.1 Mb). Stride
sampling has no such bias but converts the whole run, at 3.5 to 4 times the cost (PR #42).
A probe that recognises a sorted run before fetching would let windows stay cheap for the
rest. This script measures whether such a probe works. It does not change `popseq fetch`.

Read order belongs to an archive's copy, not to the run: NCBI can hold a run in genome order
whose ENA FASTQ is not (ERR12543635). So the probe reads both copies: the head of the ENA
mate-1 FASTQ, and --reads spots of the SRA copy at the start, 1/3 and 2/3 of the run
(`fastq-dump -N/-X`, remote, no `.sra` download; needs the sratoolkit module).

The score: stream the first --reads records of mate 1 (a few MB). In a sorted file a
read overlaps the reads just before it, because they come from neighbouring positions. In an
unsorted file neighbours come from random places and rarely share sequence. The score is the
fraction of reads (after the first --window) that share a canonical k-mer with any of the
previous --window reads. Canonical, because a FASTQ exported from a BAM can hold reads of both
strands. Mate 1 only, because the two mates of one fragment overlap whatever the order.

    module load sratoolkit/3.2.0
    uv run python scripts/popseq_sortcheck.py probe runs.txt --out sortcheck.tsv
    uv run python scripts/popseq_sortcheck.py sample --n 200 --out runs.txt   # from the panel

Network: ENA portal (one lookup per 200 runs), ENA HTTPS (one partial download per run),
NCBI runinfo (one lookup per run) and NCBI SRA (three fastq-dump slices per run).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import gzip
import io
import os
import random
import sys
import time
import urllib.request
from collections import Counter, deque

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

FIELDS = ["run", "submitted_format", "ena_reads", "read_len", "ena_score", "ena_status",
          "sra_spots", "sra_scores", "sra_score_max", "sra_status", "seconds"]
SRA_POSITIONS = (0.0, 1 / 3, 2 / 3)
_COMP = bytes.maketrans(b"ACGTN", b"TGCAN")


def canonical_kmers(seq: bytes, k: int, step: int = 1) -> set[bytes]:
    """Canonical k-mers (the smaller of a k-mer and its reverse complement), skipping N."""
    seq = seq.upper()
    rc = seq.translate(_COMP)[::-1]
    n = len(seq)
    out = set()
    for i in range(0, n - k + 1, step):
        f = seq[i:i + k]
        if b"N" in f:
            continue
        r = rc[n - i - k:n - i]
        out.add(f if f < r else r)
    return out


def neighbour_score(seqs: list[bytes], k: int = 25, window: int = 50, step: int = 5) -> float:
    """Fraction of reads, after the first `window`, that share a canonical k-mer with one of
    the `window` reads before them.

    The previous reads contribute every k-mer; the current read is tested on every `step`-th
    k-mer, which is enough to find an overlap of more than about k + step bases."""
    held: Counter = Counter()
    recent: deque = deque()
    hits = tested = 0
    for i, s in enumerate(seqs):
        if i >= window:
            tested += 1
            if any(km in held for km in canonical_kmers(s, k, step)):
                hits += 1
        ks = canonical_kmers(s, k)
        recent.append(ks)
        held.update(ks)
        if len(recent) > window:
            for km in recent.popleft():
                held[km] -= 1
                if held[km] == 0:
                    del held[km]
    return hits / tested if tested else float("nan")


class _Counted(io.RawIOBase):
    """A read-only stream that counts the compressed bytes taken from the network."""

    def __init__(self, raw):
        self.raw, self.n = raw, 0

    def readable(self):
        return True

    def readinto(self, b):
        data = self.raw.read(len(b))
        self.n += len(data)
        b[:len(data)] = data
        return len(data)


def head_reads(url: str, n: int, timeout: int = 120) -> tuple[list[bytes], int]:
    """The sequences of the first n FASTQ records at an ENA HTTPS URL, and the compressed
    bytes read. Streams, and closes the connection once n records are in."""
    req = urllib.request.Request(url, headers={"User-Agent": "PopSeq-sortcheck"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = _Counted(resp)
        return read_seqs(gzip.GzipFile(fileobj=io.BufferedReader(raw)), n), raw.n


def read_seqs(handle, n: int) -> list[bytes]:
    """Sequence lines of the first n records of a FASTQ handle (bytes lines)."""
    seqs = []
    while len(seqs) < n:
        rec = [handle.readline() for _ in range(4)]
        if not rec[3]:
            break
        seqs.append(rec[1].rstrip(b"\r\n"))
    return seqs


def ena_meta(runs: list[str]) -> dict[str, dict]:
    from popseq import lookup  # noqa: PLC0415
    import requests  # noqa: PLC0415
    out = {}
    for i in range(0, len(runs), 200):
        r = requests.post(lookup.ENA_SEARCH, data={
            "result": "read_run", "includeAccessions": ",".join(runs[i:i + 200]),
            "fields": "run_accession,fastq_ftp,submitted_format", "format": "tsv",
            "limit": 0}, timeout=120)
        r.raise_for_status()
        for row in csv.DictReader(io.StringIO(r.text), delimiter="\t"):
            out[row["run_accession"]] = row
    return out


def sra_probe(run: str, spots: int, n: int, k: int, window: int, workdir: str) -> list[float]:
    """Scores of n mate-1 spots of the SRA copy at each of SRA_POSITIONS (fastq-dump, remote)."""
    import subprocess  # noqa: PLC0415
    import tempfile  # noqa: PLC0415
    scores = []
    for frac in SRA_POSITIONS:
        start = int(spots * frac) + 1
        end = min(spots, start + n - 1)
        with tempfile.TemporaryDirectory(dir=workdir) as tmp:
            subprocess.run(["fastq-dump", "-N", str(start), "-X", str(end), "--split-files",
                            "-O", tmp, run], stdin=subprocess.DEVNULL,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=600)
            fq = sorted(f for f in os.listdir(tmp) if f.endswith(".fastq"))
            if not fq:
                scores.append(float("nan"))
                continue
            with open(os.path.join(tmp, fq[0]), "rb") as fh:   # _1, or the only file
                seqs = read_seqs(fh, n)
            scores.append(neighbour_score(seqs, k=k, window=window)
                          if len(seqs) > window else float("nan"))
    return scores


def probe_one(run: str, meta: dict, n: int, k: int, window: int, workdir: str,
              sra: bool = True) -> dict:
    row = {f: "" for f in FIELDS}
    row.update(run=run, submitted_format=meta.get("submitted_format", ""))
    t0 = time.monotonic()
    urls = [u for u in (meta.get("fastq_ftp") or "").split(";") if u]
    mate1 = [u for u in urls if u.endswith("_1.fastq.gz")] or urls[:1]
    if not mate1:
        row["ena_status"] = "no_fastq_at_ena"
    else:
        try:
            seqs, _ = head_reads("https://" + mate1[0], n)
            row.update(ena_reads=len(seqs),
                       read_len=round(sum(map(len, seqs)) / len(seqs)) if seqs else "")
            if len(seqs) > window:
                row.update(ena_score=round(neighbour_score(seqs, k=k, window=window), 4),
                           ena_status="ok")
            else:
                row["ena_status"] = "too_few_reads"
        except Exception as exc:  # noqa: BLE001  network errors are recorded, not raised
            row["ena_status"] = f"error: {type(exc).__name__}"
    if sra:
        from popseq import lookup  # noqa: PLC0415
        spots, _ = lookup.ncbi_spots(run)
        if not spots:
            row["sra_status"] = "not_at_ncbi_or_no_answer"
        else:
            try:
                sc = sra_probe(run, spots, n, k, window, workdir)
                ok = [x for x in sc if x == x]
                row.update(sra_spots=spots, sra_scores=";".join(f"{x:.4f}" for x in sc),
                           sra_score_max=round(max(ok), 4) if ok else "",
                           sra_status="ok" if ok else "no_reads")
            except Exception as exc:  # noqa: BLE001
                row["sra_status"] = f"error: {type(exc).__name__}"
    row["seconds"] = round(time.monotonic() - t0, 1)
    return row


BINS = [(0.0, 0.1), (0.1, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v == v else None


def summarise(rows: list[dict], exclude: set[str] = frozenset()) -> list[dict]:
    """Counts of runs per score bin, for the ENA head and for the SRA copy (highest of the
    three positions), over rows whose run is not in `exclude`."""
    rows = [r for r in rows if r["run"] not in exclude]
    out = []
    for lo, hi in BINS:
        out.append({"score_bin": f"{lo}-{min(hi, 1.0)}",
                    "ena": sum(1 for r in rows if (v := _f(r["ena_score"])) is not None
                               and lo <= v < hi),
                    "sra_max": sum(1 for r in rows if (v := _f(r["sra_score_max"])) is not None
                                   and lo <= v < hi)})
    out.append({"score_bin": "no_score",
                "ena": sum(1 for r in rows if _f(r["ena_score"]) is None),
                "sra_max": sum(1 for r in rows if _f(r["sra_score_max"]) is None)})
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("probe", help="score runs listed one per line")
    a.add_argument("runs")
    a.add_argument("--out", required=True)
    a.add_argument("--reads", type=int, default=4000)
    a.add_argument("--k", type=int, default=25)
    a.add_argument("--window", type=int, default=50)
    a.add_argument("--workers", type=int, default=4, help="runs probed at once")
    a.add_argument("--no-sra", action="store_true", help="probe the ENA copy only")
    a.add_argument("--workdir", default=os.environ.get("SCRATCH") or ".",
                   help="temporary fastq-dump output (default $SCRATCH)")
    c = sub.add_parser("summary", help="score bins from a probe TSV")
    c.add_argument("probe")
    c.add_argument("--exclude", help="file of runs to leave out (e.g. the pilot runs)")
    c.add_argument("--out", required=True)
    b = sub.add_parser("sample", help="random runs from a panel's metadata.csv")
    b.add_argument("--panel", default=os.path.join(HERE, "panels", "pop_sra_ge10_cap300",
                                                   "metadata.csv"))
    b.add_argument("--n", type=int, default=200)
    b.add_argument("--seed", type=int, default=0)
    b.add_argument("--out", required=True)
    args = p.parse_args()

    if args.cmd == "summary":
        with open(args.probe, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        excl = set(open(args.exclude).read().split()) if args.exclude else set()
        out = summarise(rows, excl)
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(out)
        return

    if args.cmd == "sample":
        with open(args.panel, newline="") as fh:
            runs = [r["run_accession"] for r in csv.DictReader(fh)]
        pick = random.Random(args.seed).sample(runs, args.n)
        with open(args.out, "w") as fh:
            fh.write("\n".join(pick) + "\n")
        return

    with open(args.runs) as fh:
        runs = [x.strip() for x in fh if x.strip() and not x.startswith("#")]
    meta = ena_meta(runs)
    rows = []
    with cf.ThreadPoolExecutor(args.workers) as ex:
        futs = {ex.submit(probe_one, r, meta.get(r, {}), args.reads, args.k, args.window,
                          args.workdir, not args.no_sra): r for r in runs}
        for f in cf.as_completed(futs):
            rows.append(f.result())
            r = rows[-1]
            print(f"{r['run']}\tena {r['ena_status']} {r['ena_score']}\t"
                  f"sra {r['sra_status']} {r['sra_scores']}\t{r['submitted_format']}",
                  file=sys.stderr)
    order = {r: i for i, r in enumerate(runs)}
    rows.sort(key=lambda r: order[r["run"]])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
