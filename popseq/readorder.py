"""Is a run's FASTQ in genome-position order? The read-order score behind `--sampling auto`.

Why: tier `coverage` with windows sampling keeps contiguous blocks of reads. In a run sorted by
genome position those blocks reach only `kept bases / depth` of the genome. Stride sampling has no
such bias but converts the whole run, at 3.5 to 4 times the cost (PR #42). `--sampling auto`
scores the bytes it fetched and takes stride only for a sorted run.

The score of a block of mate-1 reads is the fraction of reads, after the first `window`, that
share a canonical k-mer with one of the `window` reads before them. In a sorted file a read
overlaps its neighbours; in an unsorted file they come from random places. Canonical, because a
FASTQ exported from a BAM holds reads of both strands. Mate 1 only, because the two mates of one
fragment overlap whatever the order.

Calibration (PR #43, `benchmarks/2026-11-11_sortcheck_*`): three blocks of 4,000 spots at the
start, 1/3 and 2/3 of the run; stride when the median is `THRESHOLD` (0.8) or more. On 30 runs
fetched both ways it caught 9 of the 10 runs whose windows sample reached under 0.9 of the stride
sample's genome size; the miss (SRR32416523, median 0.40) came out 13% short. Controls scored
under 0.3. In 200 random runs of `pop_sra_ge10_cap300`, 10 had a median of 0.8 or more.
"""

from __future__ import annotations

import gzip
import statistics
from collections import Counter, deque
from typing import Iterable, Optional

K = 25
WINDOW = 50
STEP = 5
BLOCK_SPOTS = 4000
POSITIONS = (0.0, 1 / 3, 2 / 3)
THRESHOLD = 0.8

_COMP = bytes.maketrans(b"ACGTN", b"TGCAN")


def canonical_kmers(seq: bytes, k: int = K, step: int = 1) -> set:
    """Canonical k-mers (the smaller of a k-mer and its reverse complement), skipping N."""
    seq = seq.upper()
    rc = seq.translate(_COMP)[::-1]
    n = len(seq)
    out = set()
    for i in range(0, n - k + 1, step):
        f = seq[i : i + k]
        if b"N" in f:
            continue
        r = rc[n - i - k : n - i]
        out.add(f if f < r else r)
    return out


def neighbour_score(seqs: list, k: int = K, window: int = WINDOW, step: int = STEP) -> float:
    """Fraction of reads, after the first `window`, sharing a canonical k-mer with one of the
    `window` reads before them. NaN when there are not more than `window` reads.

    The previous reads contribute every k-mer; the current read is tested on every `step`-th
    k-mer, which finds any overlap longer than about k + step bases."""
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


def read_seqs(handle, n: Optional[int] = None) -> list:
    """Sequence lines (bytes) of the first n records of a binary FASTQ handle (all if n is None)."""
    seqs = []
    while n is None or len(seqs) < n:
        rec = [handle.readline() for _ in range(4)]
        if not rec[3]:
            break
        seqs.append(rec[1].rstrip(b"\r\n"))
    return seqs


def file_seqs(path: str, n: Optional[int] = None) -> list:
    """`read_seqs` of a plain or gzipped FASTQ file."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    with (gzip.open(path, "rb") if gz else open(path, "rb")) as fh:
        return read_seqs(fh, n)


def blocks_from_file(path: str, blocks: list) -> list:
    """Sequences of each `(start, count)` block (0-based spots) of one FASTQ file, in one pass.

    Stops reading after the last block, so a block at 2/3 reads 2/3 of the file."""
    want = sorted(range(len(blocks)), key=lambda i: blocks[i][0])
    out = [[] for _ in blocks]
    if not blocks:
        return out
    end = max(s + c for s, c in blocks)
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    with (gzip.open(path, "rb") if gz else open(path, "rb")) as fh:
        i = 0
        while i < end:
            rec = [fh.readline() for _ in range(4)]
            if not rec[3]:
                break
            for j in want:
                s, c = blocks[j]
                if s <= i < s + c:
                    out[j].append(rec[1].rstrip(b"\r\n"))
            i += 1
    return out


def block_starts(total: int, block: int = BLOCK_SPOTS, positions=POSITIONS) -> list:
    """`[(start, count), ...]`: one block per position, 0-based, inside a run of `total` spots."""
    out = []
    for frac in positions:
        start = min(int(total * frac), max(0, total - block))
        out.append((start, min(block, total)))
    return out


def decide(scores: Iterable[float], threshold: float = THRESHOLD) -> tuple:
    """`(sampling, median)`: `stride` when the median of the scored blocks is `threshold` or
    more, else `windows`. With no usable score the answer is `stride`: unknown order is treated
    as possibly sorted, because the cost of a wrong windows sample is a biased sample."""
    ok = [s for s in scores if s == s]
    if not ok:
        return "stride", None
    med = statistics.median(ok)
    return ("stride" if med >= threshold else "windows"), med
