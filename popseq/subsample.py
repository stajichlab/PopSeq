"""Subsampling arithmetic and readers for the depth tiers (spec 002, tier `coverage`).

Pure functions. Nothing here touches the network or an archive tool, so every rule in the
spec can be tested on small files.

Spot numbers inside this module are 0-based. The SRA tools count spots from 1: a window that
starts at spot `s` here is `-N s+1 -X s+count` for `fastq-dump`.
"""

from __future__ import annotations

import gzip
import math
from contextlib import ExitStack
from typing import Iterable, Optional


def spots_for_coverage(genome_mb: float, coverage: float, bases_per_spot: float) -> int:
    """Spots needed to reach `coverage` x of a genome. A paired spot counts both mates."""
    if genome_mb <= 0 or coverage <= 0 or bases_per_spot <= 0:
        raise ValueError("genome_mb, coverage and bases_per_spot must be positive")
    return math.ceil(genome_mb * 1e6 * coverage / bases_per_spot)


def windows(total_spots: int, n_spots: int, n_windows: int = 8) -> list[tuple[int, int]]:
    """`(start, count)` for evenly spaced windows covering `n_spots` of a run of `total_spots`.

    Each window holds `w = ceil(n_spots / n_windows)` consecutive spots. Window `i` starts at
    `floor(i * (total_spots - w) / (n_windows - 1))`, so the first window starts at spot 0 and
    the last ends at the last spot. If `n_spots >= total_spots` the whole run is one window:
    the fetcher never pads, repeats or fails a short run.
    """
    if total_spots <= 0:
        raise ValueError("the run has no spots")
    if n_spots <= 0 or n_windows <= 0:
        raise ValueError("n_spots and n_windows must be positive")
    if n_spots >= total_spots:
        return [(0, total_spots)]
    n_windows = min(n_windows, n_spots)
    w = math.ceil(n_spots / n_windows)
    if n_windows == 1:
        return [(0, w)]
    out = [(i * (total_spots - w) // (n_windows - 1), w) for i in range(n_windows)]
    # Rounding w up can make the windows overlap on a run barely longer than the target. Overlap
    # would repeat reads, so shrink such a request to non-overlapping windows instead.
    for (s0, c0), (s1, _) in zip(out, out[1:]):
        if s0 + c0 > s1:
            w = total_spots // n_windows
            return [(i * w, w) for i in range(n_windows)]
    return out


def stride(total_spots: int, target_spots: int) -> int:
    """Keep every `N`th spot, `N = max(1, total // target)`."""
    if target_spots <= 0:
        raise ValueError("target_spots must be positive")
    return max(1, total_spots // target_spots)


def _records(handle) -> Iterable[tuple[bytes, bytes, bytes, bytes]]:
    while True:
        head = handle.readline()
        if not head:
            return
        rec = (head, handle.readline(), handle.readline(), handle.readline())
        if not rec[3]:
            raise ValueError("truncated FASTQ: the last record has fewer than 4 lines")
        yield rec


def select_spots(inputs: list[str], outputs: list[str], keep) -> tuple[int, int, list[int]]:
    """Copy the spots chosen by `keep(i) -> bool` from gzip FASTQ `inputs` to `outputs`.

    All input files are read in step, so mates stay together: spot `i` is taken from every
    file or from none. Every input must hold the same number of spots. Returns
    `(spots_kept, spots_read, bases_kept_per_file)`. Reading stops early when `keep` has
    `stop_after` set (see `Window.stop_after`).
    """
    if len(inputs) != len(outputs):
        raise ValueError("one output per input")
    kept = read = 0
    bases = [0] * len(inputs)
    stop_after = getattr(keep, "stop_after", None)
    with ExitStack() as stack:
        ins = [stack.enter_context(gzip.open(p, "rb")) for p in inputs]
        outs = [stack.enter_context(gzip.open(p, "wb")) for p in outputs]
        iters = [_records(h) for h in ins]
        while True:
            if stop_after is not None and read >= stop_after:
                break
            recs = [next(it, None) for it in iters]
            if all(r is None for r in recs):
                break
            if any(r is None for r in recs):
                raise ValueError("the mates hold different numbers of spots")
            if keep(read):
                for k, rec in enumerate(recs):
                    outs[k].write(b"".join(rec))
                    bases[k] += len(rec[1].rstrip(b"\r\n"))
                kept += 1
            read += 1
    return kept, read, bases


class Window:
    """`keep` predicate for a list of `(start, count)` windows; `stop_after` ends the read."""

    def __init__(self, spans: list[tuple[int, int]]):
        self.spans = sorted(spans)
        self.stop_after = max(s + c for s, c in self.spans)

    def __call__(self, i: int) -> bool:
        return any(s <= i < s + c for s, c in self.spans)


class Every:
    """`keep` predicate for every `n`th spot, starting with spot 0."""

    def __init__(self, n: int, limit: Optional[int] = None):
        self.n = n
        self.stop_after = None if limit is None else limit * n

    def __call__(self, i: int) -> bool:
        return i % self.n == 0


class Evenly:
    """`keep` predicate for exactly `n` of `total` spots, evenly spaced, whatever the read order.

    Spot `i` is kept when `floor((i + 1) * n / total)` is larger than `floor(i * n / total)`, so the
    kept spots are spread over the whole file with gaps differing by at most one, and exactly `n`
    are kept (`n` is capped at `total`). This is the sampling that reaches every part of the genome
    even if the run is sorted by genome position: contiguous windows cannot, however many there
    are, because the genome fraction they reach is the bases kept divided by the depth (found in
    the phase 1 pilot, where two of ten runs came out as a 2.2 and a 2.3 Mb genome).
    """

    def __init__(self, total: int, n: int):
        if total <= 0 or n <= 0:
            raise ValueError("total and n must be positive")
        self.total, self.n = total, min(n, total)
        self.stop_after = None

    def __call__(self, i: int) -> bool:
        return (i + 1) * self.n // self.total > i * self.n // self.total


class First:
    """`keep` predicate for the first `n` spots (tier `head`)."""

    def __init__(self, n: int):
        self.stop_after = n

    def __call__(self, i: int) -> bool:
        return i < self.stop_after
