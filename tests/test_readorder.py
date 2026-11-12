"""Tests for popseq/readorder.py, the read-order score behind `popseq fetch --sampling auto`."""
from __future__ import annotations

import gzip

from popseq import readorder as ro


def _fq(path, seqs):
    with gzip.open(path, "wb") as fh:
        for i, s in enumerate(seqs):
            fh.write(b"@r%d\n%s\n+\n%s\n" % (i, s, b"I" * len(s)))


def test_block_starts_three_positions_inside_the_run():
    assert ro.block_starts(12_000) == [(0, 4000), (4000, 4000), (8000, 4000)]
    # a short run: blocks are clamped to fit, never past the end
    assert ro.block_starts(5000) == [(0, 4000), (1000, 4000), (1000, 4000)]
    assert ro.block_starts(100) == [(0, 100), (0, 100), (0, 100)]


def test_blocks_from_file_takes_each_block_in_one_pass(tmp_path):
    p = tmp_path / "a_1.fastq.gz"
    _fq(p, [b"A%03d" % i for i in range(30)])
    got = ro.blocks_from_file(str(p), [(0, 2), (10, 3), (28, 5)])
    assert got == [[b"A000", b"A001"], [b"A010", b"A011", b"A012"], [b"A028", b"A029"]]


def test_decide_median_threshold_and_unknown_is_stride():
    assert ro.decide([0.99, 0.02, 0.97]) == ("stride", 0.97)
    assert ro.decide([0.85, 0.02, 0.01]) == ("windows", 0.02)   # sorted at the start only
    assert ro.decide([float("nan"), 0.1, 0.2])[0] == "windows"
    assert ro.decide([float("nan")] * 3) == ("stride", None)
    assert ro.decide([0.8, 0.8, 0.1]) == ("stride", 0.8)


def test_file_seqs_reads_plain_and_gzip(tmp_path):
    g = tmp_path / "x.fastq.gz"
    _fq(g, [b"ACGT", b"TTTT"])
    t = tmp_path / "x.fastq"
    t.write_bytes(b"@a\nACGT\n+\nIIII\n@b\nTTTT\n+\nIIII\n")
    assert ro.file_seqs(str(g)) == ro.file_seqs(str(t)) == [b"ACGT", b"TTTT"]
