"""Fixture tests for scripts/popseq_sortcheck.py. No network: simulated reads only."""
from __future__ import annotations

import gzip
import io
import os
import random
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))

import popseq_sortcheck as sc  # noqa: E402

_COMP = bytes.maketrans(b"ACGT", b"TGCA")


def _reads(n=1500, genome_len=2_000_000, read_len=150, depth=20, seed=1):
    """Unsorted: n reads from anywhere in the genome. Sorted: n reads at `depth` from the
    start of the genome, in position order, as the head of a sorted file would be."""
    rng = random.Random(seed)
    genome = bytes(rng.choice(b"ACGT") for _ in range(genome_len))
    anywhere = [rng.randrange(genome_len - read_len) for _ in range(n)]
    span = n * read_len // depth
    head = sorted(rng.randrange(span) for _ in range(n))
    return genome, anywhere, head, read_len, rng


def test_canonical_kmers_are_strand_free():
    s = b"ACGTTGCAAGGCTTACGATCGATCGGA"
    assert sc.canonical_kmers(s, 7) == sc.canonical_kmers(s.translate(_COMP)[::-1], 7)
    assert sc.canonical_kmers(b"ACGNNACG", 3) == {b"ACG"} | sc.canonical_kmers(b"ACG", 3)


def test_random_order_scores_low_and_sorted_order_scores_high():
    genome, anywhere, head, L, rng = _reads()
    unsorted = [genome[s:s + L] for s in anywhere]
    srt = [genome[s:s + L] for s in head]
    assert sc.neighbour_score(unsorted) < 0.05
    assert sc.neighbour_score(srt) > 0.9


def test_sorted_even_at_low_depth_scores_high():
    genome, anywhere, head, L, rng = _reads(depth=3, seed=3)
    assert sc.neighbour_score([genome[s:s + L] for s in head]) > 0.9


def test_sorted_bam_export_with_mixed_strands_still_scores_high():
    # samtools fastq writes reverse-strand reads reverse-complemented back.
    genome, anywhere, head, L, rng = _reads(seed=2)
    srt = []
    for s in head:
        r = genome[s:s + L]
        srt.append(r.translate(_COMP)[::-1] if rng.random() < 0.5 else r)
    assert sc.neighbour_score(srt) > 0.9


def test_read_seqs_takes_whole_records_from_gzip():
    fq = b"".join(b"@r%d\nACGT%d\n+\nIIII\n" % (i, i) for i in range(10))
    h = gzip.GzipFile(fileobj=io.BytesIO(gzip.compress(fq)))
    assert sc.read_seqs(h, 3) == [b"ACGT0", b"ACGT1", b"ACGT2"]
    h = gzip.GzipFile(fileobj=io.BytesIO(gzip.compress(fq)))
    assert len(sc.read_seqs(h, 50)) == 10


def test_summarise_bins_and_excludes():
    rows = [dict(run="A", ena_score="0.05", sra_score_max="0.95"),
            dict(run="B", ena_score="0.4", sra_score_max=""),
            dict(run="P", ena_score="0.99", sra_score_max="0.99")]
    out = {r["score_bin"]: r for r in sc.summarise(rows, exclude={"P"})}
    assert (out["0.0-0.1"]["ena"], out["0.9-1.0"]["sra_max"]) == (1, 1)
    assert (out["0.3-0.5"]["ena"], out["no_score"]["sra_max"]) == (1, 1)
    assert out["0.9-1.0"]["ena"] == 0


def test_truth_pick_takes_all_high_and_some_controls():
    rows = [dict(run=f"H{i}", sra_score_max="0.5") for i in range(3)]
    rows += [dict(run=f"L{i}", sra_score_max="0.01") for i in range(10)]
    rows += [dict(run=f"M{i}", sra_score_max="0.2") for i in range(10)]
    rows += [dict(run="X", sra_score_max="")]
    pick = [r["run"] for r in sc.truth_pick(rows, n_low=2, n_mid=1)]
    assert pick[:3] == ["H0", "H1", "H2"] and len(pick) == 6 and "X" not in pick


def test_truth_spots_coverage_cap_and_default():
    # 25x of 10 Mb at 300 bases per spot = 833,333 spots, under a third of 9M spots
    assert sc.truth_spots(9_000_000, 9_000_000 * 300, 10.0) == 833_333
    assert sc.truth_spots(1_200_000, 1_200_000 * 300, 10.0) == 400_000   # a third of the run
    assert sc.truth_spots(9_000_000, 9_000_000 * 300, None) == 2_000_000


def test_truth_summary_ratio():
    truth = [dict(run="A", tier="windows", bases="1", mash_genome="2200000", mash_cov="100"),
             dict(run="A", tier="stride", bases="1", mash_genome="11000000", mash_cov="20"),
             dict(run="B", tier="windows", bases="1", mash_genome="5", mash_cov="1")]
    (a,) = sc.truth_summary(truth, {"A": dict(sra_score_max="0.99",
                                              sra_scores="0.99;0.1;0.95")})
    assert a["ratio"] == 0.2 and a["genome_stride"] == "11000000"
    assert a["sra_score_median"] == 0.95
