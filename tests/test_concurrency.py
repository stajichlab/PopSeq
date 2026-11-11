"""Fixture tests for scripts/popseq_concurrency.py. No network, no catalog contents."""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))

import popseq_concurrency as pc  # noqa: E402


def _t(level, start, end, exit=0, bytes_=100e6, sec=10, node="n1"):
    return dict(level=str(level), run="R", taxid="1", node=node, start_epoch=str(start),
                end_epoch=str(end), exit=str(exit), bytes=str(bytes_), seconds_fetch=str(sec))


def test_in_flight_counts_overlap_time_weighted():
    # A runs 0-10. B runs 5-15. C runs 20-30, alone.
    ts = [_t(5, 0, 10), _t(5, 5, 15), _t(5, 20, 30)]
    assert pc.in_flight(ts) == [1.5, 1.5, 1.0]


def test_summary_per_level_stream_rate_failures_and_aggregate():
    ts = [_t(1, 0, 10, sec=10), _t(1, 10, 20, sec=5),
          _t(2, 100, 110, sec=10), _t(2, 100, 110, exit=13)]
    lv1, lv2 = pc.summarise(ts)
    assert (lv1["level"], lv1["tasks"], lv1["ok"], lv1["failed"]) == (1, 2, 2, 0)
    assert lv1["realised_concurrency_median"] == 1.0
    assert lv1["stream_mb_s_median"] == 15.0          # 10 and 20 MB/s
    assert lv1["aggregate_mb_s"] == 10.0              # 200 MB over 20 s
    assert (lv2["ok"], lv2["failed"], lv2["realised_concurrency_median"]) == (1, 1, 2.0)


def test_plan_size_band_ncbi_only_and_species_spread():
    rows = [dict(run_accession=f"S{i}", species_taxid=str(i % 2), species_name="x",
                 base_count="1000000000", ncbi_status="in_ncbi", fastq_ftp="f")
            for i in range(6)]
    rows += [dict(run_accession="BIG", species_taxid="0", species_name="x",
                  base_count="9000000000", ncbi_status="in_ncbi", fastq_ftp="f"),
             dict(run_accession="ENA", species_taxid="0", species_name="x",
                  base_count="1000000000", ncbi_status="ena_only", fastq_ftp="f")]
    pool = pc.plan(rows, 4, 0.8, 1.6)
    accs = {r["run_accession"] for r in pool}
    assert len(pool) == 4 and not accs & {"BIG", "ENA"}
    assert sorted(r["species_taxid"] for r in pool) == ["0", "0", "1", "1"]
