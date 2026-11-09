"""Tests for the download record (issue #14, spec 002 "Performance capture").

No network. Every test writes into pytest's tmp_path, never the real
benchmarks/downloads.csv.
"""
from __future__ import annotations

import csv
import multiprocessing as mp
import subprocess
import sys

import pytest

from popseq import downloads as dl


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def test_header_is_the_spec_header_and_written_once(tmp_path):
    f = tmp_path / "downloads.csv"
    dl.record(f, run="SRR1", route="aws_odp", outcome="ok", verified=True,
              bytes_=10_000_000, seconds=2.0, node="compute")
    dl.record(f, run="SRR2", route="aws_odp", outcome="transfer_failed",
              node="compute", error="curl: (22) 503")
    lines = f.read_text().splitlines()
    assert lines[0] == ",".join(dl.FIELDS)
    assert sum(1 for x in lines if x.startswith("utc,")) == 1
    assert len(lines) == 3


def test_mb_per_s_is_computed_by_the_writer(tmp_path):
    f = tmp_path / "d.csv"
    row = dl.record(f, run="SRR1", route="ncbi_prefetch", outcome="ok",
                    verified=True, bytes_=50_000_000, seconds=4.0, node="login")
    assert row["mb_per_s"] == "12.500"


def test_a_failure_is_recorded_with_empty_unknowns_not_NA(tmp_path):
    f = tmp_path / "d.csv"
    dl.record(f, run="ERR2188722", route="aws_odp", outcome="unavailable_at_ncbi",
              node="compute", error="HTTP 404")
    (r,) = _rows(f)
    assert r["outcome"] == "unavailable_at_ncbi"
    assert r["bytes"] == "" and r["mb_per_s"] == "" and r["verified"] == ""
    assert "NA" not in f.read_text().split(",")


@pytest.mark.parametrize("kw", [
    dict(outcome="succeeded"),                      # not a spec outcome
    dict(outcome="ok"),                             # ok without verification
    dict(outcome="ok", verified=False),
    dict(outcome="verify_failed", verified=True),
    dict(outcome="ok", verified=True, route="ftp"),  # unknown route
    dict(outcome="ok", verified=True, tier="sample"),
    dict(outcome="ok", verified=True, bytes_=-1),
])
def test_rows_the_spec_forbids_are_rejected(tmp_path, kw):
    base = dict(run="SRR1", route="aws_odp", node="compute")
    base.update(kw)
    with pytest.raises(ValueError):
        dl.record(tmp_path / "d.csv", **base)
    assert not (tmp_path / "d.csv").exists()


def test_error_text_stays_on_one_line(tmp_path):
    f = tmp_path / "d.csv"
    dl.record(f, run="SRR1", route="ena_https", outcome="transfer_failed",
              node="compute", error="line one\nline two\r\n\tthree")
    assert len(f.read_text().splitlines()) == 2
    assert _rows(f)[0]["error"] == "line one line two three"


def test_node_type_follows_slurm(monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    assert dl.node_type() == "login"
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    assert dl.node_type() == "compute"


def _writer(args):
    path, i = args
    for j in range(50):
        dl.record(path, run=f"SRR{i}_{j}", route="aws_odp", outcome="transfer_failed",
                  node="compute", error="x" * 200)


def test_concurrent_writers_never_interleave(tmp_path):
    f = str(tmp_path / "d.csv")
    with mp.Pool(8) as pool:
        pool.map(_writer, [(f, i) for i in range(8)])
    rows = _rows(f)
    assert len(rows) == 400
    assert all(r["error"] == "x" * 200 for r in rows)
    assert len({r["run"] for r in rows}) == 400


def test_shell_entry_point(tmp_path):
    f = tmp_path / "d.csv"
    ok = subprocess.run(
        [sys.executable, "-m", "popseq.downloads", "append", "--file", str(f),
         "--run", "SRR1", "--route", "ena_aspera", "--outcome", "ok",
         "--verified", "true", "--bytes", "1000000", "--seconds", "1",
         "--node-type", "compute"], capture_output=True, text=True)
    assert ok.returncode == 0, ok.stderr
    bad = subprocess.run(
        [sys.executable, "-m", "popseq.downloads", "append", "--file", str(f),
         "--run", "SRR1", "--route", "aws_odp", "--outcome", "ok"],
        capture_output=True, text=True)
    assert bad.returncode == 2
    assert len(_rows(f)) == 1


def test_unpaired_is_an_outcome(tmp_path):
    f = tmp_path / "downloads.csv"
    dl.record(f, run="ERR2286963", route="ena_https", outcome="unpaired",
              error="--require-pairs: no _1/_2 at ENA, spots_with_mates 0 at NCBI")
    assert "unpaired" in f.read_text()
