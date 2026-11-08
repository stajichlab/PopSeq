"""Tests for the benchmark 7.1 driver's verification and bookkeeping (issue #12).

No network. The transfer step is not exercised; `finish` is fed files made
here, which is the part that decides whether a download counts as `ok`.
Needs `pigz` on PATH (it is on the HPCC; the fetcher uses it too).
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import os
import shutil

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "bench", os.path.join(HERE, "..", "scripts", "popseq_bench_routes.py"))
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)

pytestmark = pytest.mark.skipif(shutil.which("pigz") is None, reason="needs pigz")


def _fastq(path, n):
    with gzip.open(path, "wt") as fh:
        for i in range(n):
            fh.write(f"@r{i}\nACGT\n+\nIIII\n")
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def _plan(tmp, run, n, md5s, extra_unpaired=False):
    names = [f"{run}_1.fastq.gz", f"{run}_2.fastq.gz"]
    if extra_unpaired:
        names = [f"{run}.fastq.gz"] + names
        md5s = ["0" * 32] + md5s
    return {"run": run, "taxid": "5476", "read_count": str(n),
            "fastq_ftp": ";".join(f"ftp.example/{x}" for x in names),
            "fastq_md5": ";".join(md5s)}


def _res(wd, files):
    return {"route": "ena_https", "rc": 0, "error": "", "files": files,
            "bytes": sum(os.path.getsize(f) for f in files), "seconds": 2.0,
            "file_kind": "fastq.gz", "unavailable": False, "node_type": "compute",
            "host": "x", "slurm_job": "1", "utc_start": "2026-11-07T00:00:00Z"}


def _go(tmp_path, n_written, n_expected, extra_unpaired=False, corrupt_md5=False):
    wd = str(tmp_path)
    m = [_fastq(os.path.join(wd, f"ERR1_{k}.fastq.gz"), n_written) for k in (1, 2)]
    if corrupt_md5:
        m[1] = "f" * 32
    plan = _plan(tmp_path, "ERR1", n_expected, m, extra_unpaired)
    files = [os.path.join(wd, f"ERR1_{k}.fastq.gz") for k in (1, 2)]
    det, dls = tmp_path / "detail.csv", tmp_path / "downloads.csv"
    d = bench.finish(plan, _res(wd, files), wd, 1, {}, str(det), str(dls))
    (row,) = list(csv.DictReader(open(dls)))
    return d, row


def test_ena_files_with_matching_md5_are_ok(tmp_path):
    d, row = _go(tmp_path, 25, 25)
    assert d["outcome"] == "ok" and d["md5_ok"] == "true"
    assert row["outcome"] == "ok" and row["verified"] == "true"


def test_ena_read_count_is_not_the_target(tmp_path):
    """Adopted rule: ENA read_count counts reads for most ERR runs, so an ENA
    download is judged by MD5, not by read_count (here 2x the pairs)."""
    d, _ = _go(tmp_path, 25, 50)
    assert d["outcome"] == "ok"


@pytest.mark.parametrize("m1,m2,exp,ok", [
    (100, 100, "100", True),
    (100, 100, "200", False),     # ENA-style 2x figure: wrong for NCBI files
    (100, 99, "100", False),
    (100, 100, "", False),        # nothing stored: cannot be verified
])
def test_sra_routes_are_judged_against_ncbi_spots(m1, m2, exp, ok):
    plan = {"ncbi_spots_with_mates": exp}
    assert bench.judge("aws_odp", plan, m1, m2, "")[0] is ok
    assert bench.judge("ncbi_prefetch", plan, m1, m2, "")[0] is ok


def test_ena_mates_must_agree_with_each_other():
    assert bench.judge("ena_https", {}, 10, 10, "true")[0]
    assert not bench.judge("ena_https", {}, 10, 9, "true")[0]
    assert not bench.judge("ena_https", {}, 10, 10, "false")[0]


def test_rejudge_reclassifies_old_rows():
    plan = {"ncbi_spots_with_mates": "9139801"}
    old = {"outcome": "verify_failed", "route": "aws_odp", "reads_mate1": "9139801",
           "reads_mate2": "9139801", "md5_ok": ""}
    assert bench.rejudge(old, plan) == "ok"
    assert bench.rejudge(old, {"ncbi_spots_with_mates": "0"}) == "unavailable_at_ncbi"
    pigz = dict(old, reads_mate1="", reads_mate2="")
    assert bench.rejudge(pigz, plan) == "verify_failed"


def test_md5_mismatch_fails_even_when_counts_agree(tmp_path):
    d, _ = _go(tmp_path, 25, 25, corrupt_md5=True)
    assert d["md5_ok"] == "false" and d["outcome"] == "verify_failed"


def test_an_unpaired_file_listed_by_ena_does_not_shift_the_md5s(tmp_path):
    """ERR2188722 lists <run>.fastq.gz before the mates; pairing is by name."""
    d, _ = _go(tmp_path, 10, 10, extra_unpaired=True)
    assert d["outcome"] == "ok"


def test_transfer_failure_and_unavailable_are_recorded(tmp_path):
    wd = str(tmp_path)
    plan = _plan(tmp_path, "ERR2", 1, ["a", "b"])
    for unavailable, want in ((False, "transfer_failed"), (True, "unavailable_at_ncbi")):
        res = _res(wd, [])
        res.update(rc=22, error="AWS HTTP 404", unavailable=unavailable, route="aws_odp")
        d = bench.finish(plan, res, wd, 1, {}, str(tmp_path / "d.csv"),
                         str(tmp_path / "dl.csv"))
        assert d["outcome"] == want
    rows = list(csv.DictReader(open(tmp_path / "dl.csv")))
    assert [r["outcome"] for r in rows] == ["transfer_failed", "unavailable_at_ncbi"]


def test_quantile():
    assert bench._q([1, 2, 3, 4, 5], 0.5) == 3
    assert bench._q([10, 20], 0.1) == pytest.approx(11)
    assert bench._q([], 0.5) == ""
