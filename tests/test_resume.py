"""Resumable ENA transfer, stall timeout and lookup retries (issue #51). No network."""

from __future__ import annotations

import csv
import os

import pytest
import requests

from popseq import fetch as fx
from popseq import lookup, routes
from tests.test_fetch import Fake, env  # noqa: F401  (env is a fixture)


def fake_run(script, dest):
    """A stand-in for routes._run: each step is (rc, bytes the file holds after the call)."""
    calls = []
    it = iter(script)

    def run(cmd, timeout=None, **kw):
        calls.append(cmd)
        rc, size = next(it)
        if size is None:
            if os.path.exists(dest):
                os.remove(dest)
        else:
            with open(dest, "wb") as fh:
                fh.write(b"x" * size)
        return rc, 0.1, "", f"curl: ({rc})"

    return run, calls


def test_a_transfer_cut_at_18_is_resumed_with_c_dash_and_finishes(tmp_path, monkeypatch):
    dest = str(tmp_path / "f.fastq.gz")
    run, calls = fake_run([(18, 400), (18, 700), (0, 1000)], dest)
    monkeypatch.setattr(routes, "_run", run)
    rc, _s, _o, err = routes.curl_resume("https://h/f", dest, 1000, sleep=lambda s: None)
    assert rc == 0 and os.path.getsize(dest) == 1000
    assert len(calls) == 3
    assert all("-C" in c and "--speed-time" in c for c in calls)


def test_no_progress_in_a_row_gives_up_and_says_where(tmp_path, monkeypatch):
    dest = str(tmp_path / "f.fastq.gz")
    run, calls = fake_run([(18, 400)] + [(18, 400)] * 9, dest)
    monkeypatch.setattr(routes, "_run", run)
    rc, _s, _o, err = routes.curl_resume(
        "https://h/f", dest, 1000, max_no_progress=3, sleep=lambda s: None
    )
    assert rc == 18 and "gave up" in err and "400 bytes" in err
    assert len(calls) == 4  # one that made progress, then three that did not


def test_progress_resets_the_give_up_counter(tmp_path, monkeypatch):
    dest = str(tmp_path / "f.fastq.gz")
    script = [(18, 100), (18, 100), (18, 200), (18, 200), (18, 300), (0, 1000)]
    run, calls = fake_run(script, dest)
    monkeypatch.setattr(routes, "_run", run)
    rc, *_ = routes.curl_resume("https://h/f", dest, 1000, max_no_progress=2, sleep=lambda s: None)
    assert rc == 0 and len(calls) == 6


def test_an_http_error_is_not_retried(tmp_path, monkeypatch):
    dest = str(tmp_path / "f.fastq.gz")
    run, calls = fake_run([(22, None)], dest)
    monkeypatch.setattr(routes, "_run", run)
    rc, *_ = routes.curl_resume("https://h/f", dest, 1000, sleep=lambda s: None)
    assert rc == 22 and len(calls) == 1


def test_a_complete_file_is_not_requested_again(tmp_path, monkeypatch):
    dest = str(tmp_path / "f.fastq.gz")
    open(dest, "wb").write(b"x" * 1000)
    monkeypatch.setattr(routes, "_run", lambda *a, **k: pytest.fail("curl was called"))
    assert routes.curl_resume("https://h/f", dest, 1000)[0] == 0


def test_a_range_error_on_a_file_that_is_already_complete_is_success(tmp_path, monkeypatch):
    dest = str(tmp_path / "f.fastq.gz")
    run, calls = fake_run([(33, 1000)], dest)
    monkeypatch.setattr(routes, "_run", run)
    rc, *_ = routes.curl_resume("https://h/f", dest, 1000, sleep=lambda s: None)
    assert rc == 0 and len(calls) == 1


def test_a_server_that_will_not_resume_restarts_the_file(tmp_path, monkeypatch):
    dest = str(tmp_path / "f.fastq.gz")
    run, calls = fake_run([(18, 400), (33, None), (0, 1000)], dest)
    monkeypatch.setattr(routes, "_run", run)
    rc, *_ = routes.curl_resume("https://h/f", dest, 1000, sleep=lambda s: None)
    assert rc == 0 and os.path.getsize(dest) == 1000 and len(calls) == 3


# ------------------------------------------------------------------ lookups
class Resp:
    def __init__(self, status, text=""):
        self.status_code, self.text = status, text

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} Error", response=self)


def test_with_retries_504_then_ok():
    seq = [Resp(504), Resp(504), Resp(200, "ok")]

    def fn():
        r = seq.pop(0)
        r.raise_for_status()
        return r

    assert lookup.with_retries(fn, "x", sleep=lambda s: None).text == "ok"
    assert seq == []


def test_with_retries_404_is_not_retried():
    n = []

    def fn():
        n.append(1)
        Resp(404).raise_for_status()

    with pytest.raises(lookup.LookupFailed):
        lookup.with_retries(fn, "x", sleep=lambda s: None)
    assert len(n) == 1


def test_with_retries_gives_up_with_the_last_error():
    def fn():
        raise requests.ConnectionError("reset")

    with pytest.raises(lookup.LookupFailed, match="ConnectionError"):
        lookup.with_retries(fn, "ENA portal search", tries=3, sleep=lambda s: None)


def test_ncbi_counts_error_after_retries_is_error_not_absent(monkeypatch):
    monkeypatch.setattr(lookup.time, "sleep", lambda s: None)
    monkeypatch.setattr(lookup.requests, "get", lambda *a, **k: Resp(503))
    assert lookup.ncbi_counts("SRR1")[0] == "error"


def test_a_dead_ena_lookup_is_transient_exit_13_with_a_row_not_a_traceback(env):
    class Dead(fx.Lookups):
        def ena(self, run):
            raise lookup.LookupFailed("ENA portal search: HTTPError 504: Gateway Time-out")

    r = fx.Fetcher(Fake(), Dead()).fetch(
        "SRR1", env["opt"], crow={"species_taxid": "5141", "sample_accession": "S"}
    )
    assert (r.outcome, r.exit_code) == ("transient", 13)
    assert "504" in r.message
    rows = list(csv.DictReader(open(env["opt"].downloads_csv)))
    assert rows[0]["route"] == "ena_https" and rows[0]["outcome"] == "transfer_failed"
    assert rows[0]["error"].startswith("lookup:")
