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


# ------------------------------------------------------------------ the head tier and the CLI options
import gzip  # noqa: E402
import stat  # noqa: E402

from tests.test_fetch import write_fastq  # noqa: E402


def fake_curl(tmp_path, monkeypatch, payload: bytes, exit_code: int):
    """Put a `curl` first on PATH that prints `payload`, records its arguments and exits `exit_code`."""
    data = tmp_path / "payload.bin"
    data.write_bytes(payload)
    args = tmp_path / "args.txt"
    exe = tmp_path / "bin" / "curl"
    exe.parent.mkdir()
    exe.write_text(f'#!/bin/bash\necho "$@" > {args}\ncat {data}\nexit {exit_code}\n')
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{exe.parent}:{os.environ['PATH']}")
    return args


def gz_reads(tmp_path, n):
    p = tmp_path / "src.fastq.gz"
    write_fastq(str(p), n)
    return p.read_bytes()


def test_ena_head_has_a_stall_cutoff_on_its_curl(tmp_path, monkeypatch):
    args = fake_curl(tmp_path, monkeypatch, gz_reads(tmp_path, 20), 0)
    rc, _s, reads, err = routes.ena_head(
        "https://h/f", str(tmp_path / "o.gz"), 10, stall_seconds=45
    )
    a = args.read_text()
    assert rc == 0 and reads == 10 and err == ""
    assert "--speed-time 45" in a and "--speed-limit" in a and "--connect-timeout" in a


def test_a_stream_cut_by_curl_is_an_error_not_a_short_run(tmp_path, monkeypatch):
    whole = gz_reads(tmp_path, 200)
    fake_curl(tmp_path, monkeypatch, whole[: len(whole) // 2], 28)  # stalled part way
    rc, _s, reads, err = routes.ena_head("https://h/f", str(tmp_path / "o.gz"), 1000)
    assert rc != 0 and "curl exit 28" in err


def test_a_run_shorter_than_n_with_a_clean_curl_exit_is_still_a_short_run(tmp_path, monkeypatch):
    fake_curl(tmp_path, monkeypatch, gz_reads(tmp_path, 30), 0)
    rc, _s, reads, err = routes.ena_head("https://h/f", str(tmp_path / "o.gz"), 1000)
    assert rc == 0 and reads == 30 and err == ""


def test_the_stall_setting_reaches_the_resumable_transfer(tmp_path, monkeypatch):
    dest = str(tmp_path / "f.fastq.gz")
    run, calls = fake_run([(0, 10)], dest)
    monkeypatch.setattr(routes, "_run", run)
    routes.curl_resume("https://h/f", dest, 10 + 1, stall_seconds=33, sleep=lambda s: None)
    assert "33" in calls[0] and "--speed-time" in calls[0]


def test_fetch_cli_has_timeout_and_stall_seconds(capsys):
    from popseq import cli

    with pytest.raises(SystemExit):
        cli.main(["fetch", "--help"])
    out = capsys.readouterr().out
    assert "--timeout" in out and "--stall-seconds" in out


# ------------------------------------------------ resume trace on stderr (PopSeq #57)
def test_a_resumed_transfer_says_so_on_stderr(tmp_path, monkeypatch, capsys):
    dest = str(tmp_path / "ERR9_2.fastq.gz")
    run, _ = fake_run([(18, 400), (52, 400), (0, 1000)], dest)
    monkeypatch.setattr(routes, "_run", run)
    assert routes.curl_resume("https://h/x/ERR9_2.fastq.gz", dest, None, sleep=lambda s: None)[0] == 0
    err = capsys.readouterr().err
    assert "curl_resume ERR9_2.fastq.gz: 2 resume(s)" in err
    assert "curl 18 at 400 bytes" in err and "curl 52 at 400 bytes" in err and "done at 1000 bytes" in err


def test_a_transfer_that_needed_no_resume_is_silent(tmp_path, monkeypatch, capsys):
    dest = str(tmp_path / "f.fastq.gz")
    run, _ = fake_run([(0, 1000)], dest)
    monkeypatch.setattr(routes, "_run", run)
    routes.curl_resume("https://h/f", dest, 1000, sleep=lambda s: None)
    assert capsys.readouterr().err == ""


# ------------------------------------------------ aws_odp route: what the curl does with each outcome
def aws_curl(tmp_path, monkeypatch, code: str, exit_code: int, body: bytes = b"SRA"):
    """A `curl` on PATH that writes `body` to the -o file, prints `code` (the -w http code) and exits."""
    args = tmp_path / "args.txt"
    exe = tmp_path / "bin" / "curl"
    exe.parent.mkdir(exist_ok=True)
    exe.write_text(
        '#!/bin/bash\necho "$@" > ' + str(args) + '\nwhile [ $# -gt 0 ]; do [ "$1" = -o ] && out=$2; shift; done\n'
        f'printf %s {body.decode()!r} > "$out"\nprintf %s {code}\nexit {exit_code}\n'
    )
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{exe.parent}:{os.environ['PATH']}")
    return args


def aws_transfer(tmp_path):
    return routes.transfer({"run": "SRR9"}, "aws_odp", str(tmp_path), None, None)


def test_aws_odp_ok_returns_the_sra_file(tmp_path, monkeypatch):
    args = aws_curl(tmp_path, monkeypatch, "200", 0)
    r = aws_transfer(tmp_path)
    assert r["rc"] == 0 and r["files"] == [str(tmp_path / "SRR9.sra")] and r["file_kind"] == "sra"
    assert "--retry" in args.read_text()


@pytest.mark.parametrize("code", ["403", "404"])
def test_aws_odp_403_and_404_are_unavailable_and_leave_no_file(tmp_path, monkeypatch, code):
    aws_curl(tmp_path, monkeypatch, code, 0)
    r = aws_transfer(tmp_path)
    assert r["unavailable"] and r["rc"] != 0 and r["files"] == []
    assert not (tmp_path / "SRR9.sra").exists()


def test_aws_odp_a_cut_transfer_is_a_failure_not_unavailable_and_not_a_file(tmp_path, monkeypatch):
    # curl 18 (partial transfer): today this route has no resume (UNMEASURED how often S3 cuts); it must not look like success
    aws_curl(tmp_path, monkeypatch, "200", 18)
    r = aws_transfer(tmp_path)
    assert r["rc"] == 18 and not r["unavailable"] and r["files"] == []


# ------------------------------------------------ offline end-to-end smoke test (issue #19): real Backend, fake curl, no network
FLAKY_CURL = r"""#!/bin/bash
# serves $PAYLOAD/<basename of url>; the first request for each file is cut halfway (exit 18), a resume (-C -) sends the rest
while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift 2;; --connect-timeout|--speed-limit|--speed-time|-C) shift 2;; -*) shift;; *) url=$1; shift;; esac; done
name=$(basename "$url"); src="$PAYLOAD/$name"; total=$(wc -c < "$src")
if [ ! -e "$STATE/$name" ]; then : > "$STATE/$name"; head -c $((total / 2)) "$src" > "$out"; echo "curl: (18) transfer closed" >&2; exit 18; fi
have=0; [ -e "$out" ] && have=$(wc -c < "$out")
tail -c +$((have + 1)) "$src" >> "$out"
"""


# (the real back-off sleeps 5 s after each cut: about 10 s for the two files)
def test_offline_full_fetch_over_a_flaky_ena_resumes_verifies_and_records(env, monkeypatch, tmp_path):
    from popseq.fetch import Backend, Fetcher
    from tests.test_fetch import md5, plan_for, rows

    payload, state, bindir = tmp_path / "payload", tmp_path / "state", tmp_path / "bin"
    for d in (payload, state, bindir):
        d.mkdir()
    plan = plan_for(ncbi_status="absent", genome=0.01)
    for u in plan.fastq_ftp.split(";"):
        nm = os.path.basename(u)
        write_fastq(str(payload / nm), 500)
        plan.ena_md5[nm] = md5(str(payload / nm))
        plan.ena_bytes[nm] = (payload / nm).stat().st_size
    exe = bindir / "curl"
    exe.write_text(FLAKY_CURL)
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    monkeypatch.setenv("PAYLOAD", str(payload))
    monkeypatch.setenv("STATE", str(state))
    opt = env["opt"]
    r = Fetcher(Backend(), None).fetch(plan.run, opt, plan=plan)
    assert (r.outcome, r.exit_code, r.route, r.reads_kept) == ("ok", 0, "ena_https", 500)
    last = rows(env)[-1]  # earlier rows are the SRA routes tried first (the fake curl serves no .sra)
    assert (last["outcome"], last["route"], last["verified"]) == ("ok", "ena_https", "true")
    assert sorted(os.listdir(os.path.join(opt.cache_root, "5141", "SRR1")))[-1] == "SRR1_2.fastq.gz"
