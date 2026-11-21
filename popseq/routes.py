"""Route transfer and verification helpers shared by `popseq fetch` and benchmark 7.1.

These were written for the route benchmark (scripts/popseq_bench_routes.py, issue #12) and
moved here so that the fetcher (issue #13) reuses them and does not copy them. The benchmark
script imports them from this module, so benchmark 7.1 stays reproducible.

Rules these helpers keep (spec 002):

* stdin is closed for every subprocess: `prefetch`, `curl` and `ascp` read stdin and will
  eat a caller's loop input.
* The expected count comes from the archive that served the bytes. SRA routes: NCBI
  `spots_with_mates`, stored at plan time. ENA routes: ENA's MD5 per file, and mate 1 count
  equals mate 2 count. ENA `read_count` is never a target.
* Both mates of a run come from one source.
"""

from __future__ import annotations

import hashlib
import os
import socket
import subprocess
import sys
import time

from popseq import downloads as dl

AWS_URL = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"


def _run(cmd, timeout=None, **kw):
    """Run a command, stdin closed (prefetch/curl/ascp read stdin; spec 002)."""
    t0 = time.monotonic()
    try:
        p = subprocess.run(
            cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=timeout, **kw
        )
        rc, err = p.returncode, (p.stderr or "")[-2000:]
        out = p.stdout
    except subprocess.TimeoutExpired:
        rc, err, out = 124, f"timeout after {timeout}s", ""
    return rc, time.monotonic() - t0, out, err


# --- resumable curl (issue #51) -------------------------------------------------------------
#
# `curl --retry` does not retry a transfer that stops part way (curl exit 18, "transfer closed with N
# bytes remaining"), and without `-C -` any retry would start the file again. A stalled connection also
# has no timeout by default. `--retry-all-errors` would help but needs curl 7.71 (the HPCC has 7.61), so
# the retry loop is here.
STALL_SECONDS = 120  # a transfer below 10 kB/s for this long is cut and resumed
STALL_SPEED = 10_000
MAX_NO_PROGRESS = 5  # consecutive attempts that added no bytes before the transfer is given up
# curl exits that mean "the connection broke", not "the URL is wrong": resolve (6), connect (7),
# partial file (18), timeout (28), TLS (35), nothing received (52), send/recv (55, 56), plus our 124.
RESUMABLE_RC = {6, 7, 18, 28, 35, 52, 55, 56, 124}


def _trace(url, notes, how):
    """Say on stderr (the per-run log) how a transfer that needed resuming went; silent when it did not."""
    if notes:
        print(f"curl_resume {os.path.basename(url)}: {len(notes)} resume(s) [{'; '.join(notes)}], {how}", file=sys.stderr)


def curl_resume(
    url,
    dest,
    expected=None,
    timeout=None,
    stall_seconds=STALL_SECONDS,
    max_no_progress=MAX_NO_PROGRESS,
    sleep=time.sleep,
):
    """Download `url` to `dest`, resuming (`-C -`) after a broken or stalled connection.

    Returns (rc, seconds, out, err) like `_run`. `expected` is the byte size ENA reports; with it a
    finished file is recognised without asking the server (a range request on a complete file is an
    error for some servers). Gives up after `max_no_progress` attempts in a row that added no bytes,
    or at once on an error that is not a broken connection (HTTP 404, 403 and so on).
    """
    t0 = time.monotonic()
    size = lambda: os.path.getsize(dest) if os.path.exists(dest) else 0  # noqa: E731
    cmd = [
        "curl", "-sS", "-f", "-L", "-C", "-", "--connect-timeout", "30",
        "--speed-limit", str(STALL_SPEED), "--speed-time", str(stall_seconds),
        "-o", dest, url,
    ]  # fmt: skip
    stuck, last, notes = 0, size(), []
    while True:
        if expected and size() == expected:
            _trace(url, notes, "complete")
            return 0, time.monotonic() - t0, "", ""
        rc, _sec, out, err = _run(cmd, timeout)
        now = size()
        if rc == 0:
            _trace(url, notes, f"done at {now} bytes")
            return 0, time.monotonic() - t0, out, err
        if expected and now == expected:  # 33 or 22 on a range request for a complete file
            _trace(url, notes, f"complete at {now} bytes after rc {rc}")
            return 0, time.monotonic() - t0, out, err
        if rc == 33 or (expected and now > expected):
            # the server will not resume, or the partial file is wrong: start again, once per loop
            if os.path.exists(dest):
                os.remove(dest)
            now = 0
        elif rc not in RESUMABLE_RC:
            return rc, time.monotonic() - t0, out, err
        notes.append(f"curl {rc} at {now} bytes")
        stuck = 0 if now > last else stuck + 1
        last = now
        if stuck >= max_no_progress:
            msg = f"gave up after {len(notes)} attempts ({'; '.join(notes[-3:])}): {err.strip()}"
            return rc, time.monotonic() - t0, out, msg[-500:]
        sleep(min(5 * 2 ** (stuck), 60))


def _version(cmd):
    try:
        p = subprocess.run(
            cmd, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=30
        )
        txt = (p.stdout + p.stderr).strip().splitlines()
        return " ".join(txt[0].split()) if txt else ""
    except (OSError, subprocess.TimeoutExpired):
        return ""


def tool_versions():
    return {
        "curl": _version(["curl", "--version"]),
        "prefetch": _version(["prefetch", "--version"]),
        "parallel-fastq-dump": _version(["parallel-fastq-dump", "--version"]),
        "ascp": _version(["ascp", "--version"]),
    }


def _count_reads(path):
    """Reads in a gzipped FASTQ: lines / 4, decompressed with pigz."""
    p1 = subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE)
    lines = 0
    for chunk in iter(lambda: p1.stdout.read(1 << 22), b""):
        lines += chunk.count(b"\n")
    p1.wait()
    if p1.returncode != 0:
        raise RuntimeError(f"pigz failed on {path} (rc {p1.returncode})")
    return lines // 4


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def transfer(plan, route, wd, timeout, aspera_key):
    """Move the bytes. Returns a dict; never raises on a network failure."""
    run = plan["run"]
    res = {
        "route": route,
        "rc": 0,
        "error": "",
        "bytes": 0,
        "seconds": 0.0,
        "files": [],
        "file_kind": "",
        "unavailable": False,
        "node_type": dl.node_type(),
        "host": socket.gethostname(),
        "slurm_job": os.environ.get("SLURM_JOB_ID", ""),
        "utc_start": dl.utc_now(),
    }
    if route == "aws_odp":
        dest = os.path.join(wd, f"{run}.sra")
        cmd = [
            "curl",
            "-sS",
            "-L",
            "--retry",
            "3",
            "-o",
            dest,
            "-w",
            "%{http_code}",
            AWS_URL.format(run=run),
        ]
        rc, sec, out, err = _run(cmd, timeout)
        code = out.strip()
        if code == "404" or code == "403":
            # S3 answers 403 for a missing key in a bucket that disallows listing.
            res.update(unavailable=True, error=f"AWS HTTP {code}")
            rc = rc or 22
            if os.path.exists(dest):
                os.remove(dest)
        elif rc == 0 and code != "200":
            rc, err = 22, f"AWS HTTP {code}"
        res.update(
            rc=rc,
            seconds=sec,
            error=res["error"] or err.strip(),
            files=[dest] if rc == 0 else [],
            file_kind="sra",
        )
    elif route == "ncbi_prefetch":
        cmd = ["prefetch", run, "-O", wd, "--max-size", "30G"]
        rc, sec, out, err = _run(cmd, timeout)
        found = []
        for root, _, names in os.walk(os.path.join(wd, run)):
            found += [os.path.join(root, n) for n in names if n.endswith((".sra", ".sralite"))]
        txt = out + err
        if rc != 0 and (
            "not found" in txt.lower()
            or "cannot resolve" in txt.lower()
            or "can't find" in txt.lower()
        ):
            res["unavailable"] = True
        res.update(
            rc=rc,
            seconds=sec,
            error=" ".join(txt.split())[-500:] if rc else "",
            files=found if rc == 0 else [],
            file_kind=("sralite" if any(f.endswith(".sralite") for f in found) else "sra"),
        )
    elif route in ("ena_aspera", "ena_https"):
        srcs = [
            s
            for s in plan["fastq_aspera" if route == "ena_aspera" else "fastq_ftp"].split(";")
            if s.endswith(("_1.fastq.gz", "_2.fastq.gz"))
        ]
        total, files, rc, err = 0.0, [], 0, ""
        for src in srcs:
            dest = os.path.join(wd, os.path.basename(src))
            if route == "ena_aspera":
                host, path = src.split(":", 1)
                cmd = [
                    "ascp",
                    "-QT",
                    "-l",
                    "1000m",
                    "-P33001",
                    "-i",
                    aspera_key,
                    f"era-fasp@{host}:{path}",
                    wd,
                ]
                rc, sec, out, e = _run(cmd, timeout)
            else:
                want = (plan.get("ena_bytes") or {}).get(os.path.basename(src))
                rc, sec, out, e = curl_resume(
                    "https://" + src, dest, want, timeout, plan.get("stall_seconds") or STALL_SECONDS
                )
            total += sec
            if rc != 0:
                err = " ".join((out + e).split())[-500:]
                break
            files.append(dest)
        res.update(
            rc=rc, seconds=total, error=err, files=files if rc == 0 else [], file_kind="fastq.gz"
        )
    else:
        raise ValueError(route)
    res["bytes"] = sum(os.path.getsize(f) for f in res["files"] if os.path.exists(f))
    return res


SRA_ROUTES = ("aws_odp", "ncbi_prefetch")


def expected_pairs(route, plan):
    """The pair count a route must deliver, or "" when the check is by MD5.

    Spec 002, adopted 2026-11-07: verify against the archive that served the
    bytes. ENA `read_count` counts reads, not pairs, for most ERR runs (172 of
    200 sampled), so it is not used as a target at all.
    """
    return plan.get("ncbi_spots_with_mates", "") if route in SRA_ROUTES else ""


def judge(route, plan, m1, m2, md5_ok):
    """(verified, error) under the adopted rule.

    SRA routes: each mate equals NCBI `spots_with_mates`, stored at plan time.
    ENA routes: both files match ENA's MD5 exactly, and the mates have equal
    counts.
    """
    if route in SRA_ROUTES:
        exp = expected_pairs(route, plan)
        if not exp:
            return False, "no NCBI spots_with_mates stored in the plan"
        ok = m1 == int(exp) and m2 == int(exp)
        return ok, "" if ok else f"expected {exp} pairs (NCBI), mate1 {m1}, mate2 {m2}"
    ok = md5_ok == "true" and m1 == m2
    return ok, "" if ok else f"md5_ok {md5_ok}, mate1 {m1}, mate2 {m2}"


# ---------------------------------------------------------------------------- fetch helpers
# Added for `popseq fetch` (issue #13).


def _decompress_cmd(path):
    import shutil

    return ["pigz", "-dc", path] if shutil.which("pigz") else ["gzip", "-dc", path]


def read_stats(path):
    """`(reads, bases)` of a gzipped FASTQ in one pass: lines / 4, and the sum of sequence lengths.

    Raises RuntimeError if the file does not decompress or its line count is not a multiple
    of four (a truncated file). A count that passes this is what the verification rule needs.
    """
    p1 = subprocess.Popen(_decompress_cmd(path), stdout=subprocess.PIPE, stdin=subprocess.DEVNULL)
    p2 = subprocess.Popen(
        ["awk", "NR%4==2{b+=length($0)} END{print NR, b+0}"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        text=True,
    )
    p1.stdout.close()
    out, _ = p2.communicate()
    p1.wait()
    if p1.returncode != 0 or p2.returncode != 0:
        raise RuntimeError(
            f"cannot read {path} (decompress rc {p1.returncode}, awk rc {p2.returncode})"
        )
    lines, bases = (int(x) for x in out.split())
    if lines % 4:
        raise RuntimeError(f"{path} is truncated: {lines} lines is not a multiple of 4")
    return lines // 4, bases


def sra_window(sra, start, count, outdir, tmpdir=None, timeout=None):
    """One window of a local `.sra`: spots `start` .. `start+count-1` (0-based) with fastq-dump.

    `fastq-dump` counts spots from 1, so the call is `-N start+1 -X start+count`. Returns
    `(rc, seconds, error)`. Files land in `outdir` as `<RUN>_1.fastq.gz`, `<RUN>_2.fastq.gz`
    (paired) or one file (single-end); the caller finds them by name.
    """
    cmd = [
        "fastq-dump",
        "-N",
        str(start + 1),
        "-X",
        str(start + count),
        "--split-files",
        "--gzip",
        "-O",
        outdir,
        sra,
    ]
    env = dict(os.environ)
    if tmpdir:
        env["TMPDIR"] = tmpdir
    rc, sec, _out, err = _run(cmd, timeout, env=env)
    return rc, sec, " ".join(err.split())[-400:]


def ena_head(url, dest, n_reads, timeout=None, stall_seconds=None):
    """The first `n_reads` reads of an ENA gzip FASTQ, streamed and stopped early (tier `head`).

    ENA is the only route that can stop early: `curl` is read through a gzip decompressor and
    killed as soon as `n_reads` records are complete. Writes a gzip FASTQ to `dest`.
    Returns `(rc, seconds, reads_written, error)`; `reads_written` is below `n_reads` only if
    the run is shorter. A stalled connection (under `STALL_SPEED` bytes/s for `stall_seconds`) makes
    curl exit 28, which ends the stream; a stream that ends because curl failed is an error, not a
    short run (issue #51). The stream cannot be resumed, so the caller retries the whole fetch.
    """
    import gzip
    import zlib

    t0 = time.monotonic()
    proc = subprocess.Popen(
        ["curl", "-sS", "-f", "-L", "--retry", "3", "--connect-timeout", "30",
         "--speed-limit", str(STALL_SPEED), "--speed-time", str(stall_seconds or STALL_SECONDS),
         url],  # fmt: skip
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.DEVNULL,
    )
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    buf, reads, err, eof = b"", 0, "", False
    try:
        with gzip.open(dest, "wb") as out:
            while reads < n_reads:
                chunk = proc.stdout.read(1 << 16)
                if not chunk:
                    eof = True
                    break
                try:
                    buf += dec.decompress(chunk)
                except zlib.error as exc:
                    err = f"gzip stream error: {exc}"
                    break
                lines = buf.split(b"\n")
                buf = lines.pop()  # a partial last line stays in the buffer
                whole = (len(lines) // 4) * 4
                take = min(whole // 4, n_reads - reads)
                out.write(b"\n".join(lines[: take * 4]) + b"\n" if take else b"")
                buf = b"\n".join(lines[take * 4 :] + [buf]) if take * 4 < len(lines) else buf
                reads += take
                if time.monotonic() - t0 > (timeout or 1e9):
                    err = f"timeout after {timeout}s"
                    break
        if eof and not err:
            # the stream ended before n_reads: the run is short, or curl was cut off (stall, reset)
            try:
                crc = proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                crc = 0
            if crc != 0:
                msg = (proc.stderr.read() or b"").decode("utf-8", "replace").strip()
                err = f"curl exit {crc} after {reads} reads: {msg}"[-500:]
    finally:
        proc.kill()
        proc.wait()
    rc = (
        0
        if reads and not err
        else (proc.returncode if proc.returncode not in (None, 0, -9) else 22)
    )
    if reads == 0 and not err:
        err = "no reads read from the stream"
    return rc, time.monotonic() - t0, reads, err
