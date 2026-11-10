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
            else:
                cmd = ["curl", "-sS", "-f", "-L", "--retry", "3", "-o", dest, "https://" + src]
            rc, sec, out, e = _run(cmd, timeout)
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


def ena_head(url, dest, n_reads, timeout=None):
    """The first `n_reads` reads of an ENA gzip FASTQ, streamed and stopped early (tier `head`).

    ENA is the only route that can stop early: `curl` is read through a gzip decompressor and
    killed as soon as `n_reads` records are complete. Writes a gzip FASTQ to `dest`.
    Returns `(rc, seconds, reads_written, error)`; `reads_written` is below `n_reads` only if
    the run is shorter.
    """
    import gzip
    import zlib

    t0 = time.monotonic()
    proc = subprocess.Popen(
        ["curl", "-sS", "-f", "-L", "--retry", "3", url],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.DEVNULL,
    )
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    buf, reads, err = b"", 0, ""
    try:
        with gzip.open(dest, "wb") as out:
            while reads < n_reads:
                chunk = proc.stdout.read(1 << 16)
                if not chunk:
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
