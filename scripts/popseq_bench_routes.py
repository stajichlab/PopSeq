#!/usr/bin/env python3
"""Benchmark 7.1 (issue #12): compare the four `full`-tier download routes head to head.

Spec 002 proposes a route order that rests on single-run measurements. This
script fetches the same runs by every route and records what happened,
failures included, so the order can be decided on evidence.

Routes (names as in popseq.downloads.ROUTES):

  aws_odp        curl the .sra from AWS Open Data, then parallel-fastq-dump
  ncbi_prefetch  prefetch the .sra, then parallel-fastq-dump
  ena_aspera     ascp both FASTQ mates from ENA
  ena_https      curl both FASTQ mates from ENA over HTTPS

Subcommands:

  select   choose the benchmark runs from the catalog, look up their ENA file
           sizes, MD5s and read counts once, and write the plan CSV.
           Needs internet. Run it in a job (scripts/bench_routes.slurm select).
  attempt  one (run, route) attempt: transfer, then convert and verify.
           `--stage transfer` only moves bytes (for the login-node arm, which
           must not do CPU work); `--stage finish` later converts and verifies
           what a transfer stage left behind; `--stage all` does both.
  summary  per-route n, median, p10, p90, failures, from the detail CSV.

Every attempt appends one row to benchmarks/downloads.csv (issue #14) and one
richer row to the detail CSV. Verification is the spec 002 rule: each mate's
read count must equal ENA `read_count`. ENA files are also checked against
ENA's MD5. Nothing is reported as `ok` without the read-count check.

Timing definitions, so the numbers are not compared wrongly:

  transfer_s   wall time of the network transfer only
  convert_s    .sra to gzipped FASTQ (0 for the ENA routes, which ship FASTQ)
  to_fastq_s   transfer_s + convert_s: time to a FASTQ ready to verify.
               This is the fair end-to-end figure. MB/s is not, because an
               .sra and a pair of FASTQ.gz files are different byte counts.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import random
import shutil
import socket
import statistics
import subprocess
import sys
import time

import requests

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from popseq import downloads as dl  # noqa: E402
from popseq.query import ROOT  # noqa: E402

ROUTES = ["aws_odp", "ncbi_prefetch", "ena_aspera", "ena_https"]
ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
AWS_URL = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"
PLAN_FIELDS = ["run", "taxid", "species", "study", "read_count", "base_count",
               "fastq_bytes", "fastq_md5", "fastq_ftp", "fastq_aspera", "role",
               "ncbi_spots", "ncbi_spots_with_mates"]
DETAIL_FIELDS = ["utc_start", "run", "route", "node_type", "host", "slurm_job",
                 "transfer_bytes", "transfer_s", "transfer_mb_per_s", "convert_s",
                 "to_fastq_s", "verify_s", "file_kind", "read_count_expected",
                 "reads_mate1", "reads_mate2", "md5_ok", "exit_code", "outcome",
                 "error", "tools"]
# A known `unavailable_at_ncbi` run (spec 002, MATPredict). It is fetched as a
# deliberate failure, so the failure path is exercised and recorded (issue #14).
CONTROL_RUN = "ERR2188722"


def _open(path):
    """Read .csv, .csv.gz or .csv.zst transparently (AGENTS.md section 8)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    if path.endswith(".zst"):
        import zstandard
        return io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(
            open(path, "rb")), newline="")
    return open(path, newline="")


# --------------------------------------------------------------------------- select

def ena_lookup(accessions):
    """read_count, fastq_bytes, fastq_md5, URLs for each run, from the ENA portal."""
    out = {}
    fields = "run_accession,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera"
    for i in range(0, len(accessions), 200):
        chunk = accessions[i:i + 200]
        r = requests.post(ENA_SEARCH, data={
            "result": "read_run", "includeAccessions": ",".join(chunk),
            "fields": fields, "format": "tsv", "limit": 0}, timeout=120)
        r.raise_for_status()
        for row in csv.DictReader(io.StringIO(r.text), delimiter="\t"):
            out[row["run_accession"]] = row
    return out


def cmd_select(a):
    rng = random.Random(a.seed)
    by_bin = {}
    edges = [math.log10(a.min_bytes) + k * (math.log10(a.max_bytes) - math.log10(a.min_bytes))
             / (a.n - 1) for k in range(a.n)]
    targets = [10 ** e for e in edges]
    with _open(a.runs) as fh:
        for r in csv.DictReader(fh):
            if (r["instrument_platform"] != "ILLUMINA" or r["library_layout"] != "PAIRED"
                    or r["exclusion"] or r["fastq_aspera"].count(";") != 1
                    or not r["base_count"]):
                continue
            acc = r["run_accession"]
            prefix = acc[:3]
            if prefix not in ("SRR", "ERR"):
                continue
            est = int(r["base_count"]) / 3.5  # spec 002 estimate, only to pre-bin
            k = min(range(a.n), key=lambda i: abs(math.log10(max(est, 1)) - edges[i]))
            if abs(math.log10(max(est, 1)) - edges[k]) > 0.15:
                continue
            by_bin.setdefault((k, prefix), []).append(r)
    cand = []
    for key in sorted(by_bin):
        rows = by_bin[key]
        rng.shuffle(rows)
        cand.extend(rows[:a.per_bin])
    meta = ena_lookup([r["run_accession"] for r in cand] + [CONTROL_RUN])
    chosen, studies = [], set()
    for k, tgt in enumerate(targets):
        want = "ERR" if k % 2 else "SRR"
        best = None
        for r in cand:
            m = meta.get(r["run_accession"])
            if (not m or r["run_accession"][:3] != want or r["study_accession"] in studies
                    or m["fastq_bytes"].count(";") != 1):
                continue
            size = sum(int(x) for x in m["fastq_bytes"].split(";"))
            d = abs(math.log10(size) - math.log10(tgt))
            if best is None or d < best[0]:
                best = (d, r, m, size)
        if best is None:
            sys.exit(f"no candidate for target {tgt / 1e9:.2f} GB {want}")
        _, r, m, size = best
        studies.add(r["study_accession"])
        chosen.append({"run": r["run_accession"], "taxid": r["species_taxid"],
                       "species": r["species_name"], "study": r["study_accession"],
                       "read_count": m["read_count"], "base_count": m["base_count"],
                       "fastq_bytes": size, "fastq_md5": m["fastq_md5"],
                       "fastq_ftp": m["fastq_ftp"], "fastq_aspera": m["fastq_aspera"],
                       "role": "benchmark"})
    c = meta.get(CONTROL_RUN, {})
    chosen.append({"run": CONTROL_RUN, "taxid": "", "species": "", "study": "",
                   "read_count": c.get("read_count", ""), "base_count": c.get("base_count", ""),
                   "fastq_bytes": sum(int(x) for x in c["fastq_bytes"].split(";"))
                   if c.get("fastq_bytes") else "",
                   "fastq_md5": c.get("fastq_md5", ""), "fastq_ftp": c.get("fastq_ftp", ""),
                   "fastq_aspera": c.get("fastq_aspera", ""), "role": "control_unavailable"})
    _add_ncbi(chosen)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=PLAN_FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows(chosen)
    for r in chosen:
        print(f"{r['run']}\t{r['role']}\t{int(r['fastq_bytes'] or 0) / 1e9:.2f} GB\t"
              f"{r['species']}", file=sys.stderr)


def _add_ncbi(rows):
    """Store NCBI's counts at plan time: compute nodes may have no internet."""
    for r in rows:
        spots, mates = ncbi_spots(r["run"])
        r["ncbi_spots"] = "" if spots is None else spots
        r["ncbi_spots_with_mates"] = "" if mates is None else mates


def cmd_annotate(a):
    """Add the NCBI columns to a plan written before they existed."""
    with _open(a.plan) as fh:
        rows = list(csv.DictReader(fh))
    _add_ncbi(rows)
    tmp = a.plan + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=PLAN_FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, a.plan)


# --------------------------------------------------------------------------- attempt

def _run(cmd, timeout=None, **kw):
    """Run a command, stdin closed (prefetch/curl/ascp read stdin; spec 002)."""
    t0 = time.monotonic()
    try:
        p = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                           timeout=timeout, **kw)
        rc, err = p.returncode, (p.stderr or "")[-2000:]
        out = p.stdout
    except subprocess.TimeoutExpired:
        rc, err, out = 124, f"timeout after {timeout}s", ""
    return rc, time.monotonic() - t0, out, err


def _version(cmd):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, stdin=subprocess.DEVNULL,
                           timeout=30)
        txt = (p.stdout + p.stderr).strip().splitlines()
        return " ".join(txt[0].split()) if txt else ""
    except (OSError, subprocess.TimeoutExpired):
        return ""


def tool_versions():
    return {"curl": _version(["curl", "--version"]),
            "prefetch": _version(["prefetch", "--version"]),
            "parallel-fastq-dump": _version(["parallel-fastq-dump", "--version"]),
            "ascp": _version(["ascp", "--version"])}


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
    res = {"route": route, "rc": 0, "error": "", "bytes": 0, "seconds": 0.0,
           "files": [], "file_kind": "", "unavailable": False,
           "node_type": dl.node_type(), "host": socket.gethostname(),
           "slurm_job": os.environ.get("SLURM_JOB_ID", ""),
           "utc_start": dl.utc_now()}
    if route == "aws_odp":
        dest = os.path.join(wd, f"{run}.sra")
        cmd = ["curl", "-sS", "-L", "--retry", "3", "-o", dest, "-w", "%{http_code}",
               AWS_URL.format(run=run)]
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
        res.update(rc=rc, seconds=sec, error=res["error"] or err.strip(),
                   files=[dest] if rc == 0 else [], file_kind="sra")
    elif route == "ncbi_prefetch":
        cmd = ["prefetch", run, "-O", wd, "--max-size", "30G"]
        rc, sec, out, err = _run(cmd, timeout)
        found = []
        for root, _, names in os.walk(os.path.join(wd, run)):
            found += [os.path.join(root, n) for n in names
                      if n.endswith((".sra", ".sralite"))]
        txt = (out + err)
        if rc != 0 and ("not found" in txt.lower() or "cannot resolve" in txt.lower()
                        or "can't find" in txt.lower()):
            res["unavailable"] = True
        res.update(rc=rc, seconds=sec, error=" ".join(txt.split())[-500:] if rc else "",
                   files=found if rc == 0 else [],
                   file_kind=("sralite" if any(f.endswith(".sralite") for f in found)
                              else "sra"))
    elif route in ("ena_aspera", "ena_https"):
        srcs = [s for s in plan["fastq_aspera" if route == "ena_aspera"
                                 else "fastq_ftp"].split(";")
                if s.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
        total, files, rc, err = 0.0, [], 0, ""
        for src in srcs:
            dest = os.path.join(wd, os.path.basename(src))
            if route == "ena_aspera":
                host, path = src.split(":", 1)
                cmd = ["ascp", "-QT", "-l", "1000m", "-P33001", "-i", aspera_key,
                       f"era-fasp@{host}:{path}", wd]
            else:
                cmd = ["curl", "-sS", "-f", "-L", "--retry", "3", "-o", dest,
                       "https://" + src]
            rc, sec, out, e = _run(cmd, timeout)
            total += sec
            if rc != 0:
                err = " ".join((out + e).split())[-500:]
                break
            files.append(dest)
        res.update(rc=rc, seconds=total, error=err, files=files if rc == 0 else [],
                   file_kind="fastq.gz")
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


def finish(plan, res, wd, threads, tools, detail_csv, downloads_csv):
    """Convert and verify what a transfer left behind, then record both rows."""
    run, route = plan["run"], res["route"]
    convert_s = verify_s = 0.0
    m1 = m2 = ""
    md5_ok = ""
    outcome, error, verified = "ok", "", None
    try:
        if res["rc"] != 0:
            outcome = "unavailable_at_ncbi" if res["unavailable"] else "transfer_failed"
            error = res["error"] or f"exit {res['rc']}"
            raise StopIteration
        if res["file_kind"] in ("sra", "sralite"):
            if (plan.get("ncbi_spots_with_mates") or "") == "0":
                # NCBI holds the run unpaired: it cannot give split mates. Pairs
                # must come from ENA (spec 002, adopted rule).
                outcome = "unavailable_at_ncbi"
                error = (f"not paired at NCBI: spots {plan.get('ncbi_spots')}, "
                         "spots_with_mates 0")
                raise StopIteration
            sra = res["files"][0]
            scratch = os.environ.get("SCRATCH") or wd
            cmd = ["parallel-fastq-dump", "--sra-id", sra, "--threads", str(threads),
                   "--outdir", wd, "--split-files", "--gzip", "--tmpdir", scratch]
            rc, convert_s, out, err = _run(cmd)
            if rc != 0:
                outcome, error = "transfer_failed", (
                    f"parallel-fastq-dump exit {rc}: " + " ".join(err.split())[-400:])
                raise StopIteration
            mates = [os.path.join(wd, f"{run}_1.fastq.gz"), os.path.join(wd, f"{run}_2.fastq.gz")]
        else:
            # ENA can list an unpaired file next to the mates (ERR2188722 has
            # three). Pair URLs with MD5s by position, then pick _1 and _2 by name.
            want = dict(zip((os.path.basename(u) for u in plan["fastq_ftp"].split(";")),
                            plan["fastq_md5"].split(";")))
            mates = [os.path.join(wd, f"{run}_{k}.fastq.gz") for k in (1, 2)]
            md5_ok = "true" if all(
                os.path.exists(f) and _md5(f) == want.get(os.path.basename(f))
                for f in mates) else "false"
        t0 = time.monotonic()
        m1, m2 = (_count_reads(f) for f in mates)
        verify_s = time.monotonic() - t0
        verified, error = judge(route, plan, m1, m2, md5_ok)
        if not verified:
            outcome = "verify_failed"
    except StopIteration:
        pass
    except (RuntimeError, OSError) as exc:
        outcome, error, verified = "verify_failed", str(exc), False
    to_fastq = res["seconds"] + convert_s
    detail = {
        "utc_start": res["utc_start"], "run": run, "route": route,
        "node_type": res["node_type"], "host": res["host"], "slurm_job": res["slurm_job"],
        "transfer_bytes": res["bytes"], "transfer_s": f"{res['seconds']:.2f}",
        "transfer_mb_per_s": (f"{res['bytes'] / res['seconds'] / 1e6:.3f}"
                              if res["bytes"] and res["seconds"] else ""),
        "convert_s": f"{convert_s:.2f}", "to_fastq_s": f"{to_fastq:.2f}",
        "verify_s": f"{verify_s:.2f}", "file_kind": res["file_kind"],
        "read_count_expected": expected_pairs(route, plan), "reads_mate1": m1, "reads_mate2": m2,
        "md5_ok": md5_ok, "exit_code": res["rc"], "outcome": outcome,
        "error": " ".join(error.split())[:500], "tools": json.dumps(tools, sort_keys=True)}
    _append_detail(detail_csv, detail)
    dl.record(downloads_csv, run=run, taxid=plan["taxid"] or None, tier="full",
              route=route, node=res["node_type"], bytes_=res["bytes"] or None,
              seconds=res["seconds"] or None,
              verified=verified, outcome=outcome, error=error)
    return detail


def _append_detail(path, row):
    import fcntl
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o664)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        buf = io.StringIO()
        w = csv.DictWriter(buf, fieldnames=DETAIL_FIELDS, lineterminator="\n")
        if os.fstat(fd).st_size == 0:
            w.writeheader()
        w.writerow(row)
        os.write(fd, buf.getvalue().encode())
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


def _plan_row(plan_csv, run):
    with _open(plan_csv) as fh:
        for r in csv.DictReader(fh):
            if r["run"] == run:
                return r
    sys.exit(f"{run} is not in {plan_csv}")


def cmd_attempt(a):
    plan = _plan_row(a.plan, a.run)
    wd = os.path.join(a.workdir, a.run, a.route)
    state = os.path.join(wd, "transfer.json")
    tools = tool_versions()
    if a.stage in ("transfer", "all"):
        # A fresh start every time: never resume onto another attempt's file (spec 002).
        shutil.rmtree(wd, ignore_errors=True)
        os.makedirs(wd)
        res = transfer(plan, a.route, wd, a.timeout, a.aspera_key)
        with open(state, "w") as fh:
            json.dump(res, fh)
        print(f"{a.run} {a.route} transfer rc={res['rc']} {res['bytes']} B "
              f"{res['seconds']:.1f} s", file=sys.stderr)
    if a.stage in ("finish", "all"):
        with open(state) as fh:
            res = json.load(fh)
        d = finish(plan, res, wd, a.threads, tools, a.detail, a.downloads)
        print(f"{a.run} {a.route} {d['outcome']} to_fastq {d['to_fastq_s']} s "
              f"{d['error']}", file=sys.stderr)
        if not a.keep:
            shutil.rmtree(wd, ignore_errors=True)


# --------------------------------------------------------------------------- countcheck

NCBI_RUNINFO = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={run}"


def ncbi_spots(run):
    """(spots, spots_with_mates) from NCBI runinfo, or (None, None) if absent."""
    try:
        r = requests.get(NCBI_RUNINFO.format(run=run), timeout=60)
        rows = list(csv.DictReader(io.StringIO(r.text)))
    except (requests.RequestException, csv.Error):
        return None, None
    if not rows or not rows[0].get("spots"):
        return None, None
    return int(rows[0]["spots"] or 0), int(rows[0]["spots_with_mates"] or 0)


def cmd_countcheck(a):
    """Does ENA `read_count` count pairs or reads? Measure it on a random sample.

    Spec 002 verifies each mate against ENA `read_count` as if it counted
    pairs. Benchmark 7.1 found ERR runs where it counts reads instead. This
    samples paired Illumina runs per accession prefix and compares ENA
    `read_count` with NCBI `spots`, so the size of the problem is measured,
    not guessed.
    """
    rng = random.Random(a.seed)
    pool = {}
    with _open(a.runs) as fh:
        for r in csv.DictReader(fh):
            if (r["instrument_platform"] == "ILLUMINA" and r["library_layout"] == "PAIRED"
                    and r["read_count"]):
                pool.setdefault(r["run_accession"][:3], []).append(r)
    w = csv.DictWriter(sys.stdout, lineterminator="\n", fieldnames=[
        "run", "prefix", "center_name", "ena_read_count", "ncbi_spots",
        "ncbi_spots_with_mates", "ratio", "class"])
    w.writeheader()
    for prefix in sorted(pool):
        rows = pool[prefix]
        for r in rng.sample(rows, min(a.n, len(rows))):
            spots, mates = ncbi_spots(r["run_accession"])
            ena = int(r["read_count"])
            if spots is None or spots == 0:
                cls, ratio = "absent_at_ncbi", ""
            else:
                ratio = ena / spots
                cls = ("equal" if ena == spots else
                       "ena_is_2x" if abs(ratio - 2) < 0.01 else "other")
                if mates == 0:
                    cls += "_ncbi_unpaired"
                ratio = f"{ratio:.4f}"
            w.writerow({"run": r["run_accession"], "prefix": prefix,
                        "center_name": r["center_name"], "ena_read_count": ena,
                        "ncbi_spots": spots if spots is not None else "",
                        "ncbi_spots_with_mates": mates if mates is not None else "",
                        "ratio": ratio, "class": cls})
            sys.stdout.flush()


# --------------------------------------------------------------------------- summary

def _q(xs, p):
    xs = sorted(xs)
    if not xs:
        return ""
    k = (len(xs) - 1) * p
    lo, hi = math.floor(k), math.ceil(k)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def rejudge(row, plan):
    """The outcome of a detail row under the adopted rule.

    Rows written before 2026-11-07 ~03:00 UTC were judged against ENA
    `read_count`. Their mate counts and MD5 result are in the row, so they can
    be judged again without refetching. Rows with no counts keep their outcome.
    """
    if row["outcome"] not in ("ok", "verify_failed"):
        return row["outcome"]
    if row["route"] in SRA_ROUTES and plan.get("ncbi_spots_with_mates") == "0":
        return "unavailable_at_ncbi"
    if row["reads_mate1"] == "" or row["reads_mate2"] == "":
        return row["outcome"]
    ok, _ = judge(row["route"], plan, int(row["reads_mate1"]), int(row["reads_mate2"]),
                  row["md5_ok"])
    return "ok" if ok else "verify_failed"


def cmd_summary(a):
    rows = list(csv.DictReader(_open(a.detail)))
    if a.plan:
        with _open(a.plan) as fh:
            plans = {r["run"]: r for r in csv.DictReader(fh)}
        for r in rows:
            r["outcome"] = rejudge(r, plans[r["run"]])
    out = []
    for node in sorted({r["node_type"] for r in rows}):
        for route in ROUTES:
            rs = [r for r in rows if r["route"] == route and r["node_type"] == node]
            if not rs:
                continue
            ok = [r for r in rs if r["outcome"] == "ok"]
            mbps = [float(r["transfer_mb_per_s"]) for r in ok if r["transfer_mb_per_s"]]
            gbpm = [int(r["transfer_bytes"]) / 1e9 / (float(r["to_fastq_s"]) / 60)
                    for r in ok if float(r["to_fastq_s"]) > 0]
            out.append({
                "node_type": node, "route": route, "attempts": len(rs), "ok": len(ok),
                "verify_failed": sum(r["outcome"] == "verify_failed" for r in rs),
                "transfer_failed": sum(r["outcome"] == "transfer_failed" for r in rs),
                "unavailable_at_ncbi": sum(r["outcome"] == "unavailable_at_ncbi" for r in rs),
                "transfer_mb_per_s_median": f"{statistics.median(mbps):.2f}" if mbps else "",
                "transfer_mb_per_s_p10": f"{_q(mbps, .1):.2f}" if mbps else "",
                "transfer_mb_per_s_p90": f"{_q(mbps, .9):.2f}" if mbps else "",
                "to_fastq_gb_per_min_median": f"{statistics.median(gbpm):.3f}" if gbpm else "",
            })
    w = csv.DictWriter(sys.stdout, fieldnames=list(out[0]) if out else ["route"],
                       lineterminator="\n")
    w.writeheader()
    w.writerows(out)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--runs", default=os.path.join(ROOT, "catalog", "runs",
                                                  "fungi_all_runs.csv.gz"))
    s.add_argument("--out", required=True)
    s.add_argument("--n", type=int, default=10)
    s.add_argument("--min-bytes", type=float, default=0.5e9)
    s.add_argument("--max-bytes", type=float, default=10e9)
    s.add_argument("--per-bin", type=int, default=15)
    s.add_argument("--seed", type=int, default=71)
    t = sub.add_parser("attempt")
    t.add_argument("--plan", required=True)
    t.add_argument("--run", required=True)
    t.add_argument("--route", required=True, choices=ROUTES)
    t.add_argument("--stage", choices=["transfer", "finish", "all"], default="all")
    t.add_argument("--workdir", required=True, help="on /bigdata, not node-local")
    t.add_argument("--detail", required=True)
    t.add_argument("--downloads", default=dl.DOWNLOADS_CSV)
    t.add_argument("--threads", type=int, default=8)
    t.add_argument("--timeout", type=int, default=4 * 3600, help="per transfer, seconds")
    t.add_argument("--aspera-key", default=os.environ.get("ASPERA_KEY", ""))
    t.add_argument("--keep", action="store_true", help="keep the files (debugging)")
    m = sub.add_parser("summary")
    m.add_argument("--detail", required=True)
    m.add_argument("--plan", help="re-judge every row under the adopted rule")
    n = sub.add_parser("annotate", help="add NCBI spot counts to an existing plan")
    n.add_argument("--plan", required=True)
    c = sub.add_parser("countcheck", help="ENA read_count vs NCBI spots on a sample")
    c.add_argument("--runs", default=os.path.join(ROOT, "catalog", "runs",
                                                  "fungi_all_runs.csv.gz"))
    c.add_argument("--n", type=int, default=200, help="runs per accession prefix")
    c.add_argument("--seed", type=int, default=72)
    a = p.parse_args(argv)
    {"select": cmd_select, "attempt": cmd_attempt, "summary": cmd_summary,
     "countcheck": cmd_countcheck, "annotate": cmd_annotate}[a.cmd](a)


if __name__ == "__main__":
    main()
