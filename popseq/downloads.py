"""The download record: one row per fetch attempt in `benchmarks/downloads.csv`.

Spec 002 ("Performance capture") and issue #14 require a row for **every**
attempt, failures included. A route that works 95 percent of the time and one
that works 100 percent of the time look the same in a table of successes.

Columns, in this order (the spec fixes them; do not reorder):

    utc,run,taxid,tier,route,node_type,bytes,seconds,mb_per_s,verified,outcome,error

* `outcome` is one of OUTCOMES. `unavailable_at_ncbi` is a real state, not a
  failure: the run exists at ENA, but neither AWS Open Data nor `prefetch` has it.
  `unpaired` means pairs were required (`--require-pairs`) and no archive serves them.
* `seconds` is the transfer wall time. `mb_per_s` is `bytes / seconds / 1e6`,
  computed here so that every writer computes it the same way.
* An empty field means unknown. `NA` is never written (AGENTS.md section 8).

Many jobs on many nodes append to one file. Each row goes out as a single
`write()` on a descriptor opened with `O_APPEND`, under an exclusive `fcntl`
lock. GPFS (`/bigdata`) honours both across nodes, so rows never interleave.

Shell use, from a fetch script:

    uv run --project "$POPSEQ_ROOT" python -m popseq.downloads append \\
        --run SRR1553610 --taxid 5476 --tier full --route aws_odp \\
        --bytes 123456 --seconds 12.3 --verified true --outcome ok
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import fcntl
import io
import os
import socket
import stat
import sys
from typing import Optional

from popseq.query import ROOT

DOWNLOADS_CSV = os.path.join(ROOT, "benchmarks", "downloads.csv")

FIELDS = ["utc", "run", "taxid", "tier", "route", "node_type", "bytes", "seconds",
          "mb_per_s", "verified", "outcome", "error"]
OUTCOMES = ("ok", "verify_failed", "transfer_failed", "unavailable_at_ncbi", "unpaired")
TIERS = ("head", "windows", "stride", "full")
# Route names. The four candidate routes of spec 002 for tier `full`, plus the
# conversion-only step, which a benchmark may time on its own.
ROUTES = ("aws_odp", "ncbi_prefetch", "ena_aspera", "ena_https", "fasterq_dump")
NODE_TYPES = ("login", "compute")


def node_type() -> str:
    """`compute` inside a SLURM allocation, `login` outside one.

    A login-node session can also hold an allocation (`srun --pty`), and that
    counts as compute, which is the correct answer for a network measurement.
    """
    return "compute" if os.environ.get("SLURM_JOB_ID") else "login"


def utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _fmt_bool(v) -> str:
    if v is None or v == "":
        return ""
    if isinstance(v, str):
        s = v.strip().lower()
        if s in ("true", "1", "yes"):
            return "true"
        if s in ("false", "0", "no"):
            return "false"
        raise ValueError(f"verified must be true, false or empty, not {v!r}")
    return "true" if v else "false"


def make_row(*, run: str, outcome: str, route: str, tier: str = "full",
             taxid=None, node: Optional[str] = None, bytes_: Optional[int] = None,
             seconds: Optional[float] = None, verified=None, error: str = "",
             utc: Optional[str] = None) -> dict:
    """Validate one attempt and return it as a row of strings.

    Raises ValueError on anything the spec does not allow, so that a bad row
    is rejected at the writer, not discovered later by a reader.
    """
    if not run:
        raise ValueError("run is required")
    if outcome not in OUTCOMES:
        raise ValueError(f"outcome must be one of {OUTCOMES}, not {outcome!r}")
    if tier not in TIERS:
        raise ValueError(f"tier must be one of {TIERS}, not {tier!r}")
    if route not in ROUTES:
        raise ValueError(f"route must be one of {ROUTES}, not {route!r}")
    node = node or node_type()
    if node not in NODE_TYPES:
        raise ValueError(f"node_type must be one of {NODE_TYPES}, not {node!r}")
    ver = _fmt_bool(verified)
    if outcome == "ok" and ver != "true":
        raise ValueError("outcome ok requires verified=true: an unverified "
                         "download has not succeeded (AGENTS.md section 12)")
    if outcome == "verify_failed" and ver != "false":
        raise ValueError("outcome verify_failed requires verified=false")
    if bytes_ is not None and int(bytes_) < 0:
        raise ValueError("bytes must be >= 0")
    if seconds is not None and float(seconds) < 0:
        raise ValueError("seconds must be >= 0")
    mbps = ""
    if bytes_ and seconds:
        mbps = f"{int(bytes_) / float(seconds) / 1e6:.3f}"
    # One line per row: a newline inside a quoted CSV field is legal, but it
    # breaks every `grep` and `wc -l` a person will run on this file.
    err = " ".join(str(error or "").split())[:500]
    return {
        "utc": utc or utc_now(),
        "run": run,
        "taxid": "" if taxid in (None, "") else str(int(taxid)),
        "tier": tier,
        "route": route,
        "node_type": node,
        "bytes": "" if bytes_ is None else str(int(bytes_)),
        "seconds": "" if seconds is None else f"{float(seconds):.2f}",
        "mb_per_s": mbps,
        "verified": ver,
        "outcome": outcome,
        "error": err,
    }


def append(row: dict, path: str = DOWNLOADS_CSV) -> None:
    """Append one validated row. Writes the header if the file is new or empty."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o664)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        buf = io.StringIO()
        w = csv.DictWriter(buf, fieldnames=FIELDS, lineterminator="\n")
        if os.fstat(fd).st_size == 0:
            w.writeheader()
        w.writerow(row)
        os.write(fd, buf.getvalue().encode())
        # Others must be able to append too (AGENTS.md section 1).
        mode = os.fstat(fd).st_mode
        if not mode & stat.S_IWGRP:
            try:
                os.fchmod(fd, mode | stat.S_IWGRP)
            except PermissionError:
                pass
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


def record(path: str = DOWNLOADS_CSV, **kw) -> dict:
    """make_row + append. Returns the row written."""
    row = make_row(**kw)
    append(row, path)
    return row


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m popseq.downloads",
                                description="Append one fetch attempt to downloads.csv")
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("append")
    a.add_argument("--file", default=DOWNLOADS_CSV)
    a.add_argument("--run", required=True)
    a.add_argument("--taxid")
    a.add_argument("--tier", default="full", choices=TIERS)
    a.add_argument("--route", required=True, choices=ROUTES)
    a.add_argument("--node-type", choices=NODE_TYPES)
    a.add_argument("--bytes", type=int)
    a.add_argument("--seconds", type=float)
    a.add_argument("--verified", default="")
    a.add_argument("--outcome", required=True, choices=OUTCOMES)
    a.add_argument("--error", default="")
    args = p.parse_args(argv)
    try:
        row = record(args.file, run=args.run, taxid=args.taxid, tier=args.tier,
                     route=args.route, node=args.node_type, bytes_=args.bytes,
                     seconds=args.seconds, verified=args.verified,
                     outcome=args.outcome, error=args.error)
    except ValueError as exc:
        print(f"popseq.downloads: {exc}", file=sys.stderr)
        return 2
    print(f"recorded {row['run']} {row['route']} {row['outcome']} "
          f"on {socket.gethostname()}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
