"""`popseq fetch`: one run (or a manifest), one depth tier, verified, with a sidecar (issue #13).

The rules are spec 002; read it first. In short:

* Pairing is what the serving archive holds, not the `library_layout` label. NCBI serves
  pairs when `spots_with_mates` > 0 and equals `spots`; ENA when its file list has `_1` and `_2`.
* Verify against the archive that served the bytes (SRA: NCBI `spots_with_mates`, stored at
  plan time; ENA: MD5 per file and equal mates). Never ENA `read_count`.
* Both mates of a run come from one source.
* Temporary files live in `${SCRATCH:?}` and are removed whatever happens. Only declared
  outputs are kept. A cached run without a `<RUN>.json` sidecar is not trusted.
* Every attempt, failures included, is one row in `benchmarks/downloads.csv`.

The network and the archive tools sit behind `Backend` and `Lookups`, so the logic is tested
with fakes (tests/test_fetch.py).
"""

from __future__ import annotations

import csv
import json
import os
import shlex
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from typing import Optional

from popseq import downloads as dl
from popseq import lookup, routes, subsample
from popseq.query import ROOT, RUNS_TABLE, ena_layout, fastq_mates

TIERS = dl.TIERS
EXIT = {
    "ok": 0,
    "usage": 2,
    "unavailable_at_ncbi": 10,
    "unpaired": 11,
    "verify_failed": 12,
    "transient": 13,
    "hook_failed": 14,
}
GENOME_SIZES = os.path.join(ROOT, "catalog", "genome_sizes.csv")
CACHE_ROOT = os.path.join(ROOT, "cache")
DEFAULT_SPOTS = {"head": 500_000, "windows": 2_000_000}
DEFAULT_WINDOWS = 8
# The number of fetches a caller runs at once, recorded in every downloads.csv row. Measured
# 2026-11-11 (issue #27, spec 002 "Concurrency sweep results"): no failures at up to 40; speed per
# fetch is set by how many fetches share a node.
DEFAULT_CONCURRENCY = 20
FETCH_FIELDS = [
    "run",
    "tier",
    "route",
    "layout",
    "bases_kept",
    "reads_kept",
    "coverage_target",
    "coverage_reached",
    "seconds",
    "verified",
    "outcome",
]


class Usage(Exception):
    """Bad arguments or an unknown run. Exit code 2."""


# ----------------------------------------------------------------------------- the plan


@dataclass
class Plan:
    run: str
    taxid: Optional[int] = None
    species: str = ""
    sample: str = ""
    fastq_ftp: str = ""
    fastq_aspera: str = ""
    base_count: Optional[int] = None
    genome_mb: Optional[float] = None
    ena_md5: dict = field(default_factory=dict)  # file name -> md5
    ena_bytes: dict = field(default_factory=dict)  # file name -> bytes
    ncbi_status: str = "unknown"  # ok | absent | error | unknown
    ncbi_spots: Optional[int] = None
    ncbi_mates: Optional[int] = None

    @property
    def ena_layout(self) -> str:
        return ena_layout(self.fastq_ftp)

    @property
    def ncbi_layout(self) -> Optional[str]:
        """PAIRED, SINGLE, MIXED (some spots paired, some not) or None when NCBI has no answer."""
        if self.ncbi_status != "ok" or self.ncbi_spots is None or self.ncbi_mates is None:
            return None
        if self.ncbi_mates == 0:
            return "SINGLE"
        return "PAIRED" if self.ncbi_mates == self.ncbi_spots else "MIXED"

    def as_route_plan(self) -> dict:
        """The dict `routes.transfer` expects."""
        return {
            "run": self.run,
            "fastq_ftp": self.fastq_ftp,
            "fastq_aspera": self.fastq_aspera,
            "ncbi_spots": self.ncbi_spots,
            "ncbi_spots_with_mates": self.ncbi_mates,
        }


class Lookups:
    """Plan-time lookups. Replaced by a fake in tests."""

    def ena(self, run: str) -> dict:
        return lookup.ena_lookup([run]).get(run, {})

    def ncbi(self, run: str):
        return lookup.ncbi_counts(run)


def catalog_rows(runs) -> dict:
    """The catalog rows of `runs`, in one pass over the run table."""
    want, found = set(runs), {}
    if not want or not os.path.exists(RUNS_TABLE):
        return found
    with lookup.open_table(RUNS_TABLE) as fh:
        for row in csv.DictReader(fh):
            if row["run_accession"] in want:
                found[row["run_accession"]] = row
                if len(found) == len(want):
                    break
    return found


def genome_size(species_taxid) -> Optional[float]:
    if not species_taxid or not os.path.exists(GENOME_SIZES):
        return None
    with open(GENOME_SIZES, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["species_taxid"] == str(species_taxid):
                return float(row["genome_mb"])
    return None


def _int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def build_plan(
    run: str, crow: Optional[dict], *, lookups: Lookups, plan_row: Optional[dict] = None
) -> Plan:
    """Catalog row + plan-time values.

    A `plan_row` (benchmark-7.1 plan format) replaces the live lookups.

    Raises Usage for an unknown run. A lookup that fails leaves `ncbi_status` as `error`; the
    SRA routes then do not run (they cannot be verified) and the outcome is `transient`,
    never `unavailable_at_ncbi`.
    """
    p = Plan(run=run)
    if crow:
        p.taxid = _int(crow.get("species_taxid") or crow.get("tax_id"))
        p.species = crow.get("species_name") or crow.get("scientific_name") or ""
        p.sample = crow.get("sample_accession", "")
        p.fastq_ftp = crow.get("fastq_ftp", "")
        p.fastq_aspera = crow.get("fastq_aspera", "")
        p.base_count = _int(crow.get("base_count"))
    if plan_row:
        p.taxid = _int(plan_row.get("taxid")) or p.taxid
        p.fastq_ftp = plan_row.get("fastq_ftp") or p.fastq_ftp
        p.fastq_aspera = plan_row.get("fastq_aspera") or p.fastq_aspera
        p.base_count = _int(plan_row.get("base_count")) or p.base_count
        ena = {
            "fastq_md5": plan_row.get("fastq_md5", ""),
            "fastq_bytes": plan_row.get("fastq_bytes", ""),
        }
        spots, mates = _int(plan_row.get("ncbi_spots")), _int(plan_row.get("ncbi_spots_with_mates"))
        p.ncbi_status = "ok" if spots is not None and mates is not None else "unknown"
        p.ncbi_spots, p.ncbi_mates = spots, mates
    else:
        ena = lookups.ena(run) or {}
        p.ncbi_status, p.ncbi_spots, p.ncbi_mates = lookups.ncbi(run)
        p.fastq_ftp = ena.get("fastq_ftp") or p.fastq_ftp
        p.fastq_aspera = ena.get("fastq_aspera") or p.fastq_aspera
        p.base_count = _int(ena.get("base_count")) or p.base_count
    urls = [u for u in p.fastq_ftp.split(";") if u]
    md5s = [m for m in str(ena.get("fastq_md5", "")).split(";") if m]
    sizes = [b for b in str(ena.get("fastq_bytes", "")).split(";") if b]
    for i, u in enumerate(urls):
        name = os.path.basename(u)
        if i < len(md5s):
            p.ena_md5[name] = md5s[i]
        if i < len(sizes):
            p.ena_bytes[name] = _int(sizes[i])
    if not crow and not urls and p.ncbi_status != "ok":
        raise Usage(f"{run}: not in the catalog and no archive answered")
    p.genome_mb = genome_size(p.taxid)
    return p


# --------------------------------------------------------------------------- the backend


class Backend:
    """Everything that touches the network or an archive tool. Tests replace these methods."""

    def transfer(self, plan: dict, route: str, wd: str, timeout, aspera_key=None) -> dict:
        return routes.transfer(plan, route, wd, timeout, aspera_key)

    def convert_sra(self, sra: str, wd: str, threads: int):
        scratch = os.environ.get("SCRATCH") or wd
        cmd = [
            "parallel-fastq-dump",
            "--sra-id",
            sra,
            "--threads",
            str(threads),
            "--outdir",
            wd,
            "--split-files",
            "--gzip",
            "--tmpdir",
            scratch,
        ]
        rc, sec, _out, err = routes._run(cmd)
        return rc, sec, " ".join(err.split())[-400:]

    def window(self, sra: str, start: int, count: int, outdir: str, tmpdir: str):
        return routes.sra_window(sra, start, count, outdir, tmpdir)

    def ena_head(self, url: str, dest: str, n_reads: int, timeout=None):
        return routes.ena_head(url, dest, n_reads, timeout)

    def read_stats(self, path: str):
        return routes.read_stats(path)

    def md5(self, path: str) -> str:
        return routes._md5(path)

    def names(self, path: str):
        """(first read name, last read name) of a gzip FASTQ."""
        return first_last_names(path)

    def tool_versions(self) -> dict:
        return routes.tool_versions()


def first_last_names(path: str):
    p = subprocess.run(
        routes._decompress_cmd(path), stdin=subprocess.DEVNULL, capture_output=True, check=False
    )
    lines = p.stdout.split(b"\n")
    while lines and not lines[-1]:
        lines.pop()
    if p.returncode != 0 or len(lines) < 4:
        raise RuntimeError(f"cannot read {path}")
    return _name(lines[0]), _name(lines[-4])


def _name(header: bytes) -> str:
    """The read name without the `/1` `/2` suffix or the description after the first space."""
    n = header.decode("ascii", "replace").lstrip("@").split(" ")[0]
    return n[:-2] if n.endswith(("/1", "/2")) else n


# ------------------------------------------------------------------------------- results


@dataclass
class Attempt:
    route: str
    outcome: str  # ok | unavailable_at_ncbi | unpaired | verify_failed | transfer_failed
    error: str = ""
    layout: str = ""
    files: list = field(default_factory=list)  # work-dir paths: mate 1 [, mate 2], or the one file
    bytes: int = 0
    seconds_fetch: float = 0.0
    seconds_convert: float = 0.0
    reads_kept: int = 0
    bases_kept: int = 0
    expected: Optional[int] = None
    m1: Optional[int] = None
    m2: Optional[int] = None
    verified: Optional[bool] = None
    scope: str = ""  # full | subsample
    extra: dict = field(default_factory=dict)


@dataclass
class Result:
    run: str
    tier: str
    outcome: str
    route: str = ""
    layout: str = ""
    bases_kept: Optional[int] = None
    reads_kept: Optional[int] = None
    coverage_target: Optional[float] = None
    coverage_reached: Optional[float] = None
    seconds: float = 0.0
    verified: bool = False
    exit_code: int = 0
    message: str = ""
    sidecar: str = ""
    outputs: list = field(default_factory=list)
    resumed: bool = False

    def row(self) -> list:
        def f(x):
            return "" if x is None else (f"{x:.3f}" if isinstance(x, float) else str(x))

        return [
            self.run,
            self.tier,
            self.route,
            self.layout,
            f(self.bases_kept),
            f(self.reads_kept),
            f(self.coverage_target),
            f(self.coverage_reached),
            f"{self.seconds:.1f}",
            "true" if self.verified else "false",
            self.outcome,
        ]


# ------------------------------------------------------------------------------ the fetch


@dataclass
class Options:
    tier: str = "windows"
    out: Optional[str] = None
    coverage: Optional[float] = None
    genome_mb: Optional[float] = None
    windows: int = DEFAULT_WINDOWS
    spots: Optional[int] = None
    require_pairs: bool = False
    exec_cmd: Optional[str] = None
    keep: bool = False
    threads: int = 4
    concurrency: int = DEFAULT_CONCURRENCY
    timeout: Optional[int] = None
    cache_root: str = CACHE_ROOT
    downloads_csv: str = dl.DOWNLOADS_CSV


def _concat(parts, dest):
    """Concatenate gzip files: concatenated gzip members are one valid gzip file."""
    with open(dest, "wb") as out:
        for p in parts:
            with open(p, "rb") as fh:
                shutil.copyfileobj(fh, out, 1 << 22)


def _classify(files):
    """Mate files by name: (mate1, mate2) for `_1`/`_2`, (single, None) for one bare file."""
    m1 = m2 = other = None
    for f in sorted(files):
        b = os.path.basename(f)
        if b.endswith(("_1.fastq.gz", "_1.fq.gz")):
            m1 = f
        elif b.endswith(("_2.fastq.gz", "_2.fq.gz")):
            m2 = f
        else:
            other = f
    if m1:
        return m1, m2
    return (other, None) if other else (None, None)


class Fetcher:
    def __init__(self, backend: Optional[Backend] = None, lookups: Optional[Lookups] = None):
        self.backend = backend or Backend()
        self.lookups = lookups or Lookups()

    # -- public entry
    def fetch(
        self, run: str, opt: Options, *, plan: Optional[Plan] = None, plan_row=None, crow=None
    ) -> Result:
        t0 = time.monotonic()
        if opt.tier not in TIERS:
            raise Usage(f"tier must be one of {TIERS}")
        scratch = os.environ.get("SCRATCH")
        if not scratch:
            raise Usage(
                "SCRATCH is not set: temporary files must go to node-local $SCRATCH, not /tmp "
                "(run inside a job, or set SCRATCH)"
            )
        if opt.tier != "full" and not opt.out:
            raise Usage("--out is required for every tier except full (full goes to the cache)")
        if opt.tier == "coverage" and not opt.coverage:
            raise Usage("--coverage is required for tier coverage")
        if plan is None:
            if crow is None:
                crow = catalog_rows([run]).get(run)
            plan = build_plan(run, crow, lookups=self.lookups, plan_row=plan_row)
        if opt.tier == "coverage" and not (opt.genome_mb or plan.genome_mb):
            raise Usage(f"{run}: no genome size for species taxid {plan.taxid}; give --genome-mb")

        dest = self._dest(plan, opt)
        sidecar = os.path.join(dest, f"{run}.json")
        done = self._resumable(sidecar, plan, opt)
        if done is not None:
            return done

        wd = tempfile.mkdtemp(prefix=f"popseq_fetch_{run}_", dir=scratch)
        try:
            attempts = []
            for route in self._route_order(plan, opt):
                try:
                    a = self._attempt(plan, opt, route, wd)
                except OSError as exc:  # a missing tool (prefetch, fastq-dump) is a failed attempt
                    a = Attempt(
                        route, "transfer_failed", error=f"{type(exc).__name__}: {exc}"[:400]
                    )
                attempts.append(a)
                self._log(plan, opt, a)
                if a.outcome == "ok":
                    return self._finish(plan, opt, a, dest, sidecar, t0)
                if a.outcome == "verify_failed":
                    break  # a mismatch fails the task; the caller retries
                self._clean(wd)  # delete partials before the next route starts
            return self._failure(plan, opt, attempts, t0)
        finally:
            shutil.rmtree(wd, ignore_errors=True)

    # -- helpers
    @staticmethod
    def _clean(wd):
        for n in os.listdir(wd):
            p = os.path.join(wd, n)
            if os.path.isdir(p):
                shutil.rmtree(p, ignore_errors=True)
            else:
                os.remove(p)

    @staticmethod
    def _dest(plan, opt):
        if opt.tier == "full":
            return os.path.join(opt.cache_root, str(plan.taxid or "unknown"), plan.run)
        return opt.out

    @staticmethod
    def _params(plan, opt):
        cov = opt.tier == "coverage"
        return {
            "tier": opt.tier,
            "coverage": opt.coverage if cov else None,
            "windows": opt.windows if opt.tier in ("windows", "coverage") else None,
            "spots": opt.spots,
            "require_pairs": opt.require_pairs,
            "genome_mb": (opt.genome_mb or plan.genome_mb) if cov else None,
        }

    def _resumable(self, sidecar, plan, opt):
        """A verified sidecar for the same tier and parameters, with its files: nothing to do."""
        if not os.path.exists(sidecar):
            return None
        try:
            with open(sidecar) as fh:
                s = json.load(fh)
        except (OSError, ValueError):
            return None
        if not s.get("verified") or s.get("params") != self._params(plan, opt):
            return None
        d = os.path.dirname(sidecar)
        if s.get("reads_on_disk", True) and not all(
            os.path.exists(os.path.join(d, f)) for f in s.get("files", [])
        ):
            return None
        r = Result(
            plan.run,
            opt.tier,
            "ok",
            route=s.get("route", ""),
            layout=s.get("layout", ""),
            bases_kept=s.get("bases_kept"),
            reads_kept=s.get("reads_kept"),
            coverage_target=s.get("coverage_target"),
            coverage_reached=s.get("coverage_reached"),
            verified=True,
            message="already fetched and verified",
            sidecar=sidecar,
            resumed=True,
        )
        r.outputs = [os.path.join(d, f) for f in s.get("files", [])]
        return r

    @staticmethod
    def _route_order(plan, opt):
        ena = ["ena_https"] if plan.fastq_ftp else []
        sra = ["aws_odp", "ncbi_prefetch"]
        return ena + sra if opt.tier == "head" else sra + ena  # only ENA can stop early

    def _log(self, plan, opt, a: Attempt):
        dl.record(
            opt.downloads_csv,
            run=plan.run,
            taxid=plan.taxid,
            tier=opt.tier,
            route=a.route,
            bytes_=a.bytes or None,
            seconds=(a.seconds_fetch + a.seconds_convert) or None,
            verified=a.verified,
            outcome=a.outcome,
            error=a.error,
            concurrency=opt.concurrency,
        )

    def _failure(self, plan, opt, attempts, t0) -> Result:
        outs = [a.outcome for a in attempts]
        if "verify_failed" in outs:
            code, msg = (
                "verify_failed",
                next(a.error for a in attempts if a.outcome == "verify_failed"),
            )
        elif opt.require_pairs and "unpaired" in outs:
            code, msg = "unpaired", next(a.error for a in attempts if a.outcome == "unpaired")
        elif "transfer_failed" in outs:
            code = "transient"
            msg = "; ".join(
                f"{a.route}: {a.error}" for a in attempts if a.outcome == "transfer_failed"
            )
        else:
            code = "unavailable_at_ncbi"
            msg = "; ".join(f"{a.route}: {a.error}" for a in attempts) or "no route has the run"
        return Result(
            plan.run,
            opt.tier,
            code,
            route=attempts[-1].route if attempts else "",
            seconds=time.monotonic() - t0,
            exit_code=EXIT[code],
            message=msg[:500],
        )

    # -- one route
    def _attempt(self, plan: Plan, opt: Options, route: str, wd: str) -> Attempt:
        wd_route = os.path.join(wd, route)
        os.makedirs(wd_route, exist_ok=True)
        if route == "ena_https":
            return self._attempt_ena(plan, opt, wd_route)
        return self._attempt_sra(plan, opt, route, wd_route)

    def _attempt_sra(self, plan, opt, route, wd) -> Attempt:
        be, a = self.backend, Attempt(route=route, outcome="ok")
        if plan.ncbi_status == "absent":
            return Attempt(route, "unavailable_at_ncbi", error="NCBI has no such run (ENA only)")
        if plan.ncbi_status != "ok":
            return Attempt(
                route,
                "transfer_failed",
                error="no NCBI spot counts (lookup failed): an SRA download cannot be verified",
            )
        lay = plan.ncbi_layout
        if lay == "MIXED":
            return Attempt(
                route,
                "unavailable_at_ncbi",
                error=f"NCBI holds mixed spots ({plan.ncbi_spots}; {plan.ncbi_mates} with mates)",
            )
        if lay == "SINGLE" and opt.require_pairs:
            return Attempt(
                route,
                "unavailable_at_ncbi",
                error=f"not paired at NCBI ({plan.ncbi_spots} spots, none with mates): skipped",
            )
        res = be.transfer(plan.as_route_plan(), route, wd, opt.timeout)
        a.bytes, a.seconds_fetch = res["bytes"], res["seconds"]
        if res["rc"] != 0:
            a.outcome = "unavailable_at_ncbi" if res["unavailable"] else "transfer_failed"
            a.error = res["error"] or f"exit {res['rc']}"
            return a
        sra = res["files"][0]
        a.layout, a.scope = lay, ("full" if opt.tier in ("full", "stride") else "subsample")
        try:
            if opt.tier in ("full", "stride"):
                return self._sra_full(plan, opt, a, sra, wd)
            return self._sra_windows(plan, opt, a, sra, wd)
        except (RuntimeError, OSError, ValueError) as exc:
            a.outcome, a.error, a.verified = "verify_failed", str(exc)[:400], False
            return a

    def _sra_full(self, plan, opt, a, sra, wd):
        be = self.backend
        rc, sec, err = be.convert_sra(sra, wd, opt.threads)
        a.seconds_convert = sec
        if rc != 0:
            a.outcome, a.error = "transfer_failed", f"parallel-fastq-dump exit {rc}: {err}"
            return a
        m1, m2 = _classify([os.path.join(wd, f) for f in os.listdir(wd) if f.endswith(".fastq.gz")])
        if not m1:
            a.outcome, a.error, a.verified = "verify_failed", "conversion produced no FASTQ", False
            return a
        expected = plan.ncbi_mates if a.layout == "PAIRED" else plan.ncbi_spots
        a.expected = expected
        n1, b1 = be.read_stats(m1)
        n2, b2 = be.read_stats(m2) if m2 else (None, 0)
        a.m1, a.m2 = n1, n2
        ok = (
            (m2 is not None and n1 == expected and n2 == expected)
            if a.layout == "PAIRED"
            else (m2 is None and n1 == expected)
        )
        a.verified = ok
        if not ok:
            a.outcome = "verify_failed"
            a.error = f"expected {expected} (NCBI), mate1 {n1}, mate2 {n2}, layout {a.layout}"
            return a
        a.files, a.reads_kept, a.bases_kept = ([m1, m2] if m2 else [m1]), n1, b1 + b2
        if opt.tier == "stride":
            step = subsample.stride(a.reads_kept, opt.spots or DEFAULT_SPOTS["windows"])
            outs = [os.path.join(wd, "stride_" + os.path.basename(f)) for f in a.files]
            kept, _read, _b = subsample.select_spots(a.files, outs, subsample.Every(step))
            return self._after_cut(a, outs, kept, f"stride {step}")
        return a

    def _after_cut(self, a, outs, kept, note):
        """Check what a cutter wrote: counts equal the request, and the mates' names pair up."""
        be = self.backend
        stats = [be.read_stats(f) for f in outs]
        reads = [s[0] for s in stats]
        if any(r != kept for r in reads):
            a.outcome, a.verified, a.error = (
                "verify_failed",
                False,
                f"wrote {reads} reads, expected {kept} per mate",
            )
            return a
        if len(outs) == 2 and be.names(outs[0]) != be.names(outs[1]):
            a.outcome, a.verified = "verify_failed", False
            a.error = "read names of the mates disagree at the first or last read"
            return a
        a.files, a.reads_kept, a.bases_kept = outs, kept, sum(s[1] for s in stats)
        a.extra["cut"] = note
        return a

    def _window_spec(self, plan, opt, total, bases_per_spot):
        """`([(start, count), ...], coverage_target)` for the subsampled tiers."""
        if opt.tier == "head":
            return [(0, min(opt.spots or DEFAULT_SPOTS["head"], total))], None
        if opt.tier == "windows":
            return subsample.windows(
                total, opt.spots or DEFAULT_SPOTS["windows"], opt.windows
            ), None
        n = subsample.spots_for_coverage(
            opt.genome_mb or plan.genome_mb, opt.coverage, bases_per_spot
        )
        return subsample.windows(total, n, opt.windows), opt.coverage

    def _sra_windows(self, plan, opt, a, sra, wd):
        be = self.backend
        total = plan.ncbi_spots
        if not total:
            a.outcome, a.error, a.verified = (
                "verify_failed",
                "no NCBI spot count stored for this run",
                False,
            )
            return a
        bps = None
        if opt.tier == "coverage":
            if plan.base_count:
                bps = plan.base_count / total
            else:  # no base count in the plan: measure 1,000 spots
                probe = os.path.join(wd, "probe")
                os.makedirs(probe, exist_ok=True)
                rc, _s, err = be.window(sra, 0, min(1000, total), probe, wd)
                if rc != 0:
                    a.outcome, a.error = "transfer_failed", f"fastq-dump probe exit {rc}: {err}"
                    return a
                st = [
                    be.read_stats(os.path.join(probe, f))
                    for f in os.listdir(probe)
                    if f.endswith(".fastq.gz")
                ]
                bps = sum(s[1] for s in st) / max(1, st[0][0])
        spans, target = self._window_spec(plan, opt, total, bps)
        a.extra["windows"], a.extra["coverage_target"] = spans, target
        a.expected = sum(c for _s, c in spans)
        parts = {"m1": [], "m2": []}
        conv = 0.0
        for i, (start, count) in enumerate(spans):
            wdir = os.path.join(wd, f"win{i}")
            os.makedirs(wdir, exist_ok=True)
            rc, sec, err = be.window(sra, start, count, wdir, wd)
            conv += sec
            if rc != 0:
                a.outcome, a.error = "transfer_failed", f"fastq-dump window {i} exit {rc}: {err}"
                return a
            m1, m2 = _classify(
                [os.path.join(wdir, f) for f in os.listdir(wdir) if f.endswith(".fastq.gz")]
            )
            bad = None
            if not m1:
                bad = f"window {i} produced no FASTQ"
            elif (a.layout == "PAIRED") != (m2 is not None):
                got = "paired" if m2 else "single"
                bad = f"window {i}: {a.layout} at NCBI but the files are {got}"
            else:
                n1 = be.read_stats(m1)[0]
                n2 = be.read_stats(m2)[0] if m2 else n1
                if n1 != count or n2 != count:
                    bad = f"window {i}: expected {count} reads per mate, got {n1} and {n2}"
                elif m2 and be.names(m1) != be.names(m2):
                    bad = f"window {i}: read names of the mates disagree"
            if bad:
                a.outcome, a.error, a.verified = "verify_failed", bad, False
                return a
            parts["m1"].append(m1)
            if m2:
                parts["m2"].append(m2)
        a.seconds_convert = conv
        outs = []
        for suffix, key in zip(["_1", "_2"] if a.layout == "PAIRED" else [""], ("m1", "m2")):
            dest = os.path.join(wd, f"{plan.run}{suffix}.fastq.gz")
            _concat(parts[key], dest)
            outs.append(dest)
        stats = [be.read_stats(f) for f in outs]
        a.files, a.reads_kept, a.bases_kept = outs, stats[0][0], sum(s[1] for s in stats)
        a.m1, a.m2 = stats[0][0], (stats[1][0] if len(stats) > 1 else None)
        a.verified = a.reads_kept == a.expected and (a.m2 is None or a.m1 == a.m2)
        if not a.verified:
            a.outcome, a.error = (
                "verify_failed",
                f"concatenated {a.m1} reads, expected {a.expected}",
            )
        return a

    def _attempt_ena(self, plan, opt, wd) -> Attempt:
        be, a = self.backend, Attempt(route="ena_https", outcome="ok")
        m1u, m2u, _others = fastq_mates(plan.fastq_ftp)
        if not m1u:
            return Attempt(
                "ena_https", "unavailable_at_ncbi", error="ENA lists no FASTQ for this run"
            )
        a.layout = "PAIRED" if m2u else "SINGLE"
        if opt.require_pairs and a.layout != "PAIRED":
            return Attempt(
                "ena_https",
                "unpaired",
                layout=a.layout,
                error="pairs required, and ENA serves one file (and NCBI does not serve pairs)",
            )
        if opt.tier == "head":
            return self._ena_head(plan, opt, a, wd, [u for u in (m1u, m2u) if u])
        sub = dict(plan.as_route_plan())
        sub["fastq_ftp"] = ";".join(u for u in (m1u, m2u) if u)
        res = (
            be.transfer(sub, "ena_https", wd, opt.timeout)
            if m2u
            else self._transfer_single(sub, wd, opt)
        )
        a.bytes, a.seconds_fetch, a.scope = res["bytes"], res["seconds"], "full"
        if res["rc"] != 0:
            a.outcome, a.error = "transfer_failed", res["error"] or f"exit {res['rc']}"
            return a
        try:
            files = res["files"]
            for f in files:
                want = plan.ena_md5.get(os.path.basename(f))
                if not want or be.md5(f) != want:
                    a.outcome, a.verified = "verify_failed", False
                    a.error = (
                        f"md5 of {os.path.basename(f)} does not match ENA ({want or 'none stored'})"
                    )
                    return a
            stats = [be.read_stats(f) for f in files]
            a.m1, a.m2 = stats[0][0], (stats[1][0] if len(stats) > 1 else None)
            if a.m2 is not None and a.m1 != a.m2:
                a.outcome, a.verified, a.error = (
                    "verify_failed",
                    False,
                    f"mates differ: {a.m1} and {a.m2}",
                )
                return a
            a.verified, a.files = True, files
            a.reads_kept, a.bases_kept = a.m1, sum(s[1] for s in stats)
            return a if opt.tier == "full" else self._ena_cut(plan, opt, a, wd)
        except (RuntimeError, OSError, ValueError) as exc:
            a.outcome, a.error, a.verified = "verify_failed", str(exc)[:400], False
            return a

    def _ena_head(self, plan, opt, a, wd, urls):
        """Tier `head` from ENA: stream each mate, stop after N reads, no whole-file download."""
        be, n = self.backend, opt.spots or DEFAULT_SPOTS["head"]
        outs, reads = [], []
        for i, url in enumerate(urls):
            dest = os.path.join(wd, f"head_{plan.run}_{i + 1}.fastq.gz")
            rc, sec, got, err = be.ena_head("https://" + url, dest, n, opt.timeout)
            a.seconds_fetch += sec
            if rc != 0:
                a.outcome, a.error = "transfer_failed", err or f"exit {rc}"
                return a
            outs.append(dest)
            reads.append(got)
        a.bytes, a.scope = sum(os.path.getsize(f) for f in outs), "subsample"
        a.m1, a.m2 = reads[0], (reads[1] if len(reads) > 1 else None)
        a.expected = n
        # Equal mates, and N reads each unless the run is shorter than N.
        if len(set(reads)) != 1 or reads[0] > n:
            a.outcome, a.verified, a.error = (
                "verify_failed",
                False,
                f"mates hold {reads} reads, wanted {n}",
            )
            return a
        a.verified = True
        a.expected = reads[0]
        return self._after_cut(a, outs, reads[0], f"head {n}")

    @staticmethod
    def _transfer_single(sub, wd, opt):
        """routes.transfer keeps only `_1`/`_2` files; a single-end run is one bare file."""
        url = sub["fastq_ftp"]
        dest = os.path.join(wd, os.path.basename(url))
        rc, sec, out, e = routes._run(
            ["curl", "-sS", "-f", "-L", "--retry", "3", "-o", dest, "https://" + url], opt.timeout
        )
        files = [dest] if rc == 0 else []
        return {
            "rc": rc,
            "seconds": sec,
            "error": " ".join((out + e).split())[-500:] if rc else "",
            "files": files,
            "bytes": sum(os.path.getsize(f) for f in files if os.path.exists(f)),
        }

    def _ena_cut(self, plan, opt, a, wd):
        total = a.reads_kept
        if opt.tier == "stride":
            step = subsample.stride(total, opt.spots or DEFAULT_SPOTS["windows"])
            keep, note = subsample.Every(step), f"stride {step}"
        else:
            spans, target = self._window_spec(plan, opt, total, a.bases_kept / max(1, total))
            a.extra["windows"], a.extra["coverage_target"] = spans, target
            keep, note = subsample.Window(spans), "windows"
        names = ["_1", "_2"] if len(a.files) == 2 else [""]
        outs = [os.path.join(wd, f"cut_{plan.run}{s}.fastq.gz") for s in names]
        kept, _read, _b = subsample.select_spots(a.files, outs, keep)
        return self._after_cut(a, outs, kept, note)

    # -- success: place outputs, write the sidecar last, run the hook
    def _finish(self, plan, opt, a: Attempt, dest, sidecar, t0) -> Result:
        keep_reads = opt.tier == "full" or opt.keep or not opt.exec_cmd
        names = [
            f"{plan.run}.fastq.gz" if len(a.files) == 1 else f"{plan.run}_{i + 1}.fastq.gz"
            for i in range(len(a.files))
        ]
        genome = opt.genome_mb or plan.genome_mb
        reached = (a.bases_kept / (genome * 1e6)) if (genome and opt.tier == "coverage") else None
        target = a.extra.get("coverage_target") or (
            opt.coverage if opt.tier == "coverage" else None
        )
        side = {
            "run": plan.run,
            "taxid": plan.taxid,
            "biosample": plan.sample,
            "route": a.route,
            "tier": opt.tier,
            "params": self._params(plan, opt),
            "fetched_utc": dl.utc_now(),
            "node_type": dl.node_type(),
            "concurrency": opt.concurrency,
            "bytes_total": a.bytes,
            "seconds_fetch": round(a.seconds_fetch, 2),
            "seconds_convert": round(a.seconds_convert, 2),
            "mb_per_s": round(a.bytes / a.seconds_fetch / 1e6, 3)
            if a.bytes and a.seconds_fetch
            else None,
            "layout": a.layout,
            "ncbi_spots": plan.ncbi_spots,
            "ncbi_spots_with_mates": plan.ncbi_mates,
            "read_count_expected": a.expected,
            "read_count_mate1": a.m1,
            "read_count_mate2": a.m2,
            "ena_md5": plan.ena_md5 if a.route == "ena_https" else None,
            "reads_kept": a.reads_kept,
            "bases_kept": a.bases_kept,
            "genome_mb": genome,
            "coverage_target": target,
            "coverage_reached": reached,
            "windows": a.extra.get("windows"),
            "verified": True,
            "verify_scope": a.scope,
            "reads_on_disk": keep_reads,
            "files": names if keep_reads else [],
            "tools": self.backend.tool_versions(),
        }
        os.makedirs(dest, exist_ok=True)
        hook_rc = self._hook(opt.exec_cmd, plan.run, a.files) if opt.exec_cmd else None
        if keep_reads:
            for src, name in zip(a.files, names):
                tmp = os.path.join(dest, f".{name}.part")
                shutil.move(src, tmp)
                os.replace(tmp, os.path.join(dest, name))
        tmp = os.path.join(dest, f".{plan.run}.json.part")
        with open(tmp, "w") as fh:
            json.dump(side, fh, indent=2, default=str)
            fh.write("\n")
        os.replace(tmp, sidecar)  # the sidecar goes last: it marks a complete run
        r = Result(
            plan.run,
            opt.tier,
            "ok",
            route=a.route,
            layout=a.layout,
            bases_kept=a.bases_kept,
            reads_kept=a.reads_kept,
            coverage_target=target,
            coverage_reached=reached,
            seconds=time.monotonic() - t0,
            verified=True,
            sidecar=sidecar,
            outputs=[os.path.join(dest, n) for n in names] if keep_reads else [],
        )
        if hook_rc not in (None, 0):
            r.outcome, r.exit_code = "hook_failed", EXIT["hook_failed"]
            r.message = f"--exec returned {hook_rc}; the download is ok and verified"
        return r

    @staticmethod
    def _hook(cmd, run, files):
        """`CMD RUN FASTQ_1 [FASTQ_2]`, stdin closed, after verification and before the reads go."""
        p = subprocess.run(
            shlex.split(cmd) + [run] + list(files), stdin=subprocess.DEVNULL, check=False
        )
        return p.returncode


# --------------------------------------------------------------------------- manifests


def read_manifest(path: str) -> list[str]:
    """Run accessions from a CSV or TSV (plain, .gz or .zst) with a `run_accession` column."""
    with lookup.open_table(path) as fh:
        head = fh.readline()
    delim = "\t" if "\t" in head else ","
    with lookup.open_table(path) as fh:
        rows = csv.DictReader(fh, delimiter=delim)
        if "run_accession" not in (rows.fieldnames or []):
            raise Usage(f"{path}: no run_accession column")
        return [r["run_accession"].strip() for r in rows if r["run_accession"].strip()]
