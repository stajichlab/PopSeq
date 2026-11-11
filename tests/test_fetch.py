"""popseq fetch (issue #13), tested against a fake backend: no network, no archive tools.

The fake writes small gzip FASTQ files, so counting, MD5, pairing names, concatenation and
the subsample cutters all run for real. What is faked is only the transfer and the SRA
conversion (`transfer`, `convert_sra`, `window`, `ena_head`).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os

import pytest

from popseq import fetch as fx
from popseq import subsample as ss
from popseq.fetch import EXIT, Backend, Fetcher, Options, Plan, Usage


# --------------------------------------------------------------------------- helpers
def write_fastq(path, n, start=0, read_len=50):
    """`n` reads named by spot index `start`..`start+n-1`, so mates and windows can be checked."""
    seq, qual = ("ACGT" * 40)[:read_len], "I" * read_len
    with gzip.open(path, "wt") as fh:
        for i in range(start, start + n):
            fh.write(f"@spot{i}/1\n{seq}\n+\n{qual}\n")


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def _res(route, wd, files, rc=0, unavailable=False, error=""):
    return {
        "route": route,
        "rc": rc,
        "error": error,
        "bytes": sum(os.path.getsize(f) for f in files if os.path.exists(f)),
        "seconds": 1.0,
        "files": files if rc == 0 else [],
        "file_kind": "sra" if route != "ena_https" else "fastq.gz",
        "unavailable": unavailable,
        "node_type": "login",
        "host": "h",
        "slurm_job": "",
        "utc_start": "t",
    }


class Fake(Backend):
    """Behaviours per route: 'ok', '404' (unavailable), 'fail' (transfer fails), 'missing' (no tool)."""

    def __init__(
        self, behave=None, sra_spots=1000, convert_reads=None, ena_reads=1000, ena_files=2
    ):
        self.behave = {"aws_odp": "ok", "ncbi_prefetch": "ok", "ena_https": "ok", **(behave or {})}
        self.sra_spots, self.convert_reads = sra_spots, convert_reads
        self.ena_reads, self.ena_files = ena_reads, ena_files
        self.calls = []
        self.paired_sra = True
        self.corrupt_md5 = False
        self.hook_calls = []

    def transfer(self, plan, route, wd, timeout, aspera_key=None):
        self.calls.append(("transfer", route))
        b = self.behave[route]
        if b == "missing":
            raise FileNotFoundError(2, "No such file or directory", "prefetch")
        if b == "404":
            return _res(route, wd, [], rc=22, unavailable=True, error="HTTP 404")
        if b == "fail":
            return _res(route, wd, [], rc=7, error="connection reset")
        run = plan["run"]
        if route == "ena_https":
            files = []
            names = [os.path.basename(u) for u in plan["fastq_ftp"].split(";")]
            for nm in names:
                p = os.path.join(wd, nm)
                write_fastq(p, self.ena_reads)
                files.append(p)
            return _res(route, wd, files)
        sra = os.path.join(wd, f"{run}.sra")
        open(sra, "wb").write(b"SRA" * 100)
        return _res(route, wd, [sra])

    def convert_sra(self, sra, wd, threads):
        self.calls.append(("convert", sra))
        run = os.path.basename(sra).split(".")[0]
        n = self.convert_reads if self.convert_reads is not None else self.sra_spots
        if self.paired_sra:
            for k in (1, 2):
                write_fastq(os.path.join(wd, f"{run}_{k}.fastq.gz"), n)
        else:
            write_fastq(os.path.join(wd, f"{run}.fastq.gz"), n)
        return 0, 0.5, ""

    def window(self, sra, start, count, outdir, tmpdir):
        self.calls.append(("window", start, count))
        run = os.path.basename(sra).split(".")[0]
        if self.paired_sra:
            for k in (1, 2):
                write_fastq(os.path.join(outdir, f"{run}_{k}.fastq.gz"), count, start)
        else:
            write_fastq(os.path.join(outdir, f"{run}_1.fastq.gz"), count, start)
        return 0, 0.1, ""

    def ena_head(self, url, dest, n_reads, timeout=None):
        self.calls.append(("ena_head", url))
        n = min(n_reads, self.ena_reads)
        write_fastq(dest, n)
        return 0, 0.1, n, ""

    def md5(self, path):
        return "bad" if self.corrupt_md5 else md5(path)

    def tool_versions(self):
        return {"fake": "1"}


def plan_for(
    run="SRR1",
    layout="PAIRED",
    ncbi="PAIRED",
    spots=1000,
    ena_files=2,
    base_count=None,
    genome=None,
    ncbi_status="ok",
):
    """A plan with matching ENA MD5s computed on demand by `ena_md5_for`."""
    names = (
        [f"{run}_1.fastq.gz", f"{run}_2.fastq.gz"]
        if ena_files == 2
        else ([f"{run}.fastq.gz"] if ena_files == 1 else [])
    )
    p = Plan(
        run=run,
        taxid=5141,
        species="Neurospora crassa",
        sample="SAMN1",
        fastq_ftp=";".join("ftp.sra.ebi.ac.uk/vol1/fastq/" + n for n in names),
        base_count=base_count,
        genome_mb=genome,
    )
    p.ncbi_status = ncbi_status
    if ncbi_status == "ok":
        p.ncbi_spots = spots
        p.ncbi_mates = {"PAIRED": spots, "SINGLE": 0, "MIXED": spots // 2}[ncbi]
    return p


def with_md5(plan, reads):
    """Fill `ena_md5` with the MD5 of the files the fake writes for `reads` reads."""
    import tempfile

    d = tempfile.mkdtemp()
    for u in plan.fastq_ftp.split(";"):
        if not u:
            continue
        nm = os.path.basename(u)
        write_fastq(os.path.join(d, nm), reads)
        plan.ena_md5[nm] = md5(os.path.join(d, nm))
    return plan


@pytest.fixture
def env(tmp_path, monkeypatch):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("SCRATCH", str(scratch))
    out = tmp_path / "out"
    opt = Options(
        tier="full",
        out=str(out),
        cache_root=str(tmp_path / "cache"),
        downloads_csv=str(tmp_path / "downloads.csv"),
    )
    return {"scratch": scratch, "out": out, "opt": opt, "tmp": tmp_path}


def run_fetch(env, backend, plan, **kw):
    opt = env["opt"]
    for k, v in kw.items():
        setattr(opt, k, v)
    return Fetcher(backend, None).fetch(plan.run, opt, plan=plan)


def rows(env):
    p = env["opt"].downloads_csv
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def assert_scratch_clean(env):
    assert os.listdir(env["scratch"]) == [], "temporary files were left in $SCRATCH"


# ---------------------------------------------------------------- full tier, SRA routes
def test_full_paired_from_aws_is_verified_cached_with_a_sidecar(env):
    plan, be = plan_for(spots=1000), Fake()
    r = run_fetch(env, be, plan)
    assert (r.outcome, r.exit_code, r.route, r.layout) == ("ok", 0, "aws_odp", "PAIRED")
    d = os.path.join(env["opt"].cache_root, "5141", "SRR1")
    assert sorted(os.listdir(d)) == ["SRR1.json", "SRR1_1.fastq.gz", "SRR1_2.fastq.gz"]
    s = json.load(open(os.path.join(d, "SRR1.json")))
    assert s["verified"] and s["route"] == "aws_odp" and s["tier"] == "full"
    assert (
        s["read_count_expected"] == 1000 and s["read_count_mate1"] == s["read_count_mate2"] == 1000
    )
    assert s["layout"] == "PAIRED" and s["verify_scope"] == "full"
    assert [(x["route"], x["outcome"], x["verified"]) for x in rows(env)] == [
        ("aws_odp", "ok", "true")
    ]
    assert_scratch_clean(env)


def test_a_mate_count_that_is_not_ncbi_spots_with_mates_fails_the_task_and_deletes_the_files(env):
    r = run_fetch(env, Fake(convert_reads=999), plan_for(spots=1000))
    assert (r.outcome, r.exit_code) == ("verify_failed", EXIT["verify_failed"])
    assert not os.path.exists(os.path.join(env["opt"].cache_root, "5141", "SRR1"))
    assert [(x["outcome"], x["verified"]) for x in rows(env)] == [("verify_failed", "false")]
    assert_scratch_clean(env)


def test_a_verify_failure_is_not_retried_on_another_route(env):
    be = Fake(convert_reads=5)
    run_fetch(env, be, plan_for(spots=1000))
    assert [c for c in be.calls if c[0] == "transfer"] == [("transfer", "aws_odp")]


def test_aws_404_then_prefetch_absent_then_ena_serves_verified_pairs(env):
    plan = with_md5(plan_for(spots=1000), 1000)
    be = Fake({"aws_odp": "404", "ncbi_prefetch": "404"})
    r = run_fetch(env, be, plan)
    assert (r.outcome, r.route) == ("ok", "ena_https")
    assert [(x["route"], x["outcome"]) for x in rows(env)] == [
        ("aws_odp", "unavailable_at_ncbi"),
        ("ncbi_prefetch", "unavailable_at_ncbi"),
        ("ena_https", "ok"),
    ]
    assert json.load(open(r.sidecar))["ena_md5"]
    assert_scratch_clean(env)


def test_nothing_serves_the_run_is_exit_10_and_not_retried(env):
    plan = plan_for(spots=1000, ena_files=0)
    r = run_fetch(env, Fake({"aws_odp": "404", "ncbi_prefetch": "404"}), plan)
    assert (r.outcome, r.exit_code) == ("unavailable_at_ncbi", 10)


def test_a_run_absent_at_ncbi_goes_straight_to_ena(env):
    plan = with_md5(plan_for(ncbi_status="absent"), 1000)
    be = Fake()
    r = run_fetch(env, be, plan)
    assert (r.outcome, r.route) == ("ok", "ena_https")
    assert not any(c == ("transfer", "aws_odp") for c in be.calls), (
        "an absent run must not be tried at AWS"
    )


def test_mates_come_from_one_source_a_failed_sra_attempt_leaves_nothing_for_ena(env):
    plan = with_md5(plan_for(spots=1000), 1000)
    be = Fake({"aws_odp": "fail", "ncbi_prefetch": "fail"})
    r = run_fetch(env, be, plan)
    assert r.route == "ena_https"
    d = os.path.join(env["opt"].cache_root, "5141", "SRR1")
    assert not any(f.endswith(".sra") for f in os.listdir(d))


# ------------------------------------------------------------- pairing and single-end
def test_ncbi_unpaired_with_require_pairs_skips_sra_and_takes_pairs_from_ena(env):
    plan = with_md5(plan_for(ncbi="SINGLE", spots=1000), 1000)
    be = Fake()
    r = run_fetch(env, be, plan, require_pairs=True)
    assert (r.outcome, r.route, r.layout) == ("ok", "ena_https", "PAIRED")
    assert [c for c in be.calls if c[0] == "transfer"] == [("transfer", "ena_https")]
    sra_rows = [x for x in rows(env) if x["route"] in ("aws_odp", "ncbi_prefetch")]
    assert sra_rows and all(x["outcome"] == "unavailable_at_ncbi" for x in sra_rows)
    assert "not paired at NCBI" in sra_rows[0]["error"]


def test_require_pairs_with_no_pairs_anywhere_is_exit_11(env):
    plan = with_md5(plan_for(ncbi="SINGLE", ena_files=1, spots=1000), 1000)
    r = run_fetch(env, Fake(), plan, require_pairs=True)
    assert (r.outcome, r.exit_code) == ("unpaired", 11)
    assert rows(env)[-1]["outcome"] == "unpaired"


def test_ncbi_unpaired_without_require_pairs_is_fetched_as_one_verified_file(env):
    be = Fake()
    be.paired_sra = False
    r = run_fetch(env, be, plan_for(ncbi="SINGLE", spots=1000))
    assert (r.outcome, r.layout, r.route) == ("ok", "SINGLE", "aws_odp")
    d = os.path.join(env["opt"].cache_root, "5141", "SRR1")
    assert sorted(os.listdir(d)) == ["SRR1.fastq.gz", "SRR1.json"]
    s = json.load(open(os.path.join(d, "SRR1.json")))
    assert s["read_count_expected"] == 1000 and s["read_count_mate2"] is None


def test_a_single_end_sra_run_that_comes_back_as_pairs_fails_verification(env):
    be = Fake()  # NCBI says single-end, the converter returns mates
    r = run_fetch(env, be, plan_for(ncbi="SINGLE", spots=1000))
    assert r.outcome == "verify_failed"


def test_single_end_run_from_ena_is_verified_by_md5_as_one_file(env, monkeypatch):
    plan = with_md5(plan_for(ncbi_status="absent", ena_files=1), 1000)
    be = Fake()
    monkeypatch.setattr(
        Fetcher,
        "_transfer_single",
        staticmethod(
            lambda sub, wd, opt: be.transfer(
                {**sub, "fastq_ftp": sub["fastq_ftp"]}, "ena_https", wd, None
            )
        ),
    )
    r = run_fetch(env, be, plan)
    assert (r.outcome, r.layout, r.route) == ("ok", "SINGLE", "ena_https")
    assert json.load(open(r.sidecar))["files"] == ["SRR1.fastq.gz"]


def test_mixed_paired_and_unpaired_spots_at_ncbi_skip_the_sra_routes(env):
    plan = with_md5(plan_for(ncbi="MIXED", spots=1000), 1000)
    be = Fake()
    r = run_fetch(env, be, plan)
    assert r.route == "ena_https"
    assert "mixed" in rows(env)[0]["error"]


def test_ena_md5_mismatch_is_verify_failed_and_ena_read_count_is_not_used(env):
    plan = with_md5(plan_for(ncbi_status="absent"), 1000)
    be = Fake()
    be.corrupt_md5 = True
    r = run_fetch(env, be, plan)
    assert (r.outcome, r.exit_code) == ("verify_failed", 12)
    # no md5 stored at all is also a failure: nothing to verify against
    plan2 = plan_for(run="SRR2", ncbi_status="absent")
    r2 = Fetcher(Fake(), None).fetch("SRR2", env["opt"], plan=plan2)
    assert r2.outcome == "verify_failed" and "none stored" in r2.message


def test_ena_mates_with_different_counts_fail(env):
    plan = with_md5(plan_for(ncbi_status="absent"), 1000)
    be = Fake()
    orig = be.transfer

    def uneven(p, route, wd, t, a=None):
        res = orig(p, route, wd, t, a)
        write_fastq(res["files"][1], 990)
        plan.ena_md5[os.path.basename(res["files"][1])] = md5(res["files"][1])
        return res

    be.transfer = uneven
    assert run_fetch(env, be, plan).outcome == "verify_failed"


# ------------------------------------------------------------- failures and transients
def test_every_route_failing_to_transfer_is_transient_exit_13(env):
    plan = with_md5(plan_for(), 1000)
    r = run_fetch(
        env, Fake({"aws_odp": "fail", "ncbi_prefetch": "fail", "ena_https": "fail"}), plan
    )
    assert (r.outcome, r.exit_code) == ("transient", 13)
    assert all(x["outcome"] == "transfer_failed" for x in rows(env))


def test_a_failed_ncbi_lookup_is_transient_not_unavailable(env):
    plan = plan_for(ncbi_status="error", ena_files=0)
    r = run_fetch(env, Fake(), plan)
    assert (r.outcome, r.exit_code) == ("transient", 13)


def test_a_missing_tool_is_a_failed_attempt_not_a_crash(env):
    plan = with_md5(plan_for(), 1000)
    r = run_fetch(env, Fake({"aws_odp": "missing", "ncbi_prefetch": "missing"}), plan)
    assert r.outcome == "ok" and r.route == "ena_https"
    assert (
        rows(env)[0]["outcome"] == "transfer_failed"
        and "FileNotFoundError" in rows(env)[0]["error"]
    )


def test_usage_errors(env, monkeypatch):
    f = Fetcher(Fake(), None)
    with pytest.raises(Usage, match="--out"):
        f.fetch("SRR1", Options(tier="windows"), plan=plan_for())
    with pytest.raises(Usage, match="--coverage"):
        f.fetch("SRR1", Options(tier="coverage", out="o"), plan=plan_for())
    with pytest.raises(Usage, match="genome size"):
        f.fetch("SRR1", Options(tier="coverage", coverage=30, out="o"), plan=plan_for())
    monkeypatch.delenv("SCRATCH")
    with pytest.raises(Usage, match="SCRATCH"):
        f.fetch("SRR1", Options(tier="full"), plan=plan_for())


# ------------------------------------------------------------------- resume and sidecar
def test_a_second_call_skips_a_verified_run(env):
    plan, be = plan_for(), Fake()
    run_fetch(env, be, plan)
    n = len(be.calls)
    r = run_fetch(env, be, plan)
    assert r.resumed and r.outcome == "ok" and len(be.calls) == n


def test_a_run_without_a_sidecar_is_not_trusted(env):
    plan, be = plan_for(), Fake()
    r = run_fetch(env, be, plan)
    os.remove(r.sidecar)  # the files are there, the sidecar is not: refetch and verify
    n = len(be.calls)
    r2 = run_fetch(env, be, plan)
    assert not r2.resumed and len(be.calls) > n and os.path.exists(r2.sidecar)


def test_a_missing_output_file_forces_a_refetch(env):
    plan, be = plan_for(), Fake()
    r = run_fetch(env, be, plan)
    os.remove(r.outputs[0])
    assert not run_fetch(env, be, plan).resumed


def test_changed_parameters_are_not_resumed(env):
    plan, be = plan_for(genome=0.03, base_count=300_000), Fake()
    run_fetch(env, be, plan, tier="coverage", coverage=10.0)
    assert run_fetch(env, be, plan, tier="coverage", coverage=10.0).resumed
    assert not run_fetch(env, be, plan, tier="coverage", coverage=20.0).resumed


def test_a_killed_task_leaves_no_complete_looking_output(env):
    """The sidecar goes last: a kill between the FASTQ and the sidecar leaves no sidecar."""
    plan, be = plan_for(), Fake()
    orig = os.replace

    def dying(src, dst):
        if str(dst).endswith(".json"):
            raise KeyboardInterrupt
        return orig(src, dst)

    import popseq.fetch as mod

    mod.os.replace = dying
    try:
        with pytest.raises(KeyboardInterrupt):
            run_fetch(env, be, plan)
    finally:
        mod.os.replace = orig
    d = os.path.join(env["opt"].cache_root, "5141", "SRR1")
    assert not os.path.exists(os.path.join(d, "SRR1.json"))
    assert_scratch_clean(env)
    assert run_fetch(env, be, plan).outcome == "ok"  # and a rerun completes it


# --------------------------------------------------------------------------- the hook
def test_the_hook_runs_after_verification_and_the_reads_are_gone_afterwards(env, tmp_path):
    log = tmp_path / "hook.log"
    script = tmp_path / "hook.sh"
    script.write_text(
        f'#!/bin/sh\necho "$@" >> {log}\nfor f in "${{@:2}}"; do test -s "$f" || exit 9; done\n'.replace(
            '"${@:2}"', '"$2" "$3"'
        )
    )
    script.chmod(0o755)
    r = run_fetch(
        env, Fake(), plan_for(spots=1000), tier="windows", spots=100, exec_cmd=str(script)
    )
    assert r.outcome == "ok" and r.outputs == []
    line = log.read_text().split()
    assert line[0] == "SRR1" and line[1].endswith("_1.fastq.gz") and line[2].endswith("_2.fastq.gz")
    assert os.listdir(env["out"]) == ["SRR1.json"], "no reads are kept after a hook"
    s = json.load(open(os.path.join(env["out"], "SRR1.json")))
    assert s["reads_on_disk"] is False and s["verified"]
    assert_scratch_clean(env)


def test_a_failing_hook_is_exit_14_and_the_download_is_still_ok(env, tmp_path):
    r = run_fetch(env, Fake(), plan_for(spots=1000), tier="windows", spots=100, exec_cmd="false")
    assert (r.outcome, r.exit_code, r.verified) == ("hook_failed", 14, True)
    assert rows(env)[-1]["outcome"] == "ok"


def test_the_hook_does_not_run_on_unverified_reads(env, tmp_path):
    marker = tmp_path / "ran"
    r = run_fetch(env, Fake(convert_reads=3), plan_for(spots=1000), exec_cmd=f"touch {marker}")
    assert r.outcome == "verify_failed" and not marker.exists()


def test_keep_with_exec_keeps_the_reads(env):
    r = run_fetch(
        env, Fake(), plan_for(spots=1000), tier="windows", spots=100, exec_cmd="true", keep=True
    )
    assert len(r.outputs) == 2 and all(os.path.exists(o) for o in r.outputs)


# ---------------------------------------------------------- subsampled tiers from SRA
def test_windows_tier_reads_evenly_spaced_windows_and_checks_each(env):
    be = Fake()
    r = run_fetch(env, be, plan_for(spots=10_000), tier="windows", spots=800, windows=8)
    assert r.outcome == "ok" and r.reads_kept == 800
    assert [c[1:] for c in be.calls if c[0] == "window"] == ss.windows(10_000, 800, 8)
    s = json.load(open(r.sidecar))
    assert (
        s["verify_scope"] == "subsample" and s["read_count_mate1"] == s["read_count_mate2"] == 800
    )
    out1 = os.path.join(env["out"], "SRR1_1.fastq.gz")
    names = [l.split("/")[0][1:] for l in gzip.open(out1, "rt").read().split("\n")[::4] if l]
    assert names[0] == "spot0" and names[-1] == "spot9999"  # windows touch both ends of the run
    assert_scratch_clean(env)


def test_a_window_that_returns_the_wrong_number_of_reads_fails(env):
    be = Fake()
    orig = be.window

    def short(sra, start, count, outdir, tmpdir):
        return orig(sra, start, count - 1, outdir, tmpdir)

    be.window = short
    r = run_fetch(env, be, plan_for(spots=10_000), tier="windows", spots=800)
    assert r.outcome == "verify_failed" and "expected" in r.message


def test_head_tier_from_ena_streams_and_does_not_download_the_file(env):
    plan = with_md5(plan_for(spots=1000), 1000)
    be = Fake()
    r = run_fetch(env, be, plan, tier="head", spots=100)
    assert (r.outcome, r.route, r.reads_kept) == ("ok", "ena_https", 100)
    assert [c[0] for c in be.calls].count("ena_head") == 2 and (
        "transfer",
        "ena_https",
    ) not in be.calls


def test_head_tier_falls_back_to_the_sra_when_ena_has_no_files(env):
    r = run_fetch(env, Fake(), plan_for(spots=1000, ena_files=0), tier="head", spots=50)
    assert (r.outcome, r.route, r.reads_kept) == ("ok", "aws_odp", 50)


def test_stride_tier_verifies_the_whole_run_first_then_keeps_every_nth(env):
    be = Fake(sra_spots=1000)
    r = run_fetch(env, be, plan_for(spots=1000), tier="stride", spots=100)
    assert r.outcome == "ok" and r.reads_kept == 100
    assert json.load(open(r.sidecar))["verify_scope"] == "full"
    r2 = run_fetch(
        env, Fake(convert_reads=900), plan_for(run="SRR3", spots=1000), tier="stride", spots=100
    )
    assert r2.outcome == "verify_failed", "the full-run check comes before the subsample"


# ------------------------------------------------------------------ the coverage tier
def test_coverage_tier_keeps_bases_equal_to_genome_times_coverage(env):
    # 10,000 spots x 100 bases (2 mates x 50) = 1 Mb; genome 0.03 Mb x 10x = 0.3 Mb -> 3,000 spots
    plan = plan_for(spots=10_000, base_count=1_000_000, genome=0.03)
    be = Fake()
    r = run_fetch(env, be, plan, tier="coverage", coverage=10.0, sampling="windows")
    assert r.outcome == "ok"
    assert r.reads_kept == 3000 and r.bases_kept == 300_000
    assert r.coverage_target == 10.0 and r.coverage_reached == pytest.approx(10.0)
    s = json.load(open(r.sidecar))
    assert s["coverage_reached"] == pytest.approx(10.0) and len(s["windows"]) == 8
    assert r.row()[0:2] == ["SRR1", "coverage"]


def test_coverage_genome_size_option_overrides_a_missing_catalog_size(env):
    plan = plan_for(spots=10_000, base_count=1_000_000, genome=None)
    r = run_fetch(
        env, Fake(), plan, tier="coverage", coverage=10.0, sampling="windows", genome_mb=0.03
    )
    assert r.outcome == "ok" and r.coverage_reached == pytest.approx(10.0)


def test_a_run_shorter_than_the_target_keeps_the_whole_run_and_reports_what_it_reached(env):
    plan = plan_for(spots=1000, base_count=100_000, genome=0.03)  # needs 3,000 spots, has 1,000
    r = run_fetch(env, Fake(), plan, tier="coverage", coverage=10.0, sampling="windows")
    assert (r.outcome, r.reads_kept) == ("ok", 1000)
    assert r.coverage_reached == pytest.approx(10.0 / 3, rel=1e-6) and r.coverage_target == 10.0


def test_coverage_from_ena_is_cut_from_the_md5_verified_file(env):
    plan = with_md5(plan_for(ncbi_status="absent", genome=0.01), 2000)
    be = Fake(ena_reads=2000)
    r = run_fetch(env, be, plan, tier="coverage", coverage=10.0, sampling="windows")
    # 2,000 reads x 100 bases: 0.01 Mb x 10x = 100,000 bases = 1,000 reads
    assert (r.outcome, r.route, r.reads_kept) == (
        "ok",
        "ena_https",
        1000,
    )
    assert r.coverage_reached == pytest.approx(10.0, rel=0.05)


# -------------------------------------------------------- manifests, CLI, exit codes
def test_read_manifest_accepts_csv_tsv_and_gzip(tmp_path):
    for name, text, opener in [
        ("m.csv", "run_accession,x\nSRR1,1\nSRR2,2\n", open),
        ("m.tsv", "run_accession\tx\nSRR1\t1\nSRR2\t2\n", open),
        ("m.csv.gz", "run_accession,x\nSRR1,1\nSRR2,2\n", lambda p, m: gzip.open(p, "wt")),
    ]:
        p = tmp_path / name
        with opener(p, "w") as fh:
            fh.write(text)
        assert fx.read_manifest(str(p)) == ["SRR1", "SRR2"]
    bad = tmp_path / "bad.csv"
    bad.write_text("acc\nSRR1\n")
    with pytest.raises(Usage):
        fx.read_manifest(str(bad))


def test_the_result_row_matches_the_documented_fields(env):
    r = run_fetch(env, Fake(), plan_for(spots=1000))
    assert len(r.row()) == len(fx.FETCH_FIELDS)
    assert r.row()[fx.FETCH_FIELDS.index("outcome")] == "ok"
    assert r.row()[fx.FETCH_FIELDS.index("verified")] == "true"


def test_exit_codes_are_the_documented_values():
    assert EXIT == {
        "ok": 0,
        "usage": 2,
        "unavailable_at_ncbi": 10,
        "unpaired": 11,
        "verify_failed": 12,
        "transient": 13,
        "hook_failed": 14,
    }


def test_the_download_log_accepts_the_coverage_tier(env):
    run_fetch(
        env,
        Fake(),
        plan_for(spots=10_000, base_count=1_000_000, genome=0.03),
        tier="coverage",
        coverage=5.0,
        sampling="windows",
    )
    assert rows(env)[-1]["tier"] == "coverage"


# ---------------------------------------------------------------- subsample arithmetic
def test_windows_are_evenly_spaced_equal_and_touch_both_ends():
    w = ss.windows(1000, 100, 8)
    assert len(w) == 8 and {c for _s, c in w} == {13}
    assert w[0][0] == 0 and w[-1][0] + w[-1][1] == 1000
    gaps = [w[i + 1][0] - w[i][0] for i in range(7)]
    assert max(gaps) - min(gaps) <= 1


def test_windows_never_overlap_and_a_short_run_is_one_window():
    for total, n in [(103, 100), (100, 99), (1000, 999), (50, 100)]:
        w = ss.windows(total, n, 8)
        assert all(s0 + c0 <= s1 for (s0, c0), (s1, _c) in zip(w, w[1:]))
        assert all(0 <= s and s + c <= total for s, c in w)
    assert ss.windows(50, 100, 8) == [(0, 50)]


def test_spots_for_coverage_and_stride():
    assert ss.spots_for_coverage(40, 30, 300) == 4_000_000
    assert ss.stride(1000, 100) == 10 and ss.stride(50, 100) == 1
    with pytest.raises(ValueError):
        ss.spots_for_coverage(0, 30, 300)


def test_select_spots_keeps_mates_in_step_and_refuses_unequal_mates(tmp_path):
    a, b, c = (str(tmp_path / n) for n in ("a.gz", "b.gz", "c.gz"))
    write_fastq(a, 20)
    write_fastq(b, 20)
    write_fastq(c, 19)
    oa, ob = str(tmp_path / "oa.gz"), str(tmp_path / "ob.gz")
    kept, read, bases = ss.select_spots([a, b], [oa, ob], ss.Every(5))
    assert (kept, read, bases) == (4, 20, [200, 200])
    na = [l for l in gzip.open(oa, "rt").read().split("\n")[::4] if l]
    nb = [l for l in gzip.open(ob, "rt").read().split("\n")[::4] if l]
    assert na == nb and na[0].startswith("@spot0") and na[1].startswith("@spot5")
    with pytest.raises(ValueError, match="different numbers"):
        ss.select_spots([a, c], [oa, ob], ss.Every(1))


def test_every_log_row_and_the_sidecar_carry_the_declared_concurrency(env):
    plan = with_md5(plan_for(), 1000)
    be = Fake({"aws_odp": "404", "ncbi_prefetch": "404"})
    r = run_fetch(env, be, plan, concurrency=7)
    assert [x["concurrency"] for x in rows(env)] == ["7", "7", "7"]
    assert json.load(open(r.sidecar))["concurrency"] == 7
    assert fx.DEFAULT_CONCURRENCY == 20 and Options().concurrency == 20


def test_genome_size_empty_catalog_value_is_unknown_not_a_crash(tmp_path):
    # Issue #36 added species with no NCBI assembly as rows with an empty size.
    p = tmp_path / "genome_sizes.csv"
    p.write_text(
        "species_taxid,species_name,genome_mb,accession,assembly_level,note\n"
        "1,Aus bus,31.5,GCA_1,Contig,\n"
        "2,Aus cus,,,,no reference assembly at NCBI\n"
    )
    assert fx.genome_size(1, str(p)) == 31.5
    assert fx.genome_size(2, str(p)) is None
    assert fx.genome_size(3, str(p)) is None


# ------------------------------------------------------------- the evenly spaced sampler
def test_evenly_keeps_exactly_n_spread_over_the_whole_file():
    for total, n in [(1000, 100), (1000, 333), (4_387_476, 1_225_000), (10, 7), (7, 7), (5, 9)]:
        keep = ss.Evenly(total, n)
        kept = [i for i in range(total) if keep(i)] if total < 50_000 else None
        if kept is not None:
            assert len(kept) == min(n, total)
            gaps = [b - a for a, b in zip(kept, kept[1:])]
            if gaps:
                assert max(gaps) - min(gaps) <= 1, (total, n)
            assert (
                kept[0] <= total // n + 1 and kept[-1] >= total - total // n - 2
            )  # both ends are reached
        else:
            assert sum(1 for i in range(total) if keep(i)) == n
    with pytest.raises(ValueError):
        ss.Evenly(0, 5)


def test_evenly_reaches_every_part_of_a_position_sorted_file_where_windows_do_not():
    """A file sorted by genome position: spot i is at position i. 18% of positions are covered by 8 contiguous windows
    holding 18% of the spots; evenly spaced spots of the same number cover the whole range."""
    total, n = 100_000, 18_000
    windows_pos = {i for s, c in ss.windows(total, n, 8) for i in range(s, s + c)}
    even = ss.Evenly(total, n)
    even_pos = [i for i in range(total) if even(i)]
    bins = 100
    covered = lambda pos: len({p * bins // total for p in pos})  # noqa: E731
    assert covered(windows_pos) < 35, "8 windows touch about a quarter of the 100 regions"
    assert covered(even_pos) == bins, "evenly spaced spots touch every region"


# ------------------------------------------- tier coverage by evenly spaced reads (the default)
def hook_script(tmp_path, body, name="hook.sh"):
    p = tmp_path / name
    p.write_text("#!/bin/sh\n" + body)
    p.chmod(0o755)
    return str(p)


def spot_numbers(path):
    names = [l.split("/")[0][1:] for l in gzip.open(path, "rt").read().split("\n")[::4] if l]
    return [int(n.replace("spot", "")) for n in names]


def test_coverage_tier_is_evenly_spaced_from_the_verified_full_run_by_default(env):
    plan = plan_for(spots=10_000, genome=0.03)  # no base_count needed: bases per spot are measured
    be = Fake(sra_spots=10_000)
    r = run_fetch(env, be, plan, tier="coverage", coverage=10.0)
    assert (r.outcome, r.sampling, r.fallback_from) == ("ok", "stride", "")
    assert r.reads_kept == 3000 and r.bases_kept == 300_000
    assert r.coverage_reached == pytest.approx(10.0)
    kinds = [c[0] for c in be.calls]
    assert kinds.count("convert") == 1 and "window" not in kinds
    s = json.load(open(r.sidecar))
    assert s["verify_scope"] == "full" and s["sampling"] == "stride" and s["windows"] is None
    nums = spot_numbers(os.path.join(env["out"], "SRR1_1.fastq.gz"))
    assert len(nums) == 3000 and nums[0] <= 5 and nums[-1] >= 9990, (
        "both ends of the run are reached"
    )
    assert (
        max(b - a for a, b in zip(nums, nums[1:])) - min(b - a for a, b in zip(nums, nums[1:])) <= 1
    )
    assert_scratch_clean(env)


def test_the_full_run_check_comes_before_the_cut_for_coverage_by_stride(env):
    r = run_fetch(
        env,
        Fake(convert_reads=9_000),
        plan_for(spots=10_000, genome=0.03),
        tier="coverage",
        coverage=10.0,
    )
    assert r.outcome == "verify_failed"


def test_a_run_shorter_than_the_target_is_kept_whole_by_stride(env):
    r = run_fetch(env, Fake(), plan_for(spots=1000, genome=0.03), tier="coverage", coverage=10.0)
    assert (r.outcome, r.reads_kept) == ("ok", 1000)  # 3,000 spots needed, 1,000 held
    assert r.coverage_reached == pytest.approx(100_000 / 0.03e6) and r.coverage_target == 10.0


def test_coverage_by_stride_from_ena_is_cut_from_the_md5_verified_files(env):
    plan = with_md5(plan_for(ncbi_status="absent", genome=0.01), 2000)
    r = run_fetch(env, Fake(ena_reads=2000), plan, tier="coverage", coverage=10.0)
    assert (r.outcome, r.route, r.sampling, r.reads_kept) == ("ok", "ena_https", "stride", 1000)
    assert r.coverage_reached == pytest.approx(10.0)


def test_the_hook_gets_the_genome_size_the_coverage_and_the_sampling(env, tmp_path):
    log = tmp_path / "env.log"
    hook = hook_script(
        tmp_path,
        f'echo "$POPSEQ_RUN $POPSEQ_TIER $POPSEQ_SAMPLING $POPSEQ_GENOME_MB $POPSEQ_COVERAGE $POPSEQ_ROUTE $POPSEQ_LAYOUT" > {log}\n',
    )
    run_fetch(
        env,
        Fake(sra_spots=10_000),
        plan_for(spots=10_000, genome=0.03),
        tier="coverage",
        coverage=10.0,
        exec_cmd=hook,
    )
    assert log.read_text().split() == [
        "SRR1",
        "coverage",
        "stride",
        "0.03",
        "10.000",
        "aws_odp",
        "PAIRED",
    ]


# ----------------------------------------------- windows with a fallback to stride (hook exit 75)
FLIP = """n=$(cat {c} 2>/dev/null || echo 0); n=$((n+1)); echo $n > {c}
echo "$POPSEQ_SAMPLING $2" >> {l}
[ "$n" -eq {fail_times} ] && exit 75
[ "$n" -lt {fail_times} ] && exit 75
exit 0
"""


def flip_hook(tmp_path, fail_times=1):
    c, l = tmp_path / "count", tmp_path / "calls.log"
    return hook_script(tmp_path, FLIP.format(c=c, l=l, fail_times=fail_times)), l


def test_a_hook_that_rejects_windowed_reads_makes_fetch_resample_by_stride_without_a_new_download(
    env, tmp_path
):
    hook, log = flip_hook(tmp_path)
    plan, be = plan_for(spots=10_000, genome=0.03, base_count=1_000_000), Fake(sra_spots=10_000)
    r = run_fetch(env, be, plan, tier="coverage", coverage=10.0, sampling="windows", exec_cmd=hook)
    assert (r.outcome, r.exit_code, r.sampling, r.fallback_from) == ("ok", 0, "stride", "windows")
    assert r.reads_kept == 3000 and r.coverage_reached == pytest.approx(10.0)
    calls = [l.split()[0] for l in log.read_text().splitlines()]
    assert calls == ["windows", "stride"], "the hook ran on the windows, then on the resample"
    kinds = [c for c in be.calls]
    assert kinds.count(("transfer", "aws_odp")) == 1, "the .sra was downloaded once"
    assert (
        sum(1 for c in kinds if c[0] == "window") == 8
        and sum(1 for c in kinds if c[0] == "convert") == 1
    )
    s = json.load(open(r.sidecar))
    assert s["sampling"] == "stride" and s["fallback_from"] == "windows" and s["hook_exit"] == 0
    assert s["verify_scope"] == "full" and s["windows"] is None
    last = rows(env)[-1]
    assert last["outcome"] == "ok" and "resampled evenly spaced" in last["error"]
    assert "resampled evenly spaced" in r.message
    assert_scratch_clean(env)


def test_the_fallback_also_works_for_an_ena_transfer_and_does_not_download_again(env, tmp_path):
    hook, log = flip_hook(tmp_path)
    plan = with_md5(plan_for(ncbi_status="absent", genome=0.01), 2000)
    be = Fake(ena_reads=2000)
    r = run_fetch(env, be, plan, tier="coverage", coverage=10.0, sampling="windows", exec_cmd=hook)
    assert (r.outcome, r.route, r.sampling, r.fallback_from, r.reads_kept) == (
        "ok",
        "ena_https",
        "stride",
        "windows",
        1000,
    )
    assert be.calls.count(("transfer", "ena_https")) == 1
    assert [l.split()[0] for l in log.read_text().splitlines()] == ["windows", "stride"]


def test_tier_windows_falls_back_to_stride_with_the_same_number_of_spots(env, tmp_path):
    hook, _ = flip_hook(tmp_path)
    r = run_fetch(
        env,
        Fake(sra_spots=10_000),
        plan_for(spots=10_000),
        tier="windows",
        spots=800,
        exec_cmd=hook,
    )
    assert (r.outcome, r.sampling, r.fallback_from, r.reads_kept) == (
        "ok",
        "stride",
        "windows",
        800,
    )


def test_if_the_check_fails_even_by_stride_the_run_is_hook_failed_not_looped(env, tmp_path):
    hook, log = flip_hook(tmp_path, fail_times=99)
    r = run_fetch(
        env,
        Fake(sra_spots=10_000),
        plan_for(spots=10_000, genome=0.03),
        tier="coverage",
        coverage=10.0,
        sampling="windows",
        exec_cmd=hook,
    )
    assert (r.outcome, r.exit_code) == ("hook_failed", 14)
    assert "even with evenly spaced sampling" in r.message and r.fallback_from == "windows"
    assert len(log.read_text().splitlines()) == 2, "the hook ran twice, no more"


def test_a_stride_or_full_sample_that_fails_the_check_is_hook_failed_with_no_resampling(
    env, tmp_path
):
    for tier in ("coverage", "full"):
        d = tmp_path / tier
        d.mkdir()
        hook, log = flip_hook(d, fail_times=99)
        extra = {"coverage": 10.0} if tier == "coverage" else {}
        plan = plan_for(run=f"R{tier}", spots=10_000, genome=0.03)
        r = run_fetch(env, Fake(sra_spots=10_000), plan, tier=tier, exec_cmd=hook, **extra)
        assert (r.outcome, r.exit_code) == ("hook_failed", 14), tier
        assert len(log.read_text().splitlines()) == 1, f"{tier}: the hook ran once"
        assert r.fallback_from == "", tier


def test_other_hook_failures_do_not_trigger_a_resample(env, tmp_path):
    r = run_fetch(
        env,
        Fake(sra_spots=10_000),
        plan_for(spots=10_000, genome=0.03),
        tier="coverage",
        coverage=10.0,
        sampling="windows",
        exec_cmd="false",
    )
    assert (r.outcome, r.fallback_from) == ("hook_failed", "")


def test_a_failed_resample_reports_why_and_keeps_the_exit_code_of_the_cause(env, tmp_path):
    hook, _ = flip_hook(tmp_path)
    be = Fake(sra_spots=10_000)
    orig = be.convert_sra

    def short(sra, wd, threads):  # the conversion for the resample returns the wrong count
        be.convert_reads = 9_000
        return orig(sra, wd, threads)

    be.convert_sra = short
    r = run_fetch(
        env,
        be,
        plan_for(spots=10_000, genome=0.03),
        tier="coverage",
        coverage=10.0,
        sampling="windows",
        exec_cmd=hook,
    )
    assert (r.outcome, r.exit_code) == ("verify_failed", 12) and "resampling" in r.message


def test_changing_the_sampling_is_a_new_request_not_a_resume(env):
    plan, be = plan_for(spots=1000, genome=0.03), Fake()
    run_fetch(env, be, plan, tier="coverage", coverage=10.0)
    assert run_fetch(env, be, plan, tier="coverage", coverage=10.0).resumed
    assert not run_fetch(env, be, plan, tier="coverage", coverage=10.0, sampling="windows").resumed


def test_the_sampling_option_is_validated_and_in_the_result_row_and_the_cli(env, capsys):
    with pytest.raises(Usage, match="sampling"):
        run_fetch(env, Fake(), plan_for(), sampling="random")
    r = run_fetch(
        env,
        Fake(),
        plan_for(spots=1000, genome=0.03),
        tier="coverage",
        coverage=10.0,
        sampling="stride",
    )
    row = dict(zip(fx.FETCH_FIELDS, r.row()))
    assert row["sampling"] == "stride" and row["fallback_from"] == ""
    from popseq import cli

    with pytest.raises(SystemExit):
        cli.main(["fetch", "--help"])
    assert "--sampling" in capsys.readouterr().out
