"""popseq.profile: seconds per stage and peak memory of a fetch.

The process tree is read from a fake /proc on every platform; a real child that allocates memory is
checked on Linux only (macOS has no /proc).
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time

import pytest
from test_fetch import (  # noqa: F401  (fixtures: env)
    Fake,
    env,
    plan_for,
    with_md5,
)

from popseq import profile
from popseq.fetch import Fetcher


def fake_proc(root, procs):
    """procs: {pid: (ppid, vmrss_kb)} written as /proc/<pid>/{stat,status}."""
    for pid, (ppid, kb) in procs.items():
        d = root / str(pid)
        d.mkdir(parents=True)
        # comm with a space and a parenthesis: the parser must cut at the last ")"
        (d / "stat").write_text(f"{pid} (bad (name) x) S {ppid} 1 1 0\n")
        (d / "status").write_text(f"Name:\tx\nVmRSS:\t{kb} kB\n")
    (root / "self").mkdir()  # non-numeric entries are ignored


def test_tree_rss_adds_the_process_and_all_its_descendants(tmp_path):
    fake_proc(
        tmp_path,
        {10: (1, 1024), 11: (10, 2048), 12: (11, 4096), 13: (10, 1024), 99: (1, 999999)},
    )
    assert profile.tree_rss_mib(10, str(tmp_path)) == pytest.approx(8.0)  # 99 is not a descendant
    assert profile.tree_rss_mib(12, str(tmp_path)) == pytest.approx(4.0)


def test_tree_rss_is_none_without_a_process_table(tmp_path):
    assert profile.tree_rss_mib(1, str(tmp_path / "nope")) is None


def test_a_process_that_vanishes_is_skipped(tmp_path):
    fake_proc(tmp_path, {10: (1, 1024), 11: (10, 2048)})
    os.remove(tmp_path / "11" / "status")
    assert profile.tree_rss_mib(10, str(tmp_path)) == pytest.approx(1.0)


def test_stages_do_not_nest_and_other_is_the_rest():
    p = profile.Profiler(proc_root="/nonexistent").start()
    with p.stage("transfer"):
        time.sleep(0.05)
        with p.stage("md5"):  # inside another stage: not timed again
            time.sleep(0.05)
    time.sleep(0.05)
    s = p.stop()
    assert set(s["seconds"]) == {"transfer"}
    assert s["seconds"]["transfer"] >= 0.09
    assert s["other_seconds"] >= 0.04
    assert s["tree_peak_mib"] is None and s["tree_sampled"] is False


def test_the_peak_is_attributed_to_the_stage_it_was_seen_in(tmp_path):
    fake_proc(tmp_path, {10: (1, 1024)})
    p = profile.Profiler(interval=0.01, proc_root=str(tmp_path), root_pid=10).start()
    with p.stage("convert"):
        (tmp_path / "10" / "status").write_text("VmRSS:\t5242880 kB\n")  # 5 GiB
        time.sleep(0.1)
        (tmp_path / "10" / "status").write_text("VmRSS:\t1024 kB\n")
    s = p.stop()
    assert s["tree_peak_mib"] == pytest.approx(5120.0)
    assert s["tree_peak_stage"] == "convert"
    assert s["tree_peak_by_stage_mib"]["convert"] == pytest.approx(5120.0)


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="needs /proc")
def test_a_child_that_allocates_is_seen_in_the_tree_peak():
    p = profile.Profiler(interval=0.05).start()
    code = "b = bytearray(300 * 1024 * 1024); import time; time.sleep(0.6)"
    with p.stage("convert"):
        subprocess.run([sys.executable, "-c", code], check=True)
    s = p.stop()
    assert s["tree_peak_mib"] > 250
    assert s["tree_peak_stage"] == "convert"
    assert s["child_peak_mib"] > 250


def test_the_backend_proxy_times_calls_by_stage_and_changes_nothing(tmp_path):
    class B:
        def md5(self, path):
            time.sleep(0.03)
            return "abc"

        def tool_versions(self):
            return {"x": 1}

    p = profile.Profiler(proc_root="/nonexistent")
    be = profile.TimedBackend(B(), p)
    assert be.md5("f") == "abc"
    assert be.tool_versions() == {"x": 1}
    assert p.seconds["md5"] >= 0.03 and "tool_versions" not in p.seconds


# ------------------------------------------------------------------ through the fetcher
def test_a_fetch_reports_where_its_time_went_and_puts_it_in_the_sidecar(env):
    plan = with_md5(plan_for(spots=1000), 1000)
    r = Fetcher(Fake(), None).fetch(plan.run, env["opt"], plan=plan)
    assert r.outcome == "ok"
    sec = r.profile["seconds"]
    assert {"transfer", "convert", "count", "finish"} <= set(sec)
    assert r.profile["python_peak_mib"] > 0
    side = json.load(open(r.sidecar))
    assert set(side["profile"]) >= {"seconds", "other_seconds", "python_peak_mib", "child_peak_mib"}
    # the backend is not left wrapped
    f = Fetcher(Fake(), None)
    be = f.backend
    f.fetch(plan.run, env["opt"], plan=plan_for(run="SRR2", spots=1000))
    assert f.backend is be and f._prof is None


def test_a_failed_fetch_still_has_a_profile(env):
    plan = plan_for(spots=1000, ena_files=0)
    r = Fetcher(Fake({"aws_odp": "404", "ncbi_prefetch": "404"}), None).fetch(
        plan.run, env["opt"], plan=plan
    )
    assert r.outcome != "ok" and r.profile and "other_seconds" in r.profile


def test_the_profile_file_gets_a_header_once_and_one_line_per_fetch(env, tmp_path):
    path = str(tmp_path / "profile.tsv")
    f = Fetcher(Fake(), None)
    for run in ("SRR1", "SRR2"):
        plan = with_md5(plan_for(run=run, spots=1000), 1000)
        r = f.fetch(run, env["opt"], plan=plan)
        profile.append_profile(path, r, r.profile, threads=4)
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    assert [x["run"] for x in rows] == ["SRR1", "SRR2"]
    assert list(rows[0]) == profile.PROFILE_FIELDS
    assert rows[0]["outcome"] == "ok" and rows[0]["threads"] == "4"
    assert open(path).read().count("run\ttier") == 1
