"""Where a `popseq fetch` spends its time and memory (issue: the phase 1 species run could not say).

The GmG phase 1 run (299 runs, 4 fetches per job) found that 80 percent of each fetch was time that the
sidecar timers (transfer, conversion) do not name, and that the largest job reached 99.98 percent of a
20 GiB request with a kept output of only about 1 GB per run. Neither could be attributed. This module
records, for every fetch:

* seconds per stage: `plan` (lookups), `transfer`, `convert` (parallel-fastq-dump, fastq-dump),
  `count` (read counts of the files), `md5`, `cut` (the stride or window cut, the concatenation),
  `postcut` (the checks after a cut), `check` (read-name checks outside a cut), `hook`, `finish`
  (moving the files and writing the sidecar);
* the memory: the peak resident size of this Python process, of the largest process it started, and, on
  Linux, the peak of the whole process tree (this process and every descendant added together, as the
  scheduler counts it) and the stage during which that peak was seen.

Stages do not nest: a stage entered while another is active runs without its own timer, so no second is
counted twice and `other` (the total minus the stages) is what is left over.

The tree sampler reads `/proc` every `interval` seconds, so a spike shorter than that can be missed; the
peak is a lower bound. On a platform without `/proc` (macOS) the tree numbers are None.
"""

from __future__ import annotations

import os
import socket
import sys
import threading
import time
from contextlib import contextmanager
from typing import ClassVar, Optional

STAGES = (
    "plan",
    "transfer",
    "convert",
    "count",
    "md5",
    "cut",
    "postcut",
    "check",
    "hook",
    "finish",
)
PROFILE_FIELDS = (
    ["run", "tier", "sampling", "route", "outcome", "host", "threads", "seconds"]
    + [f"{s}_s" for s in STAGES]
    + ["other_s", "python_peak_mib", "child_peak_mib", "tree_peak_mib", "tree_peak_stage"]
)
_MIB = 1 << 20


def _rusage_mib(who: int) -> float:
    """Peak resident size in MiB of this process (RUSAGE_SELF) or of its largest finished child."""
    import resource

    kb = resource.getrusage(who).ru_maxrss
    # Linux reports kilobytes, macOS bytes
    return kb / _MIB if sys.platform == "darwin" else kb / 1024.0


def _read(path: str) -> str:
    with open(path) as fh:
        return fh.read()


def tree_rss_mib(root_pid: int, proc_root: str = "/proc") -> Optional[float]:
    """Resident size, in MiB, of `root_pid` and all its descendants added together, from `proc_root`.

    None if `proc_root` has no process table (macOS). A process that vanishes while being read is skipped.
    """
    try:
        pids = [p for p in os.listdir(proc_root) if p.isdigit()]
    except OSError:
        return None
    children: dict[int, list[int]] = {}
    for p in pids:
        try:
            stat = _read(os.path.join(proc_root, p, "stat"))
            # `pid (comm) state ppid ...`: comm may hold spaces and parentheses, so cut at the last ")"
            ppid = int(stat[stat.rindex(")") + 2 :].split()[1])
        except (OSError, ValueError, IndexError):
            continue
        children.setdefault(ppid, []).append(int(p))
    seen, stack, total_kb = set(), [root_pid], 0
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        stack.extend(children.get(pid, []))
        try:
            for line in _read(os.path.join(proc_root, str(pid), "status")).splitlines():
                if line.startswith("VmRSS:"):
                    total_kb += int(line.split()[1])
                    break
        except (OSError, ValueError, IndexError):
            continue
    return total_kb / 1024.0


class Profiler:
    def __init__(
        self, interval: float = 2.0, proc_root: str = "/proc", root_pid: Optional[int] = None
    ):
        self.seconds: dict[str, float] = {}
        self.children_peak_after: dict[str, float] = {}
        self.tree_peak_by_stage: dict[str, float] = {}
        self.tree_peak: Optional[float] = None
        self.tree_peak_stage = ""
        self._active = ""
        self._interval, self._proc_root = interval, proc_root
        self._root = root_pid if root_pid is not None else os.getpid()
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._t0 = time.monotonic()
        self.sampled = False

    def start(self) -> Profiler:
        """Start the tree sampler (only where there is a process table)."""
        if tree_rss_mib(self._root, self._proc_root) is not None:
            self.sampled = True
            self._sample()
            self._thread = threading.Thread(target=self._loop, daemon=True)
            self._thread.start()
        return self

    def _sample(self):
        try:
            v = tree_rss_mib(self._root, self._proc_root)
        except Exception:  # noqa: BLE001  a sampler must never break a fetch
            return
        if v is None:
            return
        stage = self._active or "other"
        if v > self.tree_peak_by_stage.get(stage, 0.0):
            self.tree_peak_by_stage[stage] = v
        if self.tree_peak is None or v > self.tree_peak:
            self.tree_peak, self.tree_peak_stage = v, stage

    def _loop(self):
        while not self._stop.wait(self._interval):
            self._sample()

    @contextmanager
    def stage(self, name: str):
        """Time `name`; a stage entered inside another stage is not timed again."""
        if self._active:
            yield
            return
        self._active, t0 = name, time.monotonic()
        try:
            yield
        finally:
            self.seconds[name] = self.seconds.get(name, 0.0) + time.monotonic() - t0
            self._active = ""
            self._sample()
            self.children_peak_after[name] = _rusage_mib(_children())

    def summary(self) -> dict:
        total = time.monotonic() - self._t0
        named = sum(self.seconds.values())
        return {
            "seconds": {k: round(v, 2) for k, v in self.seconds.items()},
            "other_seconds": round(max(0.0, total - named), 2),
            "python_peak_mib": round(_rusage_mib(_self()), 1),
            "child_peak_mib": round(_rusage_mib(_children()), 1),
            "child_peak_after_stage_mib": {
                k: round(v, 1) for k, v in self.children_peak_after.items()
            },
            "tree_peak_mib": None if self.tree_peak is None else round(self.tree_peak, 1),
            "tree_peak_stage": self.tree_peak_stage or None,
            "tree_peak_by_stage_mib": {k: round(v, 1) for k, v in self.tree_peak_by_stage.items()},
            "tree_sampled": self.sampled,
        }

    def stop(self) -> dict:
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=5)
        self._sample()
        return self.summary()


def _self():
    import resource

    return resource.RUSAGE_SELF


def _children():
    import resource

    return resource.RUSAGE_CHILDREN


class TimedBackend:
    """Wraps a fetch backend so that its archive and tool calls are timed by stage, with no change to what they do."""

    STAGE: ClassVar[dict] = {
        "transfer": "transfer",
        "ena_head": "transfer",
        "convert_sra": "convert",
        "window": "convert",
        "read_stats": "count",
        "md5": "md5",
        "names": "check",
    }

    def __init__(self, backend, prof: Profiler):
        self._b, self._p = backend, prof

    def __getattr__(self, name):
        attr = getattr(self._b, name)
        stage = self.STAGE.get(name)
        if not stage or not callable(attr):
            return attr

        def timed(*a, **k):
            with self._p.stage(stage):
                return attr(*a, **k)

        return timed


def profile_row(res, summary: dict, threads: Optional[int] = None) -> list:
    """One TSV line for `--profile`: the fetch's outcome, its stage seconds and its memory."""
    s = summary["seconds"]
    cells = {
        "run": res.run,
        "tier": res.tier,
        "sampling": res.sampling,
        "route": res.route,
        "outcome": res.outcome,
        "host": socket.gethostname(),
        "threads": "" if threads is None else threads,
        "seconds": f"{res.seconds:.1f}",
        "other_s": f"{summary['other_seconds']:.1f}",
        "python_peak_mib": summary["python_peak_mib"],
        "child_peak_mib": summary["child_peak_mib"],
        "tree_peak_mib": "" if summary["tree_peak_mib"] is None else summary["tree_peak_mib"],
        "tree_peak_stage": summary["tree_peak_stage"] or "",
    }
    for st in STAGES:
        cells[f"{st}_s"] = f"{s.get(st, 0.0):.1f}"
    return [str(cells[f]) for f in PROFILE_FIELDS]


def append_profile(path: str, res, summary: dict, threads: Optional[int] = None) -> None:
    """Append one line (header first when the file is new). Several fetches may append at once: one write per line."""
    line = "\t".join(profile_row(res, summary, threads)) + "\n"
    fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o664)
    try:
        if os.fstat(fd).st_size == 0:
            line = "\t".join(PROFILE_FIELDS) + "\n" + line
        os.write(fd, line.encode())
    finally:
        os.close(fd)
