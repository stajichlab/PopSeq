# D1: Benchmark 7.1: compare download routes head to head

**Track:** D   **Depends on:** nothing   **Parallel-safe:** yes

**This blocks the fetcher.** Spec 002's route order is a proposal, not a finding.

What exists:
- `prefetch` 20-23 MB/s. **One run.** TE_fungi_copy, 2026-10-24.
- ENA FTP over curl 0.4 MB/s. **One run.** Same source.
- AWS Open Data 10.3 MB/s median. **4,211 runs, 6,343 GB.** MATPredict, well evidenced.
- ENA Aspera. **Never measured here.** Module `aspera/4.2.13.820` is installed.

**Do.** 10 runs spanning 0.5 to 10 GB, each fetched by all four routes, on a compute node
and on a login node. Record MB/s, wall time and failures per route. Also find the
concurrency ceiling: GmG saw `parallel-fastq-dump` exit 64 at about 35 concurrent, but
**concurrency was not shown to be the cause** and the safe value is unknown.

**Acceptance.** Spec 002 moves from `blocked` to `stable` with a measured route order,
and a concurrency cap with evidence behind it.
