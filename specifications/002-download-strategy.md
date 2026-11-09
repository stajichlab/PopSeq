# Spec 002: download strategy

Status: **stable** (2026-11-08). Benchmark 7.1 settled the route order for tier `full` on
compute nodes. See "Benchmark 7.1 results". Two measurements are still UNMEASURED and are
tracked in issue #27: the login-node arm and the concurrency cap.

## Depth tiers

A caller states what it needs. The fetcher never guesses.

| Tier | What it gets | Typical size | Cached? | Destination |
|---|---|---|---|---|
| `head` | first N read pairs, default 500,000 | ~100 MB | no | `$SCRATCH` |
| `windows` | N spots drawn from W evenly spaced windows, default 2,000,000 in 8 | ~1 GB | no | `$SCRATCH` |
| `stride` | every Nth read pair | ~1 Gb of bases | no | `$SCRATCH` |
| `full` | the entire run | 0.5-15 GB | yes | `cache/<taxid>/<RUN>/` |

Only `full` is cached. The subsampled tiers cost less to refetch than to index and manage.

### Do not use `head` for anything population-scale

The three projects each chose a different subsample, and they are not equivalent:

| Project | Strategy | Stated reason |
|---|---|---|
| GmG `reads_sra` | first 500,000 pairs over an ENA range request | the only way to stop an HTTPS stream early |
| MATPredict | 2,000,000 spots in 8 evenly spaced windows, `fastq-dump -N <start> -X <end>` | "spreading windows avoids a head-of-file bias" |
| GmG `popsketch` | every Nth pair via `parallel-fastq-dump` | "old SRA files are position-sorted" |

Two independent projects concluded that the head of a run is not a random sample of it.
Reads in a FASTQ are ordered by flowcell position, so the first N pairs come from one
region of one tile, and anything that varies across the flowcell varies with them.

**Rule: `head` is only for a question that is insensitive to read origin**, such as asking
which organism a run contains. For any measurement compared across samples, use `windows`
or `stride`. The cost difference is small; the bias is not.

`windows` is the default subsample. It is the one with a stated reason to prefer it over
`stride`, and `stride` requires reading the whole run to produce its output, so it saves
storage but not transfer.

## Route order

### Tier `full`

0. **AWS Open Data over HTTPS**, resumable:
   ```
   curl -L -C - https://sra-pub-run-odp.s3.amazonaws.com/sra/<RUN>/<RUN> -o <RUN>.sra
   ```
   This is MATPredict's primary route and it is not in the TE_fungi_copy measurement at
   all. It is a plain resumable HTTPS GET of the `.sra`, so a broken transfer continues
   instead of restarting.

   **Measured: median 10.3 MB/s per task** (p10 3.2, p90 20.1, max 31.1), over **4,211 runs
   and 6,343 GB** moved, from `<run>.run.json` records. This is the best-evidenced transfer
   figure anyone in the lab has: the competing numbers for `prefetch` (20-23 MB/s) and ENA
   FTP (0.4 MB/s) each rest on a single run.

   Read it with MATPredict's own limits: it is one task on a compute node, not an aggregate;
   `secs_download` includes retry sleeps and any failed AWS attempt, so the true rate is a
   little higher; and it is `.sra` bytes, not FASTQ, which are several times larger.
   Concurrency was throttled by a SLURM array (for example `%4`) and was not recorded per
   run, so the rate's behaviour under load is unknown.

   **A 404 here is close to final. Do not retry it, and do not expect `prefetch` to rescue
   it.** Of 26 runs that missed on AWS and fell back to `prefetch`, **`prefetch` rescued
   none**: 25 were absent from NCBI and 1 failed to download. The fall-through to route 1
   below is therefore cheap to attempt once and should never be looped.

1. **NCBI `prefetch` then `parallel-fastq-dump`.**
   ```
   prefetch <RUN> -O <dir on /bigdata> --max-size 30G
   parallel-fastq-dump --sra-id <dir>/<RUN>/<RUN>.sra --threads 8 \
       --outdir <fastq dir> --split-files --gzip --tmpdir "${SCRATCH:?}"
   ```
   Measured on the UCR HPCC login node on 2026-10-24: about 20-23 MB/s. Verified end to end
   with `ERR3401300` (1.6 GB, 16,860,539 pairs): `prefetch` about 80 s, conversion about 3 min.

2. **ENA Aspera (`ascp`)**, when the ENA record carries a `fastq_aspera` field.
   **Not usable on this cluster as of 2026-11-07.** `ascp` fails authentication against
   `era-fasp@fasp.sra.ebi.ac.uk` with all four keys installed under `aspera/*`
   (`asperaweb_id_dsa.openssh` from 4.1.3.93 and 4.2.4.265, `aspera_id_rsa.openssh`,
   `aspera.openssh`, `aspera_tokenauth_id_rsa`, `NCBI.aspera.openssh`), with `ascp` 4.2.4 and
   4.2.13. The NCBI equivalent is gone too: `prefetch -t fasp SRR36290343` answers "cannot
   download using requested transport". Benchmark 7.1 records 13 Aspera attempts, all
   `transfer_failed` in about 1.4 s.
   `module load aspera/4.2.13.820`. **UNMEASURED on this cluster.** Aspera is often faster
   than FTP over long distances because it does not use TCP congestion control. That is a
   general claim about the protocol, not a measurement here, and must not be quoted as one.

3. **ENA HTTPS** (`https://ftp.sra.ebi.ac.uk/...`).

4. **`fasterq-dump`**, last resort.

### Tiers `head` and `sample`

**ENA HTTPS first**, because it is the only route that can stop early. A range request on a
gzipped FASTQ can be abandoned once N pairs are decompressed. `prefetch` has no partial mode:
asking it for the first 500,000 pairs costs the whole run.

Fall back to `fasterq-dump` with a read limit.

## The unresolved conflict

| Source | Measurement |
|---|---|
| `TE_fungi_copy/CLAUDE.md`, 2026-10-24, login node | `prefetch` 20-23 MB/s; ENA FTP over `curl` about 0.4 MB/s; one ENA connection refused |
| `GmG/reports/validation/reads_sra.md` | ENA HTTPS delivered 100 of 100 runs, 9.57 GB, `fastq-dump` fallback never needed |

These may not actually conflict. GmG pulls small HTTP byte ranges and closes the connection
once enough lines are decompressed, which is a different access pattern from streaming a whole
run. But **nobody has compared them directly**, and until someone does, the default above is
a proposal carried over from one project, not a finding.

Benchmark 7.1: 10 runs spanning 0.5 to 10 GB, each fetched by all three routes, on a compute
node and on the login node. Record MB/s, wall time, and failures.

### Benchmark 7.1 results (2026-11-07 to 2026-11-08)

10 runs from 0.5 to 10 GB (`benchmarks/2026-11-07_d1_runs.csv`), plus the control run
`ERR2188722`. Each was fetched by all four routes on compute nodes. Outcomes are judged under
the adopted verification rule below (`popseq_bench_routes.py summary --plan`). Source:
`benchmarks/2026-11-07_route_comparison.csv`.

| Route | Attempts | ok | verify_failed | transfer_failed | unavailable_at_ncbi | Transfer MB/s median (p10 to p90) |
|---|---|---|---|---|---|---|
| `aws_odp` | 12 | 10 | 0 | 0 | 2 | 10.3 (5.0 to 21.5) |
| `ncbi_prefetch` | 12 | 10 | 0 | 0 | 2 | 9.5 (6.1 to 16.2) |
| `ena_https` | 13 | 11 | 2 | 0 | 0 | 0.90 (0.84 to 0.94) |
| `ena_aspera` | 13 | 0 | 0 | 13 | 0 | none: authentication fails |

**Measured under load.** Every attempt overlapped about 10,000 MATPredict download jobs and
the sketch-store download arrays. Read these rates as the rates under normal lab load, not as
a baseline. They will be measured again if later rows in `benchmarks/downloads.csv` show they
vary. The attempt counts include one duplicate run (`ERR9802648`) from a cancelled first
array.

- The two `unavailable_at_ncbi` per SRA route are `ERR2188722` (absent at NCBI) and
  `ERR2286963` (held unpaired at NCBI). ENA served both as verified pairs.
- The two ENA `verify_failed` attempts (`SRR16988686`, `SRR25902591`) arrived at exactly
  ENA's byte size but failed MD5 and gunzip. Both ran on node i32. Both were verified `ok`
  when fetched again on node i29. The cause is unknown.
- The two SRA routes are equal within the spread. AWS stays first because it is a plain
  resumable HTTPS GET and MATPredict's 4,211-run record agrees (10.3 MB/s median).
- ENA HTTPS is about 10 times slower for whole runs. It stays the route for partial tiers and
  the source of pairs for runs NCBI holds unpaired.

This settles the order for tier `full`: AWS Open Data, then `prefetch`, then ENA HTTPS.
The TE_fungi_copy login-node figure for `prefetch` (20-23 MB/s, one run) was not reproduced
on compute nodes, and the login-node arm is UNMEASURED (#27).

## Verification is mandatory

There is no checksum for a converted FASTQ. The count check is the only protection.

1. Get the expected count once, at planning time, and store it:
   ```
   curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=<RUN>&result=read_run&fields=read_count,base_count,library_layout&format=tsv"
   ```
   Store it in the run table. **Compute nodes may have no internet**, so a job must check
   against a stored value, not a live query.
2. After conversion, count both mates. Each must equal `read_count`, which is read **pairs**.
   ```
   N=$(( $(zcat "<RUN>_1.fastq.gz" | wc -l) / 4 ))
   ```
3. A mismatch **fails the task**. It does not warn.
4. Delete the `.sra` only after the check passes.

### ENA `read_count` is not a pair count for most ERR runs (measured 2026-11-07)

Step 2 above assumes ENA `read_count` counts read pairs. **For most ERR runs it counts
reads.** Benchmark 7.1 found this first. A random sample then measured it: 200 paired
Illumina runs per accession prefix from `catalog/runs/fungi_all_runs.csv.gz`, comparing ENA
`read_count` with NCBI runinfo `spots`. Source:
`benchmarks/2026-11-07_read_count_convention.csv`, from `popseq_bench_routes.py countcheck`.

| Prefix | n | ENA = NCBI spots | ENA = 2 x NCBI spots | other | absent at NCBI |
|---|---|---|---|---|---|
| DRR | 200 | 200 | 0 | 0 | 0 |
| SRR | 200 | 199 | 0 | 0 | 1 |
| ERR | 200 | 22 | 172 | 4 | 2 |

Of the 172 ERR runs at 2x, 108 were submitted by the Wellcome Sanger Institute. In
benchmark 7.1, ENA's own FASTQ for `ERR17038970` held 9,139,801 pairs per mate against an
ENA `read_count` of 18,279,602, and its MD5s matched ENA's. So the halving is in the
metadata, not in the files.

Two further facts from the same sample bear on the route order:

- **27 of the 200 ERR runs (13.5%) are stored unpaired at NCBI** (`spots_with_mates = 0`),
  although ENA records them as `PAIRED`. An SRA-derived download of such a run cannot give
  properly split mates. For these runs ENA is the only source of pairs. 3 of 200 SRR runs
  show the same thing.
- **The archives can hold different reads for the same run.** For `ERR9802648` (submitted as
  CRAM), NCBI serves 5,946,824 pairs and ENA's FASTQ holds 5,933,626, with ENA MD5s
  matching. This is a further reason that both mates must come from one source.

**Rule, adopted 2026-11-07 by the owner:** verify against the archive that served the bytes.
This replaces steps 1 and 2 above for the choice of target; steps 3 and 4 stand.

| Route | Expected count | Integrity check |
|---|---|---|
| `aws_odp`, `ncbi_prefetch` | NCBI `spots_with_mates`, stored at plan time | each mate equals it; `spots_with_mates == spots`, or the run is not paired at NCBI and must come from ENA |
| `ena_https`, `ena_aspera` | none from metadata | ENA `fastq_md5` per file (exact), and mate 1 count equals mate 2 count |

#### Pairing and single-end runs (owner decision 2026-11-08)

Pairing is what the serving archive holds, not the `library_layout` label. ENA serves pairs
when its file list has `_1` and `_2` (catalog column `ena_layout`). NCBI serves pairs when
`spots_with_mates > 0`.

| Route | Paired run | Single-end run |
|---|---|---|
| `aws_odp`, `ncbi_prefetch` | each mate = NCBI `spots_with_mates` | read count = NCBI `spots`; `spots_with_mates = 0` |
| `ena_https`, `ena_aspera` | ENA MD5 per file, mate 1 count = mate 2 count | ENA MD5 of the one file |

- **`--require-pairs`.** When NCBI holds a run unpaired, the fetcher skips the SRA routes
  and fetches from ENA. If ENA does not serve `_1` and `_2` either, the run ends with the
  `unpaired` outcome. It is not retried.
- **Without `--require-pairs`**, a run NCBI holds unpaired may be fetched from an SRA route
  as single-end. It is verified as single-end and recorded with the layout obtained.
- A run labelled `PAIRED` that both archives serve as one file is single-end. Do not try to
  split it.
- Record the layout obtained in the `<RUN>.json` sidecar.

`benchmarks/downloads.csv` is an append-only log, so rows written before the rule was
adopted keep the outcome they were given then. Some of those `verify_failed` rows for ERR runs
are intact downloads. `popseq_bench_routes.py summary --plan` judges every row of the detail
CSV again under the adopted rule, from the stored mate counts and MD5 results. A run NCBI
holds unpaired is recorded as `unavailable_at_ncbi`, with the reason in `error`.

This rule exists because it has already caught a real corruption. On 2026-10-24 a leftover
ENA loop resumed a download onto an SRA-derived file and corrupted mate 1. Nothing else
detected it.

## Rules that are not negotiable

- **Both mates of a run come from one source.** Never mix an ENA file with an SRA-derived
  file for the same run.
- When moving a run to another route, first stop every process that can still write to it
  (`pgrep -af "fetch|curl|prefetch|ascp"`), then delete the partial files, then start over.
  A looping downloader can wake up much later and resume onto a file another route replaced.
- Add `</dev/null` to any loop containing `prefetch`, `curl`, or `ascp`. They read stdin and
  will consume the loop's input.
- Load each module on its own line. A `module load ... | tail` pipe runs in a subshell and
  has no effect.
- Downloads go to `/bigdata`. Never to node-local `/scratch`, which is lost when the job ends.
- **Cap concurrent downloads.** GmG/popsketch saw `parallel-fastq-dump` exit 64 for 28 of
  2,023 runs on the first pass, while running about 35 concurrently. **Concurrency was not
  shown to be the cause.** GmG corrected this on 2026-11-06: 7 of the 28 were transient and
  all 7 succeeded on rerun at a cap of 3, which is consistent with a cap helping but does not
  demonstrate it. A cap is cheap insurance, not a proven fix. The safe value is UNMEASURED;
  20 is a starting point. Benchmark 7.1 did not measure it, because it needs a quiet cluster
  (issue #27).
- **Retry before reporting a failure**, but distinguish the two causes. Of those 28 runs,
  21 were **permanent**: `sra-stat` returns 404 for every one, they are all `ERR` accessions
  from four ENA studies, and they exist at ENA but not at NCBI. That is
  `unavailable_at_ncbi`, not a transient error, and retrying at NCBI will never succeed.

## The sidecar

Every cached run carries `cache/<taxid>/<RUN>/<RUN>.json`:

```json
{
  "run": "ERR3401300",
  "taxid": 109871,
  "biosample": "SAMEA5704822",
  "route": "ncbi_prefetch",
  "tier": "full",
  "fetched_utc": "2026-11-06T15:04:11Z",
  "bytes_total": 1683114928,
  "seconds_fetch": 81.4,
  "seconds_convert": 184.2,
  "mb_per_s": 20.7,
  "read_count_expected": 16860539,
  "read_count_mate1": 16860539,
  "read_count_mate2": 16860539,
  "verified": true,
  "tools": {"sratoolkit": "3.2.0", "parallel-fastq-dump": "0.6.7"}
}
```

**A cached run with no sidecar is not trusted** and is re-verified before use. An unverified
file that looks complete is the failure mode this whole section exists to prevent.

## Size estimation

Before anything is downloaded:

```
gzipped_fastq_bytes  ~=  base_count / 3.5
coverage             =   base_count / (assumed_genome_mb * 1e6)
```

The `/3.5` divisor is **inherited from the TE_fungi_copy SRA survey and has not been checked
here.** Benchmark 7.2 checks it against the 42 FASTQ files already on disk in that project.
Every storage estimate in `PLAN.md` depends on it, so it is checked before any of them is
quoted to anyone.

## Performance capture

Every attempt, including every failure, appends one row to `benchmarks/downloads.csv`:

```csv
utc,run,taxid,tier,route,node_type,bytes,seconds,mb_per_s,verified,outcome,error
```

`outcome` is one of:

| State | Meaning |
|---|---|
| `ok` | downloaded and verified against the serving archive |
| `verify_failed` | downloaded, and the check against the serving archive failed. The file is deleted. |
| `transfer_failed` | the transfer did not complete on any route |
| `unavailable_at_ncbi` | the run exists at ENA but neither AWS Open Data nor `prefetch` has it |
| `unpaired` | pairs were required (`--require-pairs`) and neither archive serves them. Not retried. |

`unavailable_at_ncbi` is a real and separate state, not a failure. Two projects measured it
independently:

- MATPredict: **25 of 4,548 runs, 0.55 percent, about 1 in 180** (example `ERR2188722`).
  This supersedes the "about 1 in 100" in their handoff; it is the larger count.
- GmG: **21 of 300** *P. oryzae* runs, all `ERR` accessions from four ENA studies
  (`PRJEB128687`, `PRJEB14371`, `PRJEB22035`, `PRJEB45710`), each confirmed by a `sra-stat`
  404. List at `popsra/p_oryzae_missing.txt`.

GmG's rate is far above MATPredict's because the state is **not evenly distributed**: it
clusters by submitting archive and therefore by study. A species whose data came mostly from
European submitters will have many; one from NCBI submitters will have none. So a pipeline
must not treat 1 percent as a budget. It must check per run.

**This is the strongest argument for ENA remaining in the route list** even if the benchmark
shows NCBI is faster. For these runs ENA is not a fallback, it is the only source. Whether
ENA actually delivers them is untested; GmG notes this explicitly.

A pipeline that treats this state as a failure will retry it forever. Record it and move on.

Failures are recorded because a route that works 95 percent of the time and a route that
works 100 percent of the time look identical in a table of successes.

## Changes

- 2026-11-06: first draft. Route order proposed, blocked on benchmark 7.1.
- 2026-11-06: revised after reading `MATpredict_spec_handoff.md` and `GmG_spec_handoff.md`.
  Added the AWS Open Data route (route 0), which the first draft missed entirely. Added the
  `windows` and `stride` tiers and the rule against `head` for cross-sample measurement.
  Added the concurrency cap and the `unavailable_at_ncbi` state.
- 2026-11-06: corrected the concurrency claim after GmG revised their handoff. The first
  version of this spec implied ~35 concurrent downloads caused the 28 failures. GmG showed
  21 of the 28 are permanently absent from NCBI and only 7 were transient, so concurrency
  is not demonstrated as the cause. Added their measured per-study clustering of
  `unavailable_at_ncbi`.
- 2026-11-06: added MATPredict's measured AWS Open Data transfer rate (10.3 MB/s median over
  4,211 runs), their corrected `unavailable_at_ncbi` rate of 0.55 percent, and the finding
  that a `prefetch` fallback after an AWS 404 rescued 0 of 26 runs.
- 2026-11-07: measured that ENA `read_count` counts reads, not pairs, for 172 of 200 sampled
  ERR runs, and that 27 of those 200 are unpaired at NCBI. Adopted the rule "verify against
  the archive that served the bytes". ENA Aspera could not be measured: every key on the
  cluster fails authentication.
- 2026-11-08: status moved from blocked to stable. Added the benchmark 7.1 results, labelled
  as measured under load. The login-node arm and the concurrency cap moved to issue #27 and
  stay UNMEASURED.
- 2026-11-08: added the pairing rule (the archives decide, not the label), single-end
  verification, and `--require-pairs`.
