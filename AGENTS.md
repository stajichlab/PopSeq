# AGENTS.md — rules for the PopSeq repository

This file is the contract for anyone, human or agent, working in
`/bigdata/stajichlab/shared/projects/Fungi/PopSeq`.

`CLAUDE.md` is a symlink to this file. Edit this one.

Read `PLAN.md` before changing anything structural. This file says *how* to work here;
`PLAN.md` says *what* is being built and *why*.

## 1. This is a shared resource

Other projects depend on this repository's paths and file formats. A path that moves breaks
someone else's pipeline silently, usually while they are not watching.

- **Never move or rename a file under `catalog/`, `panels/`, or `cache/` without a deprecation
  period.** Add the new path, leave the old one, announce the change in `CHANGELOG.md`,
  remove the old one in a later commit.
- The group-writable setgid bit is set on these directories. Keep it. Run
  `chmod g+w` on anything you create that others must update.
- Assume a reader has no context. Every table gets a header. Every directory that is not
  self-evident gets a `README.md`.

## 2. Truth and evidence

These rules come from the lab's standing instructions and are not negotiable here.

- **Never report an inflated result.** State only what the data supports.
- **Make no assumptions.** When a decision has no data behind it, write that it has no data
  behind it, and mark the measurement that would settle it.
- Any number in a document says where it came from: a file in this repo, a named benchmark,
  or an external source with an accession. A number with no source is a bug.
- A figure inherited from another project is labelled as inherited and is re-measured before
  it is relied on. `PLAN.md` section 6.4 is the worked example.
- Use `UNMEASURED` as a literal marker in documents. It is greppable. `TODO` is not specific
  enough.

## 3. Writing style

- Short sentences. One idea per sentence. Active voice.
- Plain vocabulary. No idioms, no flourish.
- Prefer a table to a paragraph when the content is a comparison.

## 4. Python, with uv

`uv` is the only package manager used here. Do not use `pip install` into a shared
environment, and do not create a conda environment in this repository.

```bash
cd /bigdata/stajichlab/shared/projects/Fungi/PopSeq
uv sync                      # create or update .venv from pyproject.toml + uv.lock
uv run popseq --help         # run inside the environment without activating it
uv add <package>             # add a dependency; updates pyproject.toml and uv.lock
```

- `uv.lock` **is committed.** It is what makes a run reproducible.
- `.venv/` is **not** committed.
- Pin the Python version in `pyproject.toml` with `requires-python`. Do not rely on the
  system interpreter: the HPCC default `python3` is 3.9 and is too old for this code.
- In a SLURM script, call `uv run` with an absolute path to the project. Do not activate
  the venv and do not rely on the inherited `PATH`.

## 5. SLURM rules

- **Never use `$(dirname "${BASH_SOURCE[0]}")` in a script that runs under SLURM.** It
  resolves to the wrong path in a SLURM working directory. Pass the project root in an
  environment variable, `POPSEQ_ROOT`, or use an absolute path.
- `$SCRATCH` is **node-local** (`/scratch/$USER/<jobid>`), not shared. Use it for temporary
  files inside a job. Write `${SCRATCH:?}` so a missing variable fails loudly rather than
  writing to `/`.
- Copy anything that must survive to `/bigdata` before the job ends.
- Never hardcode a `/scratch/$USER/...` path.
- **Use `$SCRATCH` for temporary files, never `/tmp`.** This holds for interactive and
  login-node work too, not only batch jobs. `/tmp` on a compute node is small and shared
  with every other user of that node, so filling it breaks other people's jobs; `$SCRATCH`
  is per-job, node-local, faster, and cleaned up automatically.
- **Downloads go to `/bigdata`, never to `/scratch`.** A node-local download is lost.
- Size a scatter-gather unit toward **1 to 1.5 hours** of real runtime. Job count has its own
  cost in scheduling overhead and queue contention. If a completed run shows tasks finishing
  in minutes, bundle more work per job.
- Request what was measured, not what seems safe. An oversized request sits in `PENDING`.
  A 32-core, 96 GB, 1-day request for a job that used 3 GB and 13 minutes waited 24 hours.
- Load each module on its own line. A `module load ... | tail` pipe runs in a subshell and
  does nothing.
- Add `</dev/null` to any loop containing `prefetch`, `curl`, or `ascp`. They consume stdin
  and will eat the loop's input.

## 6. Downloading data

The full strategy is `PLAN.md` section 6 and `specifications/`. The rules that must never be
broken:

- **Verify every download against the archive that served it.** For an SRA-derived FASTQ,
  each mate's read count must equal NCBI `spots_with_mates`, stored at plan time. For an ENA
  FASTQ, each file must match ENA's MD5 and the two mates must hold equal counts. Fail the
  job on a mismatch. Do not use ENA `read_count` as the target: for most ERR runs it counts
  reads, not pairs (spec 002). There is no checksum for a converted FASTQ, so the count check
  is its only protection, and it has already caught a real corruption.
- Delete the `.sra` only **after** verification passes.
- **Take both mates of a run from one source.** Never mix an ENA file and an SRA-derived file
  for the same run.
- When switching a run to another route, stop every process that can still write to it
  (`pgrep -af "fetch|curl|prefetch|ascp"`), delete the partial files, then start over. A
  looping downloader can resume much later onto a file another route already replaced.
- Every cached run gets a `<RUN>.json` sidecar with route, bytes, seconds, read counts, the
  ENA values checked against, and tool versions. **A cached run with no sidecar is not
  trusted.** Re-verify it before use.
- Record every download in `benchmarks/downloads.csv`, including failures. The cache-versus-
  refetch decision depends on this record existing from the first download, not on someone
  starting to collect it later.

## 7. File naming and organisation

### Catalog

```
catalog/species/<taxid>.<Genus>_<species>.yaml
catalog/runs/<taxid>.<Genus>_<species>.runs.csv.gz
catalog/biosamples/<taxid>.<Genus>_<species>.biosample.csv.gz
```

- **The taxid comes first** and is the real key. Names change; taxids are stable. A species
  renamed at NCBI keeps its file and gets a `synonyms:` entry.
- Genus capitalised, species lower case, joined by `_`, matching the NCBI scientific name
  with spaces replaced. No strain, no subspecies in the filename.
- One file per species. Do not create a combined "all species" file in the catalog; it is a
  derived product and belongs in `popseq.duckdb`, which is not tracked.

### Panels

```
panels/<panel_name>/panel.yaml
panels/<panel_name>/samples.csv
panels/<panel_name>/metadata.csv
```

`<panel_name>` is lower case with underscores, and names the question, not the person:
`bd_global_334`, `cocci_popgen_554`, not `jason_panel2`.

### Cache

```
cache/<taxid>/<RUN>/<RUN>_1.fastq.gz
cache/<taxid>/<RUN>/<RUN>_2.fastq.gz
cache/<taxid>/<RUN>/<RUN>.json
```

Keyed by organism, not by project, so two projects asking for the same run get a cache hit.

### Everything else

- Scripts in `scripts/`, lower case with underscores, with the extension.
- A dated artefact uses `YYYY-MM-DD` first: `benchmarks/2026-11-06_route_comparison.csv`.
- SLURM logs to `logs/<jobname>.<jobid>.log`. `logs/` is not tracked.

## 8. Data formats

- **CSV is the default.** Comma separated, quoted where needed, header always present,
  `\n` line endings.
- **Compress anything large.** Use `.gz` for files other projects or collaborators read,
  because every bioinformatics tool reads gzip. Use `.zst` (`zstd -T0`) for internal
  intermediates only.
- **Write every new Python script to read compressed input transparently.** Accept `.csv`,
  `.csv.gz`, and `.csv.zst` for the same argument. Do not make the caller decompress first.
- If a tool genuinely cannot read compressed input, decompress into `$SCRATCH` for the life
  of the job, never into shared storage.
- **An empty field means unknown. Never write `NA`.** `NA` is a valid country code and a real
  strain name, and the ambiguity is unrecoverable once it is in the table.
- **`curated_*` columns are ours and are never overwritten by a refresh.** Everything else is
  harvested and may be replaced.

## 9. The catalog is the source of truth; DuckDB is not

`catalog/popseq.duckdb` is a **derived index**. It is regenerated from the YAML and CSV files
and is not tracked in git.

- Never edit the database directly.
- Never let a value exist only in the database.
- If regenerating the index from the catalog does not reproduce it exactly, the index is
  wrong, not the catalog.

This exists so that corrections are reviewable as diffs. A hand-correction to a species
record must be visible in a pull request.

## 10. Git

- Commit the catalog, panels, specs, scripts, benchmarks, and `uv.lock`.
- Never commit: `cache/`, `logs/`, `.venv/`, `catalog/popseq.duckdb`, any FASTQ, any `.sra`.
- A refresh commit says what changed in counts: *"Bd runs 326 to 331, 5 new from PRJNA1012345"*.
  A commit message of "update catalog" is not acceptable in a living database.
- **A refresh that removes rows is reviewed by a human before it is committed.** Public
  records are rarely withdrawn; a drop usually means the query changed.
- Commit messages end with:

```
Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
```

## 11. Secrets

- `NCBI_API_KEY` comes from the environment. It is never written into a file in this
  repository, never committed, and never printed in a log.
- The Entrez email is a required argument and may be committed. It is not a secret.

## 12. Before you claim something works

Run it and show the output. "Evidence before assertions" applies to every claim in a report,
a commit message, and a reply:

- A script "works" when it has run and its output has been looked at.
- A download "succeeded" when the read-count check passed, not when the command exited 0.
- A benchmark number is quoted with its n, and with the spread, not the mean alone.
