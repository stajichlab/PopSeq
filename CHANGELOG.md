# Changelog

## 2026-11-08

Repository rules.

- `AGENTS.md` section 7 (owner-approved): `catalog/runs/fungi_all_runs.csv.gz` is the
  canonical run table. The rule against a combined file is replaced by a rule that every
  derived table in `catalog/` must be reproducible by a named script. No file moved.

Sample sheet and catalog scope.

- `popseq runs --format samples` and `popseq.to_sample_sheet` now choose mates by file name
  (`_1`, `_2`), not by position. Before this, 15,105 runs that ENA lists with three files
  (an unpaired `<RUN>.fastq.gz` plus `_1` and `_2`) could get the unpaired file as
  `fastq_1`. The unpaired file is left out of the sheet. It stays in `fastq_ftp`, and the
  CLI prints a stderr caveat with the number of runs affected. New:
  `popseq.query.fastq_mates`, `popseq.query.sample_sheet_dropped`.
- `catalog/index.json` and `catalog/api/index.json` gain a caveat: the run table is
  WGS GENOMIC Illumina only (issue #24).

Per-run quality flags (specification 008, draft).

- `specifications/008-run-qc-flags.md`: a registry of per-run findings (contamination for now) from any project. The catalog records the finding; the consumer decides. A source replaces only its own rows; `curated_*` columns survive every import.
- `catalog/qc/sources.csv` and `catalog/qc/<species_taxid>.<Genus>_<species>.qc.csv.gz`: first source `matpredict_reads_screen_2026-11-08`, 423 flagged runs (236 mixed, 143 borderline, 44 genus_mismatch) in 39 species files; 423 of 423 runs resolved through `cache/fungi_runs_ena_2026-11-06.tsv.gz` and `catalog/taxonomy_species_map.csv.gz`.
- `scripts/popseq_qc_import.py` and `tests/test_qc_import.py`.

Benchmark 7.1 and the verification rule.

- **Issue #12 finished on compute nodes.** Spec 002 moves from blocked to stable. The route
  order for `full` is AWS Open Data, then `prefetch`, then ENA HTTPS. The rates were
  measured while other lab download arrays were running, and spec 002 says so. The
  login-node arm and the concurrency cap stay UNMEASURED (issue #27).
- `AGENTS.md` section 6 and `PLAN.md` section 6.3 now state the verification rule adopted on
  2026-11-07 (verify against the archive that served the bytes) instead of the ENA
  `read_count` check.

## 2026-11-07

Added columns only. No path moved, no column renamed or removed.

- **Issue #23.** `catalog/runs/fungi_all_runs.csv.gz` gains `country_iso`, `country_name`,
  `country_status`, `host`, `host_taxid`, `host_name`, `host_status`. The 27 existing
  columns of all 245,822 rows are byte-identical to the previous build (compared with
  pandas after a rebuild from `cache/fungi_runs_ena_2026-11-06.tsv.gz`). New directory
  `catalog/vocab/` holds the curated and generated vocabularies. `popseq runs` and
  `popseq panel` take `--country CH` or `--country Switzerland`, and a new `--host` takes a
  taxid or a name. The MCP tools take the same arguments. See spec 001.
- **Issue #14.** `popseq/downloads.py` appends one validated row per fetch attempt to
  `benchmarks/downloads.csv`, failures included. It is safe for concurrent writers on GPFS.
  It refuses an `ok` row that was not verified.
- **Issue #12.** `scripts/popseq_bench_routes.py` and `scripts/bench_routes.slurm` run
  benchmark 7.1. Its results are in `benchmarks/2026-11-07_route_comparison.csv` and
  `benchmarks/downloads.csv`. The run plan is in `benchmarks/2026-11-07_d1_runs.csv`.
- **Spec 002 verification rule changed** (adopted by the owner): verify against the archive
  that served the bytes. Use NCBI `spots_with_mates` for SRA downloads, and ENA's per-file MD5
  plus equal mate counts for ENA downloads. ENA `read_count` counts reads, not pairs, for 172
  of 200 sampled ERR runs (`benchmarks/2026-11-07_read_count_convention.csv`).
- ENA Aspera and `prefetch -t fasp` both fail on the cluster. See spec 002, route 2.
- New dependency: `pycountry`.

## 2026-11-06

First commit. Repository skeleton, plan, rules, and three draft specifications.

- `PLAN.md`: what is being built, why, and the measurements it depends on.
- `AGENTS.md` (with `CLAUDE.md` symlinked to it): the rules for working here.
- `specifications/001-catalog-schema.md`, `002-download-strategy.md`,
  `003-sample-sheets.md`: draft interface contracts.
- `scripts/popseq_biosamples.py`: batched NCBI BioSample harvester. Tested against NCBI on
  2026-11-06; fetched `SAMN12327137` and `SAMEA5704822` and reproduced the documented output
  of the script it replaces.
- `pyproject.toml` and `uv.lock`: uv-managed environment, Python >= 3.11.

Specs 001 and 002 were revised the same day after reading `MATpredict_spec_handoff.md` and
`GmG_spec_handoff.md`, which the MATPredict and GmG sessions wrote independently. That
revision added the AWS Open Data download route, the `windows` and `stride` subsample tiers,
the rule against using `head` for any cross-sample measurement, the download concurrency cap,
the `unavailable_at_ncbi` run state, and the study-level exclusion categories.

Nothing has been harvested into `catalog/` yet.
