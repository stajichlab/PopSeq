# Changelog

## 2026-11-09 (genome sizes for species with 10 or more isolates; issue #36)

No path moved, no column added, renamed or removed.

- **`catalog/genome_sizes.csv`: 115 to 459 species** (`popseq_genome_sizes.py
  --min-isolates 10`). 418 have an NCBI reference size; 41 have none and keep an empty
  size with the reason "no reference assembly at NCBI" (hybrids, unnamed `sp.` taxa and
  species with no assembly). All 41 were fetched a second time with the same result.
  The 115 existing rows are unchanged.
- **`est_depth_x` filled for 12,109 runs in 308 species** in
  `catalog/runs/fungi_all_runs.csv.gz`. Every changed cell was blank before. No other
  cell changed.
- `catalog/summary/coverage.csv.gz` grows from 111 to 418 species. The existing rows are
  unchanged. `catalog/species_layout.csv.gz` gains single-only depth values for 62
  species.
- One study to know about: *Hortaea werneckii* PRJNA428320 has 54 runs labelled PAIRED
  that ENA serves as one file, at a median 914x against a 25.33 Mb assembly. Its paired
  runs and the other studies sit at 85x to 290x. The size is not the problem.
- **New panel `panels/pop_sra_ge10_cap300/`** (`scripts/popseq_panel_set.py`), the input
  for GmG's per-species PCoA and PCA pages. It holds up to 300 isolates for every species
  with at least 10 paired Illumina isolates (owner decision 2026-11-09). Each species is
  selected by `popseq.panel(n=300)` with the default filters. The panel has 412 species
  and 34,673 isolates. 47 species reach the cap of 300; 177 have 10 to 29 isolates.
  - 1,888 isolates in 36 species have no depth estimate, because those species have no
    genome size. 549 isolates are below 10x.
  - Fetching 25x per isolate needs about 33.8 Tbases, against 135.3 Tbases for whole runs.
  - `samples.csv` holds ENA URLs. Nothing has been fetched (#13).
  - It was built from catalog commit 7b08547, which `panel.yaml` records.

## 2026-11-09 (later)

Curated species labels for assemblies (issue #33; owner answers of 2026-11-09). New paths only: no path moved, no
column of an existing table changed.

- **New spec `specifications/009-curated-species-labels.md`.** A label layer beside the NCBI taxonomy: an
  authority's species (`label_*` columns) and a curator's (`curated_*`, never overwritten), with the evidence.
- **New table `catalog/kingdom/fungi_assemblies.csv.gz`** (`scripts/popseq_assemblies.py`): 26,826 fungal GenBank
  assemblies from NCBI's assembly summary, snapshot 2026-11-09; 26,819 have a BioSample, 669 a RefSeq copy. This is
  the link from an assembly to its BioSample and runs.
- **New `catalog/labels/`** (`scripts/popseq_labels_import.py`), first source `fusariumid_gmg_f1_2026-06-23`:
  FusariumID labels for the 2,111 BFD fusarioid assemblies from GmG bench F1 in-silico barcoding. 1,916 confident,
  143 ties, 52 no call; 1,040 name a different species than NCBI (mostly *F. oxysporum* to a FOSC species).
  `run_labels.csv.gz` carries them to 1,000 runs that share a BioSample with 848 of the assemblies.
- `catalog/index.json` lists the four new tables. Tests: `tests/test_labels_import.py`.

## 2026-11-09

BGISEQ and DNBSEQ runs (owner decision 2026-11-08). Rows added only: no path moved, no
column added, renamed or removed.

- **10,267 runs added** to `catalog/runs/fungi_all_runs.csv.gz` (4,210 BGISEQ, 6,057
  DNBSEQ), from a supplementary snapshot, `cache/fungi_runs_ena_bgi_2026-11-09.tsv.gz`, made
  by the new `scripts/popseq_snapshot.py`. The table now has 256,089 rows. The 245,822
  existing rows are identical in every column but one: 78 `est_depth_x` cells for
  *Kluyveromyces marxianus* were blank and now have a value, because that species reached
  100 isolates and got a genome size.
- **134 new species**, 6,620 in the species summary. 10,094 new BioSamples harvested
  (`fungi_biosamples.csv.gz`); `host_taxa.csv.gz` rebuilt for their hosts. Genome sizes
  added for *K. marxianus* and *Laburnicola rhizohalophila*, the two species that now reach
  100 isolates, and `catalog/api/` gains their files.
- **New runs are `not_checked` against NCBI.** `popseq_runs_table.py` now takes the
  reconciled set from the snapshot the reconciliation used (`--reconciled-ena`), not from
  the species summary.
- **The default `--platform` filter stays ILLUMINA.** Panels and `catalog/api/` panels do not
  include BGISEQ or DNBSEQ unless asked. `index.json` lists both snapshots in
  `archive_snapshots`.
- **New derived table `catalog/species_layout.csv.gz`** (`scripts/popseq_species_layout.py`):
  paired versus single-end-only isolates per species, by `ena_layout`. 69 species have only
  single-end isolates (182 isolates). 161 species gain isolates from single-end-only data.
  4,776 of the 6,480 single-only isolates are estimated at 10x or more.
- **`scripts/popseq_species_summary.py`** builds `fungi_species_summary.csv.gz`, which had no
  script before. Run on the 2026-11-06 snapshot, it reproduces the committed file byte for
  byte.
- **Spot check:** 4 BGISEQ/DNBSEQ runs over AWS and ENA HTTPS
  (`benchmarks/2026-11-09_bgi_spotcheck.csv`). Results matched the served-layout prediction:
  labelled PAIRED, both routes `ok`; labelled SINGLE but served as `_1`/`_2`, NCBI holds them
  unpaired (`unavailable_at_ncbi`) and ENA delivered verified pairs.

## 2026-11-08

Pairing is what ENA serves (owner decision 2026-11-08). One column added; nothing moved.

- `catalog/runs/fungi_all_runs.csv.gz` gains `ena_layout`: `PAIRED` when ENA serves `_1`
  and `_2`, `SINGLE` for one file, empty for none. The 34 existing columns of all 245,822
  rows are identical to the previous build. `library_layout` disagrees with it for 2,274
  runs: 806 labelled PAIRED serve one file, 1,468 labelled SINGLE serve pairs.
- `--layout` (CLI, package, MCP) and the default panels now filter on `ena_layout`. 30
  species files in `catalog/api/` change: 643 more paired isolates in total. Largest moves:
  *Cryptococcus neoformans* +607, *Xylaria necrophora* 0 to 166, *C. gattii* +105,
  *Schizosaccharomyces pombe* -242, *Saccharomyces cerevisiae* -51.
- Specs 001 and 002 state the pairing rule, single-end verification and `--require-pairs`.
  `read_count` is no longer described as a verification target. `popseq/downloads.py`
  gains the outcome `unpaired`.

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
