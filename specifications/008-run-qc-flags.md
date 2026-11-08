# 008 — Run quality flags: a registry of contamination and other per-run findings

Status: **draft**

## Why

Projects that read the same public runs find the same problems in them: a run that is mostly another genus, a run
that is a mixture, a run whose reads do not match the species in its metadata. Each project found these on its own
and kept them in its own tables. A second project then downloads the same bad run and finds the problem again.

This spec gives those findings one place in the catalog, keyed like everything else, so any project can see them
before it downloads, and so a curator can review them once.

It follows the rule of spec 001: **the catalog records the finding and the consumer decides.** A flag never removes a
run from the catalog or from a panel. A consumer that wants flagged runs can still have them.

## Files

```
catalog/qc/sources.csv                                   one row per contributing source (the registry)
catalog/qc/<species_taxid>.<Genus>_<species>.qc.csv.gz   one row per run, source and check
```

The per-species files follow the catalog naming rule of `AGENTS.md` section 7 (taxid first, species taxid from
`catalog/taxonomy_species_map.csv.gz`). A run whose taxid does not resolve to a species goes to
`catalog/qc/unresolved.qc.csv.gz`, never silently dropped.

### `sources.csv`

| column | meaning |
|---|---|
| `source_id` | stable id of one producer run of one method, for example `matpredict_reads_screen_2026-11-08` |
| `project` | the project that produced it (MATPredict, GmG, TE_fungi_copy, ...) |
| `tool` | the tool and version that made the call |
| `check` | what was checked; `contamination` for now |
| `method` | one line on how the call is made |
| `detection_limit` | what the method can and cannot see, with its source |
| `reference` | where the method and its calibration are documented (repository path or URL) |
| `date` | when the rows were produced (`YYYY-MM-DD`) |
| `rows` | rows of this source in `catalog/qc/` (written by the importer) |
| `notes` | anything a reader needs, for example known false-positive patterns |

### `<species>.qc.csv.gz`

| column | meaning |
|---|---|
| `run_accession` | ENA/SRA run |
| `sample_accession` | BioSample, from the PopSeq run table |
| `tax_id` | taxid as deposited for the run |
| `species_taxid`, `species_name` | from `taxonomy_species_map.csv.gz` |
| `source_id` | a row of `sources.csv` |
| `check` | `contamination` |
| `status` | the source's own status word (MATPredict screen: `mixed`, `genus_mismatch`, `borderline`) |
| `detected_genus` | the genus the method assigned to the reads |
| `second_taxon` | the second taxon it saw, if any |
| `second_hashes` | its evidence count (method-specific units) |
| `second_share` | its share of the evidence, 0-1 |
| `deeper_check` | `yes` when the source checked again with more reads |
| `subsample` | what was examined, for example `2000000 spots, 8 windows` |
| `strain` | as the source recorded it (a convenience for review; not a key) |
| `source_note` | anything else the source attached |
| `curated_status` | **curator only; never overwritten by an import** (`confirmed`, `false_positive`, `relabel`, ...) |
| `curated_note` | curator only; never overwritten |

An empty field means unknown (`AGENTS.md` section 8). Status words are the source's own; they are not merged across
sources, because the methods differ in what they can see.

## Rules

1. **A source replaces only its own rows.** Importing `source_id` X removes the earlier rows of X and writes the new
   ones; rows of other sources are untouched. Two sources that flag the same run give two rows.
2. **`curated_*` columns survive every import.** They are carried over by `(run_accession, source_id, check)`.
3. **Only findings are required.** A source may also record `ok` rows, but need not; absence of a row means the run
   was not flagged by any source, not that it was checked and clean. `sources.csv` says what each source examined.
4. **A metadata error found this way is fixed in the catalog, not here.** A run that is the wrong species gets a
   curated correction in the species record; the QC row stays as the evidence.

## Importing

`scripts/popseq_qc_import.py` reads a producer table (CSV or TSV, `.gz`/`.zst` read transparently), resolves each run
through the PopSeq kingdom run table and the taxonomy map, and writes the per-species files and the `sources.csv` row.
The column mapping for each producer format is in the script (`FORMATS`); `matpredict_screen` is the first.

## Changes

- 2026-11-08: first draft, with the MATPredict reads screen (423 flagged runs of 38 species) as the first source.
