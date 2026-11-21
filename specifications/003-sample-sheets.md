# Spec 003: sample sheets

Status: **draft**

The two files that describe a panel. Any pipeline consuming PopSeq data reads these and
nothing else. They are joined on the `sample` column.

## `samples.csv`

Deliberately identical to the nf-core sample sheet convention, so nf-core pipelines consume
it without a conversion step.

| Column | Required | Definition |
|---|---|---|
| `sample` | yes | the sample identifier. Stable, unique within the panel, and the join key to `metadata.csv`. |
| `fastq_1` | yes | absolute path to mate 1, or to the single-end file |
| `fastq_2` | yes (column), may be empty (value) | absolute path to mate 2. Empty for single-end. |

Rules:

- **One line per run.** A sample sequenced in three runs has three lines with the same
  `sample` value. The consumer decides whether to merge them; PopSeq does not merge.
- **Paths are absolute.** A relative path breaks as soon as a pipeline changes directory.
- The `fastq_2` column is always present even when every row is single-end. A missing column
  and an empty value are different failures, and tools report the first one badly.
- Sample identifiers use letters, digits, `_`, and `-`. No spaces, no `/`, no `.` as the
  first character. A `.` inside is allowed but discouraged, because many tools split on it.

## `metadata.csv`

One line per **sample**, not per run. Deliberately wide: it carries everything public that
was recoverable, so that a question nobody anticipated can still be asked.

### Mandatory columns

| Column | Definition |
|---|---|
| `sample` | join key, matches `samples.csv` |
| `biosample` | NCBI BioSample accession, `SAMN…`/`SAMEA…`/`SAMD…` |
| `taxid` | NCBI taxonomy id. **Mandatory.** A panel without taxids cannot be regenerated. |
| `scientific_name` | as recorded at NCBI, not corrected |

### `read_state`

Whether the panel's files are the reads as sequenced or the output of a trimming step.
A panel that mixes raw and trimmed reads gives results that cannot be compared across
samples, and without this column nothing in the panel reveals it (issue #5).

| Value | Meaning |
|---|---|
| `trimmed` | every file of the sample is the output of a trimming tool |
| `raw` | the reads are known to be as sequenced. Set only from a record that says so, never inferred from a file name |
| `unknown` | anything else. The default |

Rules:

- **Inferred from the file names only.** A file is trimmed when its basename contains
  `_val_` (trim_galore paired output, `<x>_1_val_1.fq.gz`), `_trimmed` (trim_galore
  single-end and others) or `.trim.`. Matching ignores case. A directory name such as
  `trimmed/` does not count.
- **Never guess `raw`.** A file name with no trimming marker is `unknown`, because the
  absence of a marker is not evidence that nothing was done to the reads.
- **Per sample.** `metadata.csv` has one line per sample, so the value is per sample. A
  sample is `trimmed` only when every one of its files is; a sample whose files mix
  trimmed and unmarked names, or that has no files, is `unknown`. A panel that mixes raw
  and trimmed samples is not collapsed to one value: the column shows the mix sample by
  sample, and `popseq_inventory.py` prints the per-panel counts (`read_state_trimmed=`,
  `read_state_unknown=`).
- Written by `scripts/popseq_inventory.py`, the converter for local subprojects (spec
  005), on every panel it emits. Panels built from archive URLs by
  `scripts/popseq_panel_set.py` do not carry the column yet; a consumer treats a missing
  column as `unknown`.

### Harvested columns

Everything else harvested from the BioSample record keeps NCBI's `harmonized_name` where one
exists, and the raw `attribute_name` otherwise. Common examples: `strain`, `isolate`, `host`,
`geo_loc_name`, `collection_date`, `lat_lon`, `isolation_source`, `study`, `center_name`.

These are **replaced** on a catalog refresh. Do not hand-edit them.

### Curated columns

Any column whose name begins `curated_` is ours. It is entered by hand and **never
overwritten by a refresh**. Use these to record what we believe when the public record is
wrong or incomplete:

- `curated_clade`, `curated_formae_speciales`, `curated_mating_type`
- `curated_exclude` — a non-empty value means a consumer should skip this sample
- `curated_notes` — free text, including why a sample was excluded

The `curated_exclude` column exists because deleting a row loses the information that the
row was ever considered. Exclusion is a finding, not an absence.

### Empty values

**An empty field means unknown. Never write `NA`.** `NA` is a valid ISO country code and a
real strain name. Once `NA` is in the table, the ambiguity cannot be undone.

## `panel.yaml`

```yaml
name: bd_global_334
title: Global Batrachochytrium dendrobatidis isolates
taxid: 109871
created: 2026-11-06
created_by: jstajich
question: |
  What this panel was assembled to answer. One paragraph. A panel without a
  stated question gets reused for a question it does not fit.
selection:
  source_catalog_refresh: 2026-11-06
  criteria: paired Illumina WGS, >= 15x against 24 Mb, one run per biosample
  n_samples: 334
  n_runs: 341
  n_excluded: 7
  exclusion_reasons: see curated_exclude in metadata.csv
depth_tier: full
```

## Changes

- 2026-11-06: first draft.
- 2026-11-21: `read_state` column (`raw`, `trimmed`, `unknown`), inferred from file names
  by the inventory converter (issue #5).
