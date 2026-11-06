# Spec 001: catalog schema

Status: **draft**

The catalog is the source of truth. It is text, it is tracked in git, and it is regenerable
from ENA and NCBI. `catalog/popseq.duckdb` is a derived index built from these files and is
never edited.

## Why the taxid leads every filename

```
catalog/species/109871.Batrachochytrium_dendrobatidis.yaml
catalog/runs/109871.Batrachochytrium_dendrobatidis.runs.csv.gz
catalog/biosamples/109871.Batrachochytrium_dendrobatidis.biosample.csv.gz
```

Names change. Taxids do not. *Candida auris* became *Candidozyma auris*; *Botryotinia
fuckeliana* became *Botrytis cinerea*. A catalog keyed on the name forks into two records at
every rename, and the fork is silent. A catalog keyed on the taxid keeps one record and adds
a `synonyms:` entry.

When NCBI renames a species, the file is renamed and the old name goes into `synonyms`. The
taxid prefix does not change, so `git mv` shows it as a rename rather than a delete plus add.

## `catalog/species/<taxid>.<Genus>_<species>.yaml`

```yaml
taxid: 109871
scientific_name: Batrachochytrium dendrobatidis
rank: species
lineage: [Fungi, Chytridiomycota, Chytridiomycetes, Rhizophydiales]
synonyms: []
common_names: [Bd, chytrid fungus]

# Needed for every coverage calculation. Source it and say where from.
assumed_genome_mb: 24
assumed_genome_mb_source: GCA_000149865.1 (JEL423)

reference_assemblies:
  - accession: GCA_000149865.1
    strain: JEL423
    note: most widely used reference
  - accession: GCA_003595275.1
    strain: RTP6
    note: single-copy baseline donor in TE_fungi_copy

counts:
  last_refreshed: 2026-11-06
  ena_query: 'tax_tree(109871) AND library_strategy="WGS" AND library_source="GENOMIC"'
  runs_total: 0
  runs_paired_illumina_wgs: 0
  isolates_distinct: 0
  isolates_above_15x: 0
  estimated_download_tb: 0.0

projects_using:
  - TE_fungi_copy

notes: |
  Free text. Who asked for this species, what is known to be wrong in the
  public records, and anything a future reader needs in order not to repeat
  a mistake.
```

`counts` is **derived** and is rewritten by `popseq refresh-species`. Everything above it is
curated and is never overwritten by a refresh.

## `catalog/runs/<taxid>.<name>.runs.csv.gz`

One line per **run**. Harvested from the ENA Portal API, which returns every run for a clade
in one request. Columns are ENA's own field names, kept unchanged so the mapping back to the
API is obvious.

| Column | Source | Notes |
|---|---|---|
| `run_accession` | ENA | primary key, `SRR`/`ERR`/`DRR` |
| `sample_accession` | ENA | the BioSample. **Runs are not isolates; collapse on this.** |
| `secondary_sample_accession` | ENA | `SRS`/`ERS` |
| `study_accession` | ENA | `PRJ…` |
| `scientific_name` | ENA | as recorded, not corrected |
| `tax_id` | ENA | may differ from the catalog taxid for a subspecies |
| `instrument_platform` | ENA | `ILLUMINA`, `OXFORD_NANOPORE`, `PACBIO_SMRT` |
| `instrument_model` | ENA | |
| `library_layout` | ENA | `PAIRED` or `SINGLE` |
| `library_strategy` | ENA | `WGS`, `RNA-Seq`, `AMPLICON`, … |
| `library_source` | ENA | `GENOMIC`, `TRANSCRIPTOMIC`, … |
| `read_count` | ENA | read **pairs** for paired runs. The verification target. |
| `base_count` | ENA | drives the size and coverage estimates |
| `fastq_ftp` | ENA | semicolon-separated |
| `fastq_aspera` | ENA | empty when Aspera is not offered |
| `fastq_md5` | ENA | ENA's own file checksums, not valid for converted FASTQ |
| `first_public` | ENA | |
| `center_name` | ENA | often an `SUB…` submission id rather than a lab |
| `est_coverage` | derived | `base_count / (assumed_genome_mb * 1e6)` |
| `est_gz_bytes` | derived | `base_count / 3.5`. **UNMEASURED divisor**, see spec 002. |
| `eligible` | derived | `true` when it passes the filter below |
| `eligible_reason` | derived | why not, when `eligible` is false |

### The eligibility filter

One place, one definition. Three projects wrote this separately and got three slightly
different answers.

**Run level:**

```
library_strategy == "WGS"
library_source   == "GENOMIC"
instrument_platform == "ILLUMINA"
library_layout   == "PAIRED"
read_count       >= 600000
est_coverage     >= 15
two fastq files present
```

**Study level.** A run can pass every test above and still be wrong for a population
question. MATPredict and GmG both filter on the study, and both arrived at the same four
categories. They are recorded per study in an `excluded_as` column, not deleted:

| `excluded_as` | What it catches |
|---|---|
| `mutant` | knockouts, engineered strains |
| `evolution` | experimental evolution, serial passage, mutation accumulation |
| `cross` | progeny, segregants, tetrads, recombinant inbred lines |
| `not_wgs` | RNA-seq, ChIP, ATAC, bisulfite, Hi-C, metagenome, amplicon, ITS/16S, single-cell |

This is a keyword filter on study and sample titles. **It has false positives and misses**,
and MATPredict says so explicitly. A study whose title names natural or clinical sampling
(`surveillance`, `outbreak`, `GWAS`, `natural isolates`, `bloodstream`) is not excluded on a
stray keyword.

Note that an exclusion is **not** always correct for every question. TE_fungi_copy
deliberately *wants* mutant and laboratory-evolved panels, because a transposon burst in a
mutant background is the thing it is looking for. So the catalog records the category and
the consumer decides. A filter that silently deletes rows cannot be reversed by a consumer
who wants them.

Defaults, overridable per call. The thresholds are recorded in `panel.yaml` for any panel
built from them, because a panel built under a different threshold is a different panel.

### Clonal sets are not populations

*Candidozyma auris* has about 25,000 public runs, most of them outbreak surveillance of a
largely clonal population. Counting them as 25,000 independent strains is wrong for any
frequency estimate. MATPredict flags this; TE_fungi_copy hit the same problem and solved it
by collapsing calls against a distance matrix before counting events.

The catalog records the run count. **It does not claim they are independent.** Any panel that
needs independence states how it established it, in `panel.yaml`.

### Runs are not isolates

The most common error with this table is counting rows. A re-sequenced isolate has several
runs; a mutant panel has one isolate per run. **Collapse on `sample_accession` before
reporting a strain count.** When one run per isolate is needed, take the smallest run that
still passes the coverage filter, so the download is as small as the question allows.

## `catalog/biosamples/<taxid>.<name>.biosample.csv.gz`

One line per **BioSample**. Harvested from NCBI Entrez, which carries far richer sample
attributes than ENA does.

| Column | Notes |
|---|---|
| `biosample` | primary key, `SAMN`/`SAMEA`/`SAMD` |
| `taxid` | |
| `scientific_name` | |
| everything else | one column per BioSample attribute, using NCBI's `harmonized_name` where one exists and the raw `attribute_name` otherwise |

The table is **wide and sparse**, because BioSample attributes are not standardised across
submitters. That is accepted deliberately: a column present for 3 percent of samples still
answers a question for those samples, and dropping it loses information that cannot be
recovered without re-harvesting.

Harmonized names are preferred because they collapse `geographic location (country and sea)`,
`geo_loc_name`, and `country` onto one column. This is the single most valuable thing NCBI
does to this metadata.

**An empty field means unknown. Never write `NA`.**

## Regeneration

```bash
uv run popseq refresh-runs       --taxid 109871
uv run popseq refresh-biosamples --taxid 109871
uv run popseq refresh-species    --taxid 109871
```

Refreshes are incremental: a record already present is not re-queried unless `--update` is
passed. The diff is reviewed before it is committed.

**A refresh that removes rows is reviewed by a human.** Public records are withdrawn rarely.
A drop almost always means the query changed, not the archive.

## Changes

- 2026-11-06: first draft.
- 2026-11-06: added the study-level exclusion categories and the note on clonal sets, after
  reading `MATpredict_spec_handoff.md` and `GmG_spec_handoff.md`. Both projects had
  independently written the same four exclusion categories.
