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
| `library_layout` | ENA | `PAIRED` or `SINGLE`: the submitter's label. **Metadata only**; it does not decide pairing (see `ena_layout`) |
| `ena_layout` | derived | `PAIRED` when ENA serves files named `_1` and `_2`, `SINGLE` when it serves one file, empty when it serves none. **This decides pairing** |
| `library_strategy` | ENA | `WGS`, `RNA-Seq`, `AMPLICON`, … |
| `library_source` | ENA | `GENOMIC`, `TRANSCRIPTOMIC`, … |
| `read_count` | ENA | read pairs for SRR and DRR runs, but **reads** for most ERR runs (spec 002). **Not** a verification target |
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
ena_layout       == "PAIRED"     # what ENA serves, not the label
read_count       >= 600000
est_coverage     >= 15
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

### Pairing is what the archives serve (owner decision 2026-11-08)

`library_layout` is a label the submitter typed. In the catalog it disagrees with the files
ENA serves for 2,274 runs: 806 labelled `PAIRED` serve one file, and 1,468 labelled
`SINGLE` serve `_1` and `_2` (`scripts/popseq_runs_table.py`, build of 2026-11-08). In
BGISEQ/DNBSEQ data, not yet in the catalog, about 2,400 WGS runs labelled `PAIRED` serve
one file.

- `ena_layout` decides pairing in the catalog and in every `--layout` filter.
- NCBI can disagree with ENA. 27 of 200 sampled ERR runs are paired at ENA and unpaired at
  NCBI (`spots_with_mates = 0`). The catalog does not hold `spots_with_mates`. The fetcher
  reads it at plan time and, with `--require-pairs`, takes such a run from ENA (spec 002).
- A run with three files (`_1`, `_2` and an unpaired `<RUN>.fastq.gz`) is `PAIRED`. The
  unpaired file is kept in `fastq_ftp` and left out of sample sheets.

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

### Normalised country and host (issue #23)

`country` and the BioSample host are free text. *Zymoseptoria tritici* alone has 316
distinct country strings. A consumer must not re-parse them, so the kingdom run table
(`catalog/runs/fungi_all_runs.csv.gz`) carries normalised columns next to the raw ones.
**The raw `country` column is never changed.** It is the `country_raw` the issue asked for,
kept under its original name so no consumer breaks.

| Column | Notes |
|---|---|
| `country_iso` | ISO 3166-1 alpha-2 (`XK` for Kosovo). Blank unless `country_status` is `ok` |
| `country_name` | the ISO short name of `country_iso` |
| `country_status` | `ok`, `absent`, `missing`, `region`, `historical`, `conflicting`, `unresolved` |
| `host` | raw host: first non-empty of the BioSample `host`, `host_plant`, `specific hosts`, `host_common_name`. Not `lab_host` |
| `host_taxid`, `host_name` | NCBI taxid and its scientific name |
| `host_status` | `ok`, `absent`, `missing`, `not_a_host`, `same_as_organism`, `ambiguous`, `conflicting`, `unresolved` |

The statuses are what make "no value" explicit. `absent` means the field is empty: the
submitter wrote nothing. `missing` means the submitter wrote an INSDC missing term such as
`not collected`. That is a statement, not a gap, and it is kept separate. `not_a_host`
covers a lab or environment term written in the host field (`laboratory environment`,
`culture`, `soil`). `same_as_organism` covers a host field that names the sequenced
organism itself, which is usually a lab strain.

The rules live in `popseq/normalise.py`. Every curated decision is a reviewable row:

| File | Edited by | Holds |
|---|---|---|
| `catalog/vocab/country_aliases.csv` | hand | spellings pycountry does not know, with a status and a note |
| `catalog/vocab/host_overrides.csv` | hand | `not_a_host` terms, and common names mapped to a taxon *by name* (the builder resolves the taxid) |
| `catalog/vocab/host_taxa.csv.gz` | `scripts/popseq_vocab.py hosts` | every distinct host string, its taxid, and the method that resolved it |
| `catalog/vocab/country_unresolved.csv` | `scripts/popseq_vocab.py countries` | the curator's work list |

Coverage on the 2026-11-06 snapshot (245,822 runs; source: the `popseq_runs_table.py` log):

| | ok | absent | missing | other |
|---|---|---|---|---|
| country | 128,824 | 101,951 | 15,003 | region 32, historical 6, unresolved 5, conflicting 1 |
| host | 62,410 | 179,753 | 32 | not_a_host 1,665, unresolved 1,590, same_as_organism 215, ambiguous 138, conflicting 19 |

The 5 unresolved country runs say only `Korea`. They are left unresolved on purpose, because
the string does not say which Korea. Host resolution is exact-name only (scientific name
first, then synonyms and common names, then the binomial prefix or `<Genus> sp.`). It does
no fuzzy matching, so a misspelt host stays `unresolved` rather than becoming a wrong taxon.

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
- 2026-11-07: added the normalised country and host columns and their vocabularies
  (issue #23).
- 2026-11-08: added `ena_layout`, which decides pairing. `library_layout` is now metadata
  only. `read_count` is no longer described as a verification target. The eligibility
  filter uses `ena_layout == "PAIRED"` in place of the label and the two-file test.
