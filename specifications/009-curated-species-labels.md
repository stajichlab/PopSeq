# 009 — Curated species labels for genome assemblies

Status: **draft**

## Why

NCBI species names lag the taxonomy that the expert databases use. For the 2,111 fusarioid assemblies in BFD, in-silico
barcoding against FusariumID (https://www.fusarium.org/page/ID, download of 2026-06-23) gives a label for 2,059, and
1,040 of them name a different species than NCBI (`catalog/labels/assembly_labels.csv.gz`, `label_agreement`). Most
are *Fusarium oxysporum* assemblies that FusariumID places in one of the species of the *F. oxysporum* species complex
(FOSC). Every consumer that groups by the NCBI name inherits the old name. Issue #33 proposed this layer; the owner
agreed on 2026-11-09:

- PopSeq holds assembly-level metadata, because an assembly is easier to match than a run.
- FusariumID is an acceptable authority for the fusarioid genera. It is refreshed 2-3 times a year.
- Other genera can follow from other authorities, for example NCBI RefSeq ITS, LSU and SSU for *Aspergillus*.

The layer sits **beside** the NCBI taxonomy and never replaces it. `species_taxid` and `ncbi_species` stay as NCBI has
them; the authority's view is in `label_*` columns, and a curator's view in `curated_*` columns.

## Files

```
catalog/kingdom/fungi_assemblies.csv.gz   one row per fungal GenBank assembly (NCBI assembly summary)
catalog/labels/sources.csv                one row per label source (the registry)
catalog/labels/assembly_labels.csv.gz     one row per assembly and source
catalog/labels/run_labels.csv.gz          the labels carried to runs of the same BioSample
```

### `fungi_assemblies.csv.gz`

Built by `scripts/popseq_assemblies.py build` from a dated snapshot of NCBI's `assembly_summary.txt` for the fungi
group, GenBank and RefSeq (`popseq_assemblies.py fetch` writes
`cache/ncbi_assembly_summary_fungi_{genbank,refseq}_<date>.txt.gz` and a `.query.txt` with URL, bytes and sha256).
One row per GCA accession; a RefSeq copy adds `refseq_accession` and `refseq_category`. NCBI's `na` becomes an empty
field. `species_name` is the NCBI taxonomy name of `species_taxid` from `catalog/taxonomy_species_map.csv.gz`, empty
when the map does not hold it. A suppressed assembly is not in NCBI's current summary and so not in the table.

Snapshot 2026-11-09: 26,826 assemblies, 26,819 with a BioSample, 669 with a RefSeq copy.

### `sources.csv`

| column | meaning |
|---|---|
| `source_id` | stable id of one import, for example `fusariumid_gmg_f1_2026-06-23` (authority, producer, authority release) |
| `authority` | the curated reference the label comes from (FusariumID, ...) |
| `authority_release` | the release or download date of the authority's data |
| `project`, `tool` | who produced the calls, and with what |
| `method` | one line on how the call is made |
| `call_rule` | the rule that turns the evidence into a label (written by the importer) |
| `reference` | where the method and its validation are documented |
| `date` | when the rows were imported |
| `refresh` | how often the source is to be refreshed |
| `assemblies`, `runs` | rows of this source in `assembly_labels` and `run_labels` (written by the importer) |
| `notes` | what a reader needs, for example where the labels are weak |

### `assembly_labels.csv.gz`

| column | meaning |
|---|---|
| `assembly_accession` | GCA accession (a producer's GCF is resolved to its GCA) |
| `sample_accession`, `tax_id`, `species_taxid` | from `fungi_assemblies.csv.gz`; empty when the assembly is not in it |
| `ncbi_species` | NCBI taxonomy name of `species_taxid` (else the deposited organism name) |
| `source_id` | a row of `sources.csv` |
| `label_species` | the authority's species; `;`-joined for a tie; empty for no call |
| `label_complex` | the complex all names of `label_species` share; empty when they differ |
| `label_status` | `confident`, `tie` or `no_call` |
| `label_agreement` | `same`, `tie_includes_ncbi`, `differs`, `ncbi_unnamed` (NCBI name is sp./cf.) or `no_call` |
| `nearest_species` | the nearest authority species even when beyond the call rule (evidence, not a label) |
| `distance`, `differences`, `aligned_bp`, `margin`, `runner_up`, `loci` | the evidence, in the producer's units |
| `curated_species`, `curated_note` | **curator only; never overwritten by an import** |

### `run_labels.csv.gz`

`run_accession`, `sample_accession`, `assembly_accession`, `source_id`, `ncbi_species`, `label_species`,
`label_complex`, `label_status`, `label_agreement`, `curated_species`. A run gets a row when its BioSample (any
value of a `;`-joined `sample_accession`) is the BioSample of a labelled assembly. It is rebuilt from every source at
each import, so it always matches `assembly_labels.csv.gz`.

## Rules

1. **A label never replaces the NCBI fields.** Consumers choose: `curated_species`, else `label_species`, else the
   NCBI species.
2. **A source replaces only its own rows.** Importing `source_id` X removes the earlier rows of X; rows of other
   sources are untouched. Two authorities that label the same assembly give two rows.
3. **`curated_*` columns survive every import,** carried over by `(assembly_accession, source_id)`.
4. **Ties stay explicit.** A tie keeps every name. When the tied species share a complex, `label_complex` names it,
   so a consumer can fall back to the complex.
5. **An unresolved assembly is kept,** with empty NCBI fields, never dropped.
6. **Evidence stays with the label.** The distance, margin and runner-up let a reader judge a call without
   re-running the producer.

## Sources

| source | authority | assemblies | calls | refresh |
|---|---|---|---|---|
| `fusariumid_gmg_f1_2026-06-23` | FusariumID download 2026-06-23 | 2,111 BFD fusarioid | 1,916 confident, 143 tie, 52 no call | 2-3 times a year |

**FusariumID by in-silico barcoding (GmG bench F1).** The FusariumID loci (tef1, rpb1, rpb2, CaM, tub2) are found
in each assembly by megablast. Distance is (mismatches + gap openings) / aligned length over the concatenated loci to
the nearest FusariumID strain. A label is given within distance 0.01. The method and its checks are in
`github.com/stajichlab/GmG` `bench/F1_fusariumid/README.md`. Known weakness: FOSC species are 99.6-99.99% genome ANI
apart and 86% of FOSC strains are within 2 tef1 differences of another FOSC species, so a FOSC label is reliable at the
complex level and weaker at species level (same source). The complexes come from GmG
`crates/gmg-core/data/species_complexes.tsv`, made from the same FusariumID download.

### Refreshing FusariumID (2-3 times a year)

1. Download the FusariumID sequence set (https://www.fusarium.org/page/ID).
2. In GmG, rerun bench F1's genome barcoding (`bench/F1_fusariumid/run_bridge.sbatch`, then `genome_resolution.py`)
   and `scripts/make_species_complexes.py`.
3. Import under a new `source_id` with the new release date. Keep the old source until a curator has compared the
   two, then remove it by importing an empty table under the old id, or by deleting its rows in a reviewed commit.

UNMEASURED: how many labels change between FusariumID releases.

## Next authorities

- **NCBI RefSeq Targeted Loci (ITS, LSU, SSU)** for genera with no expert database, *Aspergillus* first (owner,
  issue #33). ITS rarely separates *Aspergillus* species, so the expected grain is section or series, not species.
  UNMEASURED: the resolution of each locus per genus; measure it before importing any label.

## Changes

- 2026-11-09: first draft, with the assembly table and the FusariumID labels from GmG bench F1 as the first source.
