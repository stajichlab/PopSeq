# Spec 005: local and unarchived read data

Status: **draft.** Written 2026-11-06 from a measured inventory of
`/bigdata/stajichlab/shared/projects/Population_Genomics`.

Public archives are not the only source. The lab holds read data that was never deposited,
was deposited under an accession nobody recorded, or was deposited and then the local copy
drifted from it. The catalog has to describe that data too, or every project keeps
rediscovering it.

## 1. What is actually there

Measured 2026-11-06 by walking each subproject's `config.txt` and its declared `SAMPFILE`.
29 configured subprojects, plus 5 directories with no `config.txt` at all.

| | Count |
|---|---|
| Strain rows across all sample sheets | **3,350** |
| Rows whose `FileBase` is an SRA run accession | **685** (20%) |
| Rows resolved to a local file by prefix match | 795 |

The gap is not one problem. It is four, and they need different responses.

| Cause | Example | Rows | Response |
|---|---|---|---|
| Paths contain subdirectories | `Afumigatus_Global`: `AF100.run1/lane1/AF100-2-A_S2_L001_R[12]_001.fastq.gz` | 1,152 | the resolver must walk, not list. Data is present: 1,224 FASTQ, 4.5 GB. |
| `input/` is empty, rows have accessions | `Coccidioides/C_posadasii`: 62 rows, 62 accessions, 0 files | 62 | refetchable. Catalog it as archived, cache on demand. |
| `input/` is empty, rows have **no** accession | `Rhodotorula_mucilaginosa`: 424 rows, 121 accessions, 0 files | ~303 | **this is the real problem.** See section 4. |
| Sample sheet lists far more than exists | `Rhizopus_stolonifer`: 113 rows, 2 FASTQ | 111 | a sheet is a plan, not an inventory. See section 5. |

## 2. The existing convention is good and should be kept

`PopGenomics_template` already defines one: `config.txt` declares `SAMPFILE`,
`FASTQFOLDER`, `FASTQEXT`, and `POPYAML`, and the sheet is `Strain,Filebase`. That is a real
convention and most subprojects follow it. PopSeq should **read** it, not replace it.

What it does not survive is that the convention was applied inconsistently.

### Header case varies

`Strain,Filebase` (template), `Strain,FileBase`, `STRAIN,FILEBASE`, and one file with no
header at all (`Afumigatus_Global`). A reader must match case-insensitively and must handle
a missing header row.

### `FileBase` means at least four different things

| Meaning | Example | Where |
|---|---|---|
| an SRA run accession | `SRR1636746` | `B_dendrobatidis` |
| a glob with a mate wildcard | `100D-1_S34_L001_R[12]_001.fastq.gz` | `Rhodotorula_mucilaginosa` |
| a bare prefix, no extension | `22687-1_ACTTGA_L001` | `O_ophiodiicola` |
| several relative paths, `;` separated, one per lane | `AF100.run1/lane1/...;AF100.run1/lane2/...` | `Afumigatus_Global` |

A resolver must handle all four. It must not guess: each form is unambiguous once the
`FASTQFOLDER` and `FASTQEXT` from `config.txt` are applied.

### A second schema exists

`Coccidioides/C_immitis` uses `RunAcc,Strain,BioSample,Center,Experiment,Project,Organism`.
That is **richer** than the template, and closer to what spec 001 wants. It should be
recognised rather than normalised away.

## 3. Provenance classes

Every sample in the catalog gets one. This is the field that makes local data describable.

| `provenance` | Meaning | Can be refetched? | Identity is |
|---|---|---|---|
| `archived` | in SRA/ENA, accession known | yes | the run accession |
| `archived_local` | in SRA/ENA, and a local copy exists | yes | the run accession, plus a local path |
| `local_submitted` | deposited, but the accession was never recorded here | **unknown** | a local identifier, until the accession is recovered |
| `local_only` | never deposited | **no** | a local identifier and a checksum |
| `listed_absent` | named in a sample sheet, no file found | no | nothing yet |

`local_only` data is **irreplaceable**. It deserves different handling from everything else
in this project, and the catalog should make that visible rather than letting it sit in an
`input/` directory indistinguishable from a cache.

## 4. Identity without an accession

An archived run has a run accession, and its integrity is checked against the ENA
`read_count` (spec 002). Local data has neither.

So a local sample is identified by:

```
local_id        <subproject>/<strain>          stable, human-readable, already unique
file_sha256     per file                        the only integrity check available
read_count      counted once, at registration   what a later check compares against
bytes, mtime    recorded at registration
```

**The checksum is computed once, when the file is registered, and never recomputed silently.**
A later mismatch is a finding: either the file changed or the record is wrong. Both matter.

This is the local equivalent of the ENA count check, and it exists for the same reason: there
is otherwise no way to tell a truncated file from a complete one.

### Recovering a lost accession

For `local_submitted`, the accession is often recoverable without re-deposition:

1. Search NCBI BioSample for the strain name within the expected organism.
2. Compare `read_count` of the candidate run against the local file's counted reads.
3. A match on both strain and read count is strong. **A match on strain alone is not**, since
   strain names repeat across studies and submitters.

A recovered accession is recorded with how it was recovered, in `curated_notes`. An
unverified guess is not recorded at all.

## 5. A sample sheet is a plan, not an inventory

`Rhizopus_stolonifer` lists 113 strains and has 2 FASTQ files. `Rhodotorula_mucilaginosa`
lists 424 and has none. These sheets are not wrong; they record an intent, or they record a
state that has since changed.

**PopSeq must never infer that data exists because a sheet names it.** The inventory is
produced by looking at the filesystem, and the sheet is one input to it. A row that does not
resolve becomes `listed_absent`, which is a recorded state, not an error and not a deletion.

This matters most in the opposite direction: a sheet that lists 424 strains with no files is
how a project silently loses data without anyone noticing, because the sheet still reads as
if the data is there.

## 6. What the inventory tool does

`scripts/popseq_inventory.py` (to be written):

1. For each subproject, read `config.txt` for `SAMPFILE`, `FASTQFOLDER`, `FASTQEXT`.
2. Parse the sheet, case-insensitively, tolerating a missing header and the second schema.
3. Resolve each `FileBase` by all four forms in section 2, **walking** `FASTQFOLDER`.
4. Classify each row into a provenance class from section 3.
5. Write `panels/<subproject>/samples.csv` and `metadata.csv` per spec 003, with absolute
   paths, plus a `local_manifest.csv` carrying `local_id`, size, mtime and checksum.
6. Never move, rename, or delete anything. The subprojects stay exactly as they are.

Point 6 is a hard rule. These are other people's working directories.

## 7. What this does not solve

- **Whether the local copy and the archived run are the same data.** For `archived_local`,
  nobody has compared them. A local file that was trimmed, subsampled, or re-paired is not
  the archived run, and treating them as interchangeable is how an analysis becomes
  irreproducible. Comparing read counts is cheap and should be part of registration.
- **The 5 subprojects with no `config.txt`.** They need a human to say what they are.
- **Which of the ~303 unaccessioned *Rhodotorula* rows ever had data.** The sheet records
  strain names; the files are gone. Whether they are elsewhere, archived, or lost is unknown,
  and this spec cannot determine it.

## Changes

- 2026-11-06: first draft, written from a measured inventory rather than from the template.
