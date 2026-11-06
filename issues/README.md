# PopSeq build issues

Twenty issues in six tracks. Each file states its track, what it depends on, and whether
it is safe to run in parallel with the rest of its track.

These are filed as GitHub issues on `stajichlab/PopSeq`. The files here are the
source text; the issue is the place to discuss and claim work.

Filter by track with the `track:A` to `track:F` labels. `parallel-safe` marks an issue
that can run alongside others in its track. `blocked` marks one with a dependency.

## Tracks

| Track | Theme | Can start now | Blocked by |
|---|---|---|---|
| **A** | catalog correctness | [A1](https://github.com/stajichlab/PopSeq/issues/1), [A2](https://github.com/stajichlab/PopSeq/issues/2), [A3](https://github.com/stajichlab/PopSeq/issues/3) | [A4](https://github.com/stajichlab/PopSeq/issues/4) waits on [A3](https://github.com/stajichlab/PopSeq/issues/3) |
| **B** | local and unarchived data | [B1](https://github.com/stajichlab/PopSeq/issues/5), [B2](https://github.com/stajichlab/PopSeq/issues/6), [B4](https://github.com/stajichlab/PopSeq/issues/8) | [B3](https://github.com/stajichlab/PopSeq/issues/7) waits on [B2](https://github.com/stajichlab/PopSeq/issues/6) |
| **C** | summaries and graphics | [C1](https://github.com/stajichlab/PopSeq/issues/9) | [C2](https://github.com/stajichlab/PopSeq/issues/10), [C3](https://github.com/stajichlab/PopSeq/issues/11) chain after |
| **D** | fetching | [D1](https://github.com/stajichlab/PopSeq/issues/12) | [D2](https://github.com/stajichlab/PopSeq/issues/13) waits on [D1](https://github.com/stajichlab/PopSeq/issues/12); [D3](https://github.com/stajichlab/PopSeq/issues/14) on [D2](https://github.com/stajichlab/PopSeq/issues/13) |
| **E** | other clades, metagenomes, MAGs | — | all wait on [A1](https://github.com/stajichlab/PopSeq/issues/1) |
| **F** | infrastructure | [F1](https://github.com/stajichlab/PopSeq/issues/18) running, [F2](https://github.com/stajichlab/PopSeq/issues/19) | [F3](https://github.com/stajichlab/PopSeq/issues/20) waits on [A4](https://github.com/stajichlab/PopSeq/issues/4) |

## Dependency graph

Issue numbers in brackets. Links do not render inside a code block, so they are plain here.

```
A1 [1] ──┬─> E1 [15]
         └─> E2 [16] ──> E3 [17]

A3 [3] ────> A4  [4] ──> F3 [20]

A2 [2]     independent, highest value
B2 [6] ────> B3  [7]
C1 [9] ────> C2 [10] ──> C3 [11]
D1 [12] ───> D2 [13] ──> D3 [14]
B1 [5]     independent
B4 [8]     independent
F1 [18]    running
F2 [19]    independent
```

## What to do first

**[A2](https://github.com/stajichlab/PopSeq/issues/2) is the most valuable single issue.** Every isolate count in the catalog is currently
a floor, because an ENA `tax_tree` query cannot see a run ENA holds with no taxid
assigned. For *Cercospora beticola* that is 752 against roughly 1,644. Until [A2](https://github.com/stajichlab/PopSeq/issues/2) is done,
every number the catalog reports is an undercount of unknown size.

**Four tracks can run in parallel from today:** A ([A1](https://github.com/stajichlab/PopSeq/issues/1), [A2](https://github.com/stajichlab/PopSeq/issues/2), A3), B ([B1](https://github.com/stajichlab/PopSeq/issues/5), [B2](https://github.com/stajichlab/PopSeq/issues/6), B4),
C (C1), D (D1). They touch different files and different services.

**[D1](https://github.com/stajichlab/PopSeq/issues/12) unblocks the whole of track D** and is the only remaining question in spec 002. It
needs cluster time, not design time.

**E waits on [A1](https://github.com/stajichlab/PopSeq/issues/1).** Generalising the clade parameter before the species key is fixed would
carry the rank-artifact bug into every new clade.

## Issues

| | Title | Track | Parallel |
|---|---|---|---|
| [A1](https://github.com/stajichlab/PopSeq/issues/1) | Key the catalog on resolved taxid and rank | A | yes |
| [A2](https://github.com/stajichlab/PopSeq/issues/2) | Close the ENA tax_tree blind spot | A | yes |
| [A3](https://github.com/stajichlab/PopSeq/issues/3) | sample_accession is multi-valued | A | yes |
| [A4](https://github.com/stajichlab/PopSeq/issues/4) | Split the BioSample table into core and long | A | yes |
| [B1](https://github.com/stajichlab/PopSeq/issues/5) | Record whether reads are raw or trimmed | B | yes |
| [B2](https://github.com/stajichlab/PopSeq/issues/6) | Verify archived_local copies against the archive | B | yes |
| [B3](https://github.com/stajichlab/PopSeq/issues/7) | Recover lost accessions for deposited local data | B | yes |
| [B4](https://github.com/stajichlab/PopSeq/issues/8) | Decide what the five no-config subprojects are | B | yes |
| [C1](https://github.com/stajichlab/PopSeq/issues/9) | Summary tables from the kingdom catalog | C | yes |
| [C2](https://github.com/stajichlab/PopSeq/issues/10) | Graphics from the summary tables | C | yes |
| [C3](https://github.com/stajichlab/PopSeq/issues/11) | Publish the dataset summary as a page | C | yes |
| [D1](https://github.com/stajichlab/PopSeq/issues/12) | Benchmark 7.1, compare download routes | D | yes |
| [D2](https://github.com/stajichlab/PopSeq/issues/13) | Implement popseq fetch with the four tiers | D | no |
| [D3](https://github.com/stajichlab/PopSeq/issues/14) | Record every download including failures | D | yes |
| [E1](https://github.com/stajichlab/PopSeq/issues/15) | Rename to projects/PopSeq | E | no |
| [E2](https://github.com/stajichlab/PopSeq/issues/16) | Add the metagenome entity, metadata only | E | yes |
| [E3](https://github.com/stajichlab/PopSeq/issues/17) | Add the MAG entity | E | yes |
| [F1](https://github.com/stajichlab/PopSeq/issues/18) | Finish the 114-species reconciliation | F | no |
| [F2](https://github.com/stajichlab/PopSeq/issues/19) | Tests, and a no-network smoke test | F | yes |
| [F3](https://github.com/stajichlab/PopSeq/issues/20) | Decide the DuckDB question with a measurement | F | yes |
