# PopSeq build issues

Twenty issues in six tracks. Each file states its track, what it depends on, and whether
it is safe to run in parallel with the rest of its track.

There is no git remote yet. When one exists these push with:

```bash
for f in issues/[A-F][0-9]*.md; do
  gh issue create --title "$(head -1 "$f" | sed 's/^# //')" --body-file "$f"
done
```

## Tracks

| Track | Theme | Can start now | Blocked by |
|---|---|---|---|
| **A** | catalog correctness | A1, A2, A3 | A4 waits on A3 |
| **B** | local and unarchived data | B1, B2, B4 | B3 waits on B2 |
| **C** | summaries and graphics | C1 | C2, C3 chain after |
| **D** | fetching | D1 | D2 waits on D1; D3 on D2 |
| **E** | other clades, metagenomes, MAGs | — | all wait on A1 |
| **F** | infrastructure | F1 running, F2 | F3 waits on A4 |

## Dependency graph

```
A1 ──┬─> E1
     ├─> E2 ──> E3
A3 ──> A4 ──> F3
A2   (independent, highest value)
B2 ──> B3
C1 ──> C2 ──> C3
D1 ──> D2 ──> D3
F1   (running)
F2   (independent)
```

## What to do first

**A2 is the most valuable single issue.** Every isolate count in the catalog is currently
a floor, because an ENA `tax_tree` query cannot see a run ENA holds with no taxid
assigned. For *Cercospora beticola* that is 752 against roughly 1,644. Until A2 is done,
every number the catalog reports is an undercount of unknown size.

**Four tracks can run in parallel from today:** A (A1, A2, A3), B (B1, B2, B4),
C (C1), D (D1). They touch different files and different services.

**D1 unblocks the whole of track D** and is the only remaining question in spec 002. It
needs cluster time, not design time.

**E waits on A1.** Generalising the clade parameter before the species key is fixed would
carry the rank-artifact bug into every new clade.

## Issues

| | Title | Track | Parallel |
|---|---|---|---|
| A1 | Key the catalog on resolved taxid and rank | A | yes |
| A2 | Close the ENA tax_tree blind spot | A | yes |
| A3 | sample_accession is multi-valued | A | yes |
| A4 | Split the BioSample table into core and long | A | yes |
| B1 | Record whether reads are raw or trimmed | B | yes |
| B2 | Verify archived_local copies against the archive | B | yes |
| B3 | Recover lost accessions for deposited local data | B | yes |
| B4 | Decide what the five no-config subprojects are | B | yes |
| C1 | Summary tables from the kingdom catalog | C | yes |
| C2 | Graphics from the summary tables | C | yes |
| C3 | Publish the dataset summary as a page | C | yes |
| D1 | Benchmark 7.1, compare download routes | D | yes |
| D2 | Implement popseq fetch with the four tiers | D | no |
| D3 | Record every download including failures | D | yes |
| E1 | Rename to projects/PopSeq | E | no |
| E2 | Add the metagenome entity, metadata only | E | yes |
| E3 | Add the MAG entity | E | yes |
| F1 | Finish the 114-species reconciliation | F | no |
| F2 | Tests, and a no-network smoke test | F | yes |
| F3 | Decide the DuckDB question with a measurement | F | yes |
