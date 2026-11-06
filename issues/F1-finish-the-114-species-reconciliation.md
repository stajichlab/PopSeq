# F1: Finish the 114-species reconciliation and publish both difference lists

**Track:** F   **Depends on:** nothing   **Parallel-safe:** no, serialise

**In progress.** Two workers, 81 of 114 species done at the time of writing.

Three-species trial result: ENA only 14 runs (0.46 percent, close to MATPredict's
independently measured 0.55 percent); NCBI only 914 runs, 23 percent of its set, which
led directly to issue A2.

**Acceptance.** `catalog/kingdom/fungi_byspecies.ena_only_runs.csv.gz` and
`.ncbi_only_runs.csv.gz` committed, with per-species counts. The `ena_only` list feeds
`unavailable_at_ncbi` in spec 002.
