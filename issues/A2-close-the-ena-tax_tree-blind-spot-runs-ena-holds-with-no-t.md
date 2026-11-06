# A2: Close the ENA tax_tree blind spot: runs ENA holds with no taxid

**Track:** A   **Depends on:** nothing   **Parallel-safe:** yes

**Problem.** `tax_tree(N)` cannot return a run whose taxonomy ENA has not populated.
Measured 2026-11-06: of 40 runs NCBI had and ENA's query did not, **39 were present in
ENA with an empty `tax_id` and `scientific_name`**, and 1 was absent entirely.

For *Cercospora beticola*, ENA's tax_tree reports 752 runs and NCBI reports about 1,644.

**Every isolate count in the catalog today is a floor, not a total.** This is a
systematic blind spot, not an occasional miss, and it is invisible from inside ENA: the
query returns a clean answer that is simply short.

**Do.** Add a second discovery route that does not go through ENA taxonomy. Options to
evaluate, in order of expected cost:
1. Use the NCBI run list (issue F1 output) as the authoritative accession set and pull
   ENA file URLs per accession via `filereport`.
2. Query ENA by `study_accession` for every study the NCBI side names.
3. Report untaxed runs to ENA. Slow, but it fixes the source.

**Acceptance.** For three test species including *C. beticola*, the catalog's run count
matches NCBI's within 1 percent, and the residual difference is characterised.

**Note.** This invalidates the "ENA is simply better" conclusion in an earlier commit.
ENA is faster. It is not more complete. Spec 002 needs that correction.
