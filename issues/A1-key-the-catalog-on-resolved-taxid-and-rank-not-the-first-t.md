# A1: Key the catalog on resolved taxid and rank, not the first two words of scientific_name

**Track:** A   **Depends on:** nothing   **Parallel-safe:** yes

**Problem.** The `species` column in `catalog/kingdom/fungi_species_summary.csv.gz` is
built by taking the first two words of the deposited `scientific_name`. A run deposited
at a rank above species produces a row that looks like a species but is not.

Measured 2026-11-06, among the 114 species with 100+ isolates, five such rows exist:
`Leucocoprinus sp.` (240), `Fungi` (217), `Saccharomycotina` (145), `Coemansia sp.`
(128), `Histoplasma` (105).

This is not cosmetic. `popseq_reconcile.py` resolved `Fungi` to taxid 4751 and spent
17 minutes subdividing the entire kingdom inside a job meant to survey one species. A
guard now catches it (`RankArtifact`), but the guard treats a symptom.

**Do.** Resolve every deposited name to a taxid, look up its rank, and key the catalog on
the taxid. A name at genus rank or above gets its own row class, not a species row.
`Leucocoprinus sp.` and `Ascomycota sp.` are buckets of unrelated organisms and must be
labelled as such (see spec 006).

**Acceptance.** No row in the species summary has a rank above species. The five known
artifacts are each classified correctly. `popseq_reconcile.py --species-list` runs with
zero `RankArtifact` skips.
