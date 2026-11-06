# E3: Add the MAG entity keyed on derived_from_run

**Track:** E   **Depends on:** E2   **Parallel-safe:** yes

A MAG gives a reference; the reads it came from give the within-sample variation. The
link is the point.

Measured 2026-11-06 via ENA `result=assembly` for fungi: 490 flagged
`Metagenome-Assembled Genome (MAG)`, 2,584 `clone or isolate`, and **25,227 with no
`assembly_type` at all**. So the flag finds MAGs that were labelled and misses the rest.
Cross-check NCBI BioSample (`derived from`) as a second source; do not assume they agree.

Uncultured taxa resolve: `Ascomycota sp.` is taxid 1709939, rank species. Storage is
fine. What breaks is that such a taxid is a bucket of unrelated organisms, and the real
information in a name like `Ascomycota sp. Antarctolithica obscura` is in the suffix,
which nothing automatic preserves. Put it in `curated_*` fields.

**Acceptance.** `catalog/mags/` with `derived_from_sample` and `derived_from_run`
populated where recoverable, and an explicit statement of how incomplete the census is.
