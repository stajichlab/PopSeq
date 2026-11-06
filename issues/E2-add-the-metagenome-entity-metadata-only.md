# E2: Add the metagenome entity, metadata only

**Track:** E   **Depends on:** A1   **Parallel-safe:** yes

Measured 2026-11-06: `tax_tree(4751) AND library_source="METAGENOMIC"` returns **22,555
runs** that the current `GENOMIC` filter discards. That is about 9 percent again on top
of the isolate catalog, and it is currently invisible.

**Do.** Add `catalog/metagenomes/` with the run table only. **Make no eligibility claim.**
An isolate run can be judged from metadata; a metagenome cannot, because the question is
"does this sample contain enough of the target taxon", which needs a classifier run
against the reads.

**Acceptance.** 22,555 runs catalogued with study, environment and sample metadata. No
`eligible` column.
