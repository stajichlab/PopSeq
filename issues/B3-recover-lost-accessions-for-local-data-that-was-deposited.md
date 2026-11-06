# B3: Recover lost accessions for local data that was deposited

**Track:** B   **Depends on:** B2   **Parallel-safe:** yes

**Problem.** 686 samples across four subprojects are `local_only`: no accession recorded.
Some were almost certainly deposited and the accession was simply never written down.
Those are `local_submitted`, not `local_only`, and the distinction decides whether the
data is irreplaceable.

**Do.** Implement the recovery procedure in spec 005 section 4: search NCBI BioSample for
the strain name within the expected organism, then compare the candidate run's
`read_count` with the local file's counted reads. **A match on strain name alone is not
sufficient** — strain names repeat across studies.

Record how each accession was recovered in `curated_notes`. An unverified guess is not
recorded at all.

**Acceptance.** Every `local_only` sample is either reclassified with a verified
accession, or confirmed as having no public counterpart found.
