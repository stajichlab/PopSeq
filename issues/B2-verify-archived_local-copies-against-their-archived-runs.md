# B2: Verify archived_local copies against their archived runs

**Track:** B   **Depends on:** nothing   **Parallel-safe:** yes

**Problem.** Spec 005 section 7 states this and nothing does it. For a sample classed
`archived_local`, nobody has compared the local file with the archived run. A local copy
that was trimmed, subsampled or re-paired is not the archived run, and treating them as
interchangeable makes an analysis irreproducible.

314 samples in `B_dendrobatidis` alone are in this class.

**Do.** For each `archived_local` sample, count reads in the local file and compare with
the ENA `read_count` for its accession. Record `local_matches_archive` as true, false, or
unknown, plus both counts.

**Acceptance.** B_dendrobatidis reports a verdict for all 314. Any mismatch is listed
with both counts, not silently reclassified.

**Cost.** Counting reads means decompressing. Use `${SCRATCH:?}`, and size the work into
SLURM array tasks of about an hour.
