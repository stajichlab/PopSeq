# D3: Record every download in benchmarks/downloads.csv, including failures

**Track:** D   **Depends on:** D2   **Parallel-safe:** yes

Spec 002 requires it from the first download, not from whenever someone starts
collecting. A route that works 95 percent of the time and one that works 100 percent
look identical in a table of successes.

Columns: `utc,run,taxid,tier,route,node_type,bytes,seconds,mb_per_s,verified,outcome,error`.
`outcome` is one of `ok`, `verify_failed`, `transfer_failed`, `unavailable_at_ncbi`.

**Acceptance.** Every fetch appends a row. A deliberately failed fetch appends one too.
