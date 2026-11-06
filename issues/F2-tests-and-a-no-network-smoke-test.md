# F2: Tests, and a smoke test that runs without network

**Track:** F   **Depends on:** nothing   **Parallel-safe:** yes

There are no tests. Every bug found so far was found by running the thing at scale, and
three of them produced plausible wrong answers rather than failures:

- `write_table` wrote plain text into a `.gz` filename, because the temp name hid the
  extension from the gzip chooser. Only visible at 205 MB.
- A process substitution read twice gave an empty species scope, and the script reported
  a clean reconciliation of nothing against nothing.
- `done < file </dev/null` made a loop read `/dev/null` and report 0 of everything.

Each of those printed a plausible number. That is the failure mode worth testing against.

**Do.** Unit tests over recorded fixtures for: FileBase resolution in all four forms,
mate pairing, gzip round trip, shard resume, the rank-artifact guard, and multi-valued
`sample_accession`.

**Acceptance.** `uv run pytest` passes with no network access.
