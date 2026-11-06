# D2: Implement popseq fetch with the four depth tiers

**Track:** D   **Depends on:** D1   **Parallel-safe:** no, serialise

Spec 002 defines `head`, `windows`, `stride` and `full`. None is implemented.

**Do.** One fetcher, taking a tier and a run accession. Enforce the spec's hard rules:
- verify both mates against the stored ENA `read_count`; fail the task on mismatch
- delete the `.sra` only after verification passes
- write the `<RUN>.json` sidecar; a cached run without one is not trusted
- both mates from one source, never mixed
- `head` is restricted to questions insensitive to read origin; `windows` is the default

An AWS 404 is close to final: of 26 runs that fell back to `prefetch`, **`prefetch`
rescued none**. Attempt the fall-through once, never loop it.

**Acceptance.** Fetches one run at each tier, verified, with a sidecar. Resumes correctly
after a kill.
