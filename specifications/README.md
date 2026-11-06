# Specifications

One file per topic. These are the interface contracts other projects code against.

| Spec | Topic | Status |
|---|---|---|
| `001-catalog-schema.md` | species YAML, run CSV, biosample CSV field definitions | draft |
| `002-download-strategy.md` | routes, depth tiers, verification, sidecars | draft, blocked on benchmark 7.1 |
| `003-sample-sheets.md` | `samples.csv` and `metadata.csv` contract | draft |
| `004-response-to-handoffs.md` | what was adopted from the handoffs, and what is still needed from each project | open |
| `005-local-and-unarchived-data.md` | local read data with no accession; provenance classes; reading the Population_Genomics convention | draft |
| `006-beyond-fungal-isolates.md` | other organism groups, metagenomes, MAGs, uncultured taxa | draft |

## Handoffs

Written by each project's own session, describing what that project actually did. They are
**evidence, not specification**: they record measurements and decisions made elsewhere, and
nothing in them is approved for PopSeq.

| File | From | What it covers |
|---|---|---|
| `MATpredict_spec_handoff.md` | MATPredict | fungal-wide ENA survey, AWS Open Data fetch route, window subsampling, per-run state records, contamination screening |
| `GmG_spec_handoff.md` | GmG / popsketch | 2,002 runs already fetched across 7 species, concurrency failures, what centralizing needs |
| `MATpredict_reply_to_004.md` | MATPredict | measured AWS transfer rate over 4,211 runs, corrected `unavailable_at_ncbi` rate, a correction to my reading of their filters |

**Naming:** numbered files are specifications. Handoffs and replies use
`<Project>_<kind>.md` so the two schemes cannot collide. `MATpredict_reply_to_004.md` was
renamed from `005-matpredict-reply.md` on 2026-11-06 for this reason.

Specs 001 and 002 were revised after reading these. See the Changes section of each.

## How a spec changes

A spec is an interface other projects depend on. Changing one breaks their code.

1. **Additive changes** (a new optional column, a new enum value) can be made in a commit that
   says what was added. Readers that ignore unknown columns keep working.
2. **Breaking changes** (renaming a column, changing a key, changing a path) need the old
   form to keep working for a deprecation period. Add the new form, write both, announce it,
   remove the old form in a later commit.
3. Every spec carries a **Changes** section at the bottom with dates. A spec with no history
   cannot be trusted by someone writing against it.

## Status values

- **draft** — being written, do not code against it yet.
- **blocked** — written, but depends on a measurement that has not been made. The blocking
  item is named.
- **stable** — other projects may code against it. Breaking changes follow rule 2 above.
