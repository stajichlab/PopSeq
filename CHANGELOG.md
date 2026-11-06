# Changelog

## 2026-11-06

First commit. Repository skeleton, plan, rules, and three draft specifications.

- `PLAN.md`: what is being built, why, and the measurements it depends on.
- `AGENTS.md` (with `CLAUDE.md` symlinked to it): the rules for working here.
- `specifications/001-catalog-schema.md`, `002-download-strategy.md`,
  `003-sample-sheets.md`: draft interface contracts.
- `scripts/popseq_biosamples.py`: batched NCBI BioSample harvester. Tested against NCBI on
  2026-11-06; fetched `SAMN12327137` and `SAMEA5704822` and reproduced the documented output
  of the script it replaces.
- `pyproject.toml` and `uv.lock`: uv-managed environment, Python >= 3.11.

Specs 001 and 002 were revised the same day after reading `MATpredict_spec_handoff.md` and
`GmG_spec_handoff.md`, which the MATPredict and GmG sessions wrote independently. That
revision added the AWS Open Data download route, the `windows` and `stride` subsample tiers,
the rule against using `head` for any cross-sample measurement, the download concurrency cap,
the `unavailable_at_ncbi` run state, and the study-level exclusion categories.

Nothing has been harvested into `catalog/` yet.
