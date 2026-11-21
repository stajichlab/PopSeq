# catalog/qc — per-run quality flags

Findings about individual runs (contamination, mixtures, reads that do not match the species in the metadata), from any
project that screens public runs. The format and the rules are in `specifications/008-run-qc-flags.md`.

- `sources.csv`: the registry. One row per contributing source: who made the calls, with which tool, what the method
  can see, and where it is documented.
- `<species_taxid>.<Genus>_<species>.qc.csv.gz`: one row per flagged run, source and check.
- `unresolved.qc.csv.gz`: rows whose run did not resolve to a species in the PopSeq run table.

- `findings/`: hand-entered inputs (findings reported in an issue, archive sizes measured by a consumer); see its README.

A flag is evidence. Nothing here removes a run from the catalog; `popseq runs` and `popseq panel` skip, by default, the
runs whose status the source lists in `exclude_statuses`, and `--include-flagged` keeps them (spec 008). The
`curated_status` and `curated_note` columns are for the curator and are never overwritten by an import.

Add or refresh a source with `scripts/popseq_qc_import.py` (see its `--help`). Checks PopSeq evaluates from catalog
data (`low_coverage`, `oversized_archive_object`) run through `scripts/popseq_qc_checks.py`; `low_coverage` covers the runs named in issues #37 and #49 only (`findings/issue_named_runs.csv`).
