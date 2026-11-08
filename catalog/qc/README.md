# catalog/qc — per-run quality flags

Findings about individual runs (contamination, mixtures, reads that do not match the species in the metadata), from any
project that screens public runs. The format and the rules are in `specifications/008-run-qc-flags.md`.

- `sources.csv`: the registry. One row per contributing source: who made the calls, with which tool, what the method
  can see, and where it is documented.
- `<species_taxid>.<Genus>_<species>.qc.csv.gz`: one row per flagged run, source and check.
- `unresolved.qc.csv.gz`: rows whose run did not resolve to a species in the PopSeq run table.

A flag is evidence, not an exclusion. Nothing here removes a run from the catalog or a panel. The `curated_status`
and `curated_note` columns are for the curator and are never overwritten by an import.

Add or refresh a source with `scripts/popseq_qc_import.py` (see its `--help`).
