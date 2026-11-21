# catalog/qc/findings — hand-entered QC inputs

Inputs for run-QC sources that have no producer file of their own: findings reported in an issue, and measurements a
consumer made. They are edited by hand and imported into `catalog/qc/` (spec 008). The per-species `*.qc.csv.gz` files
are the output; edit these inputs, not the outputs, and import again.

| file | used by | source_id | what |
|---|---|---|---|
| `matpredict_afum_species_id_2026-11-10.csv` | `popseq_qc_import.py --format popseq_findings` | same as the file name | 6 *A. fumigatus* runs whose reads are another species or unresolved (issue #37, MATPredict) |
| `matpredict_afum_read_sets_2026-11-10.csv` | `popseq_qc_import.py --format popseq_findings` | same as the file name | 54 tiny read sets filed as WGS (issue #37, MATPredict) |
| `archive_sizes.csv` | `popseq_qc_checks.py archive-size --sizes` | `gmg_archive_size_2026-11-21` | archive object sizes measured by a consumer (issue #49, GmG P6) |

The `popseq_findings` columns are the QC column names (`run_accession`, `status`, and any of `filed_name`,
`detected_genus`, `detected_species`, `subsample`, `strain`, `source_note`). Every row's `source_note` names its issue.
The exact import commands are in `specifications/008-run-qc-flags.md` (section "Checks and sources added 2026-11-21").
