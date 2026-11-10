# catalog/labels — curated species labels (specification 009)

Species names for genome assemblies from a curated authority, beside the NCBI names. Nothing here replaces an NCBI
field. Read `specifications/009-curated-species-labels.md` for the columns and rules.

| file | grain |
|---|---|
| `sources.csv` | one row per label source: authority, release, method, call rule, refresh |
| `assembly_labels.csv.gz` | one row per assembly and source |
| `run_labels.csv.gz` | the same labels on every run whose BioSample is the assembly's |

Which name to use: `curated_species` when set, else `label_species`, else the NCBI species. A tie (`label_status`
`tie`) keeps every name, `;`-joined; `label_complex` names the complex they share.

`curated_species` and `curated_note` are for the curator and are never overwritten by an import. Everything else is
rewritten by `scripts/popseq_labels_import.py`.
