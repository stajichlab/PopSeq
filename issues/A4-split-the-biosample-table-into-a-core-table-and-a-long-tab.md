# A4: Split the BioSample table into a core table and a long table

**Track:** A   **Depends on:** A3   **Parallel-safe:** yes

**Problem.** The harvested table is 218,608 rows by 783 columns and **1.38 percent of
cells are filled**. 5 columns appear in at least half the rows, 46 in at least 1 percent,
and 611 of 783 in fewer than 1 in 1,000.

Gzip absorbs this well (196 MB to 4.1 MB), so storage is not the issue. Usability is: a
783-column CSV is awkward to read and is no longer reviewable as a diff, which is the
argument spec 001's format choice rests on.

**Do.** Write two files per clade:
- `*.biosample_core.csv.gz` — the columns present in at least 1 percent of rows.
- `*.biosample_attr.csv.gz` — long form, `biosample,attribute,value`, for the tail.

**Acceptance.** Both files regenerate from the shards. Their union reproduces the wide
table exactly. The core table is under 60 columns.

**Note.** This is the cheaper answer to the DuckDB question spec 001 left open. Do this
first and re-measure before considering DuckDB (issue F3).
