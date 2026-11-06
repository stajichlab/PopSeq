# F3: Decide the DuckDB question with a measurement

**Track:** F   **Depends on:** A4   **Parallel-safe:** yes

Spec 001 left it open deliberately. The question is not abstract: time the query
"every paired Illumina WGS run above 20x for a given taxid, across all species" over the
CSVs with pandas, and over a DuckDB index built from them.

Do this **after** A4, because the core and long split may remove the need entirely.

**Acceptance.** Spec 001 records either a measured decision to adopt DuckDB, or a
measured statement that it was not needed, with the timings behind it.
