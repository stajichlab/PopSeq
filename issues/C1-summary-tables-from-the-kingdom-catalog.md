# C1: Summary tables from the kingdom catalog

**Track:** C   **Depends on:** nothing   **Parallel-safe:** yes

All input exists: `fungi_species_summary.csv.gz` and the 218,608-row BioSample table.
No fetching required.

**Do.** Produce, as committed CSVs:
- runs and isolates by species by year of first public release
- isolates per study per species, as a clonality proxy
- country coverage per species, from the harmonized `geo_loc_name` (53.1 percent filled)
- estimated coverage distribution per species, from `base_count` over genome size
- host and isolation_source breakdown where populated (28.2 and 36.5 percent)

**Acceptance.** Each table regenerates from a single script and carries the date and the
query it came from.

**Depends on nothing. Start any time.**
