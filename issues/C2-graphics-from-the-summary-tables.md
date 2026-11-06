# C2: Graphics from the summary tables

**Track:** C   **Depends on:** C1   **Parallel-safe:** yes

**Do.**
- cumulative runs by year, whole kingdom and per candidate species
- isolates-per-species rank curve on log axes; 2 species exceed 10,000 isolates and
  5,435 of 5,549 have fewer than 100, so the shape matters
- coverage histograms for the candidate species
- a world map of `geo_loc_name`

Load the `dataviz` skill before writing any chart code.

**Acceptance.** Figures regenerate from a script, render legibly in light and dark, and
each states its n and its date.
