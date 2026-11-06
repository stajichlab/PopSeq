# PopSeq

A shared, regenerable catalog of public fungal population sequencing data, and one verified
way to fetch it.

Built so that any project asking a question across many strains of a fungus does not have to
write its own ENA query, its own eligibility filter, and its own download loop. Three
projects in this lab already did.

**Status: the catalog is harvested and queryable. Nothing has been downloaded yet.**

245,822 runs / 217,911 isolates / 6,486 species taxids, from the ENA snapshot of
2026-11-06, reconciled against NCBI for the 113 species with 100 or more isolates.

## Using the catalog from another tool

Read `catalog/index.json` first. It names every table, its grain, its key, its
columns, its SHA-256, and the caveats that apply. Nothing else needs reading.

Then ask it questions, four ways. One implementation behind all of them, so
they cannot drift apart. No network, no database, no credentials.

**Shell**

```bash
uv run popseq species --min-isolates 100 --format tsv

# 100 Aspergillus fumigatus strains at 30x or better, spread across studies,
# as an nf-core sample sheet
uv run popseq panel --taxid 746128 --n 100 --min-depth 30 --format samples > samples.csv

uv run popseq describe --taxid 1047171
```

Exit 0 on success, 2 when the species is not found or the name is ambiguous,
3 when the filters leave nothing, so a caller branches without parsing.

**Python**

```python
import popseq
popseq.species(min_isolates=100)
popseq.panel(taxid=746128, n=100, min_depth=30)      # dict, JSON-ready
popseq.to_sample_sheet(popseq.panel(taxid=746128, n=100))
popseq.preload(taxids)      # asking about many species? read the table once
```

**As an agent tool (MCP)**

```bash
claude mcp add popseq -- uv --directory /bigdata/stajichlab/shared/projects/Fungi/PopSeq \
    run --extra mcp python -m popseq.mcp_server
```

Five tools: `popseq_index`, `popseq_species`, `popseq_describe`, `popseq_panel`,
`popseq_runs`.

**Over HTTP, with no clone**

```
https://raw.githubusercontent.com/stajichlab/PopSeq/main/catalog/api/index.json
https://raw.githubusercontent.com/stajichlab/PopSeq/main/catalog/api/746128.json
```

One file per species with 100+ isolates: the full description plus a
pre-selected 100-strain panel.

Verb names, flag names, function names, JSON keys and exit codes are the
contract; table paths and columns are not. See
`specifications/007-query-interface.md`.

Three caveats the interface will not hide from you: every count is a floor
(issue #2), a run count is not a strain count (52 of the 113 largest species
take half or more of their isolates from one study), and `est_depth_x` comes
from metadata rather than from the reads.

[What is in the archive](figures/README.md) has the figures.

## Start here

| File | What it is |
|---|---|
| `PLAN.md` | what is being built, why, and the measurements it depends on |
| `AGENTS.md` | the rules for working in this repo (`CLAUDE.md` symlinks to it) |
| `specifications/` | the detailed specs, one per topic |

## Layout

```
catalog/     the living database: species YAML, run and biosample CSV. Tracked.
panels/      curated strain sets: samples.csv + metadata.csv. Tracked.
cache/       FASTQ, keyed by taxid. Not tracked.
scripts/     harvesters and the fetcher.
benchmarks/  measured download and storage performance. Tracked.
```

## Setup

```bash
cd /bigdata/stajichlab/shared/projects/Fungi/PopSeq
uv sync
uv run python scripts/popseq_biosamples.py --help
```

Set `NCBI_API_KEY` in the environment to raise the Entrez rate limit from 3 to 10 requests
per second. Never commit it.

## Consuming a panel

Any pipeline reads two files, joined on `sample`:

- `panels/<name>/samples.csv` — `sample,fastq_1,fastq_2`, nf-core compatible, one line per run
- `panels/<name>/metadata.csv` — one line per sample, everything known about the strain

## Open questions

The design is not settled. `PLAN.md` section 10 lists what is unresolved, and section 7 lists
the benchmarks that resolve it. The largest open question is which download route is actually
fastest on this cluster: two projects in this lab measured contradictory-looking results, and
nobody has compared them directly.
