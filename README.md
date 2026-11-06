# PopSeq

A shared, regenerable catalog of public fungal population sequencing data, and one verified
way to fetch it.

Built so that any project asking a question across many strains of a fungus does not have to
write its own ENA query, its own eligibility filter, and its own download loop. Three
projects in this lab already did.

**Status: skeleton. Nothing has been harvested yet.**

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
