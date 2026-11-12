# 007: The query interface

**Status:** implemented, schema 1.1.0
**Supersedes:** nothing. Extends 001 (catalog schema) and 003 (sample sheets).

## Why this exists

Three projects need the same thing from this catalog: *which strains of species X
can I profile, and how do I get the reads.* Before this, each one would have had
to learn the table layout, the taxid-not-name rule, the exclusion labels and the
clonality caveat on its own. Each would have learned a different subset.

So the catalog advertises one interface and keeps the tables private behind it.

## The contract

Four verbs. **Verb names, flag names, JSON key names and exit codes are the
contract.** Table paths and columns are not; they may
change under a consumer without notice, which is the point.

| verb | answers |
|---|---|
| `species` | what is in the catalog at all |
| `runs` | every run of one species, filtered |
| `panel` | N strains to profile, one per isolate, spread across studies |
| `describe` | everything known about one species, as JSON |

| exit | meaning |
|---|---|
| 0 | ok |
| 2 | species not found, or the name matched more than one |
| 3 | filters left nothing |

A caller branches on the exit code without parsing anything.

### Discovery

`catalog/index.json` is the entry point for a tool that knows nothing. It carries
the schema version, the snapshot date, every table with its grain, key, columns,
size and SHA-256, the column notes that matter, and the caveats. A consumer reads
that file first and needs no other documentation.

### Schema version

`popseq_schema` (in `catalog/index.json`, every `catalog/api/` file and every JSON
answer) is `popseq.SCHEMA_VERSION`. It versions the contract above, not the data.

| change | bump |
|---|---|
| a verb, flag, function, JSON key or exit code removed, renamed or changed in meaning | major |
| a table, column, JSON key or verb added | minor |
| a fix that leaves every name and meaning alone | patch |

A refresh that only changes counts does not change the version. A consumer that
checks the major number keeps working across minor bumps.

| version | date | what changed |
|---|---|---|
| 1.0.0 | 2026-11-06 | first contract |
| 1.1.0 | 2026-11-12 | added only. Run table: `country_iso`, `country_name`, `country_status`, `host`, `host_taxid`, `host_name`, `host_status`, `ena_layout`. New tables: `catalog/species_layout.csv.gz`, `catalog/vocab/{country_aliases,host_overrides}.csv`, `catalog/vocab/host_taxa.csv.gz`, `catalog/kingdom/fungi_assemblies.csv.gz`, `catalog/labels/{assembly_labels,run_labels}.csv.gz`, `catalog/labels/sources.csv`. Index key `archive_snapshots`. 15 new BioSample attribute columns |

The package version (`popseq.__version__`, `pyproject.toml`) is separate. It tracks
the code, and a git tag `v<version>` marks each release.

### Formats

`--format` takes `csv`, `tsv`, `json`, `accessions`, `samples`.

`samples` is the nf-core `sample,fastq_1,fastq_2` sheet from spec 003. Its FASTQ
columns are **ENA URLs, not local paths**, because a run that has not been fetched
has no local path and inventing one would be a lie. A consumer either streams from
those URLs or runs the fetcher and rewrites the sheet against local files.

## Four ways to reach it

One implementation, `popseq/query.py`. Everything below is a wrapper over it, so
none of them can drift from the others.

| reach | for | how |
|---|---|---|
| **CLI** | shell, Makefile, SLURM script | `uv run popseq panel --taxid 746128 --n 100` |
| **Python package** | a pipeline that runs in-process | `import popseq; popseq.panel(taxid=746128, n=100)` |
| **MCP server** | an agent session, as a tool | `uv run --extra mcp python -m popseq.mcp_server` |
| **Static JSON** | anything that can fetch a URL | `catalog/api/<taxid>.json` on raw.githubusercontent.com |

`scripts/popseq_query.py` still works. It is a shim onto `popseq.cli`, kept
because this spec advertised that path before the package existed.

**The Python package** raises `SpeciesNotFound` and `EmptyResult`, which are the
exceptions behind exit codes 2 and 3. A batch caller asking about many species
should call `popseq.preload(taxids)` first: the run table is read once instead of
once per species. Measured: 2.0 s to hold 220,917 rows, after which 20 panels take
0.58 s rather than 44 s. The answers are identical either way.

**The MCP server** exposes `popseq_index`, `popseq_species`, `popseq_describe`,
`popseq_panel` and `popseq_runs`. Register it with:

```bash
claude mcp add popseq -- uv --directory /bigdata/stajichlab/shared/projects/Fungi/PopSeq \
    run --extra mcp python -m popseq.mcp_server
```

It needs `mcp>=2.0`, which renamed `FastMCP` to `MCPServer`. It opens no socket
and holds no credentials.

**The static JSON** covers the 113 species with 100 or more isolates, 5.9 MB in
all. Each file holds the full `describe` payload plus a pre-selected panel of up
to 100 strains at the default filters. It answers the common request in advance;
anything else needs one of the other three, because the filter space is too large
to pre-render. Regenerate with `scripts/popseq_api.py` after the catalog changes,
or the files go stale silently.

## Design decisions, and why

**Taxid, never name.** `--species` is a convenience that resolves through the
catalog and *fails loudly* on an ambiguous match rather than guessing. Programs
should pass `--taxid`. NCBI labels one taxon with several names; a name-keyed
query silently splits it (issue #1).

**Labels, not deletions.** The `exclusion` column marks mutant, evolved, cross and
non-WGS material. `runs` and `panel` drop those by default and `--include-excluded`
keeps them. `--only-exclusion mutant` selects *just* that material, because
TE_fungi_copy wants the mutant and evolved panels that the other two projects
discard. One table serves both.

**The panel rule is explicit and stated in the output.** One run per isolate
(deepest wins), then round-robin across studies. The metadata block that comes
back names the rule and says how many studies the panel touched. It does not claim
the panel is independent: 52 of the 113 largest species take half or more of their
isolates from a single study, so study spread is the best available proxy and
nothing more.

**No network.** Every verb reads committed files. A consumer gets the same answer
offline, on a compute node, and in six months from a pinned commit.

**No database.** The run table is a 14 MB gzipped CSV and a full scan of it takes
about two seconds, which is below the threshold where a consumer notices. DuckDB
stays a derived index behind this interface if benchmark 7.5 ever justifies it
(issue #20); the interface does not change if it arrives.

## What this does not do

It does not fetch anything. It selects. The fetcher is track D (issues #12-#14),
and it will take an `accessions` list from here as its input.

It does not know about local, already-downloaded reads. Spec 005 covers the
Population_Genomics panels and their provenance classes. Merging those into
`panel` output is the obvious next step and is not built.

## For other organism groups

Nothing above is fungal. The scope line in `catalog/index.json` is the only place
the kingdom appears. Spec 006 covers what changes for bacteria, metagenomes and
MAGs; the verbs do not.
