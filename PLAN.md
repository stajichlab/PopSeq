# PopSeq: a shared catalog of fungal population sequencing data

Status: **plan, not yet built.** Nothing in this document is a measurement unless it says
where the measurement came from. Several numbers that matter are marked UNMEASURED and
are the first work items.

Location: `/bigdata/stajichlab/shared/projects/Fungi/PopSeq`

## 1. Why this exists

Three projects in the lab independently wrote code to find, filter, and download public
fungal short-read data:

| Project | What it fetches per run | Keeps the reads? |
|---|---|---|
| `TE_fungi_copy` | the whole run, full depth | yes, while k-mer counting runs |
| `MATPredict` | 2,000,000 spots from 8 spaced windows, about 1 GB | no, deleted on exit |
| `GmG` `reads_sra` | the first 500,000 pairs, about 96 MB | no, deleted after scoring |
| `GmG` `popsketch` | every Nth pair, about 1 Gb of bases | no, only a 0.1-0.4 MB sketch kept |

Each wrote its own ENA query, its own eligibility filter, its own "one run per BioSample"
rule, and its own download loop. The work that is genuinely shared is **finding and
describing runs**, not storing them. Three of the four deliberately never keep the bytes.

They also surveyed the same archive independently and got nearly the same answer, which is
a good sign for the query but a waste of four efforts:

| Survey | Date | Fungal Illumina WGS runs found |
|---|---|---|
| MATPredict, `tax_tree(4751)` | 2026-10-26 | 245,243 |
| GmG, per `first_public` year | earlier | 236,453 |

So PopSeq is a **catalog and a fetcher**, and only incidentally a cache.

### What PopSeq is

- A living, regenerable catalog of public fungal WGS runs, organised by species and NCBI taxid.
- A single verified fetcher with a depth argument, so a caller asks for what it needs.
- A convention for describing a strain panel (`samples.csv` + `metadata.csv`) that any
  pipeline can consume.
- A measured record of download performance, so the cache question is settled with data.

### What PopSeq is not

- Not an analysis pipeline. It produces inputs; it does not interpret them.
- Not a mirror of the SRA. The cache holds what a project is actively using.
- Not a replacement for ENA or NCBI as the source of truth. It is a derived index.

## 2. Design decisions, and why

### 2.1 The catalog is text; the index is derived

The catalog is the source of truth and lives in git as YAML and CSV. DuckDB, if it earns
its place, is a **derived index regenerated from those files** and is never edited directly.

The reason is review. This is a living database that will be corrected by hand when NCBI
metadata is wrong, and a wrong correction must be visible in a diff. A binary database
cannot be reviewed in a pull request. A 300,000-row CSV can be, one line at a time.

DuckDB gets evaluated on a real question, not in the abstract: *how long does it take to
select every paired Illumina WGS run above 20x for a given taxid, across all species?*
If pandas over the CSVs answers that in under two seconds, DuckDB is not needed yet.

### 2.2 YAML for species, CSV for runs

- **Species-level records are YAML**, one file per species. There are hundreds of them, they
  are hand-edited, and they carry nested fields (synonyms, reference assemblies, notes).
  YAML survives hand-editing and reviews well.
- **Run-level records are CSV**, one file per species. There are hundreds of thousands of
  rows, they are machine-generated, and they are flat. CSV is the right shape and every
  tool reads it.

JSON was considered and rejected for the species files: it has no comments, and a catalog
that records *why* a strain was excluded needs comments.

### 2.3 Catalog and cache are separate trees

The catalog is small (tens of MB) and tracked. The cache is large (TB) and ignored. They
must not be interleaved, or `git status` becomes unusable and a stray `git add` becomes a
disaster.

## 3. Directory layout

```
PopSeq/
├── AGENTS.md              rules for this repo (CLAUDE.md symlinks here)
├── PLAN.md                this file
├── README.md              short orientation
├── pyproject.toml         uv-managed environment
├── specifications/        numbered specs, one per topic
├── catalog/               THE LIVING DATABASE (tracked)
│   ├── species/           <taxid>.<genus>_<species>.yaml
│   ├── runs/              <taxid>.<genus>_<species>.runs.csv.gz
│   ├── biosamples/        <taxid>.<genus>_<species>.biosample.csv.gz
│   └── popseq.duckdb      derived index, regenerable, NOT tracked
├── panels/                curated subsets (tracked)
│   └── <panel_name>/
│       ├── panel.yaml     what this panel is and who made it
│       ├── samples.csv    sample,fastq_1,fastq_2
│       └── metadata.csv   sample,<everything known about the strain>
├── cache/                 FASTQ (NOT tracked)
│   └── <taxid>/<RUN>/<RUN>_1.fastq.gz, <RUN>_2.fastq.gz, <RUN>.json
├── scripts/
├── benchmarks/            measured download and storage performance (tracked)
├── logs/                  SLURM logs (NOT tracked)
└── docs/
```

### 3.1 Why the cache is keyed by taxid, not by project

A run belongs to an organism, not to a question. Keying by project guarantees the same run
is downloaded twice when two projects ask for it. Keying by taxid makes the second request
a cache hit. Projects reference the cache by path from their own `samples.csv`; they do not
copy the reads.

## 4. The sample sheet contract

Any pipeline that wants to consume a PopSeq panel reads two files. The join key is `sample`.

**`samples.csv`** follows the nf-core convention, so nf-core pipelines consume it unchanged.
A sample may appear on more than one line when it has several runs.

```csv
sample,fastq_1,fastq_2
Bd_AVS4,/bigdata/.../cache/109871/SRR6252319/SRR6252319_1.fastq.gz,/bigdata/.../cache/109871/SRR6252319/SRR6252319_2.fastq.gz
Bd_JAM81,/bigdata/.../cache/109871/SRR630155/SRR630155_1.fastq.gz,/bigdata/.../cache/109871/SRR630155/SRR630155_2.fastq.gz
Bd_JAM81,/bigdata/.../cache/109871/SRR630566/SRR630566_1.fastq.gz,/bigdata/.../cache/109871/SRR630566/SRR630566_2.fastq.gz
```

Rules:
- `fastq_2` is empty for single-end data. The column is always present.
- Paths are absolute. A relative path breaks the moment a pipeline changes directory.
- A sample with two lines has two runs. The consumer decides whether to merge them.

**`metadata.csv`** is one line per sample, and is deliberately wide. It carries everything
public that was recoverable, plus curated fields.

```csv
sample,biosample,taxid,scientific_name,strain,isolate,host,geo_loc_name,collection_date,sample_type,curated_clade,curated_notes
Bd_JAM81,SAMN12327137,109871,Batrachochytrium dendrobatidis,GPL,Bd_JAM81,Rana sierrae,USA:Sierra Nevada Mtns,2006,isolate,GPL,
```

Rules:
- Columns beginning `curated_` are **ours**, entered by hand, and are never overwritten by a
  catalog refresh. Everything else is harvested and may be replaced.
- An absent value is an empty field. The string `NA` is never written, because `NA` is a real
  country code and a real strain name.
- `taxid` is mandatory. A panel without taxids cannot be regenerated.

## 5. The living database: how it regenerates

Two sources, used for what each is good at.

| Source | Good at | Used for |
|---|---|---|
| **ENA Portal API** | bulk listing, one request returns every run for a clade, with `base_count`, `read_count`, `library_*`, `fastq_ftp`, `fastq_aspera` | the run table |
| **NCBI Entrez** (`esearch`/`epost`/`efetch`) | rich, harmonized BioSample attributes: host, geography, collection date, isolation source | the biosample table |

ENA alone is not enough, because ENA's sample metadata is thin. NCBI alone is not enough,
because walking it run by run at 3 requests per second does not scale to 236,000 runs.

### 5.1 Refresh procedure

1. `popseq refresh-runs --taxid <N>` queries the ENA portal with
   `query=tax_tree(<N>) AND library_strategy="WGS" AND library_source="GENOMIC"` and writes
   `catalog/runs/<taxid>.<name>.runs.csv.gz`.
2. `popseq refresh-biosamples --taxid <N>` takes the distinct `sample_accession` values from
   that table, posts them to NCBI in batches with `epost`, fetches with `efetch`, and writes
   `catalog/biosamples/<taxid>.<name>.biosample.csv.gz`.
3. `popseq refresh-species --taxid <N>` recomputes the summary in
   `catalog/species/<taxid>.<name>.yaml`: run count, distinct isolate count, coverage
   distribution, study count, estimated download size.
4. The diff is reviewed and committed. **A refresh that drops rows is a red flag**, not a
   routine update; records disappear from public archives rarely, and usually it means the
   query changed, not the archive.

Refreshes are incremental. A run already in the table is not re-queried unless `--update`
is passed. This is the behaviour of `biosample2table.py` in
`~/projects/biosample_metadata/scripts/` and it is the right behaviour; see section 9.

### 5.2 Species record

```yaml
taxid: 109871
scientific_name: Batrachochytrium dendrobatidis
common_names: [Bd, chytrid fungus]
assumed_genome_mb: 24
reference_assemblies:
  - accession: GCA_000149865.1
    strain: JEL423
    note: most used reference
  - accession: GCA_003595275.1
    strain: RTP6
    note: single-copy baseline donor in TE_fungi_copy
counts:
  last_refreshed: 2026-11-06
  runs_total: 0            # UNMEASURED until first refresh
  runs_paired_illumina_wgs: 0
  isolates_distinct: 0
  isolates_above_15x: 0
  estimated_download_tb: 0.0
notes: |
  Why this species is in the catalog, who asked for it, and anything known
  about the archive records being wrong.
```

## 6. Download strategy

### 6.1 Depth tiers

A caller states what it needs. The fetcher does not guess.

| Tier | What it gets | Typical size | Who asks for it |
|---|---|---|---|
| `head` | first N read pairs (default 500,000) | ~100 MB | GmG-style classification |
| `sample` | first N bases, default 2 Gbp | 1-2 GB | MATPredict-style typing |
| `full` | the entire run | 0.5-15 GB | TE_fungi_copy-style k-mer counting |

Only `full` is worth caching. `head` and `sample` are cheap enough to refetch and are
written to `$SCRATCH` by default.

### 6.2 Route order, and the conflict that must be resolved first

**This is the most important open question in the plan.**

`TE_fungi_copy/CLAUDE.md` records a measurement from 2026-10-24 on the UCR HPCC login node:
NCBI `prefetch` ran at 20-23 MB/s; ENA FTP over `curl` ran at about 0.4 MB/s, and one ENA
connection was refused.

`GmG/reports/validation/reads_sra.md` records that ENA HTTPS delivered all 100 of 100 runs,
9.57 GB total, with the `fastq-dump` fallback never needed.

These are not necessarily in conflict. GmG pulls small HTTP byte ranges and closes the
connection early, which is a different access pattern from streaming a whole run. But
**nobody has measured them against each other**, and a shared fetcher cannot pick a default
until someone does.

Proposed order, to be confirmed by benchmark 7.1:

1. **NCBI `prefetch` + `parallel-fastq-dump`** for `full`. This is the measured-fastest
   route for whole runs on this cluster.
2. **ENA Aspera (`ascp`)** where `fastq_aspera` is present in the ENA record.
   Module `aspera/4.2.13.820` is installed. **UNMEASURED on this cluster.** Aspera is
   commonly much faster than FTP over long distances because it does not use TCP congestion
   control, but that claim is not tested here and must not be repeated as fact until it is.
3. **ENA HTTPS** for `head` and `sample`, because a range request on a gzipped FASTQ is the
   only route that can stop early. `prefetch` cannot return a partial run.
4. **`fasterq-dump`** as a last resort.

### 6.3 Verification is mandatory

There is no checksum for a converted FASTQ. Every download records the read count of both
mates and compares it with the ENA `read_count`. A mismatch fails the job. This rule exists
because it has already caught a real corruption: on 2026-10-24 a leftover ENA loop resumed
onto an SRA-derived file and corrupted mate 1, and only the count check found it.

After verification passes, the `.sra` is deleted. It is a large duplicate.

Each cached run carries a sidecar `<RUN>.json` recording route, bytes, seconds, read counts,
the ENA values they were checked against, and the tool versions. A cached run with no
sidecar is **not trusted** and is re-verified before use.

### 6.4 Size estimation

For planning before anything is downloaded:

```
gzipped_fastq_bytes  ~=  base_count / 3.5
coverage             =   base_count / (assumed_genome_mb * 1e6)
```

The `/3.5` figure is inherited from the `TE_fungi_copy` SRA survey. It is an approximation
that has not been checked against the 42 FASTQ files now on disk. **Benchmark 7.2 checks it**,
because every storage decision in this document rests on it.

## 7. Benchmarks to run before the design is fixed

These are the measurements the plan depends on. Each writes a row to
`benchmarks/<name>.csv`. None has been run.

**7.1 Route comparison.** 10 runs spanning 0.5-10 GB, each fetched by `prefetch`, ENA
Aspera, and ENA HTTPS, on a compute node and on the login node. Record MB/s, wall time,
failures. Settles section 6.2.

**7.2 Size model.** For every run already on disk in `TE_fungi_copy`, compare actual
gzipped bytes against `base_count / 3.5`. Report the ratio distribution, not the mean alone.
Settles section 6.4.

**7.3 Cache versus refetch.** The real question is not "is storage cheaper than bandwidth",
it is "how often is a run requested twice?" Measure: across the three projects' existing
panels, how many run accessions appear in more than one? If the answer is near zero, the
cache should be small and short-lived. This is a catalog query, not a download, and should
be done first because it is free.

**7.4 Filesystem choice.** Write and read throughput for `/bigdata` against node-local
`/scratch` for the one access pattern that matters: streaming a 5 GB gzipped FASTQ through
`jellyfish count`. Decides whether the cache should ever be staged.

**7.5 DuckDB.** Time the section 2.1 query over the CSVs with pandas, and over a DuckDB
index built from them. Adopt DuckDB only if the gap is large on a realistic catalog size.

## 7a. What already exists, and should be backfilled first

This was found by reading the handoffs in `specifications/`, not by planning. It is the
cheapest material in the project, because it is already downloaded and already curated.

### GmG / popsketch: 2,002 runs across 7 species, already fetched

`/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra/<slug>/`, **23 GB**, verified
2026-11-06. Each slug holds `runs.tsv`, `runs_ftp.tsv` and `sketches/`.

| Slug | Runs selected | Sketches |
|---|---|---|
| `z_tritici`, `c_neoformans`, `c_beticola`, `p_nodorum`, `b_cinerea` | 300 each | complete |
| `n_discreta` | 223 (all that pass) | complete |
| `p_oryzae` | 300 | 279, the other 21 failed twice |

The `runs.tsv` schema is already close to spec 001:
`run_accession, study_accession, sample_accession, species, strain, base_count, read_count,
coverage, center_name, country, genome_mb`.

It carries a **per-species `genome_mb`**, which is better than the single assumed size the
first draft of this plan proposed, and it is adopted.

Selection used: paired Illumina, coverage >= 15x, one run per biosample, spread across
studies, mutant/evolved/cross/knockout titles removed; *C. auris* and *S. cerevisiae*
skipped as too large.

### MATPredict: a fungal-wide survey and a per-run state record

`datasets/reads/sra_survey/` in the MATPredict repo. 649 species with at least 4 strains
after exclusions, 74,905 strains. Per-run `<run>.run.json` records state, source, byte count
and timings, which is the same idea as the sidecar in spec 002 and validates it.

### What backfilling gets us

Nothing needs to be downloaded to populate the catalog for 7 species. The run tables exist,
the selection criteria are recorded, and the coverage and country fields are already there.
This is phase 1 work and it is free.

## 8. Phases

| Phase | Deliverable | Depends on |
|---|---|---|
| 0 | repo skeleton, `AGENTS.md`, `pyproject.toml`, uv environment | - |
| 1 | backfill the 7 popsketch species from section 7a, and benchmark 7.3. Both free. | - |
| 2 | `popseq refresh-runs` and one species end to end | 1 |
| 3 | `popseq refresh-biosamples` using batched Entrez | 2 |
| 4 | benchmarks 7.1, 7.2, 7.4 | 2 |
| 5 | `popseq fetch` with the three depth tiers and the route order 4 settled | 4 |
| 6 | `popseq panel` to emit `samples.csv` + `metadata.csv` | 3, 5 |
| 7 | backfill the remaining projects' panels into the catalog | 6 |
| 8 | benchmark 7.5, adopt DuckDB or record that it was not needed | 6 |

Phase 1 comes before any download because it can change the design of the cache for free.

## 9. Reusing the biosample scripts

`~/projects/biosample_metadata/scripts/biosample2table.py` already does the hard part:
it resolves an SRA run or a BioSample accession to NCBI's internal id, fetches the sample,
reads every `Attributes/Attribute`, prefers `harmonized_name` over `attribute_name`, and
unions the columns into a growing wide table. It also loads an existing output table and
queries only what is new.

Three changes are needed for catalog scale, and they are the reason this is a rewrite and
not an import:

1. **Batch the requests.** The script calls `esearch` then `efetch` once per accession.
   At the unauthenticated limit of 3 requests per second, 10,000 BioSamples take about two
   hours. `epost` plus a single `efetch` handles 200 accessions per round trip.
2. **Support an API key.** `Entrez.api_key` raises the limit from 3 to 10 requests per second.
   The key comes from the `NCBI_API_KEY` environment variable and is never written to a file
   in this repo.
3. **Retry and resume.** Entrez fails intermittently. A harvest of 10,000 samples that dies
   at 9,000 and starts over is unusable. Write results in batches; resume from what is on disk.

The attribute-flattening logic is kept as is, and the incremental-update behaviour is kept
as is. Both are correct.

## 10. Open questions

1. **Route order** (section 6.2). Blocks the fetcher default. Benchmark 7.1.
2. **Is Aspera actually faster here?** Untested. Benchmark 7.1.
3. **How much cross-project reuse is there?** Benchmark 7.3. If it is near zero, the cache
   in section 3 is over-designed and should shrink.
4. **Who owns a correction?** When NCBI says a strain is *Fusarium oxysporum* and we believe
   it is something else, the `curated_` column records our view. There is no process yet for
   two projects disagreeing about a curated value.
5. **Retention.** Nothing here says when a cached run is deleted. It needs a rule before the
   cache reaches a terabyte, not after.
6. **Where the data lives.** GmG asks this directly: the popsketch data sits under a personal
   projects folder "by instruction". Moving it into shared space is a curator decision, not
   an agent one. Until it is made, the catalog records the current path.
7. **Should the fetcher always sketch?** GmG proposes writing a popsketch sketch during the
   fetch step, before the reads are deleted. A sketch is 0.1-0.4 MB and cannot be recovered
   without re-downloading the run. The cost is near zero and the option expires when the
   reads do. This looks right, but it couples PopSeq to a `popsketch` binary, so it is a
   decision, not an obvious yes.
8. **What is the real concurrency limit?** Failures appeared at about 35 concurrent
   `parallel-fastq-dump` downloads. The safe ceiling is unmeasured. Benchmark 7.1 should
   find it, because every array job in every project needs the number.
9. **Canonical taxonomy.** All three projects hit species synonyms (*Magnaporthe* /
   *Pyricularia*, *Candida glabrata* / *Nakaseomyces glabratus*, *Candida* / *Candidozyma
   auris*). The taxid-first naming in spec 001 handles storage, but nothing yet decides which
   name is displayed, or which NCBI taxdump version is canonical.


## 11. Provenance of this plan

The first draft was written from `TE_fungi_copy` alone. It was then revised after reading
`specifications/MATpredict_spec_handoff.md` and `specifications/GmG_spec_handoff.md`, which
those projects wrote independently on 2026-11-06.

Four things the first draft got wrong or missed, all corrected above:

1. It had three download routes and missed **AWS Open Data**, which is MATPredict's primary.
2. It proposed a `head` subsample as a tier. Two projects had already concluded the head of a
   run is a biased sample. Spec 002 now restricts `head` and defaults to `windows`.
3. It had no concurrency limit. Failures were already observed at about 35 concurrent.
4. It assumed nothing had been fetched. 2,002 runs across 7 species already have been.

Every accession used as an example in this plan and in the specs was fetched from NCBI or ENA
on 2026-11-06 and checked. An earlier draft used invented accessions; they were replaced.
