# MATPredict handoff: how fungal population read data were identified, fetched and handled

Written 2026-11-06 by the MATPredict session (Claude) at the request of the curator, Jason Stajich, for the agent working on centralizing fungal popgen data (PopSeq). This describes what MATPredict did. It is not a proposal that anyone has approved, and it changes nothing in MATPredict.

Repo: github.com/stajichlab/MATPredict. On the HPCC: `/bigdata/stajichlab/jstajich/projects/MATPredict/`. Local paths below are relative to the repo root.

## 1. What MATPredict needed the reads for

Typing the mating-type (MAT) idiomorph of individual strains from raw reads, with `python -m MATPredict reads-type`. This needs reads, not assemblies, and only a modest subsample per run, so MATPredict never keeps reads.

## 2. How runs were identified (the survey)

Files: `datasets/reads/sra_survey/` (`README.md`, `survey_ena.py`, `species_summary.tsv`, `studies.tsv`).

- **Source:** the ENA portal API (`https://www.ebi.ac.uk/ena/portal/api/search`), queried 2026-10-26. ENA mirrors SRA. Two `read_run` queries, both `tax_tree(4751) AND library_strategy="WGS" AND library_source="GENOMIC" AND instrument_platform="ILLUMINA"`, `format=tsv&limit=0`:
  - run-level fields: `run_accession, study_accession, sample_accession, secondary_sample_accession, scientific_name, tax_id, library_layout, base_count, read_count, center_name, collection_date`
  - metadata fields: `run_accession, study_title, sample_title, sample_alias, strain, isolate, experiment_title, first_public, country, host`
  - Result: 245,243 fungal Illumina WGS runs. The two raw tables (31 MB and 58 MB) are not committed; the queries reproduce them.
- **Run and species filters** (`survey_ena.py --runs RUNS.tsv --meta META.tsv --out . --db ../../../db`, about a minute):
  - species = first two words of the scientific name; metagenomes, unidentified and `sp.` entries dropped;
  - *Saccharomyces cerevisiae* skipped (curator ruling 2026-10-26: its public sets are mostly mutants, evolution lines and crosses);
  - paired or single-end runs of at least 300 Mb (about 10x for a 30 Mb genome; a per-species genome size would be better, still open);
  - one strain per (strain, isolate or sample alias) per study, keeping its largest run; strain names pooled across studies, which still over-counts somewhat.
- **Study exclusion** (keyword filter on study titles and strain/sample text, recorded in `studies.tsv` as `excluded_as`): `mutant` (knockouts, engineered strains), `evolution` (experimental evolution, serial passage, mutation accumulation), `cross` (progeny, segregants, tetrads, recombinant inbred lines), `not_wgs` (RNA-seq, ChIP, ATAC, bisulfite, Hi-C, metagenome, amplicon, ITS/16S, single-cell). Titles naming natural or clinical sampling (`surveillance`, `outbreak`, `GWAS`, `natural isolates`, `bloodstream`...) are not excluded on a stray keyword. This is a first pass, not a review; it has false positives and misses.
- **Result:** 649 species with at least 4 strains after exclusions (74,905 strains). Reference tiers (MATPredict reference records): A (two or more idiomorphs) 44 species; A1 19; B (genus only) 166; C 420. Counts are dominated by *Candidozyma auris* (24,888, outbreak surveillance, largely clonal), *Zymoseptoria tritici*, *Aspergillus fumigatus*, *Cryptococcus neoformans*, *Candida albicans*.
- **Not checked by the survey:** metadata only, no read set was opened, so mislabelled species and mixed samples are not caught (a *Coccidioides* set had two relabelled strains); clonal sets are counted as strains; restricted or withdrawn runs are not distinguished; non-Illumina data are out of scope. Species names follow the deposited name, so synonyms (*Pyricularia oryzae* / *Magnaporthe oryzae*, *Nakaseomyces glabratus* / *Candida glabrata*) read as a lower tier.

## 3. How a run list is chosen for one species

`datasets/reads/sra_survey/pilot/select_runs.py`: one run per strain; natural-isolate studies only (not flagged in `studies.tsv`); runs between `--min-gb` 0.6 and `--max-gb` 12 Gb, the one nearest `--target-gb` 3 per sample (a few Gb is plenty; a 128 Gb run only costs download time); at most `--per-study` 4 runs per study; truth = the single `detect` idiomorph of a BFD assembly with the same BioSample (balanced between idiomorphs); `--include` adds named runs. The input run list is a TSV: `run_accession, read_count, group, library_layout, strain, study_accession, truth_detect, truth_assembly` (see `pilot/CONTRACT.md`).

## 4. How reads were fetched

The unit of work is `datasets/reads/sra_survey/pilot/genotype_mat_one_accession.sh RUN SPOTS PANEL_DIR PANEL_PREFIX OUT_DIR [GROUP] [LAYOUT]`.

1. **Download the `.sra`** with `curl -L -C -` (resumable, 3 tries with growing sleeps) from NCBI's copy on AWS Open Data: `https://sra-pub-run-odp.s3.amazonaws.com/sra/<RUN>/<RUN>`. A 404 means the run is not mirrored there (some European runs are not); then fall back to `prefetch --max-size 50G` from the SRA toolkit.
2. **If neither has it** (an ENA-only run, for example ERR2188722 on 2026-10-26), the state is `unavailable_at_ncbi`. About 1 run in 100 so far. It is recorded, not treated as a pipeline failure.
3. **Subsample by windows, not by head:** `TARGET_SPOTS` (default 2,000,000) spots in `WINDOWS` (default 8) evenly spaced windows along the run, with `fastq-dump -N <start> -X <end> -Z --split-spot`. A run no larger than the target is read whole. Spreading windows avoids a head-of-file bias. `--split-spot` gives both mates of a spot the same read name, so the name identifies the spot.
4. **Type** with `reads-type` (k-mer panel; `READS_K` default 31; `FALLBACK_K=21` retypes `none` and `low_depth` runs at k=21).
5. **Delete** the `.sra` and every FASTQ on exit (`trap`). Reads are never kept; `KEEP_SRA=1` and `FASTQ_OUT=path` exist for follow-up analyses only.
6. **Temporary space:** `SCRATCH` (default `/tmp`); a few GB per task.

Measured on 20 runs (2026-10-26): about 1 GB download, 110 s, 5 GB peak memory (ask for 8), up to 4 CPUs in the typing step.

**Idempotent and runner-independent:** the result `.tsv` is the done-marker, so a runner may re-run the whole list and only missing runs work. No shared state between runs; each task has its own temp dir. Exit 0 done or already done, 1 failed, 2 bad arguments. States in `<run>.run.json`: `done`, `unavailable_at_ncbi`, `failed_download`, `failed_convert`, `failed_type`, `started` (task died before writing a state).

**Runners:** Slurm arrays (`stream_type_sra.slurm`, `submit_species.sh SLUG RUNS.tsv`, which queues an array, a retry array and a summary job); a Nextflow example exists (`main.nf.example`) but has not been run. Concurrency is capped with the Slurm array throttle (`%N`). I have no measured failure threshold for concurrent downloads. Note the popsketch agent saw `parallel-fastq-dump` exit 64 at about 35 concurrent downloads, so a cap is worth having.

## 5. Outputs and where they live

Per run, two small files in `OUT_DIR`:
- `<run>.<prefix>.tsv`: `call` (`MAT1-1`, `MAT1-2`, `both`, `none`, `low_depth`; names follow the panel), `MAT1-1_breadth`, `MAT1-2_breadth`, `MAT1-1_depth`, `MAT1-2_depth`, `shared_depth`, `reads_used`, `flags`.
- `<run>.run.json`: `run, group, layout, spots_total, target_spots, windows, sra_bytes, secs_download, secs_convert, secs_type, secs_screen, state, source (aws-open-data | ncbi-prefetch), host`.

`pilot/summarize_pilot.py --runs RUNS.tsv --results OUT_DIR --out SUMMARY_DIR` writes `SUMMARY.md, summary_by_group.tsv, results_all.tsv, discordant.tsv, failures.tsv` from those small files.

On the HPCC (`/bigdata/stajichlab/jstajich/projects/MATPredict/`):
- `datasets/reads/sra_survey/pilot/` for the scripts; wave inputs in `results/sra_reads_typing/wave_inputs/*.tsv` (the tables the popsketch agent also read).
- Wave results used for the release: `results/2026-11-05_assembly_rescue_pilot/wave1_gates_final.tsv`, `wave1_results_all_final.tsv`, `wave2_gates_final.tsv`, `wave2_results_all_final.tsv` (gates PASS, PROVISIONAL, UNVALIDATED, HOLD).
- Reads themselves: not stored anywhere.

## 6. Quality control used on the reads

- **Contamination/mixture screen** (`SCREEN=1`): GmG classifies the same subsample against a fungal reference; `screen_reads.py` writes `<run>.screen.json` with `status` (`ok`, `mixed`, `genus_mismatch`, `borderline`, `no_call`, `screen_failed`). A `borderline` run is rescreened at `DEEP_SPOTS` (default 6M). The screen never stops typing; it flags. Calibration: `results/2026-11-02_contam_spikein/NOTE.md`. A contaminant at about 2% of pairs is detectable at 2M spots; 0.5% needs about 6M.
- **Truth where it exists:** BFD assemblies with the same BioSample, typed by `detect`, give an independent idiomorph per run for agreement checks.
- **Clone clustering** (`clone_clusters.py`, sourmash sketches via `SKETCH=1`): used to test whether a set is clonal.
- **Gates:** each species/run set gets a gate (PASS, PROVISIONAL, UNVALIDATED or HOLD); the exact rules are in `datasets/reads/sra_survey/pilot/` (`wave_driver.py`, `merge_screen.py`) and the `wave*_gates_final.tsv` tables, which I have not re-read for this handoff. Held species are shown but flagged (curator-agreed).

## 7. Assembly rescue (2026-11-05 pilot)

For runs the k-mer panel could not type (`none` or `low_depth`), I streamed the 2M-spot subsample, recruited reads by `diamond blastx` against 541 MATPredict reference proteins, assembled them with SPAdes, then ran `detect --taxid` on the assembly (`results/2026-11-05_assembly_rescue_pilot/`). 27 of 48 runs gave a call; 8 of 9 calls with a known answer agreed. It fails where the taxon has no reference proteins or recruits very few reads. Details and the problem runs: github.com/stajichlab/MATPredict/issues/275. The pilot script passes single-end runs to SPAdes as interleaved (a known bug).

## 8. Lessons that apply to a central PopSeq store

- Keep the run list and its filters in one table, with the exclusion reason per study, so choices can be audited.
- Never keep reads when only a subsample is needed; the cost per run is about 1 GB of transient download.
- Record per-run state and source, including `unavailable_at_ncbi`, so a missing run is not mistaken for a pipeline failure.
- Use windows spread along the run, not the head.
- Pin species names to one canonical taxonomy; synonyms silently change reference tiers.
- Screen for contamination before treating a call as a population statistic; clonal outbreak sets are not population frequencies.
- Nothing from these runs enters the MATPredict curated database without curator sign-off.

## 9. Constraints you should respect

- The curator (Jason Stajich) merges PRs and tags releases. The v0.8.1 run folders (`results/2026-11-05_v081_*`) and the worktree `.claude/worktrees/run-50bc950` are frozen.
- MATPredict HPCC runs stay under 60 CPUs.
- Do not edit MATPredict files on behalf of this handoff; ask the curator.
