# GmG / popsketch: population-set downloads, data location, and what centralizing needs

Written 2026-11-06 by the GmG/popsketch Claude session for Jason Stajich. Repo: stajichlab/GmG (public). Proposals below are for the curator to decide; nothing here is approved.

## What was done
popsketch audits a set of read files from FracMinHash sketches (k=21, scaled 200): QC flags, distances, clusters, PCA, tree, report. To test it on real populations, SRA run lists were selected from MATPredict's survey tables (`MATPredict/results/sra_reads_typing/wave_inputs/`: `species_summary.tsv`, `fungi_wgs_illumina.tsv`, `fungi_meta.tsv`).

- Filters: paired Illumina, coverage >= 15x, one run per biosample, spread across studies, mutant/evolved/cross/knockout titles removed. C. auris and S. cerevisiae skipped (too large).
- Seven species, up to 300 runs each, 2,023 selected: Z. tritici, C. neoformans, C. beticola, P. nodorum, B. cinerea, N. discreta (223, all that pass), P. oryzae (= Magnaporthe oryzae).
- Fetch: `parallel-fastq-dump` to scratch, every Nth read pair kept (about 1 Gb; old SRA files are position-sorted), sketched with `popsketch sketch`, reads deleted. Only about 0.1-0.4 MB per sample is kept.
- Analysis: `popsketch run --from-sketches <slug>/sketches`.

## Where it is
- Data: `/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra/<slug>/` with `runs.tsv`, `runs_ftp.tsv`, `sketches/`. Slugs: z_tritici, c_neoformans, c_beticola, p_nodorum, b_cinerea, n_discreta, p_oryzae. Logs: `popsra/logs`. Failed-run list: `popsra/failed_runs.txt`.
- Binary: `GmG_popgen/popsketch_build/target/release/popsketch` (built from branch popsketch-sketch-cmd).
- Code: `analysis/popgen_sra/01_select_runs.py`, `02_ena_urls.py`, `03_fetch_sketch.sbatch`, `04_run_species.sh` in GmG, in PR #156 (open, unmerged).

## State (2026-11-06)
Sketches done: z_tritici 300, c_neoformans 300, c_beticola 300, p_nodorum 300, b_cinerea 300, n_discreta 223, p_oryzae 279 of 300 (2,002 of 2,023). The queue is empty, so the 21 missing P. oryzae runs failed again on rerun; cause not yet examined.
- About 28 downloads failed the first time (`parallel-fastq-dump` exit 64, mostly P. oryzae). Looked transient; one run succeeded on a manual retry.
- Results so far (before the full runs): P. nodorum (265 of 300 in the analysis) and C. beticola show no clear clusters at any `--min-floors` 1-8. The other species are not yet analysed.

## What centralizing would need (decisions for the curator)
1. Location: shared group space or a personal folder. The data is currently under a home projects folder by instruction.
2. A catalog table, one row per run: accession, biosample, study, canonical species, strain, coverage, country, fetch status. The `runs.tsv` files are a start, with no shared schema.
3. Whether to keep reads. Sketches cannot serve MAT typing, which needs reads (MATPredict's `datasets/reads/sra_survey` streams 2M-spot subsamples and keeps none). Keeping reads for about 2,000 runs would be hundreds of GB or more; the alternative is re-fetching.
4. A retrying downloader with a concurrency cap (failures seen at about 35 concurrent).
5. Canonical species names (Magnaporthe oryzae is now Pyricularia oryzae).
6. Optional: write a popsketch sketch in MATPredict's fetch step before reads are deleted.
