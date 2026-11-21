#!/usr/bin/env bash
# Import the run-QC sources added 2026-11-21 (issues #37 and #49; spec 008). Re-running it reproduces their rows in
# catalog/qc/ and their rows of catalog/qc/sources.csv; the rows of other sources and every curated_* column are kept.
#
#   POPSEQ_ROOT=/path/to/PopSeq bash scripts/popseq_qc_sources_2026-11-21.sh
set -euo pipefail
cd "${POPSEQ_ROOT:?set POPSEQ_ROOT to the repository root}"
RUNS=catalog/runs/fungi_all_runs.csv.gz
MP_REF="MATPredict results/2026-11-10_afum_bafA_yunnan/ on the HPCC (species_id.tsv, species_id.py); stajichlab/PopSeq#37"
MP_TOOL="MATPredict species_id.py (solid-k-mer containment, k=21, of a 30x subsample sketch against 18 Aspergillus genomes + Af293)"

# 1. Issue #37 (a): A. fumigatus runs whose reads are another species (mismatch) or could not be placed (unverified_species).
uv run python scripts/popseq_qc_import.py --source-id matpredict_afum_species_id_2026-11-10 --format popseq_findings \
    --input catalog/qc/findings/matpredict_afum_species_id_2026-11-10.csv --runs "$RUNS" \
    --project MATPredict --tool "$MP_TOOL" --check species_match \
    --method "runs with mean all-against-all similarity below 0.95 in the MATPredict A. fumigatus store (4,849 sketches; store median 0.995) were compared with 18 Aspergillus reference genomes plus Af293; the fraction of the run's solid k-mers in each genome is converted to an identity; mismatch = nearest reference is another species; unverified_species = no reference placed the run" \
    --detection-limit "only the 18 Aspergillus references and Af293; a species outside that set (e.g. A. hiratsukae) cannot be named; found by MATPredict and transcribed from the issue, not re-measured by PopSeq" \
    --reference "$MP_REF" --date 2026-11-10 --exclude-statuses "mismatch;unverified_species" \
    --notes "transcribed from PopSeq#37; the runs are not relabelled (a run_curation.csv relabel is the curator's decision, spec 009); identity is in source_note"

# 2. Issue #37 (b): tiny read sets filed as WGS.
uv run python scripts/popseq_qc_import.py --source-id matpredict_afum_read_sets_2026-11-10 --format popseq_findings \
    --input catalog/qc/findings/matpredict_afum_read_sets_2026-11-10.csv --runs "$RUNS" \
    --project MATPredict --tool "$MP_TOOL" --check low_coverage \
    --method "runs of the MATPredict A. fumigatus store with too few reads for a genome-wide call (hundreds to 136k reads per run); study titles suggest targeted sequencing registered as library_strategy WGS" \
    --detection-limit "the A. fumigatus store only (one run per BioSample); found by MATPredict and transcribed from the issue, not re-measured by PopSeq" \
    --reference "$MP_REF" --date 2026-11-10 --exclude-statuses "tiny_read_set" \
    --notes "transcribed from PopSeq#37 (54 runs of 5 studies); 53 of them are also below_1x in popseq_low_coverage (run list of the issues), SRR33550994 is 1.016x there"

# 3. Issues #37 and #49: base_count / genome size below 1x, for the runs those issues name ONLY (owner decision: no catalog-wide
#    pass; add --catalog-wide instead of --run-list for one, which is not committed).
uv run python scripts/popseq_qc_checks.py low-coverage --runs "$RUNS" --run-list catalog/qc/findings/issue_named_runs.csv --below 1 --date 2026-11-21 \
    --notes "evaluated by PopSeq from catalog data for the runs named in stajichlab/PopSeq#37 and #49 only (no catalog-wide pass); est_depth_x and the runs/panel --min-depth filter use the same ratio"

# 4. Issue #49: archive bytes per base. The threshold 2.0 is a placeholder; only SRR065689 is examined (no catalog-wide archive sizes).
uv run python scripts/popseq_qc_checks.py archive-size --sizes catalog/qc/findings/archive_sizes.csv \
    --source-id gmg_archive_size_2026-11-21 --project GmG --max-bytes-per-base 2 --date 2026-11-21 \
    --reference "stajichlab/PopSeq#49; GmG P6 transfer log" \
    --notes "UNMEASURED: whether other runs exceed the threshold (no catalog-wide pass over archive sizes); the threshold 2.0 bytes per base is not calibrated. Owner judgement (PopSeq#49): old record, very short reads (about 76 bp per mate), exclude"
