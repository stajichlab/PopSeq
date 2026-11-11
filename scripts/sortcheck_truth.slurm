#!/bin/bash
# Ground truth for the read-order probe (scripts/popseq_sortcheck.py, PR #43): does a windows
# sample of this run reach less of the genome than a stride sample of the same size?
# One array task = one run: `popseq fetch --tier windows --spots N`, then `--tier stride
# --spots N`, the same N, each with scripts/sortcheck_mash_hook.sh as --exec.
#
#   D=~/projects/PopSeq/bench/2026-11-11_sortcheck   # truth_plan.csv here; results to truth.tsv
#   sbatch --array=0-$((K-1))%8 --export=ALL,D=$D,POPSEQ_ROOT=$D/wt scripts/sortcheck_truth.slurm
#
# Resources: tier stride converts the whole run (PR #42 measured a peak of 3.1 GB); 4 cores.
#SBATCH -J sortgt
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -t 02:00:00

set -uo pipefail
: "${POPSEQ_ROOT:?}"
: "${D:?}"
: "${SCRATCH:?SCRATCH is not set; run under SLURM}"
module load sratoolkit/3.2.0
module load parallel-fastq-dump/0.6.7
module load mash/2.3
UV="$HOME/.local/bin/uv"
export SORTCHECK_OUT="$D/truth.tsv"

i=${SLURM_ARRAY_TASK_ID:?}
line=$(tail -n +2 "$D/truth_plan.csv" | sed -n "$((i + 1))p")
RUN=$(echo "$line" | cut -d, -f1)
SPOTS=$(echo "$line" | cut -d, -f2)
[ -n "$RUN" ] || { echo "no run for task $i"; exit 2; }
for tier in windows stride; do
    "$UV" run --project "$POPSEQ_ROOT" popseq fetch "$RUN" --tier "$tier" --spots "$SPOTS" \
        --out "$SCRATCH/out.$tier" --threads 4 --concurrency 8 \
        --exec "bash $POPSEQ_ROOT/scripts/sortcheck_mash_hook.sh" \
        --results "$D/truth_fetch.tsv" </dev/null
    echo "task $i $RUN $tier exit $?"
done
