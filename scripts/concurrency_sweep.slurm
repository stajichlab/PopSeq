#!/bin/bash
# Download concurrency sweep, one level (issue #27). One array task = one `popseq fetch --tier
# full` of one run from the pool. The array's `%N` throttle sets the level N; the summary measures
# the concurrency each fetch actually saw (scripts/popseq_concurrency.py).
#
# POPSEQ_ROOT is a throwaway checkout: its cache/ and benchmarks/downloads.csv receive the
# downloads, so the shared cache is never touched. Each task deletes its cached run after reading
# the sidecar, so the next level downloads it again instead of skipping it as already verified.
#
#   W=~/projects/PopSeq/bench/2026-11-10_concurrency   # work dir: pool.csv, tasks.tsv, logs/
#   POPSEQ_ROOT=$W/wt
#   for N in 1 5 10 20 40; do ... sbatch --array=0-$((K-1))%$N --export=ALL,LEVEL=$N ...; done
#   (scripts/concurrency_sweep.sh does this, one level after another)
#
# Resources: 4 cores for the conversion; 8 GB (the acceptance run peaked well below that for a
# 0.54 GB run; peak memory for full-tier runs of 0.8-1.6 Gbases is UNMEASURED).
#SBATCH -J conc
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -t 01:30:00

set -uo pipefail
: "${POPSEQ_ROOT:?set POPSEQ_ROOT to the throwaway checkout}"
: "${W:?set W to the sweep work directory}"
: "${LEVEL:?set LEVEL to the array throttle N}"
: "${SCRATCH:?SCRATCH is not set; run under SLURM}"

module load sratoolkit/3.2.0
module load parallel-fastq-dump/0.6.7
UV="$HOME/.local/bin/uv"

i=${SLURM_ARRAY_TASK_ID:?}
line=$(tail -n +2 "$W/pool.csv" | sed -n "$((i + 1))p")
RUN=${line%%,*}
TAXID=$(echo "$line" | cut -d, -f2)
[ -n "$RUN" ] || { echo "no run for task $i"; exit 2; }

start=$(date +%s.%N)
"$UV" run --project "$POPSEQ_ROOT" popseq fetch "$RUN" --tier full --threads 4 \
    --concurrency "$LEVEL" --plan "$W/plan.csv" --results "$W/results.tsv" </dev/null
rc=$?
end=$(date +%s.%N)

DIR="$POPSEQ_ROOT/cache/$TAXID/$RUN"
read -r bytes sec < <("$UV" run --project "$POPSEQ_ROOT" python -c '
import json, sys
try:
    d = json.load(open(sys.argv[1]))
    print(d.get("bytes_total") or "", d.get("seconds_fetch") or "")
except OSError:
    print(" ")
' "$DIR/$RUN.json" 2>/dev/null) || true
[ -f "$DIR/$RUN.json" ] && cp "$DIR/$RUN.json" "$W/sidecars/L${LEVEL}_$RUN.json"
rm -rf "${DIR:?}"

# One line per task. flock: tasks of one level append at the same time.
(
    flock 9
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$LEVEL" "$RUN" "$TAXID" "$(hostname -s)" \
        "$start" "$end" "$rc" "${bytes:-}" "${sec:-}" >> "$W/tasks.tsv"
) 9>"$W/tasks.lock"
echo "level $LEVEL task $i $RUN exit $rc bytes ${bytes:-} fetch_s ${sec:-}"
