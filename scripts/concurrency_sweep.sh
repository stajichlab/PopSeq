#!/bin/bash
# Submit the download concurrency sweep (issue #27): one SLURM array per level, chained so that
# levels never overlap. Each level fetches the first K runs of the pool with at most N at once.
#
#   W=~/projects/PopSeq/bench/2026-11-10_concurrency POPSEQ_ROOT=$W/wt \
#       bash scripts/concurrency_sweep.sh
#
# Levels and task counts: N=1 runs 8 tasks (serial, the baseline); every other level runs
# max(2N, 20) tasks, capped at the pool size, so each level has at least two waves at its cap.
# It also records how busy the queue was when each level starts (squeue counts), because the
# rates of benchmark 7.1 were measured under load and this sweep should say whether it was.
set -euo pipefail
: "${POPSEQ_ROOT:?}"
: "${W:?}"
LEVELS=${LEVELS:-"1 5 10 20 40"}
POOL=$(($(wc -l < "$W/pool.csv") - 1))
mkdir -p "$W/logs" "$W/sidecars"
[ -f "$W/tasks.tsv" ] || printf 'level\trun\ttaxid\tnode\tstart_epoch\tend_epoch\texit\tbytes\tseconds_fetch\n' > "$W/tasks.tsv"

dep=""
for N in $LEVELS; do
    if [ "$N" -eq 1 ]; then K=8; else K=$((2 * N > 20 ? 2 * N : 20)); fi
    [ "$K" -le "$POOL" ] || K=$POOL
    # The load snapshot runs as the first thing of each level, inside a tiny job.
    snap=$(sbatch --parsable -p short -c 1 --mem 1G -t 5 ${dep:+--dependency=afterany:$dep} \
        -o "$W/logs/load.L$N.log" \
        --wrap "date -u +%FT%TZ; echo running_all \$(squeue -h -t R | wc -l); echo running_mine \$(squeue -h -t R -u $USER | wc -l); echo pending_all \$(squeue -h -t PD | wc -l)")
    job=$(sbatch --parsable --array=0-$((K - 1))%$N --dependency=afterany:$snap \
        --export=ALL,LEVEL=$N,W=$W,POPSEQ_ROOT=$POPSEQ_ROOT \
        -o "$W/logs/conc.L$N.%A_%a.log" "$POPSEQ_ROOT/scripts/concurrency_sweep.slurm")
    echo "level $N: $K tasks, array $job (after $snap)"
    dep=$job
done
