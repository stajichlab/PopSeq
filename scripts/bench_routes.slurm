#!/bin/bash
# Benchmark 7.1 (issue #12): one (run, route) attempt per array task.
#
#   export POPSEQ_ROOT=/path/to/PopSeq BENCH=/bigdata/.../bench/2026-11-07_d1_routes
#   sbatch --export=ALL --array=0-43%4 scripts/bench_routes.slurm
#
# Task i takes run i/4 of the plan and route (i + i/4) mod 4, so the four
# tasks that run together (%4) are usually one of each route, and every route
# sees a mix of times of day. The 11th plan row is the control run, a known
# unavailable_at_ncbi accession, which exercises the failure path.
#
# Resources: 8 cores for parallel-fastq-dump. Time is set by the slowest
# route: ENA over HTTPS was 0.4 MB/s in the only earlier measurement, and the
# largest run is about 10 GB. This is the first measurement, so the request
# was sized from the smoke test (tasks 0-3: peak RSS 1.4 GB, 10.4 min for 0.54 GB
# over ENA HTTPS at 0.88 MB/s). The 10 GB run over ENA HTTPS needs about 3.2 h.
#SBATCH -J d1routes
#SBATCH -p intel,batch,epyc
#SBATCH -c 8
#SBATCH --mem 4G
#SBATCH -t 06:00:00
#SBATCH -o logs/d1routes.%A_%a.log

set -euo pipefail
: "${POPSEQ_ROOT:?set POPSEQ_ROOT to the PopSeq checkout}"
: "${BENCH:?set BENCH to the benchmark work directory on /bigdata}"
: "${SCRATCH:?SCRATCH is not set; run under SLURM}"

module load sratoolkit/3.2.0
module load parallel-fastq-dump/0.6.7
module load aspera/4.2.4.265
export ASPERA_KEY=/opt/linux/rocky/8.x/x86_64/pkgs/aspera/4.2.4.265/etc/asperaweb_id_dsa.openssh

PLAN="$POPSEQ_ROOT/benchmarks/2026-11-07_d1_runs.csv"
ROUTES=(aws_odp ncbi_prefetch ena_aspera ena_https)
i=${SLURM_ARRAY_TASK_ID:?}
mapfile -t RUNS < <(tail -n +2 "$PLAN" | cut -d, -f1)
r=$(( i / 4 ))
RUN=${RUNS[$r]}
ROUTE=${ROUTES[$(( (i + r) % 4 ))]}
echo "task $i: $RUN $ROUTE on $(hostname)"

"$HOME/.local/bin/uv" run --project "$POPSEQ_ROOT" python \
    "$POPSEQ_ROOT/scripts/popseq_bench_routes.py" attempt \
    --plan "$PLAN" --run "$RUN" --route "$ROUTE" --stage all \
    --workdir "$BENCH/work" \
    --detail "$POPSEQ_ROOT/benchmarks/2026-11-07_route_comparison.csv" \
    --downloads "$POPSEQ_ROOT/benchmarks/downloads.csv" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --timeout "${BENCH_TIMEOUT:-8400}" </dev/null
# --timeout is per transfer (per mate for ENA). 2 x 8400 s stays inside the 6 h
# job limit, so a slow transfer is recorded as transfer_failed ("timeout"), not
# lost to a SLURM kill that would leave no row in downloads.csv.
