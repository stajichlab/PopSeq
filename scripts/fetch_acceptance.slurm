#!/bin/bash
# Acceptance run for `popseq fetch` (issue #13). One job, sequential steps, each step checks its
# own exit code and its own postcondition. Runs from $SCRATCH on a compute node; the checkout it
# runs from keeps its own cache/ and benchmarks/downloads.csv, so the shared cache is not touched.
#
#   export POPSEQ_ROOT=~/projects/PopSeq_fetch_build ACC=$POPSEQ_ROOT/acceptance
#   sbatch --export=ALL scripts/fetch_acceptance.slurm
#
# Steps map to the acceptance list of issue #13:
#   tiers      one run fetched at each tier (full, head, windows, stride, coverage), verified, sidecar
#   unpaired   an NCBI-unpaired ERR run with --require-pairs arrives as ENA pairs (exit 0, layout PAIRED)
#   single     a single-end run is verified as one file
#   hook       --exec runs on the verified reads and the reads are gone afterwards
#   resume     a killed task leaves nothing complete-looking and a rerun finishes it
#   exitcodes  usage (2); exits 10 to 14 are covered by tests/test_fetch.py
# Resources: 4 cores for parallel-fastq-dump; peak memory UNMEASURED for fetch (the benchmark smoke
# test peaked at 1.4 GB for a 0.54 GB run). Time: the slowest step is one ENA HTTPS transfer of about
# 0.2 GB at about 0.9 MB/s (benchmark 7.1, measured under load).
#SBATCH -J fetchacc
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -t 02:00:00
#SBATCH -o logs/fetchacc.%j.log

set -uo pipefail
: "${POPSEQ_ROOT:?set POPSEQ_ROOT to the PopSeq checkout}"
: "${ACC:?set ACC to the acceptance work directory}"
: "${SCRATCH:?SCRATCH is not set; run under SLURM}"

module load sratoolkit/3.2.0
module load parallel-fastq-dump/0.6.7
UV="$HOME/.local/bin/uv"
# One fetch at a time, so the declared concurrency is 1 (the default 20 would be untrue here).
fetch() { "$UV" run --project "$POPSEQ_ROOT" popseq fetch "$@" --concurrency 1 </dev/null; }

RESULTS="$ACC/acceptance_results.tsv"
SUMMARY="$ACC/acceptance_summary.tsv"
rm -rf "$ACC/out" "$ACC/hook" "$ACC/logs"
mkdir -p "$ACC/out" "$ACC/hook" "$ACC/logs"
: > "$RESULTS"
printf "step\texpect_exit\tgot_exit\tcheck\tstatus\n" > "$SUMMARY"
FAIL=0

# step NAME EXPECT_EXIT CHECK_CMD -- popseq-fetch-args...
step() {
  local name=$1 expect=$2 check=$3; shift 4
  local t0=$SECONDS out rc status="ok" why=""
  out=$(fetch "$@" --results "$RESULTS" 2> "$ACC/logs/$name.err"); rc=$?
  echo "$out" > "$ACC/logs/$name.out"
  [ "$rc" = "$expect" ] || { status="FAIL"; why="exit $rc, wanted $expect"; }
  if [ "$status" = ok ] && [ -n "$check" ]; then
    why=$(eval "$check" 2>&1) || status="FAIL"
  fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$name" "$expect" "$rc" "$why" "$status" >> "$SUMMARY"
  [ "$status" = ok ] || FAIL=1
  echo "== $name: $status ($((SECONDS - t0)) s) $why"
}

AURIS=SRR36290343          # C. auris, 4.39 M spots, 1.3 Gb, 535 MB of FASTQ (benchmark 7.1)
side() { echo "$ACC/out/$1/$2.json"; }

# ---- tiers: one run at each tier, with a sidecar and a verified read count
step full 0 "test -s $POPSEQ_ROOT/cache/498019/$AURIS/$AURIS.json && grep -q '\"verified\": true' $POPSEQ_ROOT/cache/498019/$AURIS/$AURIS.json && echo sidecar verified" -- $AURIS --tier full
step head 0 "grep -q '\"reads_kept\": 100000' $(side head $AURIS) && echo 100000 pairs" -- $AURIS --tier head --spots 100000 --out "$ACC/out/head"
step windows 0 "grep -q '\"verify_scope\": \"subsample\"' $(side windows $AURIS) && echo 8 windows" -- $AURIS --tier windows --spots 400000 --out "$ACC/out/windows"
step stride 0 "grep -q '\"verify_scope\": \"full\"' $(side stride $AURIS) && echo full check before stride" -- $AURIS --tier stride --spots 400000 --out "$ACC/out/stride"
step coverage 0 "python3 - <<PY
import json
s = json.load(open('$(side coverage $AURIS)'))
r = s['coverage_reached']; t = s['coverage_target']
assert abs(r - t) / t < 0.05, (r, t)
print('coverage target', t, 'reached', round(r, 2), 'bases', s['bases_kept'])
PY" -- $AURIS --tier coverage --coverage 30 --out "$ACC/out/coverage"

# ---- unpaired at NCBI: pairs must come from ENA
step unpaired 0 "grep -q '\"route\": \"ena_https\"' $POPSEQ_ROOT/cache/4932/ERR1069479/ERR1069479.json && grep -q '\"layout\": \"PAIRED\"' $POPSEQ_ROOT/cache/4932/ERR1069479/ERR1069479.json && test -s $POPSEQ_ROOT/cache/4932/ERR1069479/ERR1069479_2.fastq.gz && echo pairs from ENA" -- ERR1069479 --tier full --require-pairs

# ---- single-end: verified as one file
step single 0 "grep -q '\"layout\": \"SINGLE\"' $POPSEQ_ROOT/cache/4932/SRR5628799/SRR5628799.json && test ! -e $POPSEQ_ROOT/cache/4932/SRR5628799/SRR5628799_2.fastq.gz && echo one file" -- SRR5628799 --tier full

# ---- hook: runs on verified reads; the reads are gone afterwards
cat > "$ACC/hook.sh" <<'HOOK'
#!/bin/bash
# $1 run, then the FASTQ paths. Count reads and write one line, as a sketcher would read them.
n=$(zcat "$2" | awk 'END{print NR/4}')
echo "$1 $# files, $n reads in the first" >> "$ACC/hook.log"
HOOK
chmod +x "$ACC/hook.sh"; : > "$ACC/hook.log"
step hook 0 "grep -q 'SRR36290343 3 files, 400000 reads' $ACC/hook.log && test \"\$(ls $ACC/hook | grep -c fastq)\" = 0 && echo reads gone, hook saw 400000" -- $AURIS --tier windows --spots 400000 --out "$ACC/hook" --exec "$ACC/hook.sh"

# ---- resume: kill a full fetch mid-transfer; nothing complete is left; a rerun finishes it
RES=SRR16988686            # 5.5 M spots, 963 MB of FASTQ
rm -rf "$POPSEQ_ROOT/cache/51922/$RES"
setsid "$UV" run --project "$POPSEQ_ROOT" popseq fetch $RES --tier full --concurrency 1 --results "$RESULTS" </dev/null > "$ACC/logs/resume_killed.out" 2>&1 &
PGID=$!
sleep 25
kill -9 -- -"$PGID" 2>/dev/null; wait "$PGID" 2>/dev/null
left=$(ls "$POPSEQ_ROOT/cache/51922/$RES" 2>/dev/null | grep -c -E 'json$|fastq.gz$')
printf "resume_killed\t-\t-\t%s complete-looking files after the kill\t%s\n" "$left" "$([ "$left" = 0 ] && echo ok || echo FAIL)" >> "$SUMMARY"
[ "$left" = 0 ] || FAIL=1
step resume 0 "grep -q '\"verified\": true' $POPSEQ_ROOT/cache/51922/$RES/$RES.json && echo finished after the kill" -- $RES --tier full

# ---- exit codes: usage, and a run that nobody has
step usage 2 "" -- $AURIS --tier coverage --out "$ACC/out/usage"
# Exit 10 (a run neither archive has) is covered by tests/test_fetch.py only: no real accession is
# known that NCBI and ENA both lack, and one must not be invented.

# The kill -9 above cannot clean up after itself, so one popseq_fetch_* directory is expected here.
echo "scratch dirs left over (1 expected, from the kill): $(ls "$SCRATCH" | grep -c popseq_fetch_)"
echo "summary:"; cat "$SUMMARY"
exit $FAIL
