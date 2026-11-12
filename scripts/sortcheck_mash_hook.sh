#!/bin/bash
# `popseq fetch --exec` hook for the read-order ground truth (scripts/sortcheck_truth.slurm).
# Called as: HOOK RUN FASTQ_1 [FASTQ_2], with POPSEQ_TIER set by popseq fetch.
# Writes one line to $SORTCHECK_OUT: run, tier, bases in the sample, and the genome size and
# k-mer coverage that `mash sketch -r -m 2` estimates from the reads. A windows sample of a
# sorted run reaches only part of the genome, so its estimated size falls well below the
# stride sample's. The reads are deleted by popseq fetch after this returns.
set -euo pipefail
: "${SORTCHECK_OUT:?}"
run=$1; shift
# The FASTQ may be gzipped: zcat -f passes plain files through unchanged.
bases=$(zcat -f "$@" | LC_ALL=C awk 'NR % 4 == 2 { n += length($0) } END { print n + 0 }')
est=$(zcat -f "$@" | mash sketch -r -m 2 -k 21 -s 1000 -o "${SCRATCH:?}/$run.$POPSEQ_TIER" - 2>&1 >/dev/null \
      | awk -F': ' '/Estimated genome size/ {g = $2} /Estimated coverage/ {c = $2} END {print g "\t" c}')
(
    flock 9
    printf '%s\t%s\t%s\t%s\n' "$run" "$POPSEQ_TIER" "$bases" "$est" >> "$SORTCHECK_OUT"
) 9>"$SORTCHECK_OUT.lock"
