# GmG reply to spec 004 and the revised specs

Written 2026-11-06 by the GmG/popsketch Claude session, at the curator's request, after reading 001-006, `004-response-to-handoffs.md` and `MATpredict_feedback_on_002_005.md`. Nothing here is approved; the curator decides.

The revised 002 and 004 carry the corrections to my first handoff correctly: the 21 permanent 404s, the 7 transient failures, and that concurrency was not shown to be the cause. No objection to anything in them. What follows is additions and two small fixes.

## 1. Spec 001: the `/3.5` size divisor is probably too high for older runs (measured, one run)

Spec 001 gives `est_gz_bytes = base_count / 3.5` and says the divisor is unmeasured. I have one data point from a test fetch of SRR1713175 (*Cercospora beticola*, an old run), converted with `parallel-fastq-dump --split-files --gzip`:

| Quantity | Value |
|---|---|
| `base_count` (ENA) | 1,110,000,000 |
| `SRR1713175_1.fastq.gz` | 547,286,252 bytes |
| `SRR1713175_2.fastq.gz` | 542,897,857 bytes |
| gz bytes per base | about 0.98 (spec assumes about 0.29) |

For this run the spec would underestimate the download by about 3.4 times. This is one run. Newer runs with binned quality scores probably compress better, so I am not proposing a replacement value. I am saying benchmark 7.2 should include old runs, and that the divisor should not be treated as one constant for all years. The 42 files already on disk in TE_fungi_copy may all be recent.

## 2. README row for `GmG_spec_handoff.md` is stale

It says "concurrency failures". It should say 21 permanent 404s at NCBI (*P. oryzae*, ERR runs from four ENA studies) and 7 transient failures that succeeded on rerun.

## 3. Spec 005 section 7: sketches as a second check that a local file is the archived run

The spec proposes comparing read counts for `archived_local`, and strain name plus read count to recover a lost accession. A popsketch sketch of the local file compared with a sketch of the archived run is an additional check. It would catch a swapped or contaminated file that happens to have the same read count. It would not catch a trimmed or subsampled copy that is otherwise the same data, so it complements the read-count check and does not replace it. Not tested; I would want to try it on a few known pairs before the spec relies on it.

## 4. Spec 006: popsketch cannot be applied to metagenomes as it stands

Spec 006 names GmG's tooling as the way to produce `target_taxon_fraction`. Part of the machinery exists: the index screen can name a second taxon and give its fraction. But popsketch's QC is cohort-relative and assumes one isolate per run. On metagenomes it would flag nearly every sample as mixed or unrelated. Using it there needs new design, not a different input. Separately, the statement that GmG "already pages by `first_public` year" is not something this session did; it may refer to earlier GmG work I have not read.

## 5. Correction to the size of a sketch (quoted in spec 004 question 6)

My handoff said a sketch is 0.1-0.4 MB. That was wrong. Measured on the 300 *C. neoformans* sketches: 2.7 GB in 1,201 files, about 8 MB per sample: the whole sketch about 3 MB (3.4 MB for DRR829879) and each of the two split halves about 2 MB. The "sketch during fetch" option in spec 004 is still cheap compared with the reads, but at about 8 MB per sample, 2,000 runs is about 16 GB, not about 0.5 GB. Only the whole sketch and the halves of the 20 samples the noise floor uses are needed for an analysis, which is about 3 MB per sample.

## 6. Agreement with the MATPredict feedback on 002

`stride` and `windows` both transfer the whole run. My pipeline downloads the full run through `parallel-fastq-dump` and keeps about 1 Gb of bases. Only the stored sketch is small. Spec 002 already says this for `stride`.

## Current state of the GmG data

- Sketches: 2,002 of 2,023 selected runs. The 21 missing are the *P. oryzae* 404 runs, **skipped by decision of the curator** on 2026-11-06 and not to be retried here.
- Analysis jobs for six species are running (the seventh, *P. nodorum*, is done). Results will be added to the GmG repo documentation, not here.
- Nothing about where the data lives has changed: `/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsra/`.

## Changes

- 2026-11-06: first version.
