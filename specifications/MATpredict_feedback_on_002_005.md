# MATPredict feedback on specs 002, 004 and 005

Written 2026-11-06 by the MATPredict session (Claude), at the curator's request, after reading the revised 002, 004 and the new 005 and 006. Nothing here is approved; the curator decides. Points are in order of how much they could mislead a reader writing code against the spec.

## 1. Spec 002: the `windows` tier has no defined route (likely contradiction)

The tier table says `windows` is the default subsample and costs about 1 GB. Section "Tiers `head` and `sample`" then says ENA HTTPS comes first "because it is the only route that can stop early", that `prefetch` "has no partial mode", and falls back to `fasterq-dump` with a read limit. Those routes serve `head`. They do not serve `windows`.

What MATPredict does for `windows` (`genotype_mat_one_accession.sh`): download the whole `.sra` archive (AWS Open Data, then `prefetch`), then run `fastq-dump -N <start> -X <end> --split-spot` once per window on the local `.sra`. Windows are spots numbered in the archive, so the method needs the `.sra`; I did not do windows over ENA gzipped FASTQ, and a gzip stream has no random access, so I do not believe it can be done there without reading the whole file.

Consequences for the spec:
- `windows` transfers the **whole run**, not about 1 GB. The ~1 GB figure is the extracted subsample (2M spots), or roughly MATPredict's own per-task transfer. Over the 4,211 completed AWS runs in my records the mean `.sra` was 1.5 GB (6,343 GB total), so budget transfer by run size, not by tier.
- Therefore `windows` and `full` share a transfer cost; `windows` saves only conversion time and storage. That belongs in the tier table next to "Cached?".
- The section name says `sample` where the table says `stride`/`windows`; one name should be used.

## 2. Spec 002: "Downloads go to `/bigdata`. Never to node-local `/scratch`"

This is right for the cached `full` tier and for a resumable download. It is not how MATPredict's uncached tiers work: their `.sra` and subsample FASTQ go to `$SCRATCH` and are deleted when the task ends (trap on exit), which the tier table also says. The rule, as worded, reads as if it contradicts the table. Suggest scoping it: "`full` and any resumable download go to `/bigdata`; transient tiers that are deleted on exit may use `$SCRATCH`." I do not know whether `$SCRATCH` on the HPCC is node-local or shared (Slurm sets it); MATPredict's tasks require it to be set and do not depend on it surviving the job.

## 3. Spec 004: stale rate in the adoption table

The first table row for `unavailable_at_ncbi` still says "about 1 run in 100". Spec 002 already carries the corrected 0.55 percent (about 1 in 180, 25 of 4,548). Please align 004's row so one number is quoted.

## 4. Spec 005: one fact that may help section 6

MATPredict's `reads-type` takes one read file per strain, and single-end reads work, so an inventory entry that resolves a strain to a single FASTQ (or a mate pair) is enough to run a MAT profile on local data. Nothing in MATPredict has been run on the Population_Genomics subprojects, and I make no claim about those data.

## Not addressed
I have no measurement for the concurrency cap (the `20` in spec 002 is the PopSeq agent's number). I have not reviewed spec 006's metagenome design beyond noting that its 0.5 percent contaminant remark matches my screening calibration (about 2 percent detectable at 2M spots, 0.5 percent needs about 6M).
