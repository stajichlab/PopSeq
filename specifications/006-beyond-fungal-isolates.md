# Spec 006: other organism groups, metagenomes, and MAGs

Status: **draft.** Design note, not an implementation. Written 2026-11-06 with the ENA
numbers below measured the same day.

Two extensions are wanted. They are not the same difficulty, and conflating them would
produce a schema that serves neither.

1. **Other organism groups** (bacteria, others). Mostly a scale problem.
2. **Metagenomes and MAGs**, including uncultured organisms. A **unit of analysis** problem.

## Part 1: other organism groups

### What already generalises

Almost everything. The only fungus-specific thing in spec 001 is the number `4751` in the
top-level query. Replace it and the same machinery runs:

- `tax_tree(<taxid>)` is a general ENA construct.
- The run table fields are organism-independent.
- The NCBI datasets reference-genome lookup is organism-independent.
- Taxid-first naming matters **more** for bacteria, where reclassification is routine.
- The four study-level exclusion categories hold.

So PopSeq should not be renamed or forked for bacteria. The clade is a parameter, and
`projects/Fungi/PopSeq` should eventually become `projects/PopSeq` with a `Fungi` catalog
inside it. **This is a rename, so spec 005 section 6 and AGENTS.md section 1 apply:**
announce, deprecate, do not move paths under people's feet.

### What does not generalise: scale

Fungi have 245,822 Illumina WGS runs in total (measured 2026-11-06). Single bacterial
species exceed that. *Salmonella enterica*, *Escherichia coli* and *Mycobacterium
tuberculosis* each have run counts in the hundreds of thousands, dominated by public-health
surveillance of near-clonal populations.

Three consequences:

- **A single `limit=0` ENA query stops being reasonable.** Page by `first_public` year, as
  GmG already does for the fungal-wide survey.
- **One CSV per species stops being reasonable.** A 500,000-row gzipped CSV is still
  readable, but it is no longer reviewable as a diff, which is the argument spec 001 rests
  on. At that size the catalog should shard by year or by study, or that species' table
  should be declared derived and the diff requirement dropped for it explicitly.
- **The clonality problem becomes the normal case, not the exception.** Fungal *Candidozyma
  auris* already shows it at about 25,000 runs. For bacteria it is most of the archive.
  Spec 001 already says the catalog records the run count and does not claim independence.
  That rule is what makes bacterial scale survivable.

### Genome size for taxa with no reference

`reference_genome()` fetches the NCBI reference assembly. Many taxa have none, and this is
the common case for the uncultured taxa in part 2. The fallback must be recorded, never
silent:

```yaml
assumed_genome_mb: 38.0
assumed_genome_mb_source: INHERITED from taxid 5141 (Neurospora crassa), nearest relative
                          with a reference assembly
assumed_genome_mb_confidence: low
```

A coverage number computed from an inherited size is not the same measurement as one from
the organism's own assembly, and any consumer filtering on coverage must be able to see
which it has.

## Part 2: metagenomes and MAGs

### The current filter throws this away on purpose

Spec 001's eligibility filter requires `library_source="GENOMIC"`. Measured 2026-11-06:

| Query | Runs |
|---|---|
| `tax_tree(4751) AND library_source="GENOMIC"` AND Illumina WGS | 245,822 |
| `tax_tree(4751) AND library_source="METAGENOMIC"` | **22,555** |

That second pool is about 9 percent again on top of the first, and it is currently invisible.

### A metagenome run is not a strain

This is the whole difficulty. Every entity in specs 001 and 003 assumes one run maps to one
isolate. A metagenome run contains an unknown number of organisms at unknown abundance, and
the strain-level variation worth studying is **inside** the run, not between runs.

So a metagenome cannot be a row in the runs table with a `sample_accession` standing for a
strain. It needs its own entity and its own eligibility question, which is not "is this
isolate well covered?" but "does this sample contain enough of the target taxon to say
anything?"

### Proposed: three entities, not one

| Entity | Unit | Key | Lives in |
|---|---|---|---|
| run | one sequencing run of one isolate | `run_accession` | `catalog/runs/` (exists) |
| metagenome | one environmental sample | `run_accession` | `catalog/metagenomes/` (new) |
| MAG | one genome binned from one metagenome | assembly accession | `catalog/mags/` (new) |

A MAG row carries `derived_from_sample` and `derived_from_run`, which is what makes the link
back to raw reads possible. That link is the point: a MAG gives a reference, and the reads it
came from give the within-sample variation.

### How to find MAGs at ENA, and why the flag is not enough

ENA's `result=assembly` exposes `assembly_type`. Measured for fungi on 2026-11-06:

| `assembly_type` | Count |
|---|---|
| *(empty)* | 25,227 |
| `clone or isolate` | 2,584 |
| `Metagenome-Assembled Genome (MAG)` | **490** |

The flag works, but **89 percent of fungal assemblies carry no `assembly_type` at all**, so
filtering on it finds MAGs that were labelled and misses everything else. It is a starting
point, not a census. NCBI BioSample attributes (`metagenome-assembled genome`,
`derived from`) are a second source and should be cross-checked, not assumed to agree.

### Uncultured taxa do have taxids

Checked 2026-11-06: `Ascomycota sp.` resolves to **taxid 1709939, rank species**. So the
taxid-first naming in spec 001 already handles these, and nothing new is needed for storage.

What breaks is everything downstream that assumes a taxid means one organism:

- **No reference assembly**, so `assumed_genome_mb` must be inherited and flagged (above).
- **`Ascomycota sp.` is not a species in any useful sense.** It is a bucket. Two samples
  under it may be unrelated at phylum level. Any panel built across one of these taxids must
  establish relatedness from the data, not from the name.
- A name like `Ascomycota sp. Antarctolithica obscura` carries real information in its
  **suffix**, placed there by the submitter, not in the taxid. That suffix belongs in
  `curated_*` fields, because nothing automatic will preserve it.

### What "eligible" means for a metagenome

A different question needs a different field. Proposed, and **entirely unvalidated**:

```
target_taxon_fraction   estimated share of reads from the target clade
target_est_coverage     target_taxon_fraction * base_count / assumed_genome_mb
```

Both need a classifier run against the reads, which means a metagenome cannot be judged from
metadata alone the way an isolate run can. That is a real cost: the catalog can list
isolate runs without fetching anything, but it cannot do the same for metagenomes.

GmG's existing tooling is the obvious way to produce these numbers, and GmG has already
measured classification accuracy against coverage on 100 runs
(`reports/validation/reads_sra.md`). **That measurement is for whole genomes of known
isolates, not for a target inside a mixture, so it does not transfer directly.** Whether a
0.5 percent member of a community can be called at all is open, and MATPredict's
contamination spike-in calibration (a contaminant at about 2 percent detectable at 2M spots,
0.5 percent needing about 6M) is the closest existing evidence.

## Recommended order

1. Finish the isolate catalog. It works and nothing above is needed for it.
2. Generalise the clade parameter and rename to `projects/PopSeq`. Cheap, and it stops the
   fungal assumption hardening.
3. Add the `metagenome` entity as **metadata only**, with no eligibility claim. Listing
   22,555 fungal metagenome runs with their study and environment is useful on its own and
   costs one query.
4. Add MAGs, keyed on `derived_from_run`, starting from the 490 flagged and accepting that
   the census is incomplete.
5. Only then attempt `target_taxon_fraction`, because it requires fetching and classifying
   reads and every number in it is currently a guess.

## Open questions

1. At what row count does a species table stop being reviewable as a diff? Spec 001's whole
   argument for CSV depends on an answer, and nobody has picked one.
2. Does the lab actually want bacterial data, or is this future-proofing? The schema work is
   cheap; the archive-scale work is not, and the two should not be bundled.
3. For an uncultured taxon, what is the unit a panel is built over, if not a strain?

## Changes

- 2026-11-06: first draft. ENA counts in this document were measured the same day.
