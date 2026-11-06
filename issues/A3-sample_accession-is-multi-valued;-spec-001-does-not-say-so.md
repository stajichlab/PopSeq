# A3: sample_accession is multi-valued; spec 001 does not say so

**Track:** A   **Depends on:** nothing   **Parallel-safe:** yes

**Problem.** 97 rows of the ENA run table carry semicolon-joined `sample_accession`
values. The BioSample harvest treated each joined string as one accession and NCBI
rejected all 15 such batches. Those lines held 131 real accessions.

**Do.** State in spec 001 that `sample_accession` is multi-valued and semicolon
separated. Make every consumer split on `;`. Audit the other ENA fields for the same
property: `fastq_ftp`, `fastq_aspera` and `fastq_md5` are known to be multi-valued and
are already handled; others may not be.

**Acceptance.** Spec 001 documents it, a test covers a joined value, and re-running the
harvest over the full BioSample list produces zero format failures.
