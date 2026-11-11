"""Fixture tests for scripts/popseq_panel_set.py. No network, no catalog contents."""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))

import popseq_panel_set as ps  # noqa: E402


def test_candidates_keep_species_at_threshold_largest_first():
    rows = [dict(species_taxid="1", isolates_paired="9"),
            dict(species_taxid="2", isolates_paired="10"),
            dict(species_taxid="3", isolates_paired="400"),
            dict(species_taxid="4", isolates_paired="")]
    assert ps.candidates(rows, 10) == [3, 2]


def test_sample_id_falls_back_to_run_for_joined_or_empty_samples():
    assert ps.sample_id(dict(run_accession="SRR1", sample_accession="SAMN1")) == "SAMN1"
    assert ps.sample_id(dict(run_accession="SRR2", sample_accession="SAMN1;SAMN2")) == "SRR2"
    assert ps.sample_id(dict(run_accession="SRR3", sample_accession="")) == "SRR3"


def test_rows_for_joins_full_row_and_keeps_mates_only():
    run = dict(run_accession="SRR9", sample_accession="SAMN9", study_accession="P1",
               strain="X1", est_depth_x="31.5", base_count="900",
               fastq_ftp="ftp.sra.ebi.ac.uk/a/SRR9.fastq.gz;"
                         "ftp.sra.ebi.ac.uk/a/SRR9_1.fastq.gz;"
                         "ftp.sra.ebi.ac.uk/a/SRR9_2.fastq.gz")
    panel = dict(species_taxid=10, species_name="Aus bus", runs=[run])
    full = {"SRR9": dict(tax_id="11", scientific_name="Aus bus var. c", isolate="i9",
                         instrument_model="Illumina NovaSeq 6000")}
    (s,), (m,) = ps.rows_for(panel, full)
    assert s == {"sample": "SAMN9",
                 "fastq_1": "https://ftp.sra.ebi.ac.uk/a/SRR9_1.fastq.gz",
                 "fastq_2": "https://ftp.sra.ebi.ac.uk/a/SRR9_2.fastq.gz"}
    assert (m["sample"], m["biosample"], m["taxid"], m["species_taxid"]) == \
        ("SAMN9", "SAMN9", "11", 10)
    assert (m["scientific_name"], m["isolate"], m["strain"]) == ("Aus bus var. c", "i9", "X1")
    assert list(m) == ps.METADATA
