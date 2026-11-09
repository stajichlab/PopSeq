"""Fixture tests for the snapshot, species summary and species layout scripts. No network."""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))

import popseq_snapshot as snap  # noqa: E402
import popseq_species_layout as lay  # noqa: E402
import popseq_species_summary as ss  # noqa: E402


def test_snapshot_query_is_the_catalog_scope():
    q = snap.query(4751, "DNBSEQ")
    assert q == ('tax_tree(4751) AND library_strategy="WGS" AND library_source="GENOMIC"'
                 ' AND instrument_platform="DNBSEQ"')
    assert snap.FIELDS[0] == "run_accession" and len(snap.FIELDS) == 22


def test_species_summary_counts_split_samples_and_skips_duplicate_runs():
    taxmap = {"1": ("10", "Aus bus"), "2": ("10", "Aus bus")}
    rows = [
        dict(run_accession="R1", tax_id="1", sample_accession="S1;S2", study_accession="P1",
             base_count="3500000000", library_layout="PAIRED", first_public="2015-01-01"),
        dict(run_accession="R2", tax_id="2", sample_accession="S2", study_accession="P2",
             base_count="0", library_layout="SINGLE", first_public="2020-01-01"),
        # The same run again, as when two snapshots overlap: counted once.
        dict(run_accession="R2", tax_id="2", sample_accession="S2", study_accession="P2",
             base_count="0", library_layout="SINGLE", first_public="2020-01-01"),
        dict(run_accession="R3", tax_id="999", sample_accession="S9"),
    ]
    (d,) = ss.summarise(rows, taxmap)
    assert (d["runs"], d["isolates"], d["studies"]) == (2, 2, 2)
    assert (d["paired_runs"], d["single_runs"]) == (1, 1)
    assert d["total_gbase"] == 3.5 and d["est_download_tb"] == 0.001
    assert (d["first_year"], d["last_year"]) == ("2015", "2020")


def test_species_layout_single_only_isolates():
    base = dict(species_taxid="10", species_name="Aus bus", exclusion="")
    rows = [
        dict(base, sample_accession="S1", library_layout="PAIRED", ena_layout="PAIRED",
             est_depth_x="30"),
        # S2: labelled PAIRED but served as one file, so single-only, and a mismatch.
        dict(base, sample_accession="S2", library_layout="PAIRED", ena_layout="SINGLE",
             est_depth_x="12"),
        dict(base, sample_accession="S2", library_layout="SINGLE", ena_layout="SINGLE",
             est_depth_x="5"),
        # S3: a SINGLE run and a PAIRED run, so it counts as paired.
        dict(base, sample_accession="S3", library_layout="SINGLE", ena_layout="SINGLE"),
        dict(base, sample_accession="S3", library_layout="PAIRED", ena_layout="PAIRED"),
        dict(base, sample_accession="S4", library_layout="PAIRED", ena_layout=""),
        dict(base, sample_accession="S5", library_layout="SINGLE", ena_layout="SINGLE",
             exclusion="mutant"),
    ]
    (d,) = lay.build(rows)
    assert d["isolates"] == 4
    assert (d["isolates_paired"], d["isolates_single_only"], d["isolates_no_fastq"]) == (2, 1, 1)
    assert d["single_only_median_depth_x"] == "12.0" and d["single_only_ge_10x"] == 1
    assert d["runs_label_mismatch"] == 1
