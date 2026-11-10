"""Specification 009: the NCBI assembly table (scripts/popseq_assemblies.py) and the curated label import
(scripts/popseq_labels_import.py). No network."""
from __future__ import annotations

import csv
import gzip
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))

import popseq_assemblies as asm
import popseq_labels_import as li

SUMMARY_COLS = ["assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid", "species_taxid",
                "organism_name", "infraspecific_name", "isolate", "version_status", "assembly_level", "release_type",
                "genome_rep", "seq_rel_date", "asm_name", "asm_submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path",
                "excluded_from_refseq", "relation_to_type_material", "asm_not_live_date", "assembly_type", "group",
                "genome_size"]


def _summary(path, rows):
    with gzip.open(path, "wt", newline="") as fh:
        fh.write("##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the columns\n")
        fh.write("#" + "\t".join(SUMMARY_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "na") for c in SUMMARY_COLS) + "\n")


def _csv_gz(path, header, rows):
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(header); w.writerows(rows)


def _read_gz(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def _assemblies(tmp_path):
    _summary(tmp_path / "gb.txt.gz", [
        {"assembly_accession": "GCA_1.1", "biosample": "SAMN1", "bioproject": "PRJ1", "taxid": "5507", "species_taxid": "5507",
             "organism_name": "Fusarium oxysporum f. sp. lycopersici 4287", "infraspecific_name": "strain=4287",
             "assembly_level": "Chromosome", "genome_rep": "Full", "seq_rel_date": "2012/05/01", "asm_name": "FO2",
             "gbrs_paired_asm": "GCF_1.1", "refseq_category": "na"},
        {"assembly_accession": "GCA_2.1", "biosample": "SAMN2", "taxid": "5128", "species_taxid": "5128", "organism_name": "Fusarium proliferatum"},
        {"assembly_accession": "GCA_3.1", "biosample": "SAMN3", "taxid": "5507", "species_taxid": "5507", "organism_name": "Fusarium oxysporum"},
        {"assembly_accession": "GCA_4.1", "biosample": "na", "taxid": "1", "species_taxid": "1", "organism_name": "Fusarium sp. X"},
    ])
    _summary(tmp_path / "rs.txt.gz", [{"assembly_accession": "GCF_1.1", "gbrs_paired_asm": "GCA_1.1", "refseq_category": "reference genome"}])
    _csv_gz(tmp_path / "tax.csv.gz", ["taxid", "name", "rank", "species_taxid", "species_name", "resolution"],
            [["5507", "Fusarium oxysporum", "species", "5507", "Fusarium oxysporum", "species"],
             ["5128", "Fusarium proliferatum", "species", "5128", "Fusarium proliferatum", "species"]])
    return asm.main(["build", "--genbank", str(tmp_path / "gb.txt.gz"), "--refseq", str(tmp_path / "rs.txt.gz"),
                     "--taxonomy", str(tmp_path / "tax.csv.gz"), "--out", str(tmp_path / "assemblies.csv.gz")])


def test_assembly_table_pairs_refseq_and_blanks_na(tmp_path):
    assert _assemblies(tmp_path) == 0
    rows = {r["assembly_accession"]: r for r in _read_gz(tmp_path / "assemblies.csv.gz")}
    a = rows["GCA_1.1"]
    assert (a["refseq_accession"], a["refseq_category"], a["sample_accession"]) == ("GCF_1.1", "reference genome", "SAMN1")
    assert (a["strain"], a["species_name"], a["seq_rel_date"]) == ("4287", "Fusarium oxysporum", "2012-05-01")
    assert rows["GCA_4.1"]["sample_accession"] == "" and rows["GCA_4.1"]["species_name"] == ""   # 'na' and unmapped -> empty
    assert list(rows) == sorted(rows)


LABEL_HDR = ["genome", "ncbi_species", "fid_label", "fid_call", "fid_dist", "fid_diffs", "fid_len", "fid_margin", "fid_other",
             "n_loci_hit", "agreement"]


def _producer(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(LABEL_HDR); w.writerows(rows)


def _setup(tmp_path):
    _assemblies(tmp_path)
    _producer(tmp_path / "labels.tsv", [
        ["GCF_1.1_FO2", "Fusarium oxysporum", "", "Fusarium lycopersici", "0.002", "5", "2500", "0.003", "Fusarium curvatum", "5", ""],
        ["GCA_2.1_x", "Fusarium proliferatum", "", "Fusarium annulatum", "0.000", "0", "2600", "0.004", "Fusarium fujikuroi", "5", ""],
        ["GCA_3.1_y", "Fusarium oxysporum", "", "Fusarium fabacearum;Fusarium wimaladesilvae", "0.0", "0", "3000", "0.001", "", "5", ""],
        ["GCA_4.1_z", "Fusarium sp. X", "", "Fusarium avenaceum", "0.02", "50", "2500", "0.001", "", "5", ""],
        ["GCA_9.1_gone", "Fusarium solani", "", "Neocosmospora solani", "0.001", "2", "2500", "0.003", "", "5", ""],
    ])
    with open(tmp_path / "complexes.tsv", "w") as fh:
        fh.write("# comment line\nspecies\tcomplex\tcollapse\n"
                 "Fusarium lycopersici\tFOSC\t1\nFusarium fabacearum\tFOSC\t1\nFusarium wimaladesilvae\tFOSC\t1\n"
                 "Fusarium annulatum\tFFSC\t0\n")
    _csv_gz(tmp_path / "runs.csv.gz", ["run_accession", "sample_accession", "tax_id"],
            [["SRR1", "SAMN1", "5507"], ["SRR2", "SAMN1;SAMN7", "5507"], ["SRR3", "SAMN2", "5128"], ["SRR9", "SAMN9", "5507"]])


def _import(tmp_path, source="fid_a", producer="labels.tsv"):
    return li.main(["--source-id", source, "--format", "gmg_fusariumid", "--input", str(tmp_path / producer),
                    "--complexes", str(tmp_path / "complexes.tsv"), "--assemblies", str(tmp_path / "assemblies.csv.gz"),
                    "--runs", str(tmp_path / "runs.csv.gz"), "--labels-dir", str(tmp_path / "labels"),
                    "--authority", "FusariumID", "--authority-release", "2026-06-23", "--project", "GmG", "--tool", "t",
                    "--date", "2026-11-09"])


def test_labels_resolve_assemblies_and_reach_runs_by_biosample(tmp_path):
    _setup(tmp_path)
    assert _import(tmp_path) == 0
    rows = {r["assembly_accession"]: r for r in _read_gz(tmp_path / "labels" / "assembly_labels.csv.gz")}
    a = rows["GCA_1.1"]                                                  # GCF id in the producer -> its GenBank twin
    assert (a["sample_accession"], a["label_species"], a["label_complex"], a["label_status"], a["label_agreement"]) == \
        ("SAMN1", "Fusarium lycopersici", "FOSC", "confident", "differs")
    assert (rows["GCA_2.1"]["label_complex"], rows["GCA_2.1"]["label_agreement"]) == ("FFSC", "differs")
    t = rows["GCA_3.1"]                                                  # a tie inside one complex keeps both names and the complex
    assert (t["label_status"], t["label_complex"]) == ("tie", "FOSC")
    z = rows["GCA_4.1"]                                                  # beyond --max-dist: no label, nearest kept as evidence
    assert (z["label_species"], z["label_status"], z["nearest_species"], z["label_agreement"]) == ("", "no_call", "Fusarium avenaceum", "no_call")
    assert rows["GCA_9.1"]["sample_accession"] == ""                     # not in the assembly table: kept, not dropped
    runs = {(r["run_accession"], r["assembly_accession"]) for r in _read_gz(tmp_path / "labels" / "run_labels.csv.gz")}
    assert runs == {("SRR1", "GCA_1.1"), ("SRR2", "GCA_1.1"), ("SRR3", "GCA_2.1")}
    with open(tmp_path / "labels" / "sources.csv") as fh:
        src = list(csv.DictReader(fh))
    assert (src[0]["assemblies"], src[0]["runs"]) == ("5", "3") and "0.01" in src[0]["call_rule"]


def test_reimport_replaces_only_its_own_rows_and_keeps_curated_columns(tmp_path):
    _setup(tmp_path); _import(tmp_path); _import(tmp_path, source="fid_b")
    path = tmp_path / "labels" / "assembly_labels.csv.gz"
    rows = _read_gz(path)
    for r in rows:
        if r["assembly_accession"] == "GCA_2.1" and r["source_id"] == "fid_a":
            r["curated_species"], r["curated_note"] = "Fusarium annulatum", "checked by hand"
    li.write_gz(str(path), rows, li.LABEL_COLS)
    _producer(tmp_path / "labels2.tsv", [
        ["GCA_2.1_x", "Fusarium proliferatum", "", "Fusarium proliferatum", "0.001", "1", "2600", "0.004", "", "5", ""]])
    _import(tmp_path, source="fid_a", producer="labels2.tsv")            # fid_a refreshed with one assembly only
    rows = {(r["assembly_accession"], r["source_id"]): r for r in _read_gz(path)}
    a = rows[("GCA_2.1", "fid_a")]
    assert (a["label_species"], a["label_agreement"], a["curated_species"], a["curated_note"]) == \
        ("Fusarium proliferatum", "same", "Fusarium annulatum", "checked by hand")
    assert ("GCA_1.1", "fid_a") not in rows and ("GCA_1.1", "fid_b") in rows   # the other source is untouched
    with open(tmp_path / "labels" / "sources.csv") as fh:
        src = {s["source_id"]: (s["assemblies"], s["runs"]) for s in csv.DictReader(fh)}
    assert src == {"fid_a": ("1", "1"), "fid_b": ("5", "3")}
    runs = _read_gz(tmp_path / "labels" / "run_labels.csv.gz")
    assert sorted((r["run_accession"], r["source_id"]) for r in runs) == \
        [("SRR1", "fid_b"), ("SRR2", "fid_b"), ("SRR3", "fid_a"), ("SRR3", "fid_b")]


def test_duplicate_assemblies_in_one_import_fail(tmp_path):
    _setup(tmp_path)
    _producer(tmp_path / "dup.tsv", [["GCA_2.1_x", "", "", "Fusarium annulatum", "0", "0", "1", "0", "", "5", ""],
                                     ["GCA_2.1_y", "", "", "Fusarium annulatum", "0", "0", "1", "0", "", "5", ""]])
    try:
        _import(tmp_path, producer="dup.tsv")
    except SystemExit as e:
        assert "more than once" in str(e)
    else:
        raise AssertionError("a duplicated assembly must stop the import")
