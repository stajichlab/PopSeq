"""Spec 008 checks PopSeq evaluates itself (scripts/popseq_qc_checks.py): low_coverage (issue #37) and
oversized_archive_object (issue #49); and the hand-entered `popseq_findings` import format. Small synthetic fixtures."""
from __future__ import annotations

import csv
import gzip
import importlib.util
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


qc = _load("popseq_qc_checks")
qi = qc.qi


def _write(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(header); w.writerows(rows)


def _qc(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def _sources(tmp_path):
    return {r["source_id"]: r for r in csv.DictReader(open(tmp_path / "qc" / "sources.csv"))}


def _world(tmp_path):
    # Af genome 29.38 Mb; R_low 0.0059x, R_edge 1.016x, R_ok 60x; R_nobase has base_count 0; R_nosize is a species with no size
    _write(tmp_path / "runs.csv", ["run_accession", "sample_accession", "tax_id", "species_taxid", "base_count"],
           [["R_low", "S1", "746128", "746128", "173443"], ["R_edge", "S2", "746128", "746128", "29851492"],
            ["R_ok", "S3", "746128", "746128", "1762800000"], ["R_nobase", "S4", "746128", "746128", "0"],
            ["R_nosize", "S5", "5007", "5007", "1000"], ["R_big", "S6", "5007", "5007", "7589121360"]])
    _write(tmp_path / "sizes.csv", ["species_taxid", "species_name", "genome_mb"], [["746128", "Aspergillus fumigatus", "29.38"], ["5007", "Brettanomyces bruxellensis", ""]])
    with gzip.open(tmp_path / "tax.csv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(["taxid", "name", "rank", "species_taxid", "species_name", "resolution"])
        w.writerows([["746128", "Aspergillus fumigatus", "species", "746128", "Aspergillus fumigatus", "species"],
                     ["5007", "Brettanomyces bruxellensis", "species", "5007", "Brettanomyces bruxellensis", "species"]])


def _common(tmp_path):
    return ["--runs", str(tmp_path / "runs.csv"), "--taxonomy", str(tmp_path / "tax.csv.gz"), "--qc-dir", str(tmp_path / "qc"), "--date", "2026-11-21"]


def test_low_coverage_flags_below_the_threshold_and_skips_runs_it_cannot_examine(tmp_path):
    _world(tmp_path)
    rows, n = qc.low_coverage(qc.iter_runs(str(tmp_path / "runs.csv")), qc.genome_sizes(str(tmp_path / "sizes.csv")), 1.0)
    assert [r["run_accession"] for r in rows] == ["R_low"]               # R_edge is 1.016x: not below 1x
    assert rows[0]["status"] == "below_1x" and rows[0]["source_note"].startswith("depth 0.0059x")
    assert n == {"examined": 3, "no_genome_size": 2, "no_base_count": 1, "flagged": 1}


def test_low_coverage_default_flags_only_runs_in_the_list(tmp_path):
    _world(tmp_path)
    (tmp_path / "list.csv").write_text("run_accession,issue\nR_ok,PopSeq#37\n")        # R_low is below 1x but is not in the list
    qc.main(["low-coverage", *_common(tmp_path), "--genome-sizes", str(tmp_path / "sizes.csv"), "--run-list", str(tmp_path / "list.csv")])
    assert not (tmp_path / "qc" / "746128.Aspergillus_fumigatus.qc.csv.gz").exists()    # R_low not flagged, R_ok is 60x
    (tmp_path / "list.csv").write_text("run_accession,issue\nR_ok,PopSeq#37\nR_low,PopSeq#37\n")
    qc.main(["low-coverage", *_common(tmp_path), "--genome-sizes", str(tmp_path / "sizes.csv"), "--run-list", str(tmp_path / "list.csv")])
    got = _qc(tmp_path / "qc" / "746128.Aspergillus_fumigatus.qc.csv.gz")
    assert [(r["run_accession"], r["status"]) for r in got] == [("R_low", "below_1x")]  # R_edge (1.016x) and R_nosize stay out regardless


def test_low_coverage_requires_a_run_list_or_the_catalog_wide_flag(tmp_path):
    import pytest
    _world(tmp_path)
    with pytest.raises(SystemExit):
        qc.main(["low-coverage", *_common(tmp_path), "--genome-sizes", str(tmp_path / "sizes.csv")])
    assert not (tmp_path / "qc").exists()


def test_run_list_reads_a_table_or_plain_accessions(tmp_path):
    (tmp_path / "a.csv").write_text("run_accession,issue\nR1,x\nR2,y\n"); (tmp_path / "b.txt").write_text("# c\nR3\nR4\n")
    assert qc.read_run_list(str(tmp_path / "a.csv")) == {"R1", "R2"} and qc.read_run_list(str(tmp_path / "b.txt")) == {"R3", "R4"}


def test_low_coverage_import_writes_rows_and_an_excluding_source(tmp_path):
    _world(tmp_path)
    qc.main(["low-coverage", *_common(tmp_path), "--genome-sizes", str(tmp_path / "sizes.csv"), "--below", "1.5", "--catalog-wide"])
    got = _qc(tmp_path / "qc" / "746128.Aspergillus_fumigatus.qc.csv.gz")
    assert [(r["run_accession"], r["check"], r["status"]) for r in got] == [("R_edge", "low_coverage", "below_1.5x"), ("R_low", "low_coverage", "below_1.5x")]
    src = _sources(tmp_path)["popseq_low_coverage"]
    assert src["exclude_statuses"] == "below_1.5x" and src["rows"] == "2" and src["check"] == "low_coverage"


def test_no_exclude_records_findings_that_skip_nothing(tmp_path):
    _world(tmp_path)
    qc.main(["low-coverage", *_common(tmp_path), "--genome-sizes", str(tmp_path / "sizes.csv"), "--no-exclude", "--catalog-wide"])
    assert _sources(tmp_path)["popseq_low_coverage"]["exclude_statuses"] == ""


def test_archive_size_flags_bytes_per_base_above_the_threshold(tmp_path):
    _world(tmp_path)
    _write(tmp_path / "archive.csv", ["run_accession", "archive_bytes", "archive_note"],
           [["R_big", "51460000000", "GmG P6"], ["R_ok", "705000000", ""], ["R_absent", "10", ""]])
    qc.main(["archive-size", *_common(tmp_path), "--sizes", str(tmp_path / "archive.csv"), "--source-id", "arch", "--project", "GmG"])
    got = _qc(tmp_path / "qc" / "5007.Brettanomyces_bruxellensis.qc.csv.gz")
    assert [(r["run_accession"], r["check"], r["status"]) for r in got] == [("R_big", "oversized_archive_object", "oversized")]
    assert got[0]["source_note"].startswith("6.78 archive bytes per base") and got[0]["source_note"].endswith("GmG P6")
    assert not (tmp_path / "qc" / "746128.Aspergillus_fumigatus.qc.csv.gz").exists()   # R_ok is 0.40 bytes per base
    assert _sources(tmp_path)["arch"]["exclude_statuses"] == "oversized"


def test_popseq_findings_format_keeps_the_hand_entered_columns(tmp_path):
    _world(tmp_path)
    _write(tmp_path / "f.csv", ["run_accession", "status", "filed_name", "detected_genus", "detected_species", "source_note"],
           [["R_ok", "mismatch", "Aspergillus fumigatus", "Aspergillus", "Aspergillus nidulans", "identity 0.998"]])
    qi.main(["--source-id", "hand", "--format", "popseq_findings", "--input", str(tmp_path / "f.csv"), "--runs", str(tmp_path / "runs.csv"),
             "--taxonomy", str(tmp_path / "tax.csv.gz"), "--qc-dir", str(tmp_path / "qc"), "--project", "MATPredict", "--tool", "t",
             "--check", "species_match", "--date", "2026-11-10", "--exclude-statuses", "mismatch;unverified_species"])
    r = _qc(tmp_path / "qc" / "746128.Aspergillus_fumigatus.qc.csv.gz")[0]
    assert (r["status"], r["filed_name"], r["detected_species"], r["source_note"], r["sample_accession"]) == \
        ("mismatch", "Aspergillus fumigatus", "Aspergillus nidulans", "identity 0.998", "S3")
    assert _sources(tmp_path)["hand"]["exclude_statuses"] == "mismatch;unverified_species"


def test_an_import_that_changes_nothing_leaves_the_file_bytes_alone(tmp_path):
    _world(tmp_path)
    args = ["low-coverage", *_common(tmp_path), "--genome-sizes", str(tmp_path / "sizes.csv"), "--catalog-wide"]
    qc.main(args)
    path = tmp_path / "qc" / "746128.Aspergillus_fumigatus.qc.csv.gz"
    before = path.read_bytes(); os.utime(path, (0, 0))
    qc.main(args)
    assert path.read_bytes() == before and path.stat().st_mtime == 0      # not rewritten
    assert before[4:8] == b"\x00\x00\x00\x00"                             # gzip header mtime 0: same rows, same bytes
