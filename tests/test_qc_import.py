"""Specification 008 rules for scripts/popseq_qc_import.py: a source replaces only its own rows, curated_* columns survive
an import, an unresolved run is kept (not dropped), and sources.csv counts each source's rows."""
from __future__ import annotations

import csv
import gzip
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("popseq_qc_import", ROOT / "scripts" / "popseq_qc_import.py")
qi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qi)


def _write(path, header, rows, delim=","):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=delim, lineterminator="\n"); w.writerow(header); w.writerows(rows)


def _qc(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def _setup(tmp_path):
    _write(tmp_path / "runs.tsv", ["run_accession", "sample_accession", "tax_id"], [["R1", "S1", "5507"], ["R2", "S2", "5507"], ["R3", "S3", "746128"]], "\t")
    with gzip.open(tmp_path / "tax.csv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(["taxid", "name", "rank", "species_taxid", "species_name", "resolution"])
        w.writerows([["5507", "Fusarium oxysporum", "species", "5507", "Fusarium oxysporum", "species"],
                     ["746128", "Aspergillus fumigatus", "species", "746128", "Aspergillus fumigatus", "species"]])
    hdr = ["run", "species", "strain", "study", "screen", "screen_genus", "second_genus", "second_hashes", "second_share", "deeper_screen", "call", "truth", "wave"]
    _write(tmp_path / "flags.tsv", hdr, [["R1", "Fusarium oxysporum", "s1", "P", "mixed", "Fusarium", "Meyerozyma", "40", "0.05", "", "MAT1-1", "", "wave1"],
                                         ["R3", "Aspergillus fumigatus", "s3", "P", "borderline", "Aspergillus", "Penicillium", "6", "0.8", "yes", "MAT1-2", "", "wave1"],
                                         ["R9", "Unknown", "s9", "P", "genus_mismatch", "Candida", "", "", "", "", "none", "", "wave1"]], "\t")


def _run(tmp_path, source="src_a", input_name="flags.tsv"):
    return qi.main(["--source-id", source, "--format", "matpredict_screen", "--input", str(tmp_path / input_name), "--runs", str(tmp_path / "runs.tsv"),
                    "--taxonomy", str(tmp_path / "tax.csv.gz"), "--qc-dir", str(tmp_path / "qc"), "--project", "MATPredict", "--tool", "GmG",
                    "--date", "2026-11-08"])


def test_rows_land_in_species_files_and_unresolved_runs_are_kept(tmp_path):
    _setup(tmp_path); _run(tmp_path)
    fo = _qc(tmp_path / "qc" / "5507.Fusarium_oxysporum.qc.csv.gz")
    assert [(r["run_accession"], r["sample_accession"], r["status"], r["second_taxon"]) for r in fo] == [("R1", "S1", "mixed", "Meyerozyma")]
    assert _qc(tmp_path / "qc" / "unresolved.qc.csv.gz")[0]["run_accession"] == "R9"
    src = list(csv.DictReader(open(tmp_path / "qc" / "sources.csv")))
    assert src[0]["source_id"] == "src_a" and src[0]["rows"] == "3"


def test_reimport_replaces_only_its_own_rows_and_keeps_curated_columns(tmp_path):
    _setup(tmp_path); _run(tmp_path)
    _run(tmp_path, source="src_b")                                                      # a second source flags the same runs
    path = tmp_path / "qc" / "5507.Fusarium_oxysporum.qc.csv.gz"
    rows = _qc(path)
    for r in rows:
        if r["source_id"] == "src_a": r["curated_status"] = "confirmed"                 # a curator reviews src_a's flag
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.DictWriter(fh, qi.QC_COLS, lineterminator="\n"); w.writeheader(); w.writerows(rows)
    _write(tmp_path / "flags2.tsv", ["run", "species", "strain", "study", "screen", "screen_genus", "second_genus", "second_hashes", "second_share", "deeper_screen", "call", "truth", "wave"],
           [["R1", "Fusarium oxysporum", "s1", "P", "genus_mismatch", "Meyerozyma", "", "", "", "", "none", "", "wave1"]], "\t")
    _run(tmp_path, source="src_a", input_name="flags2.tsv")                            # src_a refreshed: R3 and R9 no longer flagged by it
    rows = {(r["run_accession"], r["source_id"]): r for r in _qc(path)}
    assert rows[("R1", "src_a")]["status"] == "genus_mismatch" and rows[("R1", "src_a")]["curated_status"] == "confirmed"
    assert rows[("R1", "src_b")]["status"] == "mixed"                                   # the other source is untouched
    af = _qc(tmp_path / "qc" / "746128.Aspergillus_fumigatus.qc.csv.gz")
    assert [r["source_id"] for r in af] == ["src_b"]                                     # src_a's old row for R3 is gone
    src = {s["source_id"]: s["rows"] for s in csv.DictReader(open(tmp_path / "qc" / "sources.csv"))}
    assert src == {"src_a": "1", "src_b": "3"}
