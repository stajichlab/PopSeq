"""Issue #22: every run with no species taxid is listed in catalog/summary/unplaced_runs.csv.gz.

A fixture test pins the classification rules. The catalog test pins the property the
issue asks for: the table holds exactly the runs of the run table whose species_taxid is
empty, so none is dropped silently. The issue's figures (648 above species, 98 unresolved)
are the ILLUMINA runs of the 2026-11-18 run table; with BGISEQ and DNBSEQ that table has
822 and 100 (922 in all). They are asserted only while the run table is that build, so a
refresh does not fail the suite for a count that was meant to change.
"""
from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import hashlib
import pathlib

import pytest

import popseq
from popseq import cli
from popseq import query as q

ROOT = pathlib.Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("popseq_unplaced",
                                               ROOT / "scripts" / "popseq_unplaced.py")
up = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(up)

TAXMAP = {
    "5052": dict(taxid="5052", name="Aspergillus", rank="genus", species_taxid="",
                 species_name="", resolution="above_species:genus"),
    "4751": dict(taxid="4751", name="Fungi", rank="kingdom", species_taxid="",
                 species_name="", resolution="above_species:kingdom"),
    "9001": dict(taxid="9001", name="some clade", rank="clade", species_taxid="",
                 species_name="", resolution="above_species:order"),
    "9002": dict(taxid="9002", name="orphan", rank="no rank", species_taxid="",
                 species_name="", resolution="unresolved"),
    "746128": dict(taxid="746128", name="Aspergillus fumigatus", rank="species",
                   species_taxid="746128", species_name="Aspergillus fumigatus",
                   resolution="self"),
    "5085": dict(taxid="5085", name="A. fumigatus strain", rank="strain",
                 species_taxid="746128", species_name="Aspergillus fumigatus",
                 resolution="rolled_up"),
}


def _run(acc, tax_id, species_taxid=""):
    return dict(run_accession=acc, study_accession="PRJ" + acc, sample_accession="SAM" + acc,
                scientific_name="x", tax_id=tax_id, instrument_platform="ILLUMINA",
                first_public="2020-01-01", species_taxid=species_taxid)


def test_classify_each_placement():
    c = up.classify
    assert c("5052", TAXMAP) == dict(placement="above_species", resolved_rank="genus",
                                     ranked_taxid="5052", ranked_name="Aspergillus",
                                     part_species_taxids="")
    assert c("4751", TAXMAP)["resolved_rank"] == "kingdom"
    # An unranked taxid whose walk stopped at an ancestor: the map does not name it.
    d = c("9001", TAXMAP)
    assert (d["placement"], d["resolved_rank"], d["ranked_taxid"]) == ("above_species", "order", "")
    assert c("9002", TAXMAP)["placement"] == "unresolved_in_map"
    assert c("123456789", TAXMAP)["placement"] == "not_in_taxonomy_map"
    assert c("", TAXMAP)["placement"] == "no_taxid"
    m = c("5085;746128;5052;77", TAXMAP)
    assert (m["placement"], m["resolved_rank"], m["part_species_taxids"]) == \
        ("multiple_taxids", "unresolved", "746128")
    assert c("746128", TAXMAP)["placement"] == "species_in_map"


def test_table_keeps_every_empty_species_run_and_nothing_else(tmp_path):
    runs = [_run("R1", "746128", "746128"), _run("R2", "5052"), _run("R3", ""),
            _run("R4", "5085;746128"), _run("R5", "4751"), _run("R6", "31337")]
    rows = up.unplaced(runs, TAXMAP)
    assert sorted(r["run_accession"] for r in rows) == ["R2", "R3", "R4", "R5", "R6"]
    # Round trip through the writer and the reader the query layer uses.
    out = tmp_path / "unplaced.csv.gz"
    up.write(rows, str(out), "# fixture")
    again = list(q._rows(str(out)))
    assert [r["run_accession"] for r in again] == [r["run_accession"] for r in rows]
    assert list(again[0]) == up.FIELDS
    # Same inputs, same bytes (gzip mtime 0).
    first = out.read_bytes()
    up.write(rows, str(out), "# fixture")
    assert out.read_bytes() == first


def test_committed_table_accounts_for_every_run_with_no_species_taxid():
    with gzip.open(q.RUNS_TABLE, "rt", newline="") as fh:
        empty = {r["run_accession"] for r in csv.DictReader(fh) if not r["species_taxid"].strip()}
    table = [r["run_accession"] for r in q._rows(q.UNPLACED_TABLE)]
    assert len(table) == len(empty), "unplaced_runs.csv.gz is stale: rerun popseq_unplaced.py"
    assert len(set(table)) == len(table)
    assert set(table) == empty
    rows = popseq.unplaced()
    assert {r["placement"] for r in rows} <= set(up.PLACEMENTS)
    assert not [r for r in rows if r["placement"] == "species_in_map"], \
        "the run table is older than the taxonomy map"
    assert all(r["resolved_rank"] not in ("", "unresolved", "species")
               for r in rows if r["placement"] == "above_species")


# sha256 of catalog/runs/fungi_all_runs.csv.gz in the 2026-11-18 build (index.json).
ISSUE_22_RUN_TABLE = "6c11a1d61299268d221ae52818a123336662be8b037c6bac9e560042a348e5c8"


def test_issue_22_figures_on_the_table_they_were_measured_on():
    h = hashlib.sha256(pathlib.Path(q.RUNS_TABLE).read_bytes()).hexdigest()
    if h != ISSUE_22_RUN_TABLE:
        pytest.skip("run table rebuilt since issue #22; its figures no longer apply")
    rows = popseq.unplaced()
    ill = [r for r in rows if r["instrument_platform"] == "ILLUMINA"]
    assert sum(r["placement"] == "above_species" for r in ill) == 648
    assert sum(r["resolved_rank"] == "unresolved" for r in ill) == 98
    assert sum(r["placement"] == "above_species" for r in rows) == 822
    assert sum(r["resolved_rank"] == "unresolved" for r in rows) == 100
    assert len(rows) == 922


def test_committed_table_is_registered_in_the_index():
    t = popseq.index()["tables"]["unplaced_runs"]
    assert t["path"] == "catalog/summary/unplaced_runs.csv.gz"
    assert t["built_by"] == "scripts/popseq_unplaced.py"
    assert t["columns"] == up.FIELDS == q.UNPLACED_FIELDS


def test_cli_unplaced(capsys):
    assert cli.main(["unplaced", "--placement", "above_species", "--rank", "genus",
                     "--format", "json"]) == 0
    doc = json.loads(capsys.readouterr().out)
    assert doc["returned"] == len(doc["runs"]) > 0
    assert all(r["resolved_rank"] == "genus" for r in doc["runs"])
    assert cli.main(["unplaced", "--rank", "no_such_rank"]) == 3
    capsys.readouterr()
    assert cli.main(["species", "--limit", "1", "--format", "json"]) == 0
    assert json.loads(capsys.readouterr().out)["unplaced_runs"] == q.unplaced_count() > 0
