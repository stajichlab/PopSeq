"""Tests for the advertised query contract.

These check the contract other tools depend on, not the catalog's contents:
the shapes, the exit codes, the failure modes, and the rules that are easy to
break silently. A test that asserted an exact run count would fail on every
archive refresh for no reason.
"""
from __future__ import annotations

import json
import subprocess
import sys

import pytest

import popseq
from popseq import cli

ZT = 1047171          # Zymoseptoria tritici, a mid-sized species
PARAPSILOSIS = 5480   # carries two names in NCBI; the issue #1 regression


def test_index_declares_its_schema_and_tables():
    idx = popseq.index()
    assert idx["popseq_schema"] == popseq.SCHEMA_VERSION
    assert "runs" in idx["tables"]
    assert idx["tables"]["runs"]["key"] == "run_accession"
    assert idx["caveats"], "the index must carry its caveats"


def test_species_sorted_and_filtered():
    rows = popseq.species(min_isolates=100)
    assert rows, "the catalog should hold species with 100+ isolates"
    iso = [int(r["isolates"]) for r in rows]
    assert iso == sorted(iso, reverse=True)
    assert min(iso) >= 100


def test_resolve_rejects_an_ambiguous_name_rather_than_guessing():
    with pytest.raises(popseq.SpeciesNotFound) as exc:
        popseq.resolve_species(name="Candida")
    assert "matches" in str(exc.value)


def test_resolve_rejects_an_unknown_taxid():
    with pytest.raises(popseq.SpeciesNotFound):
        popseq.resolve_species(taxid=1)


def test_panel_returns_one_run_per_isolate():
    res = popseq.panel(taxid=ZT, n=50)
    samples = [r["sample_accession"] for r in res["runs"] if r["sample_accession"]]
    assert len(samples) == len(set(samples))


def test_panel_spreads_across_studies():
    """With far more studies than slots, every pick should be a different study."""
    res = popseq.panel(taxid=ZT, n=20)
    assert res["studies_available"] > 20
    assert res["studies_in_panel"] == 20


def test_panel_is_deterministic_for_a_seed():
    a = popseq.panel(taxid=ZT, n=25, seed=7)
    b = popseq.panel(taxid=ZT, n=25, seed=7)
    assert [r["run_accession"] for r in a["runs"]] == \
           [r["run_accession"] for r in b["runs"]]


def test_min_depth_drops_runs_without_an_estimate():
    res = popseq.panel(taxid=ZT, n=30, min_depth=40)
    for r in res["runs"]:
        assert r["est_depth_x"] != ""
        assert float(r["est_depth_x"]) >= 40


def test_exclusions_are_labelled_not_deleted():
    """The mutant/evolved material must still be reachable, for TE_fungi_copy."""
    clean = popseq.runs(taxid=ZT, limit=None)
    assert all(r["exclusion"] == "" for r in clean["runs"])
    evolved = popseq.runs(taxid=ZT, only_exclusion="evolution")
    assert evolved["runs"] and all(r["exclusion"] == "evolution"
                                   for r in evolved["runs"])


def test_impossible_filter_raises_empty_not_an_empty_list():
    with pytest.raises(popseq.EmptyResult):
        popseq.panel(taxid=ZT, min_depth=1e9)


def test_a_synonym_pair_resolves_to_one_taxon():
    """Issue #1. NCBI labels this taxon both Candida parapsilosis and
    Lodderomyces parapsilosis. One taxid, so one species row, so one set of runs."""
    rows = [r for r in popseq.species(min_isolates=1)
            if int(r["species_taxid"]) == PARAPSILOSIS]
    assert len(rows) == 1
    res = popseq.runs(taxid=PARAPSILOSIS, limit=5)
    assert res["species_taxid"] == PARAPSILOSIS


def test_sample_sheet_is_nf_core_shaped_with_real_urls():
    sheet = popseq.to_sample_sheet(popseq.panel(taxid=ZT, n=5))
    assert list(sheet[0]) == ["sample", "fastq_1", "fastq_2"]
    for row in sheet:
        assert row["fastq_1"].startswith("https://")


def test_preload_does_not_change_the_answer():
    before = popseq.panel(taxid=ZT, n=15, seed=3)
    popseq.preload([ZT])
    try:
        after = popseq.panel(taxid=ZT, n=15, seed=3)
    finally:
        popseq.preload_clear()
    assert before["runs"] == after["runs"]


@pytest.mark.parametrize("args,code", [
    (["species", "--min-isolates", "100", "--limit", "2"], 0),
    (["describe", "--species", "Candida"], 2),
    (["describe", "--taxid", "1"], 2),
    (["panel", "--taxid", str(ZT), "--min-depth", "1000000"], 3),
])
def test_cli_exit_codes(args, code, capsys):
    assert cli.main(args) == code
    capsys.readouterr()


def test_cli_json_is_parseable(capsys):
    assert cli.main(["panel", "--taxid", str(ZT), "--n", "3",
                     "--format", "json"]) == 0
    doc = json.loads(capsys.readouterr().out)
    assert doc["popseq_schema"] == popseq.SCHEMA_VERSION
    assert len(doc["runs"]) == 3


def test_script_shim_still_works():
    """specification 007 advertised scripts/popseq_query.py before the package."""
    r = subprocess.run([sys.executable, "scripts/popseq_query.py", "species",
                        "--limit", "1", "--format", "tsv"],
                       capture_output=True, text=True)
    assert r.returncode == 0
    assert r.stdout.startswith("species_taxid\t")


def test_country_filter_accepts_iso_code_name_and_raw():
    """Issue #23: "CH", "Switzerland" and the raw string select the same runs."""
    by_iso = popseq.runs(taxid=ZT, country="CH", include_excluded=True)["runs"]
    by_name = popseq.runs(taxid=ZT, country="switzerland", include_excluded=True)["runs"]
    assert by_iso and {r["run_accession"] for r in by_iso} == \
        {r["run_accession"] for r in by_name}
    assert all(r["country_iso"] == "CH" for r in by_iso)
    # The raw values behind them differ, which is the point of the field.
    assert len({r["country"] for r in by_iso}) > 1


def test_host_filter_by_taxid_and_name():
    a = popseq.runs(taxid=ZT, host="4565", include_excluded=True)["runs"]
    b = popseq.runs(taxid=ZT, host="Triticum aestivum", include_excluded=True)["runs"]
    assert a and {r["run_accession"] for r in a} <= {r["run_accession"] for r in b}
    assert all(r["host_taxid"] == "4565" for r in a)
