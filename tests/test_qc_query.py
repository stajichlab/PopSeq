"""Spec 008 and 007: run QC in `runs` and `panel`. A source's `exclude_statuses` skip a run by default, `include_flagged` keeps it,
a curated false_positive and a stale verdict skip nothing, and `effective_species` follows spec 009. Synthetic tables, no catalog."""
from __future__ import annotations

import csv
import gzip

import pytest

from popseq import query as q

TX = 5507


def _run(acc, sample, study="S1", depth="40", **kw):
    r = {"run_accession": acc, "sample_accession": sample, "study_accession": study, "species_taxid": str(TX), "species_name": "Fusarium oxysporum",
         "instrument_platform": "ILLUMINA", "ena_layout": "PAIRED", "library_strategy": "WGS", "exclusion": "", "est_depth_x": depth,
         "base_count": "2000000000", "first_public": "2024-01-01", "fastq_ftp": f"ftp/{acc}_1.fastq.gz;ftp/{acc}_2.fastq.gz"}
    r.update(kw)
    return r


def _gz(path, cols, rows):
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.DictWriter(fh, cols, lineterminator="\n"); w.writeheader(); w.writerows(rows)


QC_COLS = ["run_accession", "source_id", "check", "status", "filed_name", "curated_status"]


@pytest.fixture
def world(tmp_path, monkeypatch):
    qc = tmp_path / "qc"; qc.mkdir()
    _gz(qc / f"{TX}.Fusarium_oxysporum.qc.csv.gz", QC_COLS, [
        {"run_accession": "R_bad", "source_id": "gm", "check": "species_match", "status": "mismatch", "filed_name": "Fusarium oxysporum"},
        {"run_accession": "R_mixed", "source_id": "gm", "check": "species_match", "status": "mixed", "filed_name": "Fusarium oxysporum"},
        {"run_accession": "R_fp", "source_id": "gm", "check": "species_match", "status": "mismatch", "filed_name": "Fusarium oxysporum", "curated_status": "false_positive"},
        {"run_accession": "R_stale", "source_id": "gm", "check": "species_match", "status": "mismatch", "filed_name": "Fusarium oxysporum"},
        {"run_accession": "R_other", "source_id": "screen", "check": "contamination", "status": "mixed", "filed_name": ""}])
    with open(qc / "sources.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, ["source_id", "exclude_statuses"], lineterminator="\n"); w.writeheader()
        w.writerows([{"source_id": "gm", "exclude_statuses": "mismatch;host_dominated"}, {"source_id": "screen", "exclude_statuses": ""}])
    labels = tmp_path / "run_labels.csv.gz"
    _gz(labels, ["run_accession", "label_species", "curated_species"], [
        {"run_accession": "R_stale", "label_species": "", "curated_species": "Fusarium oxysporum f. sp. relabelled"}])
    rows = [_run("R_ok", "SAMa"), _run("R_bad", "SAMb"), _run("R_mixed", "SAMc"), _run("R_fp", "SAMd"), _run("R_stale", "SAMe"), _run("R_other", "SAMf")]
    monkeypatch.setattr(q, "QC_DIR", str(qc)); monkeypatch.setattr(q, "RUN_LABELS_TABLE", str(labels)); q._LABELS_CACHE.clear()
    monkeypatch.setattr(q, "RUN_CURATION_TABLE", str(tmp_path / "run_curation.csv")); q._CURATION_CACHE.clear()      # no curation file
    monkeypatch.setattr(q, "resolve_species", lambda name=None, taxid=None: {"species_taxid": str(TX), "species_name": "Fusarium oxysporum"})
    monkeypatch.setattr(q, "_species_runs", lambda tx: iter(rows))
    return rows


def _ids(result):
    return [r["run_accession"] for r in result["runs"]]


def test_runs_skips_what_the_source_excludes_and_says_how_many(world):
    res = q.runs(taxid=TX)
    # R_bad is skipped; R_mixed is not in gm's exclude_statuses; R_fp is a curated false positive; R_stale was relabelled since the verdict
    assert sorted(_ids(res)) == ["R_fp", "R_mixed", "R_ok", "R_other", "R_stale"]
    assert res["qc_skipped"] == 1 and res["qc_skipped_runs"] == ["R_bad"]


def test_include_flagged_brings_the_run_back_with_its_status(world):
    res = q.runs(taxid=TX, include_flagged=True)
    assert "R_bad" in _ids(res) and res["qc_skipped"] == 0
    by = {r["run_accession"]: r for r in res["runs"]}
    assert by["R_bad"]["qc_status"] == "species_match:mismatch"
    assert by["R_ok"]["qc_status"] == ""
    assert by["R_other"]["qc_status"] == "contamination:mixed"          # a finding that does not skip is still shown


def test_a_stale_verdict_is_marked_and_skips_nothing(world):
    by = {r["run_accession"]: r for r in q.runs(taxid=TX)["runs"]}
    assert by["R_stale"]["qc_status"] == "species_match:mismatch:stale"
    assert by["R_stale"]["effective_species"] == "Fusarium oxysporum f. sp. relabelled"
    assert by["R_ok"]["effective_species"] == "Fusarium oxysporum"       # no label: the NCBI name


def test_a_curated_false_positive_lifts_the_skip(world):
    assert "R_fp" in _ids(q.runs(taxid=TX))


def test_panel_applies_qc_before_it_picks_isolates_and_reports_the_skip(world):
    res = q.panel(taxid=TX, n=10)
    assert "R_bad" not in _ids(res) and res["qc_skipped"] == 1
    assert "R_bad" in _ids(q.panel(taxid=TX, n=10, include_flagged=True))


def test_everything_flagged_is_an_empty_result_that_names_the_flag(world, monkeypatch):
    monkeypatch.setattr(q, "_species_runs", lambda tx: iter([_run("R_bad", "SAMb")]))
    with pytest.raises(q.EmptyResult, match="include-flagged"):
        q.runs(taxid=TX)


def test_no_qc_dir_and_no_labels_change_nothing(world, monkeypatch, tmp_path):
    monkeypatch.setattr(q, "QC_DIR", str(tmp_path / "nowhere")); monkeypatch.setattr(q, "RUN_LABELS_TABLE", str(tmp_path / "none.csv.gz")); q._LABELS_CACHE.clear()
    res = q.runs(taxid=TX)
    assert len(res["runs"]) == 6 and res["qc_skipped"] == 0
    assert all(r["qc_status"] == "" for r in res["runs"])


def test_the_new_fields_are_in_the_run_contract():
    assert q.RUN_FIELDS[-2:] == ["effective_species", "qc_status"]


def test_a_run_level_curation_wins_and_makes_an_older_verdict_stale(world, monkeypatch, tmp_path):
    cur = tmp_path / "run_curation.csv"
    cur.write_text("run_accession,curated_species,curated_note,curator,date\nR_bad,Aspergillus fumigatus,reads are Aspergillus,owner,2026-11-19\nR_stale,Fusarium oxysporum f. sp. y,newer,owner,2026-11-19\n")
    monkeypatch.setattr(q, "RUN_CURATION_TABLE", str(cur)); q._CURATION_CACHE.clear()
    res = q.runs(taxid=TX)
    by = {r["run_accession"]: r for r in res["runs"]}
    # R_bad was filed Fusarium oxysporum and flagged mismatch; the curator has corrected the species, so the verdict is stale and skips nothing
    assert "R_bad" in by and by["R_bad"]["effective_species"] == "Aspergillus fumigatus"
    assert by["R_bad"]["qc_status"] == "species_match:mismatch:stale"
    # the run-level curation beats the assembly-level one (R_stale has both)
    assert by["R_stale"]["effective_species"] == "Fusarium oxysporum f. sp. y"
    assert res["qc_skipped"] == 0


def test_a_missing_curation_file_changes_nothing(world):
    assert q._curation() == {}
