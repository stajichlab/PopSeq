"""Fixture tests for the harvesters and the inventory (issues #19 and #3). No network.

Issue #19 lists the bugs that printed a plausible number instead of failing. Each test
here pins one of them, or one of the pieces the issue names:

  FileBase resolution, four forms      scripts/popseq_inventory.py candidates_for
  mate pairing                         scripts/popseq_inventory.py pair_up
  gzip round trip, the .gz temp name   scripts/popseq_biosamples.py write_table
  shard resume                         scripts/popseq_biosamples.py load_shards, main
  the rank-artifact guard              scripts/popseq_reconcile.py fetch_window, main
  a species list read twice            scripts/popseq_reconcile.py main, from a pipe
  multi-valued sample_accession        biosamples, survey, popseq.query (issue #3)
  `done < file </dev/null`             a static check over scripts/*.sh and *.slurm
"""
from __future__ import annotations

import csv
import glob
import gzip
import os
import re
import socket
import sys

import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))

import popseq_biosamples as bs  # noqa: E402
import popseq_inventory as inv  # noqa: E402
import popseq_reconcile as rec  # noqa: E402
import popseq_survey as survey  # noqa: E402

from popseq import query as q  # noqa: E402


def _touch(root, *rels):
    for rel in rels:
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"")


def _names(hits):
    return [h[0] for h in hits]


# --- FileBase resolution: the four forms of spec 005 section 2 -----------------------

@pytest.fixture
def fq(tmp_path):
    _touch(tmp_path,
           "SRR1636746_1.fastq.gz", "SRR1636746_2.fastq.gz",
           "SRR16367460_1.fastq.gz", "SRR16367460_2.fastq.gz",      # a longer accession
           "100D-1_S34_L001_R1_001.fastq.gz", "100D-1_S34_L001_R2_001.fastq.gz",
           "100D-10_S35_L001_R1_001.fastq.gz",                       # another sample
           "22687-1_ACTTGA_L001_R1_001.fq.gz", "22687-1_ACTTGA_L001_R2_001.fq.gz",
           "AF100.run1/lane1/x_R1_001.fastq.gz", "AF100.run1/lane1/x_R2_001.fastq.gz",
           "AF100.run1/lane2/x_R1_001.fastq.gz", "AF100.run1/lane2/x_R2_001.fastq.gz",
           "AF100.run1/lane3/x_R1_001.fastq.gz",                     # a lane not listed
           "notes.txt")
    return inv.index_fastq(str(tmp_path))


def test_index_holds_read_files_only(fq):
    assert len(fq) == 14 and not any(bn == "notes.txt" for _, _, bn in fq)


def test_filebase_an_sra_accession_does_not_take_a_longer_accession(fq):
    assert _names(inv.candidates_for("SRR1636746", fq)) == [
        "SRR1636746_1.fastq.gz", "SRR1636746_2.fastq.gz"]


def test_filebase_a_glob_with_a_mate_wildcard(fq):
    assert _names(inv.candidates_for("100D-1_S34_L001_R[12]_001.fastq.gz", fq)) == [
        "100D-1_S34_L001_R1_001.fastq.gz", "100D-1_S34_L001_R2_001.fastq.gz"]


def test_filebase_a_bare_prefix(fq):
    assert _names(inv.candidates_for("22687-1_ACTTGA_L001", fq)) == [
        "22687-1_ACTTGA_L001_R1_001.fq.gz", "22687-1_ACTTGA_L001_R2_001.fq.gz"]


def test_filebase_several_relative_paths_one_per_lane(fq):
    fb = "AF100.run1/lane1/x_R[12]_001.fastq.gz;AF100.run1/lane2/x_R[12]_001.fastq.gz"
    hits = inv.candidates_for(fb, fq)
    assert _names(hits) == [f"AF100.run1/lane{n}/x_R{m}_001.fastq.gz"
                            for n in (1, 2) for m in (1, 2)]
    assert len(inv.pair_up(hits)) == 2


def test_filebase_that_matches_nothing_is_empty_not_everything(fq):
    assert inv.candidates_for("ZZZ999", fq) == []
    assert inv.candidates_for(";", fq) == []


# --- mate pairing -------------------------------------------------------------------

@pytest.mark.parametrize("name, stem, mate", [
    ("s_R1_001.fastq.gz", "s", "1"),        # _R1_ before _1: not split on the 1 in 001
    ("s_R2_001.fastq.gz", "s", "2"),
    ("SRR1_1.fastq.gz", "SRR1", "1"),
    ("SRR1_2.fq", "SRR1", "2"),
    ("SRR1.fastq.gz", "SRR1", None),
])
def test_split_mate(name, stem, mate):
    assert inv.split_mate(name) == (stem, mate)


def test_pair_up_pairs_mates_and_keeps_a_lone_mate_2_visible():
    hits = [(f"d/{b}", f"/abs/d/{b}", b) for b in
            ("a_R1_001.fastq.gz", "a_R2_001.fastq.gz", "b.fastq.gz", "c_R2_001.fastq.gz")]
    assert inv.pair_up(hits) == [
        ("/abs/d/a_R1_001.fastq.gz", "/abs/d/a_R2_001.fastq.gz"),
        ("/abs/d/b.fastq.gz", ""),
        ("/abs/d/c_R2_001.fastq.gz", ""),
    ]


# --- gzip round trip; the .gz temp-name regression ----------------------------------

def test_write_table_to_a_gz_name_writes_gzip_and_reads_back(tmp_path):
    rows = {"SAMN2": {"biosample": "SAMN2", "host": "Homo sapiens"},
            "SAMN1": {"biosample": "SAMN1", "strain": "A, B"}}
    path = tmp_path / "t.csv.gz"
    bs.write_table(str(path), rows, ["strain", "host"])
    # The regression: the temp name hid .gz, so plain text went into a .gz file.
    assert path.read_bytes()[:2] == b"\x1f\x8b"
    assert not list(tmp_path.glob("*.tmp*"))
    back, cols = bs.load_existing(str(path))
    assert cols == ["biosample", "strain", "host"]
    assert back["SAMN1"]["strain"] == "A, B" and back["SAMN2"]["strain"] == ""
    assert list(back) == ["SAMN1", "SAMN2"]


def test_write_table_to_a_plain_name_writes_plain_text(tmp_path):
    path = tmp_path / "t.csv"
    bs.write_table(str(path), {"S1": {"biosample": "S1"}}, [])
    assert path.read_text().splitlines() == ["biosample", "S1"]


# --- shard resume -------------------------------------------------------------------

def test_load_shards_resumes_by_accession_and_drops_an_unreadable_shard(tmp_path, capsys):
    d = tmp_path / "shards"
    d.mkdir()
    bs.write_shard(str(d), 0, {"S1": {"biosample": "S1", "host": "x"}})
    bs.write_shard(str(d), 3, {"S2": {"biosample": "S2", "strain": "y"}})
    (d / "batch_000004.csv").write_bytes(b"biosample\n\xff\xfe\xfa")  # not UTF-8
    (d / "batch_000005.csv.tmp").write_text("biosample\nS9\n")         # never renamed
    rows, seen, nxt = bs.load_shards(str(d))
    assert set(rows) == {"S1", "S2"}
    assert {"biosample", "host", "strain"} <= seen
    assert nxt == 5                     # after the highest shard number, even an unreadable one
    assert "batch_000004.csv unreadable" in capsys.readouterr().err


def test_load_shards_on_a_missing_directory_starts_at_zero(tmp_path):
    assert bs.load_shards(str(tmp_path / "nope")) == ({}, set(), 0)


@pytest.fixture
def fake_ncbi(monkeypatch):
    """Replace the Entrez round trip. Records every batch it was asked for."""
    asked = []

    def fetch_batch(accessions, debug=False):
        asked.append(list(accessions))
        return {a: {"biosample": a, "host": f"h{a}"} for a in accessions if a.startswith("SAM")}

    monkeypatch.setattr(bs, "fetch_batch", fetch_batch)
    monkeypatch.setattr(bs.time, "sleep", lambda s: None)
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    return asked


def _harvest(monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", ["popseq_biosamples.py", "-e", "test@example.org", *args])
    return bs.main()


def test_a_resumed_harvest_asks_only_for_what_no_shard_holds(tmp_path, monkeypatch, fake_ncbi):
    d, out = tmp_path / "shards", tmp_path / "bs.csv.gz"
    d.mkdir()
    bs.write_shard(str(d), 0, {"SAMN1": {"biosample": "SAMN1", "host": "old"}})
    assert _harvest(monkeypatch, "-s", "SAMN1", "SAMN2", "SAMN3", "--batch", "1",
                    "--shard-dir", str(d), "-o", str(out)) == 0
    assert fake_ncbi == [["SAMN2"], ["SAMN3"]]
    assert sorted(os.listdir(d)) == ["batch_000000.csv", "batch_000001.csv", "batch_000002.csv"]
    rows, _ = bs.load_existing(str(out))
    assert sorted(rows) == ["SAMN1", "SAMN2", "SAMN3"] and rows["SAMN1"]["host"] == "old"


def test_a_failed_accession_is_listed_and_the_exit_code_says_so(tmp_path, monkeypatch, fake_ncbi):
    out = tmp_path / "bs.csv"
    assert _harvest(monkeypatch, "-s", "SAMN1", "XYZ1", "-o", str(out)) == 1
    assert (tmp_path / "bs.csv.failed.txt").read_text() == "XYZ1\n"


# --- multi-valued sample_accession (issue #3) ----------------------------------------

def test_split_accessions():
    assert bs.split_accessions("SAMN1;SAMN2 ; ;SAMN3\n") == ["SAMN1", "SAMN2", "SAMN3"]
    assert bs.split_accessions("") == [] and bs.split_accessions(None) == []


def test_from_runs_sends_each_biosample_of_a_joined_value(tmp_path, monkeypatch, fake_ncbi):
    runs = tmp_path / "runs.csv.gz"
    with gzip.open(runs, "wt", newline="") as fh:
        fh.write("# provenance line\n")
        w = csv.writer(fh)
        w.writerow(["run_accession", "sample_accession"])
        w.writerows([["R1", "SAMN1;SAMN2"], ["R2", "SAMN2"], ["R3", ""], ["R4", "SAMN3 ; SAMN1"]])
    out = tmp_path / "bs.csv"
    assert _harvest(monkeypatch, "--from-runs", str(runs), "-o", str(out)) == 0
    assert fake_ncbi == [["SAMN1", "SAMN2", "SAMN3"]]
    assert not any(";" in a for batch in fake_ncbi for a in batch)


def test_input_file_and_command_line_are_split_too(tmp_path, monkeypatch, fake_ncbi):
    lst = tmp_path / "acc.txt"
    lst.write_text("SAMN1;SAMN2\n\nSAMN2\n")
    _harvest(monkeypatch, "-i", str(lst), "-o", str(tmp_path / "a.csv"))
    _harvest(monkeypatch, "-s", "SAMN4;SAMN5", "-o", str(tmp_path / "b.csv"))
    assert fake_ncbi == [["SAMN1", "SAMN2"], ["SAMN4", "SAMN5"]]


def test_survey_counts_each_biosample_of_a_joined_value():
    assert survey.split_samples("S1;S2") == ["S1", "S2"]
    assert survey.isolate_key({"sample_accession": "S2;S1"}) == "S1;S2"
    assert survey.isolate_key({"sample_accession": ""}) == ""


def test_query_isolate_key_and_sample_sheet_with_a_joined_value():
    assert q.split_samples(" S1 ; S2;") == ["S1", "S2"]
    assert q.isolate_key({"run_accession": "R1", "sample_accession": "S2;S1;S2"}) == "S1;S2"
    assert q.isolate_key({"run_accession": "R1", "sample_accession": ""}) == "R1"
    ftp = "ftp.sra.ebi.ac.uk/x/{0}_1.fastq.gz;ftp.sra.ebi.ac.uk/x/{0}_2.fastq.gz"
    res = {"runs": [{"run_accession": "R1", "sample_accession": "S1;S2", "fastq_ftp": ftp.format("R1")},
                    {"run_accession": "R2", "sample_accession": "S3", "fastq_ftp": ftp.format("R2")},
                    {"run_accession": "R3", "sample_accession": "", "fastq_ftp": ftp.format("R3")}]}
    assert [r["sample"] for r in q.to_sample_sheet(res)] == ["R1", "S3", "R3"]


def test_panel_takes_one_run_of_a_joined_sample_set_whatever_its_order(monkeypatch):
    base = dict(species_taxid="1", species_name="Aus bus", study_accession="P1",
                fastq_ftp="a_1.fastq.gz;a_2.fastq.gz")
    runs = [dict(base, run_accession="R1", sample_accession="S1;S2", est_depth_x="10"),
            dict(base, run_accession="R2", sample_accession="S2;S1", est_depth_x="30"),
            dict(base, run_accession="R3", sample_accession="S1", est_depth_x="20")]
    monkeypatch.setattr(q, "resolve_species", lambda name=None, taxid=None:
                        {"species_taxid": "1", "species_name": "Aus bus"})
    monkeypatch.setattr(q, "_species_runs", lambda tx: runs)
    monkeypatch.setattr(q, "_filter", lambda rows, **kw: list(rows))
    monkeypatch.setattr(q, "_apply_qc", lambda rows, tx, inc: (rows, []))
    res = q.panel(taxid=1, n=10)
    assert sorted(r["run_accession"] for r in res["runs"]) == ["R2", "R3"]
    assert res["isolates_available"] == 2


# --- the rank-artifact guard ----------------------------------------------------------

def test_a_name_above_species_rank_raises_before_any_page_is_fetched(tmp_path, monkeypatch):
    monkeypatch.setattr(rec, "esearch", lambda taxid, extra="": (255437, "W", "1"))
    monkeypatch.setattr(rec, "get", lambda *a, **k: pytest.fail("fetched a page"))
    with pytest.raises(rec.RankArtifact, match="above species rank"):
        rec.fetch_window(4751, "", "Fungi", str(tmp_path), 500, 5000, expected=40)
    assert os.listdir(tmp_path) == []


def test_a_count_near_the_expected_one_is_fetched_and_resumes_by_page(tmp_path, monkeypatch):
    pages = []
    monkeypatch.setattr(rec, "esearch", lambda taxid, extra="": (1200, "W", "1"))
    monkeypatch.setattr(rec, "get", lambda url, params: pages.append(params["retstart"]) or "Run\n")
    monkeypatch.setattr(rec.time, "sleep", lambda s: None)
    (tmp_path / "Aus_bus.p0001.csv").write_text("Run\n")  # page 1 survived an earlier run
    got = rec.fetch_window(5, "", "Aus_bus", str(tmp_path), 500, 5000, expected=1000)
    assert got == 1200 and pages == [0, 1000]
    assert sorted(os.listdir(tmp_path)) == [f"Aus_bus.p000{i}.csv" for i in range(3)]


def test_the_guard_only_applies_at_the_top_and_with_an_expectation(tmp_path, monkeypatch):
    monkeypatch.setattr(rec, "esearch", lambda taxid, extra="": (300, "W", "1"))
    monkeypatch.setattr(rec, "get", lambda url, params: "Run\n")
    monkeypatch.setattr(rec.time, "sleep", lambda s: None)
    assert rec.fetch_window(5, "", "x", str(tmp_path), 500, 5000, expected=None) == 300
    assert rec.fetch_window(5, "", "y", str(tmp_path), 500, 5000, depth=1, expected=1) == 300


# --- reconcile: a species list read once, from a pipe ----------------------------------

def _reconcile_fixture(tmp_path):
    shards = tmp_path / "shards"
    shards.mkdir()
    (shards / "Aus_bus.p0000.csv").write_text(
        "Run,TaxID\nSRR1,11\nSRR2,10\nSRR9,99\n")       # SRR9 is another species
    with gzip.open(tmp_path / "taxmap.csv.gz", "wt") as fh:
        fh.write("taxid,species_taxid,species_name\n10,10,Aus bus\n11,10,Aus bus\n"
                 "99,99,Cus dus\n")
    with gzip.open(tmp_path / "ena.tsv.gz", "wt") as fh:
        fh.write("run_accession\ttax_id\tsample_accession\n"
                 "SRR1\t11\tS1\nERR3\t10\tS3;S4\nERR8\t99\tS8\n")
    return shards


def test_a_species_list_from_a_pipe_is_read_once_and_scopes_the_comparison(
        tmp_path, monkeypatch, capsys):
    shards = _reconcile_fixture(tmp_path)
    r, w = os.pipe()                       # what bash <(...) hands a script
    os.write(w, b"species,isolates,runs\nAus bus,150,3\nTiny one,2,2\n")
    os.close(w)
    monkeypatch.setattr(sys, "argv", [
        "popseq_reconcile.py", "--species-list", f"/dev/fd/{r}", "--compare-only",
        "--ena", str(tmp_path / "ena.tsv.gz"), "--taxmap", str(tmp_path / "taxmap.csv.gz"),
        "--shard-dir", str(shards), "--out-dir", str(tmp_path / "out")])
    try:
        assert rec.main() == 0
    finally:
        os.close(r)
    cap = capsys.readouterr()
    assert "Scoping to 1 species taxids (from 1 names)" in cap.err
    assert "ENA distinct runs:  2" in cap.err      # not 0: the list was not read empty
    assert "ENA only : 1 runs" in cap.out and "NCBI only: 1 runs" in cap.out
    assert "In both  : 1" in cap.out
    with gzip.open(tmp_path / "out" / "fungi_byspecies.ena_only_runs.csv.gz", "rt") as fh:
        assert [row["run_accession"] for row in csv.DictReader(fh)] == ["ERR3"]


def test_reconcile_species_list_skips_a_rank_artifact_and_marks_it_done(tmp_path, monkeypatch):
    shards = _reconcile_fixture(tmp_path)
    sl = tmp_path / "species.csv"
    sl.write_text("species,isolates,runs\nFungi,500,40\n")
    monkeypatch.setattr(rec, "resolve_taxid", lambda name: 4751)
    monkeypatch.setattr(rec, "esearch", lambda taxid, extra="": (255437, "W", "1"))
    monkeypatch.setattr(sys, "argv", [
        "popseq_reconcile.py", "--species-list", str(sl),
        "--ena", str(tmp_path / "ena.tsv.gz"), "--taxmap", str(tmp_path / "taxmap.csv.gz"),
        "--shard-dir", str(shards), "--out-dir", str(tmp_path / "out")])
    assert rec.main() == 0
    assert (shards / "Fungi.DONE").read_text().startswith("skipped: rank artifact")


# --- bash: `done < file </dev/null` -----------------------------------------------------

def test_no_shell_loop_reads_its_input_and_dev_null_at_once():
    """Static check. `done < list </dev/null`: the last redirect wins, the loop reads
    /dev/null and reports 0 of everything. The shell scripts have no Python entry point
    to test, so this only stops the pattern from coming back."""
    bad = re.compile(r"\bdone\b[^#\n]*<\s*[^<\s][^#\n]*<\s*/dev/null")
    hits = []
    for path in sorted(glob.glob(os.path.join(HERE, "scripts", "*.sh"))
                       + glob.glob(os.path.join(HERE, "scripts", "*.slurm"))):
        with open(path) as fh:
            hits += [f"{os.path.basename(path)}:{i}" for i, line in enumerate(fh, 1)
                     if bad.search(line)]
    assert not hits, hits
    assert bad.search("done < runs.txt </dev/null")
    assert not bad.search("done < runs.txt")


# --- the network guard in tests/conftest.py ------------------------------------------

def test_the_suite_cannot_reach_the_network():
    if os.environ.get("POPSEQ_ALLOW_NETWORK") == "1":
        pytest.skip("network guard switched off")
    # conftest.NetworkBlocked; matched by message because pytest loads conftest itself.
    with pytest.raises(RuntimeError, match="test tried to resolve 'www.ebi.ac.uk'"):
        socket.create_connection(("www.ebi.ac.uk", 443), timeout=1)
    s = socket.socket()
    try:
        with pytest.raises(RuntimeError, match="test tried to connect"):
            s.connect(("192.0.2.1", 443))
    finally:
        s.close()


def test_loopback_is_still_allowed():
    srv = socket.socket()
    srv.bind(("127.0.0.1", 0))
    srv.listen(1)
    try:
        c = socket.create_connection(srv.getsockname(), timeout=1)
        c.close()
    finally:
        srv.close()
