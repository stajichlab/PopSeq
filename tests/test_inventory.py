"""Fixture tests for scripts/popseq_inventory.py: read_state (spec 003, issue #5).

The fixture subproject is built in tmp_path. It mirrors the file naming of the real
Afumigatus_Global subproject (trim_galore output, <strain>_1_val_1.fq.gz) without
reading it: that subproject lives on the HPCC only.
"""
from __future__ import annotations

import csv
import os
import sys

import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))

import popseq_inventory as inv  # noqa: E402


@pytest.mark.parametrize("name,state", [
    ("AF100_1_val_1.fq.gz", "trimmed"),        # trim_galore, paired
    ("AF100_2_val_2.fq.gz", "trimmed"),
    ("AF100_trimmed.fq.gz", "trimmed"),        # trim_galore, single-end
    ("AF100_R1.trim.fastq.gz", "trimmed"),     # Trimmomatic-style
    ("AF100_R1_001.fastq.gz", "unknown"),      # no marker: never guessed raw
    ("SRR1636746_1.fastq.gz", "unknown"),
    ("/data/trimmed/AF100_1.fq.gz", "unknown"),  # a directory name does not count
])
def test_file_read_state(name, state):
    assert inv.file_read_state(name) == state


def test_sample_read_state_needs_every_file_trimmed():
    assert inv.sample_read_state([("/a/x_1_val_1.fq.gz", "/a/x_2_val_2.fq.gz")]) == "trimmed"
    assert inv.sample_read_state([("/a/x_trimmed.fq.gz", "")]) == "trimmed"
    # Mixed within one sample: not all trimmed, so unknown.
    assert inv.sample_read_state([("/a/x_1_val_1.fq.gz", "/a/x_2_val_2.fq.gz"),
                                  ("/a/x_L2_1.fq.gz", "/a/x_L2_2.fq.gz")]) == "unknown"
    assert inv.sample_read_state([]) == "unknown"
    assert set(inv.READ_STATES) == {"raw", "trimmed", "unknown"}


def _subproject(root, sheet_rows, files):
    sub = root / "Afumigatus_Global"
    (sub / "input").mkdir(parents=True)
    (sub / "config.txt").write_text("SAMPFILE=samples.csv\nFASTQFOLDER=input\n")
    with open(sub / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Strain", "FileBase"])
        w.writerows(sheet_rows)
    for f in files:
        (sub / "input" / f).write_bytes(b"@r\nACGT\n+\nIIII\n")
    return sub


def _metadata(out_root, name="Afumigatus_Global"):
    with open(out_root / "panels" / name / "metadata.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def test_converter_writes_read_state_trimmed_for_val_files(tmp_path):
    sub = _subproject(tmp_path / "pg",
                      [["AF100", "AF100"], ["AF293", "AF293"]],
                      ["AF100_1_val_1.fq.gz", "AF100_2_val_2.fq.gz",
                       "AF293_1_val_1.fq.gz", "AF293_2_val_2.fq.gz"])
    res = inv.process("Afumigatus_Global", str(sub), str(tmp_path / "out"), False)
    rows = _metadata(tmp_path / "out")
    assert [r["read_state"] for r in rows] == ["trimmed", "trimmed"]
    assert res["read_state_trimmed"] == 2
    with open(tmp_path / "out" / "panels" / "Afumigatus_Global" / "samples.csv") as fh:
        s = list(csv.DictReader(fh))
    assert s[0]["fastq_1"].endswith("AF100_1_val_1.fq.gz")
    assert s[0]["fastq_2"].endswith("AF100_2_val_2.fq.gz")


def test_converter_mixed_panel_is_reported_per_sample(tmp_path):
    # One trimmed sample, one unmarked sample, one listed with no files.
    sub = _subproject(tmp_path / "pg",
                      [["AF100", "AF100"], ["AF200", "AF200"], ["AF300", "AF300"]],
                      ["AF100_1_val_1.fq.gz", "AF100_2_val_2.fq.gz",
                       "AF200_R1_001.fastq.gz", "AF200_R2_001.fastq.gz"])
    res = inv.process("Afumigatus_Global", str(sub), str(tmp_path / "out"), False)
    rows = {r["sample"]: r for r in _metadata(tmp_path / "out")}
    assert rows["AF100"]["read_state"] == "trimmed"
    assert rows["AF200"]["read_state"] == "unknown"
    assert rows["AF300"]["read_state"] == "unknown"
    assert rows["AF300"]["provenance"] == "listed_absent"
    assert (res["read_state_trimmed"], res["read_state_unknown"]) == (1, 2)
    assert "raw" not in {r["read_state"] for r in rows.values()}
