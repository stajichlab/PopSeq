"""Tests for country and host normalisation (issue #23). No network.

Country tests use the committed aliases plus pycountry. Host tests use a tiny
vocabulary written to tmp_path, so they do not depend on the generated table.
"""
from __future__ import annotations

import csv

import pytest

from popseq import normalise as nz


@pytest.mark.parametrize("raw,iso,status", [
    ("Switzerland", "CH", "ok"),
    ("Switzerland: Eschikon", "CH", "ok"),
    ("Switzerland:Eschikon", "CH", "ok"),
    ("USA: California, Riverside", "US", "ok"),
    ("UK (Manchester)", "GB", "ok"),
    ("Kerala [India]", "IN", "ok"),
    ("Cerkvenjak, Slovenia", "SI", "ok"),
    ("Iran", "IR", "ok"),                 # pycountry "Iran, Islamic Republic of"
    ("Côte d'Ivoire", "CI", "ok"),        # accent folded
    ("Phillipines", "PH", "ok"),          # curated misspelling
    ("", "", "absent"),
    ("missing", "", "missing"),
    ("Not Collected", "", "missing"),
    ("missing: control sample", "", "missing"),
    ("Indian Ocean", "", "region"),
    ("Yugoslavia", "", "historical"),
    ("Korea", "", "unresolved"),          # two countries; not guessed
    ("missing;USA;USA", "US", "ok"),      # multi-sample run, parts agree
    ("USA;China", "", "conflicting"),
])
def test_country(raw, iso, status):
    got = nz.country(raw)
    assert (got["country_iso"], got["country_status"]) == (iso, status)
    assert got["country_status"] in nz.COUNTRY_STATUSES


def test_country_never_writes_NA():
    for raw in ["", "NA", "N/A", "not applicable"]:
        got = nz.country(raw)
        assert "NA" not in got.values()


def test_country_aliases_are_valid():
    """Every curated alias carries a known status and a real ISO code."""
    with open(nz.COUNTRY_ALIASES, newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert rows
    for r in rows:
        assert r["status"] in nz.COUNTRY_STATUSES, r
        if r["status"] == "ok":
            assert r["iso2"], r
        else:
            assert not r["iso2"], r


@pytest.fixture
def vocab(tmp_path):
    p = tmp_path / "host_taxa.csv"
    rows = [
        ("triticum aestivum", "4565", "Triticum aestivum", "ok"),
        ("wheat", "4565", "Triticum aestivum", "ok"),
        ("saccharomyces cerevisiae", "4932", "Saccharomyces cerevisiae", "ok"),
        ("soil", "", "", "not_a_host"),
        ("sugarcane", "", "", "ambiguous"),
        ("hordeum vulgare", "4513", "Hordeum vulgare", "ok"),
    ]
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["host_key", "host_taxid", "host_name", "host_status"])
        w.writerows(rows)
    nz._host_index.cache_clear()
    yield str(p)
    nz._host_index.cache_clear()


@pytest.mark.parametrize("raw,taxid,status", [
    ("Triticum aestivum", "4565", "ok"),
    ("  WHEAT ", "4565", "ok"),
    ("", "", "absent"),
    ("not collected", "", "missing"),
    ("Soil", "", "not_a_host"),
    ("Sugarcane", "", "ambiguous"),
    ("Something odd", "", "unresolved"),
    ("wheat;Triticum aestivum", "4565", "ok"),
    ("wheat;Hordeum vulgare", "", "conflicting"),
])
def test_host(vocab, raw, taxid, status):
    got = nz.host(raw, taxa_path=vocab)
    assert (got["host_taxid"], got["host_status"]) == (taxid, status)


def test_host_naming_the_organism_itself_is_flagged(vocab):
    got = nz.host("Saccharomyces cerevisiae", organism_taxids=(4932, 4932),
                  taxa_path=vocab)
    assert got["host_status"] == "same_as_organism"
    assert got["host_taxid"] == "4932"
    # The same string on a run of another organism is a real host.
    got = nz.host("Saccharomyces cerevisiae", organism_taxids=(5476,), taxa_path=vocab)
    assert got["host_status"] == "ok"


def test_key_folds_case_accents_and_space():
    assert nz.key("  Côte  d'Ivoire. ") == "cote d'ivoire"
