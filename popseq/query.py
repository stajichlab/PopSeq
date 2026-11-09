"""The query layer. One implementation; the CLI, the MCP server and the static
JSON renderer are all thin wrappers over this.

Design notes that matter to a caller:

* Species are identified by **taxid**. A name is a convenience that resolves
  through the catalog and raises on an ambiguous match rather than guessing.
  NCBI labels one taxon with more than one name (issue #1).
* Nothing is deleted from the tables. `exclusion` labels mutant, evolved, cross
  and non-WGS material; the default filters drop it, `include_excluded=True`
  keeps it, and `only_exclusion=...` selects just that material.
* Everything returns plain dicts and lists, JSON-serialisable as they stand.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import random
import re
from collections import defaultdict
from typing import Any, Iterable, Iterator, Optional

SCHEMA_VERSION = "1.0.0"

ROOT = os.environ.get("POPSEQ_ROOT") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))

RUNS_TABLE = os.path.join(ROOT, "catalog", "runs", "fungi_all_runs.csv.gz")
SUMMARY_TABLE = os.path.join(ROOT, "catalog", "kingdom", "fungi_species_summary.csv.gz")
CLONALITY_TABLE = os.path.join(ROOT, "catalog", "summary", "clonality.csv.gz")
COVERAGE_TABLE = os.path.join(ROOT, "catalog", "summary", "coverage.csv.gz")
INDEX_FILE = os.path.join(ROOT, "catalog", "index.json")

SPECIES_FIELDS = ["species_taxid", "species_name", "runs", "isolates", "studies",
                  "est_download_tb", "first_year", "last_year"]
RUN_FIELDS = ["run_accession", "sample_accession", "study_accession",
              "species_taxid", "species_name", "strain", "library_layout", "ena_layout",
              "read_count", "base_count", "est_depth_x", "country",
              "country_iso", "country_status", "host", "host_taxid", "host_name",
              "host_status", "collection_date", "first_public", "exclusion", "ncbi_status",
              "fastq_ftp", "fastq_aspera"]

CAVEATS = [
    "Every count is a floor. ENA holds runs with an empty tax_id that no taxonomy "
    "query reaches, so true per-species totals are higher (issue #2).",
    "A run count is not a strain count. 52 of the 113 largest species take half or "
    "more of their isolates from a single study.",
    "est_depth_x is base_count over a reference assembly size: metadata, not a "
    "measurement on the reads.",
    "Only the 113 species with 100 or more isolates were reconciled against NCBI. "
    "Elsewhere ncbi_status is not_checked, which is not the same as absent.",
]


class CatalogError(Exception):
    """Base class. Exit code 1."""


class SpeciesNotFound(CatalogError):
    """No species matched, or a name matched more than one. Exit code 2."""


class EmptyResult(CatalogError):
    """The filters left nothing. Exit code 3."""


def _open(path: str):
    return gzip.open(path, "rt", newline="") if path.endswith(".gz") else open(path, newline="")


def _rows(path: str) -> Iterator[dict]:
    """Read a catalog CSV, skipping the optional provenance comment line."""
    with _open(path) as fh:
        if not fh.readline().startswith("#"):
            fh.seek(0)
        yield from csv.DictReader(fh)


def index() -> dict:
    """The catalog manifest: tables, columns, checksums, caveats."""
    with open(INDEX_FILE) as fh:
        return json.load(fh)


def species(min_isolates: int = 1, contains: Optional[str] = None,
            limit: Optional[int] = None) -> list[dict]:
    """Species in the catalog, most isolates first."""
    out = [r for r in _rows(SUMMARY_TABLE) if int(r["isolates"]) >= min_isolates]
    if contains:
        out = [r for r in out if contains.lower() in r["species_name"].lower()]
    out.sort(key=lambda r: -int(r["isolates"]))
    if limit:
        out = out[:limit]
    return [{k: r[k] for k in SPECIES_FIELDS} for r in out]


def resolve_species(name: Optional[str] = None,
                    taxid: Optional[int] = None) -> dict:
    """Resolve a species to its catalog summary row.

    Raises SpeciesNotFound on no match and on an ambiguous name. It does not
    guess: guessing is how a query silently answers about the wrong organism.
    """
    rows = list(_rows(SUMMARY_TABLE))
    if taxid is not None:
        hit = [r for r in rows if int(r["species_taxid"]) == int(taxid)]
        if not hit:
            raise SpeciesNotFound(f"taxid {taxid} is not in the catalog")
        return hit[0]
    if not name:
        raise SpeciesNotFound("give a name or a taxid")
    want = name.strip().lower()
    exact = [r for r in rows if r["species_name"].lower() == want]
    if len(exact) == 1:
        return exact[0]
    near = [r for r in rows if want in r["species_name"].lower()]
    if len(near) == 1:
        return near[0]
    if not near:
        raise SpeciesNotFound(f'no species matching "{name}"')
    names = sorted(r["species_name"] for r in near)
    raise SpeciesNotFound(
        f'"{name}" matches {len(near)} species: {", ".join(names[:8])}'
        f'{" ..." if len(names) > 8 else ""}. Pass a taxid.')


# Optional in-memory cache. A single query scans the run table once, which takes
# about two seconds and is not worth caching. A batch job that asks about many
# species pays that per species, so it can preload instead. This is an
# optimisation only: the answers are identical either way.
_CACHE: Optional[dict[int, list[dict]]] = None


def preload(taxids: Optional[Iterable[int]] = None) -> int:
    """Read the run table once and hold the rows for these species in memory.

    Pass None to hold every species. Returns the number of rows held. Call
    `preload_clear()` to drop it.
    """
    global _CACHE
    want = set(int(t) for t in taxids) if taxids is not None else None
    cache: dict[int, list[dict]] = defaultdict(list)
    n = 0
    for r in _rows(RUNS_TABLE):
        if not r["species_taxid"]:
            continue
        tx = int(r["species_taxid"])
        if want is None or tx in want:
            cache[tx].append(r)
            n += 1
    _CACHE = dict(cache)
    return n


def preload_clear() -> None:
    global _CACHE
    _CACHE = None


def _species_runs(taxid: int) -> Iterator[dict]:
    if _CACHE is not None:
        yield from _CACHE.get(taxid, [])
        return
    for r in _rows(RUNS_TABLE):
        if r["species_taxid"] and int(r["species_taxid"]) == taxid:
            yield r


def _depth(r: dict) -> float:
    try:
        return float(r["est_depth_x"])
    except (TypeError, ValueError):
        return -1.0


def _country_matches(r: dict, want: str) -> bool:
    """An ISO code, a country name, or the raw string, case-insensitive (issue #23).

    Matches the normalised fields first, so "CH", "Switzerland" and
    "Switzerland: Eschikon" all select the same runs.
    """
    w = want.strip().lower()
    if r.get("country_status") == "ok" and w in (
            (r.get("country_iso") or "").lower(), (r.get("country_name") or "").lower()):
        return True
    return (r.get("country") or "").strip().lower() == w


def _host_matches(r: dict, want) -> bool:
    """A host taxid (exact) or a host name (case-insensitive, resolved or raw)."""
    w = str(want).strip().lower()
    if w.isdigit():
        return r.get("host_taxid") == w and r.get("host_status") == "ok"
    return (r.get("host_status") == "ok" and (r.get("host_name") or "").lower() == w) \
        or (r.get("host") or "").strip().lower() == w


def _filter(rows: Iterable[dict], *, platform="ILLUMINA", layout="PAIRED",
            strategy="WGS", min_depth=None, min_gbases=None, country=None,
            host=None, since=None, include_excluded=False, only_exclusion=None,
            require_fastq=True) -> list[dict]:
    out = []
    for r in rows:
        if platform and (r["instrument_platform"] or "").upper() != platform.upper():
            continue
        # Pairing is what ENA serves, not the submitter's label (spec 001, ena_layout).
        if layout and (r.get("ena_layout") or "").upper() != layout.upper():
            continue
        if strategy and (r["library_strategy"] or "").upper() != strategy.upper():
            continue
        if only_exclusion is not None:
            if r["exclusion"] != only_exclusion:
                continue
        elif not include_excluded and r["exclusion"]:
            continue
        if min_depth is not None and _depth(r) < min_depth:
            # A run with no estimate cannot be shown to pass, so it does not.
            continue
        if min_gbases is not None:
            try:
                if float(r["base_count"] or 0) < min_gbases * 1e9:
                    continue
            except ValueError:
                continue
        if country and not _country_matches(r, country):
            continue
        if host and not _host_matches(r, host):
            continue
        if since and (r["first_public"] or "") < since:
            continue
        if require_fastq and not r["fastq_ftp"]:
            continue
        out.append(r)
    return out


def runs(taxid: Optional[int] = None, name: Optional[str] = None,
         limit: Optional[int] = None, **filters: Any) -> dict:
    """Every run of one species that passes the filters."""
    row = resolve_species(name=name, taxid=taxid)
    tx = int(row["species_taxid"])
    hits = _filter(_species_runs(tx), **filters)
    if not hits:
        raise EmptyResult(f'{row["species_name"]}: no runs pass the filters')
    if limit:
        hits = hits[:limit]
    return {
        "popseq_schema": SCHEMA_VERSION,
        "species_taxid": tx, "species_name": row["species_name"],
        "selected": len(hits),
        "runs": [{k: r.get(k, "") for k in RUN_FIELDS} for r in hits],
        "caveats": CAVEATS,
    }


def panel(taxid: Optional[int] = None, name: Optional[str] = None, n: int = 100,
          max_per_study: Optional[int] = None, seed: int = 0,
          **filters: Any) -> dict:
    """Select n strains: one run per isolate, spread across studies.

    Within an isolate the deepest run wins. Across studies the selection takes
    one run from each study in turn, so a species dominated by one outbreak study
    does not fill the panel with that outbreak unless there is nothing else.

    This is a sampling rule. It is not a guarantee of independence, and the
    returned metadata says so.
    """
    row = resolve_species(name=name, taxid=taxid)
    tx = int(row["species_taxid"])
    hits = _filter(_species_runs(tx), **filters)
    if not hits:
        raise EmptyResult(f'{row["species_name"]}: no runs pass the filters')

    best: dict[str, dict] = {}
    for r in hits:
        key = r["sample_accession"] or r["run_accession"]
        if key not in best or _depth(r) > _depth(best[key]):
            best[key] = r

    by_study: dict[str, list[dict]] = defaultdict(list)
    for r in best.values():
        by_study[r["study_accession"]].append(r)
    for v in by_study.values():
        v.sort(key=lambda r: (-_depth(r), r["run_accession"]))

    studies = sorted(by_study, key=lambda s: (-len(by_study[s]), s))
    random.Random(seed).shuffle(studies)

    picked: list[dict] = []
    depth_i = 0
    while len(picked) < n:
        added = False
        for s in studies:
            if depth_i < len(by_study[s]) and (max_per_study is None
                                               or depth_i < max_per_study):
                picked.append(by_study[s][depth_i])
                added = True
                if len(picked) >= n:
                    break
        if not added:
            break
        depth_i += 1

    if not picked:
        raise EmptyResult(f'{row["species_name"]}: nothing survived selection')

    return {
        "popseq_schema": SCHEMA_VERSION,
        "species_taxid": tx, "species_name": row["species_name"],
        "requested": n, "selected": len(picked),
        "runs_passing_filters": len(hits),
        "isolates_available": len(best), "studies_available": len(by_study),
        "studies_in_panel": len({r["study_accession"] for r in picked}),
        "selection_rule": ("one run per isolate, deepest wins; then round-robin "
                           "across studies. A sampling rule, not a guarantee of "
                           "independence."),
        "runs": [{k: r.get(k, "") for k in RUN_FIELDS} for r in picked],
        "caveats": CAVEATS,
    }


def describe(taxid: Optional[int] = None, name: Optional[str] = None) -> dict:
    """Everything the catalog knows about one species."""
    row = resolve_species(name=name, taxid=taxid)
    tx = int(row["species_taxid"])
    clon = next((r for r in _rows(CLONALITY_TABLE)
                 if int(r["species_taxid"]) == tx), None)
    cov = next((r for r in _rows(COVERAGE_TABLE)
                if int(r["species_taxid"]) == tx), None)
    all_runs = list(_species_runs(tx))
    excl: dict[str, int] = defaultdict(int)
    for r in all_runs:
        excl[r["exclusion"] or "none"] += 1
    return {
        "popseq_schema": SCHEMA_VERSION,
        "species_taxid": tx, "species_name": row["species_name"],
        "summary": row, "clonality": clon, "coverage": cov,
        "runs_by_exclusion_label": dict(excl),
        "ncbi_status": {k: sum(1 for r in all_runs if r["ncbi_status"] == k)
                        for k in ("in_ncbi", "ena_only", "not_checked")},
        "caveats": CAVEATS,
    }


_MATE = re.compile(r"_([12])\.f(?:ast)?q\.gz$")


def fastq_mates(fastq_ftp: str) -> tuple[str, str, list[str]]:
    """Split an ENA `fastq_ftp` value into (mate 1, mate 2, other files).

    Mates are chosen by file name (`_1` / `_2`), never by position. ENA lists some paired
    runs as three files, an unpaired `<RUN>.fastq.gz` plus `_1` and `_2`, and the unpaired
    file is not always last. A run with one file is single-end: it is returned as mate 1,
    with mate 2 empty. Any file that is neither mate is returned in the third element.
    """
    files = [f for f in (fastq_ftp or "").split(";") if f]
    mates = {}
    for f in files:
        m = _MATE.search(f)
        if m and m.group(1) not in mates:
            mates[m.group(1)] = f
    if "1" in mates and "2" in mates:
        return mates["1"], mates["2"], [f for f in files if f not in mates.values()]
    if len(files) == 1:
        return files[0], "", []
    return "", "", files


def ena_layout(fastq_ftp: str) -> str:
    """What ENA actually serves for a run: PAIRED, SINGLE, or "" (no FASTQ).

    PAIRED means files named `_1` and `_2` exist. SINGLE means exactly one file. This, not
    the submitter's `library_layout`, decides pairing (owner decision 2026-11-08): the
    label is wrong for 806 catalog runs labelled PAIRED that serve one file and 1,468
    labelled SINGLE that serve `_1` and `_2`.
    """
    m1, m2, _ = fastq_mates(fastq_ftp)
    if m1 and m2:
        return "PAIRED"
    if m1:
        return "SINGLE"
    return ""


def to_sample_sheet(result: dict) -> list[dict]:
    """nf-core sample,fastq_1,fastq_2 rows from a runs() or panel() result.

    The FASTQ columns are ENA URLs. A run that has not been fetched has no local
    path, and inventing one would be a lie. Rewrite them after fetching.

    Mates come from `fastq_mates`. A file that is neither mate (the unpaired file of a
    three-file run) is left out of the sheet. It is still in the run's `fastq_ftp`, and
    `sample_sheet_dropped` counts the runs affected.
    """
    out = []
    for r in result["runs"]:
        m1, m2, _ = fastq_mates(r.get("fastq_ftp") or "")
        out.append({
            "sample": r.get("sample_accession") or r["run_accession"],
            "fastq_1": "https://" + m1 if m1 else "",
            "fastq_2": "https://" + m2 if m2 else "",
        })
    return out


def sample_sheet_dropped(result: dict) -> int:
    """How many runs in `result` have a FASTQ file that `to_sample_sheet` leaves out."""
    return sum(1 for r in result["runs"] if fastq_mates(r.get("fastq_ftp") or "")[2])
