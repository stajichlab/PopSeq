"""The query layer. One implementation; the CLI, the MCP server and the static
JSON renderer are all thin wrappers over this.

Design notes that matter to a caller:

* Species are identified by **taxid**. A name is a convenience that resolves
  through the catalog and raises on an ambiguous match rather than guessing.
  NCBI labels one taxon with more than one name (issue #1).
* Nothing is deleted from the tables. `exclusion` labels mutant, evolved, cross
  and non-WGS material; the default filters drop it, `include_excluded=True`
  keeps it, and `only_exclusion=...` selects just that material.
* The same holds for run QC (spec 008): a finding a source marks as excluding by
  default (`exclude_statuses` in `catalog/qc/sources.csv`) makes `runs` and `panel`
  skip the run, say how many they skipped, and `include_flagged=True` keeps them.
  A curator's `curated_status=false_positive` lifts the skip, and a verdict made for
  another species name than the run's effective one (a relabel since) is stale and skips nothing.
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

SCHEMA_VERSION = "1.3.0"

ROOT = os.environ.get("POPSEQ_ROOT") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))

RUNS_TABLE = os.path.join(ROOT, "catalog", "runs", "fungi_all_runs.csv.gz")
SUMMARY_TABLE = os.path.join(ROOT, "catalog", "kingdom", "fungi_species_summary.csv.gz")
CLONALITY_TABLE = os.path.join(ROOT, "catalog", "summary", "clonality.csv.gz")
COVERAGE_TABLE = os.path.join(ROOT, "catalog", "summary", "coverage.csv.gz")
UNPLACED_TABLE = os.path.join(ROOT, "catalog", "summary", "unplaced_runs.csv.gz")
INDEX_FILE = os.path.join(ROOT, "catalog", "index.json")
QC_DIR = os.path.join(ROOT, "catalog", "qc")
RUN_LABELS_TABLE = os.path.join(ROOT, "catalog", "labels", "run_labels.csv.gz")
RUN_CURATION_TABLE = os.path.join(ROOT, "catalog", "labels", "run_curation.csv")

SPECIES_FIELDS = ["species_taxid", "species_name", "runs", "isolates", "studies",
                  "est_download_tb", "first_year", "last_year"]
UNPLACED_FIELDS = ["run_accession", "study_accession", "sample_accession", "scientific_name",
                   "tax_id", "instrument_platform", "first_public", "placement",
                   "resolved_rank", "ranked_taxid", "ranked_name", "part_species_taxids"]
RUN_FIELDS = ["run_accession", "sample_accession", "study_accession",
              "species_taxid", "species_name", "strain", "library_layout", "ena_layout",
              "read_count", "base_count", "est_depth_x", "country",
              "country_iso", "country_status", "host", "host_taxid", "host_name",
              "host_status", "collection_date", "first_public", "exclusion", "ncbi_status",
              "fastq_ftp", "fastq_aspera", "effective_species", "qc_status"]

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


def unplaced(placement: Optional[str] = None, rank: Optional[str] = None,
             contains: Optional[str] = None, limit: Optional[int] = None) -> list[dict]:
    """Runs with no species taxid, which `species`, `runs` and `panel` cannot reach.

    One row per run whose tax_id resolves above species rank, to several taxids, or
    to nothing (issue #22). `placement` and `rank` filter on the columns of the same
    name; `contains` matches the deposited scientific_name, case-insensitively. These
    runs are listed, never rolled up to an ancestor: a run filed as a genus is not
    evidence for any one species of it.
    """
    out = list(_rows(UNPLACED_TABLE))
    if placement:
        out = [r for r in out if r["placement"] == placement]
    if rank:
        out = [r for r in out if r["resolved_rank"] == rank]
    if contains:
        out = [r for r in out if contains.lower() in r["scientific_name"].lower()]
    if limit:
        out = out[:limit]
    return [{k: r.get(k, "") for k in UNPLACED_FIELDS} for r in out]


def unplaced_count() -> int:
    """How many runs have no species taxid; 0 when the table is absent."""
    if not os.path.exists(UNPLACED_TABLE):
        return 0
    return sum(1 for _ in _rows(UNPLACED_TABLE))


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


_LABELS_CACHE: dict[str, dict[str, dict]] = {}


def _labels() -> dict[str, dict]:
    """Run labels (spec 009) by run accession; empty when the table is absent."""
    path = RUN_LABELS_TABLE
    if path not in _LABELS_CACHE:
        _LABELS_CACHE[path] = {r["run_accession"]: r for r in _rows(path)} if os.path.exists(path) else {}
    return _LABELS_CACHE[path]


_CURATION_CACHE: dict[str, dict[str, dict]] = {}


def _curation() -> dict[str, dict]:
    """Run-level curator corrections (spec 009, run_curation.csv) by run accession; empty when the file is absent."""
    path = RUN_CURATION_TABLE
    if path not in _CURATION_CACHE:
        _CURATION_CACHE[path] = {r["run_accession"]: r for r in _rows(path)} if os.path.exists(path) else {}
    return _CURATION_CACHE[path]


def effective_species(r: dict, labels: Optional[dict] = None) -> str:
    """The species name to use for a run: the curator's run-level correction, else the curator's assembly-level one, else the
    authority's label, else NCBI's (spec 009).

    A label tie keeps every name, ';'-joined, as the labels table does.
    """
    run = r["run_accession"]
    cur = _curation().get(run, {}).get("curated_species", "")
    lab = (labels if labels is not None else _labels()).get(run, {})
    return cur or lab.get("curated_species") or lab.get("label_species") or r.get("species_name", "")


def _qc_sources() -> dict[str, set[str]]:
    """source_id -> the status words that skip a run by default (`exclude_statuses`, ';'-separated)."""
    path = os.path.join(QC_DIR, "sources.csv")
    if not os.path.exists(path):
        return {}
    return {r["source_id"]: {w for w in (r.get("exclude_statuses") or "").split(";") if w} for r in _rows(path)}


def _qc_rows(tx: int) -> dict[str, list[dict]]:
    """The QC rows of one species by run accession (resolved runs only)."""
    out: dict[str, list[dict]] = defaultdict(list)
    if not os.path.isdir(QC_DIR):
        return out
    for f in sorted(os.listdir(QC_DIR)):
        if f.startswith(f"{tx}.") and f.endswith(".qc.csv.gz"):
            for q in _rows(os.path.join(QC_DIR, f)):
                out[q["run_accession"]].append(q)
    return out


def _qc_view(qc_rows: list[dict], sources: dict[str, set[str]], eff: str) -> tuple[str, bool]:
    """(`check:status` words, skipped by default) for one run.

    A row skips the run when its status is one its source lists in `exclude_statuses`, unless a curator marked it
    `false_positive` or it is stale: it was made for a filed name (`filed_name`) other than the run's effective species now.
    """
    words, skip = [], False
    for q in qc_rows:
        stale = bool(q.get("filed_name")) and q["filed_name"] != eff
        words.append(f'{q.get("check", "")}:{q.get("status", "")}' + (":stale" if stale else ""))
        if (not stale and q.get("curated_status") != "false_positive"
                and q.get("status", "") in sources.get(q.get("source_id", ""), set())):
            skip = True
    return ";".join(words), skip


def _apply_qc(hits: list[dict], tx: int, include_flagged: bool) -> tuple[list[dict], list[dict]]:
    """Annotate each run with `effective_species` and `qc_status`; split off the runs QC skips by default."""
    qc, sources, labels = _qc_rows(tx), _qc_sources(), _labels()
    kept, skipped = [], []
    for r in hits:
        eff = effective_species(r, labels)
        words, skip = _qc_view(qc.get(r["run_accession"], []), sources, eff)
        r = dict(r, effective_species=eff, qc_status=words)
        (skipped if skip and not include_flagged else kept).append(r)
    return kept, skipped


def _qc_note(n: int) -> str:
    return f"{n} run(s) that pass the filters are flagged by run QC (spec 008) and skipped; --include-flagged keeps them"


def runs(taxid: Optional[int] = None, name: Optional[str] = None,
         limit: Optional[int] = None, include_flagged: bool = False, **filters: Any) -> dict:
    """Every run of one species that passes the filters (and, by default, is not skipped by run QC)."""
    row = resolve_species(name=name, taxid=taxid)
    tx = int(row["species_taxid"])
    hits = _filter(_species_runs(tx), **filters)
    if not hits:
        raise EmptyResult(f'{row["species_name"]}: no runs pass the filters')
    hits, skipped = _apply_qc(hits, tx, include_flagged)
    if not hits:
        raise EmptyResult(f'{row["species_name"]}: {_qc_note(len(skipped))}, and nothing else passes')
    if limit:
        hits = hits[:limit]
    return {
        "popseq_schema": SCHEMA_VERSION,
        "species_taxid": tx, "species_name": row["species_name"],
        "selected": len(hits),
        "qc_skipped": len(skipped),
        "qc_skipped_runs": [r["run_accession"] for r in skipped],
        "runs": [{k: r.get(k, "") for k in RUN_FIELDS} for r in hits],
        "caveats": CAVEATS,
    }


def panel(taxid: Optional[int] = None, name: Optional[str] = None, n: int = 100,
          max_per_study: Optional[int] = None, seed: int = 0,
          include_flagged: bool = False, **filters: Any) -> dict:
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
    # run QC first, so a skipped run cannot take the place of an isolate's other run or a slot in the panel
    hits, skipped = _apply_qc(hits, tx, include_flagged)
    if not hits:
        raise EmptyResult(f'{row["species_name"]}: {_qc_note(len(skipped))}, and nothing else passes')

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
        "qc_skipped": len(skipped),
        "qc_skipped_runs": [r["run_accession"] for r in skipped],
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
