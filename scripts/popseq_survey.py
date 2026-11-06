#!/usr/bin/env python3
"""Survey public short-read runs for a clade and write a PopSeq catalog entry.

Registers the queries so they are reproducible rather than retyped. Every query this
script issues is printed and is stored in the YAML it writes, so a later refresh can be
compared against the exact query that produced the previous numbers.

Two modes:

  --top-level           count runs for a whole clade, no per-run table. Cheap.
  --taxid <N>           full survey of one taxon: run table + species YAML.

Nothing is downloaded. Cataloguing a species is not a commitment to fetch it, which is
why a 25,000-isolate species can be described without any transfer.

Examples:
    uv run python scripts/popseq_survey.py --top-level --taxid 4751
    uv run python scripts/popseq_survey.py --taxid 5141 --out-root .
    uv run python scripts/popseq_survey.py --taxid 5141 --include-mutants
"""

import argparse
import csv
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_TAXNAME = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"
NCBI_DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon"

# ---------------------------------------------------------------------------
# The registered queries. Change these here, not at a call site, so that every
# species in the catalog is surveyed the same way.
# ---------------------------------------------------------------------------

# Clade-wide. tax_tree() includes every descendant taxon, which is what we want:
# a run deposited under a strain-level taxid still belongs to its species.
QUERY_RUNS = (
    'tax_tree({taxid})'
    ' AND library_strategy="WGS"'
    ' AND library_source="GENOMIC"'
    ' AND instrument_platform="ILLUMINA"'
)

RUN_FIELDS = [
    "run_accession", "study_accession", "sample_accession", "secondary_sample_accession",
    "scientific_name", "tax_id", "instrument_platform", "instrument_model",
    "library_layout", "library_strategy", "library_source", "library_selection",
    "read_count", "base_count", "first_public", "center_name", "country",
    "collection_date", "strain", "isolate", "sample_title", "study_title",
    "fastq_ftp", "fastq_aspera", "fastq_md5",
]

# Study-level exclusion. Both MATPredict and GmG arrived at these four categories
# independently (see specifications/). Spec 001: rows are LABELLED, never deleted,
# because TE_fungi_copy deliberately wants the mutant and evolved panels.
EXCLUDE_PATTERNS = {
    "mutant": r"\b(mutant|mutants|knock-?out|knockout|deletion strain|gene deletion|"
              r"transform(ant|ants)|engineered|CRISPR|RNAi|overexpress)\b",
    "evolution": r"\b(experimental evolution|evolved|serial (passage|transfer)|"
                 r"mutation accumulation|adaptive laboratory|ALE\b|passaged)\b",
    "cross": r"\b(progeny|segregant|segregants|tetrad|tetrads|recombinant inbred|"
             r"\bRIL\b|backcross|\bcross\b|offspring|meiotic)\b",
    "not_wgs": r"\b(RNA-?seq|transcriptom|ChIP|ATAC|bisulfite|methylom|Hi-?C|"
               r"metagenom|amplicon|ITS\b|16S|18S|single-?cell|exome|RAD-?seq)\b",
}

# Titles naming natural or clinical sampling are not excluded on a stray keyword.
# MATPredict's rule, kept.
RESCUE_PATTERN = (
    r"\b(surveillance|outbreak|GWAS|natural isolate|natural population|clinical isolate|"
    r"bloodstream|epidemiolog|population genom|biogeograph|wild isolate)\b"
)


def get(url, params=None, tries=3):
    """GET with retries. ENA and NCBI both fail intermittently."""
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PopSeq/0.1"})
            with urllib.request.urlopen(req, timeout=300) as fh:
                return fh.read().decode("utf-8", errors="replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"GET failed after {tries} tries: {url}\n{last}")


ENA_TAXID = "https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/"


def resolve_taxon(name):
    data = json.loads(get(ENA_TAXNAME + urllib.parse.quote(name)))
    if not data:
        raise SystemExit(f"No ENA taxon for {name!r}")
    return int(data[0]["taxId"]), data[0]["scientificName"], data[0].get("rank", "")


def canonical_name(taxid):
    """The taxon's own name, rank and lineage.

    Needed because neither source gives a clean species name on its own: the ENA run
    table reports whatever name the submitter deposited, and the NCBI reference report
    gives the strain-level organism of the assembly ("Aspergillus fumigatus Af293").
    Keying the catalog on the taxid is only useful if the name beside it is the taxon's.
    """
    try:
        d = json.loads(get(ENA_TAXID + str(taxid), tries=2))
        return d.get("scientificName"), d.get("rank", ""), d.get("lineage", "")
    except Exception:  # noqa: BLE001
        return None, "", ""


def reference_genome(taxid):
    """Reference assembly accession and size from NCBI. Returns (acc, mb, name, level).

    Measured, not assumed. Every coverage number in the catalog depends on this, so a
    species whose size cannot be fetched is recorded with a null and a reason rather
    than a guess.
    """
    for flt in ("filters.reference_only=true", "filters.assembly_level=chromosome", ""):
        url = f"{NCBI_DATASETS}/{taxid}/dataset_report?page_size=1"
        if flt:
            url += "&" + flt
        try:
            reports = json.loads(get(url, tries=2)).get("reports") or []
        except Exception:  # noqa: BLE001
            continue
        if not reports:
            continue
        r = reports[0]
        size = r.get("assembly_stats", {}).get("total_sequence_length")
        if not size:
            continue
        return (r.get("accession"),
                round(int(size) / 1e6, 2),
                r.get("organism", {}).get("organism_name"),
                r.get("assembly_info", {}).get("assembly_level"))
    return (None, None, None, None)


def ena_runs(taxid, fields=RUN_FIELDS):
    query = QUERY_RUNS.format(taxid=taxid)
    body = get(ENA_SEARCH, {
        "result": "read_run", "query": query,
        "fields": ",".join(fields), "format": "tsv", "limit": 0,
    })
    reader = csv.DictReader(io.StringIO(body), delimiter="\t")
    return query, list(reader)


def classify_exclusion(row):
    """Return an excluded_as category, or '' when the row is kept.

    Spec 001: this LABELS, it does not delete. A consumer that wants mutant panels
    filters on the column; a consumer that does not, filters the other way.
    """
    text = " ".join(filter(None, (
        row.get("study_title", ""), row.get("sample_title", ""),
        row.get("strain", ""), row.get("isolate", ""),
    )))
    if not text:
        return ""
    if re.search(RESCUE_PATTERN, text, re.I):
        return ""
    for category, pattern in EXCLUDE_PATTERNS.items():
        if re.search(pattern, text, re.I):
            return category
    return ""


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--taxid", type=int, help="NCBI/ENA taxid")
    p.add_argument("--name", help="scientific name, resolved to a taxid via ENA")
    p.add_argument("--top-level", action="store_true",
                   help="count runs only; no per-run table, no YAML")
    p.add_argument("--out-root", default=".", help="PopSeq root (default: cwd)")
    p.add_argument("--min-coverage", type=float, default=15.0)
    p.add_argument("--min-reads", type=int, default=600000)
    p.add_argument("--include-mutants", action="store_true",
                   help="do not mark mutant/evolution/cross studies as excluded")
    p.add_argument("--genome-mb", type=float,
                   help="override the genome size instead of fetching it from NCBI")
    p.add_argument("--note", default="", help="free text stored in the species YAML")
    args = p.parse_args()

    if args.name and not args.taxid:
        args.taxid, resolved, rank = resolve_taxon(args.name)
        print(f"{args.name!r} -> taxid {args.taxid} ({rank}) {resolved}", file=sys.stderr)
    if not args.taxid:
        p.error("one of --taxid or --name is required")

    if args.top_level:
        query, rows = ena_runs(args.taxid, fields=["run_accession", "library_layout"])
        paired = sum(1 for r in rows if r.get("library_layout") == "PAIRED")
        print(f"query   : {query}")
        print(f"runs    : {len(rows)}")
        print(f"paired  : {paired}")
        print(f"as of   : {time.strftime('%Y-%m-%d')}")
        return 0

    canon, rank, lineage = canonical_name(args.taxid)
    acc, mb, orgname, level = reference_genome(args.taxid)
    if args.genome_mb:
        mb, acc, level = args.genome_mb, acc or "user-supplied", level or "user-supplied"
    print(f"reference: {acc} {orgname} {level} {mb} Mb", file=sys.stderr)

    query, rows = ena_runs(args.taxid)
    print(f"ENA returned {len(rows)} runs", file=sys.stderr)
    if not rows:
        print("No runs. Nothing written.", file=sys.stderr)
        return 1

    sci = canon or rows[0].get("scientific_name", f"taxid{args.taxid}")
    base_name = sci.split()
    slug = "_".join(base_name[:2]) if len(base_name) >= 2 else sci.replace(" ", "_")
    slug = re.sub(r"[^A-Za-z0-9_]", "", slug)
    stem = f"{args.taxid}.{slug}"

    for r in rows:
        bc = int(r["base_count"]) if r.get("base_count", "").isdigit() else 0
        rc = int(r["read_count"]) if r.get("read_count", "").isdigit() else 0
        r["est_coverage"] = round(bc / (mb * 1e6), 2) if mb else ""
        r["est_gz_bytes"] = int(bc / 3.5) if bc else 0   # UNMEASURED divisor, spec 002
        r["excluded_as"] = "" if args.include_mutants else classify_exclusion(r)

        reasons = []
        if r.get("library_layout") != "PAIRED":
            reasons.append("not_paired")
        if rc < args.min_reads:
            reasons.append("too_few_reads")
        if mb and r["est_coverage"] < args.min_coverage:
            reasons.append("low_coverage")
        if not r.get("fastq_ftp"):
            reasons.append("no_fastq")
        if r["excluded_as"]:
            reasons.append(r["excluded_as"])
        r["eligible"] = "false" if reasons else "true"
        r["eligible_reason"] = ";".join(reasons)

    eligible = [r for r in rows if r["eligible"] == "true"]
    isolates = {r["sample_accession"] for r in eligible if r.get("sample_accession")}
    excluded_counts = {}
    for r in rows:
        if r["excluded_as"]:
            excluded_counts[r["excluded_as"]] = excluded_counts.get(r["excluded_as"], 0) + 1

    # One run per isolate: the smallest that still passes, so the download is as small
    # as the question allows.
    per_isolate = {}
    for r in eligible:
        s = r["sample_accession"]
        if s not in per_isolate or r["est_gz_bytes"] < per_isolate[s]["est_gz_bytes"]:
            per_isolate[s] = r
    dl_tb = sum(r["est_gz_bytes"] for r in per_isolate.values()) / 1e12

    runs_dir = os.path.join(args.out_root, "catalog", "runs")
    spec_dir = os.path.join(args.out_root, "catalog", "species")
    os.makedirs(runs_dir, exist_ok=True)
    os.makedirs(spec_dir, exist_ok=True)

    runs_path = os.path.join(runs_dir, f"{stem}.runs.csv.gz")
    cols = RUN_FIELDS + ["est_coverage", "est_gz_bytes", "excluded_as",
                         "eligible", "eligible_reason"]
    with gzip.open(runs_path, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["run_accession"]):
            w.writerow(r)

    yaml_path = os.path.join(spec_dir, f"{stem}.yaml")
    excl_yaml = "\n".join(f"    {k}: {v}" for k, v in sorted(excluded_counts.items())) or "    {}"
    with open(yaml_path, "wt") as fh:
        fh.write(f"""taxid: {args.taxid}
scientific_name: {sci}
rank: {rank or 'null'}
lineage: {lineage.strip() if lineage else 'null'}
slug: {slug}

assumed_genome_mb: {mb if mb is not None else 'null'}
assumed_genome_mb_source: {acc or 'UNRESOLVED - no reference assembly found at NCBI'}
assumed_genome_mb_assembly_level: {level or 'null'}

reference_assemblies:
  - accession: {acc or 'null'}
    organism: {orgname or 'null'}
    level: {level or 'null'}
    note: NCBI datasets reference, fetched by popseq_survey.py

query:
  source: ENA Portal API read_run
  ena_query: '{query}'
  min_coverage: {args.min_coverage}
  min_reads: {args.min_reads}
  mutants_excluded: {str(not args.include_mutants).lower()}

counts:
  last_refreshed: {time.strftime('%Y-%m-%d')}
  runs_total: {len(rows)}
  runs_eligible: {len(eligible)}
  isolates_distinct_eligible: {len(isolates)}
  estimated_download_tb_one_run_per_isolate: {round(dl_tb, 4)}
  excluded_by_category:
{excl_yaml}

projects_using: []

notes: |
  {args.note or 'Generated by popseq_survey.py. No reads have been fetched.'}
  Cataloguing is not a commitment to download. est_gz_bytes uses the base_count/3.5
  divisor inherited from the TE_fungi_copy survey, which is UNMEASURED here (spec 002).
""")

    print(f"\n{sci} (taxid {args.taxid})")
    print(f"  runs total            {len(rows)}")
    print(f"  runs eligible         {len(eligible)}")
    print(f"  distinct isolates     {len(isolates)}")
    print(f"  est download (1/iso)  {dl_tb:.3f} TB")
    if excluded_counts:
        print(f"  excluded              {excluded_counts}")
    print(f"  wrote {yaml_path}")
    print(f"  wrote {runs_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
