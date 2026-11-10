"""Plan-time lookups shared by `popseq fetch` and benchmark 7.1.

A fetch needs numbers the catalog does not hold: ENA file sizes and MD5s, and NCBI's `spots`
and `spots_with_mates`. They are looked up once, before the transfer, and stored (spec 002:
a compute node may have no internet, so the job checks against a stored value).
"""

from __future__ import annotations

import csv
import gzip
import io

import requests

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
NCBI_RUNINFO = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={run}"


def open_table(path):
    """Read .csv, .csv.gz or .csv.zst transparently (AGENTS.md section 8)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    if path.endswith(".zst"):
        import zstandard

        return io.TextIOWrapper(
            zstandard.ZstdDecompressor().stream_reader(open(path, "rb")), newline=""
        )
    return open(path, newline="")


def ena_lookup(accessions):
    """read_count, fastq_bytes, fastq_md5, URLs for each run, from the ENA portal."""
    out = {}
    fields = "run_accession,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera"
    for i in range(0, len(accessions), 200):
        chunk = accessions[i : i + 200]
        r = requests.post(
            ENA_SEARCH,
            data={
                "result": "read_run",
                "includeAccessions": ",".join(chunk),
                "fields": fields,
                "format": "tsv",
                "limit": 0,
            },
            timeout=120,
        )
        r.raise_for_status()
        for row in csv.DictReader(io.StringIO(r.text), delimiter="\t"):
            out[row["run_accession"]] = row
    return out


def ncbi_counts(run):
    """(status, spots, spots_with_mates) from NCBI runinfo.

    `status` is `ok`, `absent` (NCBI answered and has no such run: the run exists only at ENA)
    or `error` (no answer: the caller must not conclude anything). Keeping `absent` apart from
    `error` matters because `unavailable_at_ncbi` is permanent and a network error is not.
    """
    try:
        r = requests.get(NCBI_RUNINFO.format(run=run), timeout=60)
        r.raise_for_status()
        rows = list(csv.DictReader(io.StringIO(r.text)))
    except (requests.RequestException, csv.Error):
        return "error", None, None
    if not rows or not rows[0].get("spots"):
        return "absent", None, None
    return "ok", int(rows[0]["spots"] or 0), int(rows[0]["spots_with_mates"] or 0)


def ncbi_spots(run):
    """(spots, spots_with_mates), or (None, None) when NCBI has no answer for the run."""
    status, spots, mates = ncbi_counts(run)
    return (spots, mates) if status == "ok" else (None, None)
