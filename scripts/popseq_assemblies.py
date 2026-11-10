#!/usr/bin/env python3
"""Fungal genome assemblies from NCBI: a dated snapshot and the catalog table built from it (specification 009).

    # 1. snapshot (network; the login node can reach NCBI): cache/ncbi_assembly_summary_fungi_{genbank,refseq}_<date>.txt.gz
    uv run python scripts/popseq_assemblies.py fetch --date 2026-11-09
    # 2. table: catalog/kingdom/fungi_assemblies.csv.gz, one row per GenBank assembly (GCA)
    uv run python scripts/popseq_assemblies.py build \\
        --genbank cache/ncbi_assembly_summary_fungi_genbank_2026-11-09.txt.gz \\
        --refseq cache/ncbi_assembly_summary_fungi_refseq_2026-11-09.txt.gz

The source is NCBI's own assembly_summary.txt for the fungi group (genomes/genbank/fungi and genomes/refseq/fungi).
It is the link from an assembly to its BioSample, so an assembly-level finding (a curated species label, spec 009)
reaches the runs of the same isolate through `sample_accession`. Assemblies that NCBI has suppressed are not in the
current summary and so are not in the table. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import os
import re
import sys
import urllib.request

URLS = {
    "genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt",
    "refseq": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt",
}
OUT_COLS = ["assembly_accession", "refseq_accession", "asm_name", "sample_accession", "study_accession", "tax_id",
            "species_taxid", "species_name", "organism_name", "strain", "isolate", "assembly_level", "genome_rep",
            "refseq_category", "relation_to_type_material", "seq_rel_date", "asm_submitter", "genome_size",
            "excluded_from_refseq"]


def snapshot_path(cache: str, db: str, date: str) -> str:
    return os.path.join(cache, f"ncbi_assembly_summary_fungi_{db}_{date}.txt.gz")


def fetch(url: str, out: str, timeout: int = 600) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "PopSeq-assemblies"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = resp.read()
    if b"assembly_accession" not in data[:2000]:
        raise SystemExit(f"{url}: not an assembly summary")
    tmp = out + ".tmp"
    with gzip.open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, out)
    rows = sum(1 for line in data.splitlines() if line and not line.startswith(b"#"))
    return {"url": url, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), "rows": rows}


def read_summary(path: str) -> list[dict]:
    """Rows of an NCBI assembly_summary.txt (plain or .gz). The header is the '#assembly_accession' line."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        lines = fh.read().split("\n")
    hdr = next(i for i, line in enumerate(lines) if line.startswith(("#assembly_accession", "# assembly_accession")))
    cols = lines[hdr].lstrip("#").strip().split("\t")
    body = [line for line in lines[hdr + 1:] if line and not line.startswith("#")]
    return list(csv.DictReader(io.StringIO("\n".join(body)), fieldnames=cols, delimiter="\t", quoting=csv.QUOTE_NONE))


def blank(v: str | None) -> str:
    """NCBI writes 'na' for a missing value; the catalog writes an empty field (AGENTS.md section 8)."""
    v = (v or "").strip()
    return "" if v.lower() == "na" else v


def strain_of(infraspecific: str) -> str:
    m = re.match(r"strain=(.*)", blank(infraspecific))
    return m.group(1).strip() if m else ""


def read_taxonomy(path: str) -> dict[str, str]:
    """species_taxid -> NCBI species name, from catalog/taxonomy_species_map.csv.gz."""
    out = {}
    with gzip.open(path, "rt", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["species_taxid"]:
                out[r["species_taxid"]] = r["species_name"]
    return out


def build(genbank: list[dict], refseq: list[dict], species_names: dict[str, str]) -> list[dict]:
    """One row per GenBank assembly; a RefSeq copy adds its GCF accession and refseq_category."""
    rs = {}
    for r in refseq:
        gca = blank(r.get("gbrs_paired_asm"))
        if gca:
            rs[gca] = r
    out = []
    for r in genbank:
        acc = r["assembly_accession"]
        f = rs.get(acc, {})
        sp = blank(r.get("species_taxid"))
        out.append({
            "assembly_accession": acc,
            "refseq_accession": f.get("assembly_accession", ""),
            "asm_name": blank(r.get("asm_name")),
            "sample_accession": blank(r.get("biosample")),
            "study_accession": blank(r.get("bioproject")),
            "tax_id": blank(r.get("taxid")),
            "species_taxid": sp,
            "species_name": species_names.get(sp, ""),
            "organism_name": blank(r.get("organism_name")),
            "strain": strain_of(r.get("infraspecific_name", "")),
            "isolate": blank(r.get("isolate")),
            "assembly_level": blank(r.get("assembly_level")),
            "genome_rep": blank(r.get("genome_rep")),
            "refseq_category": blank(f.get("refseq_category")) or blank(r.get("refseq_category")),
            "relation_to_type_material": blank(r.get("relation_to_type_material")),
            "seq_rel_date": blank(r.get("seq_rel_date")).replace("/", "-"),
            "asm_submitter": blank(r.get("asm_submitter")),
            "genome_size": blank(r.get("genome_size")),
            "excluded_from_refseq": blank(r.get("excluded_from_refseq")),
        })
    return sorted(out, key=lambda x: x["assembly_accession"])


def write_csv(path: str, rows: list[dict], cols: list[str]) -> None:
    tmp = path + ".tmp"
    with open(tmp, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as gz:
        fh = io.TextIOWrapper(gz, newline="")
        w = csv.DictWriter(fh, cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
        fh.flush()
        fh.detach()
    os.replace(tmp, path)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch", help="download the two NCBI fungi assembly summaries into the cache")
    f.add_argument("--date", required=True, help="YYYY-MM-DD, the snapshot label")
    f.add_argument("--cache", default="cache")
    b = sub.add_parser("build", help="write catalog/kingdom/fungi_assemblies.csv.gz from a snapshot")
    b.add_argument("--genbank", required=True)
    b.add_argument("--refseq", required=True)
    b.add_argument("--taxonomy", default="catalog/taxonomy_species_map.csv.gz")
    b.add_argument("--out", default="catalog/kingdom/fungi_assemblies.csv.gz")
    a = ap.parse_args(argv)

    if a.cmd == "fetch":
        os.makedirs(a.cache, exist_ok=True)
        for db, url in URLS.items():
            out = snapshot_path(a.cache, db, a.date)
            info = fetch(url, out)
            with open(out + ".query.txt", "w") as fh:
                fh.write("".join(f"{k}\t{v}\n" for k, v in dict(info, date=a.date).items()))
            print(f"{db}: {info['rows']} assemblies, {info['bytes']} bytes -> {out}")
        return 0

    rows = build(read_summary(a.genbank), read_summary(a.refseq), read_taxonomy(a.taxonomy))
    write_csv(a.out, rows, OUT_COLS)
    n_bs = sum(1 for r in rows if r["sample_accession"])
    print(f"{a.out}: {len(rows)} assemblies, {n_bs} with a BioSample, "
          f"{sum(1 for r in rows if r['refseq_accession'])} with a RefSeq copy")
    return 0


if __name__ == "__main__":
    sys.exit(main())
