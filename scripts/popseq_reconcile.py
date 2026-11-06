#!/usr/bin/env python3
"""Pull NCBI's run list for a clade and reconcile it against ENA's.

Why this exists. On 2026-11-06 the same question returned 245,822 runs from the ENA
portal and 255,437 from NCBI, a 3.9 percent difference. Neither archive is known to be
a superset: GmG found 21 runs that exist at ENA and not at NCBI, and the NCBI count here
is the larger of the two. A catalog that uses one archive is therefore incomplete by an
unknown amount, in both directions.

This script pulls the NCBI side and writes the two difference lists, so the size and the
direction of the gap become measurements instead of an open question.

What each difference means, and what to do with it:

  ENA only   -> candidates for `unavailable_at_ncbi` (spec 002). These runs must be
                fetched over ENA; an NCBI route will never find them.
  NCBI only  -> runs the ENA query did not return. Could be a genuine holding
                difference, or different library_strategy/source labelling between the
                archives for the same run. NOT assumed to be either.

Checkpointed: efetch pages are written to --shard-dir as they arrive, so an interrupted
pull resumes rather than restarting. NCBI history servers expire, so a resumed run
re-issues esearch and continues by page offset.

Usage:
    uv run python scripts/popseq_reconcile.py --taxid 4751 \
        --ena cache/fungi_runs_ena_2026-11-06.tsv.gz \
        --out-dir catalog/kingdom --shard-dir cache/_runinfo_shards
"""

import argparse
import csv
import gzip
import io
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PAGE = 500

# The NCBI equivalent of the registered ENA query in popseq_survey.py.
# [Organism:exp] explodes the taxonomy the way ENA's tax_tree() does.
NCBI_TERM = 'txid{taxid}[Organism:exp] AND "wgs"[Strategy] AND "genomic"[Source] AND "illumina"[Platform]'


def get(url, params, tries=4):
    key = os.environ.get("NCBI_API_KEY")
    if key:
        params = dict(params, api_key=key)
    full = f"{url}?{urllib.parse.urlencode(params)}"
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(full, headers={"User-Agent": "PopSeq/0.1"})
            with urllib.request.urlopen(req, timeout=300) as fh:
                return fh.read().decode("utf-8", errors="replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"eutils failed after {tries} tries: {full}\n{last}")


def esearch(taxid):
    xml = get(f"{EUTILS}/esearch.fcgi", {
        "db": "sra", "term": NCBI_TERM.format(taxid=taxid),
        "retmax": 0, "usehistory": "y",
    })
    def tag(t):
        a = xml.find(f"<{t}>")
        return xml[a + len(t) + 2: xml.find(f"</{t}>")] if a >= 0 else None
    count = int(tag("Count") or 0)
    return count, tag("WebEnv"), tag("QueryKey")


def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else open(path, mode, newline="")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--taxid", type=int, required=True)
    p.add_argument("--ena", required=True, help="ENA run table (tsv or tsv.gz) to compare against")
    p.add_argument("--out-dir", default=".")
    p.add_argument("--shard-dir", required=True, help="checkpoint directory for efetch pages")
    p.add_argument("--page", type=int, default=PAGE)
    args = p.parse_args()

    os.makedirs(args.shard_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    count, webenv, qkey = esearch(args.taxid)
    print(f"NCBI reports {count} runs for taxid {args.taxid}", file=sys.stderr)
    if not count:
        return 1

    done = {int(n[5:-4]) for n in os.listdir(args.shard_dir)
            if n.startswith("page_") and n.endswith(".csv")}
    if done:
        print(f"Resuming: {len(done)} pages already on disk", file=sys.stderr)

    start = time.time()
    for i, offset in enumerate(range(0, count, args.page)):
        if i in done:
            continue
        body = get(f"{EUTILS}/efetch.fcgi", {
            "db": "sra", "query_key": qkey, "WebEnv": webenv,
            "rettype": "runinfo", "retmode": "text",
            "retstart": offset, "retmax": args.page,
        })
        tmp = os.path.join(args.shard_dir, f"page_{i:05d}.csv.tmp")
        with open(tmp, "wt") as fh:
            fh.write(body)
        os.replace(tmp, tmp[:-4])
        if i % 20 == 0:
            el = time.time() - start
            pages = count // args.page + 1
            print(f"  page {i}/{pages}  {offset + args.page}/{count}  {el:.0f}s", file=sys.stderr)
        time.sleep(0.11)

    # Collect NCBI run accessions.
    ncbi = {}
    for name in sorted(os.listdir(args.shard_dir)):
        if not (name.startswith("page_") and name.endswith(".csv")):
            continue
        with open(os.path.join(args.shard_dir, name), "rt", newline="") as fh:
            for r in csv.DictReader(fh):
                run = (r.get("Run") or "").strip()
                if run:
                    ncbi[run] = r
    print(f"NCBI distinct runs: {len(ncbi)}", file=sys.stderr)

    ena = {}
    with open_maybe_gzip(args.ena, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run = (r.get("run_accession") or "").strip()
            if run:
                ena[run] = r
    print(f"ENA distinct runs:  {len(ena)}", file=sys.stderr)

    ena_only = sorted(set(ena) - set(ncbi))
    ncbi_only = sorted(set(ncbi) - set(ena))

    def dump(path, accs, source, fields):
        with gzip.open(path, "wt", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(fields)
            for a in accs:
                row = source[a]
                w.writerow([row.get(f, "") for f in fields])

    dump(os.path.join(args.out_dir, f"{args.taxid}.ena_only_runs.csv.gz"), ena_only, ena,
         ["run_accession", "sample_accession", "study_accession", "scientific_name",
          "read_count", "base_count", "center_name", "first_public"])
    dump(os.path.join(args.out_dir, f"{args.taxid}.ncbi_only_runs.csv.gz"), ncbi_only, ncbi,
         ["Run", "BioSample", "SRAStudy", "ScientificName", "spots", "bases",
          "CenterName", "ReleaseDate", "LibraryStrategy", "LibrarySource", "Platform"])

    pct = 100 * len(ena_only) / max(len(ena), 1)
    print(f"\nENA only : {len(ena_only)} runs ({pct:.2f}% of ENA) -> unavailable_at_ncbi candidates")
    print(f"NCBI only: {len(ncbi_only)} runs")
    print(f"In both  : {len(set(ena) & set(ncbi))}")
    print(f"\nWrote both lists to {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
