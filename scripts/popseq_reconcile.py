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
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PAGE = 500

# The NCBI equivalent of the registered ENA query in popseq_survey.py.
# [Organism:exp] explodes the taxonomy the way ENA's tax_tree() does.
NCBI_TERM = 'txid{taxid}[Organism:exp] AND "wgs"[Strategy] AND "genomic"[Source] AND "illumina"[Platform]'


def get(url, params, tries=6):
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
            # NCBI returns 500 under load often enough that a short retry ladder is
            # not enough. Back off further on a server error than on a timeout.
            time.sleep(min(60, 5 * (attempt + 1) ** 2))
    # AGENTS.md section 11: the API key is never printed. urlencode puts it in the
    # URL, so the URL must be redacted before it reaches a log or a traceback.
    raise RuntimeError(f"eutils failed after {tries} tries: {redact(full)}\n{last}")


def redact(text):
    return re.sub(r"api_key=[^&\s]+", "api_key=REDACTED", str(text))


ENA_TAXNAME = "https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/"


class RankArtifact(Exception):
    """The name resolved to a taxon far above species rank."""


def resolve_taxid(name):
    """Species name -> taxid, via ENA. Returns None when the name does not resolve."""
    import json
    try:
        d = json.loads(get_plain(ENA_TAXNAME + urllib.parse.quote(name), tries=2))
        return int(d[0]["taxId"]) if d else None
    except Exception:  # noqa: BLE001
        return None


def get_plain(url, tries=3):
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PopSeq/0.1"})
            with urllib.request.urlopen(req, timeout=120) as fh:
                return fh.read().decode("utf-8", errors="replace")
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"GET failed: {url}\n{last}")


def esearch(taxid, extra=""):
    term = NCBI_TERM.format(taxid=taxid) + extra
    xml = get(f"{EUTILS}/esearch.fcgi", {
        "db": "sra", "term": term,
        "retmax": 0, "usehistory": "y",
    })
    def tag(t):
        a = xml.find(f"<{t}>")
        return xml[a + len(t) + 2: xml.find(f"</{t}>")] if a >= 0 else None
    count = int(tag("Count") or 0)
    return count, tag("WebEnv"), tag("QueryKey")


def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode, newline="") if str(path).endswith(".gz") else open(path, mode, newline="")


def fetch_window(taxid, extra, label, shard_dir, page, max_depth, depth=0, expected=None):
    """Fetch one query window, subdividing when retstart would go deep.

    Deep retstart paging into an NCBI history server degrades badly. Measured
    2026-11-06 on the whole fungal kingdom: 4.3 s per page over the first 80 pages,
    34 s per page by page 140, heading for hours. Keeping every window shallow is
    what makes this tractable, so a window over max_depth is split by publication
    year, and a year still over the limit is split by month.
    """
    count, webenv, qkey = esearch(taxid, extra)
    if count == 0:
        return 0

    # Guard against a rank artifact. The species column of the summary table is the
    # first two words of the deposited scientific_name, so a run deposited at a rank
    # above species produces a row like "Fungi", "Saccharomycotina" or "Histoplasma".
    # Resolving that name gives the taxid of the whole clade: "Fungi" resolves to 4751
    # and pulls all 255,437 fungal runs. Measured 2026-11-06, after it ran for 17
    # minutes subdividing the kingdom by year inside a job meant to do one species.
    if depth == 0 and expected and count > max(10 * expected, expected + 1000):
        raise RankArtifact(
            f"NCBI returns {count} runs for taxid {taxid} but the ENA summary expects "
            f"about {expected}. This name is almost certainly above species rank."
        )

    if count > max_depth and depth < 2:
        if depth == 0:
            parts = [(f" AND {y}[PDAT]", f"{label}.{y}") for y in range(2008, 2028)]
        else:
            year = label.rsplit(".", 1)[-1]
            parts = [(f" AND {year}/{m:02d}[PDAT]", f"{label}.{m:02d}") for m in range(1, 13)]
        total = 0
        for sub_extra, sub_label in parts:
            total += fetch_window(taxid, extra + sub_extra, sub_label,
                                  shard_dir, page, max_depth, depth + 1, expected)
        return total

    if count > max_depth:
        print(f"    {label}: {count} runs still over {max_depth} after month split; "
              f"paging deep anyway", file=sys.stderr)

    fetched = 0
    for i, offset in enumerate(range(0, count, page)):
        out = os.path.join(shard_dir, f"{label}.p{i:04d}.csv")
        if os.path.exists(out):
            fetched += page
            continue
        body = get(f"{EUTILS}/efetch.fcgi", {
            "db": "sra", "query_key": qkey, "WebEnv": webenv,
            "rettype": "runinfo", "retmode": "text",
            "retstart": offset, "retmax": page,
        })
        tmp = out + ".tmp"
        with open(tmp, "wt") as fh:
            fh.write(body)
        os.replace(tmp, out)
        fetched += page
        time.sleep(0.11)
    return min(fetched, count)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--taxid", type=int, help="single taxon to reconcile")
    p.add_argument("--species-list",
                   help="species summary CSV (gz ok) with a 'species' and 'isolates' "
                        "column; reconciles each species separately so retstart stays "
                        "shallow")
    p.add_argument("--min-isolates", type=int, default=100,
                   help="with --species-list, skip species below this (default 100)")
    p.add_argument("--max-depth", type=int, default=5000,
                   help="split a query by year, then month, above this many runs")
    p.add_argument("--ena", required=True, help="ENA run table (tsv or tsv.gz) to compare against")
    p.add_argument("--out-dir", default=".")
    p.add_argument("--shard-dir", required=True, help="checkpoint directory for efetch pages")
    p.add_argument("--page", type=int, default=PAGE)
    args = p.parse_args()

    os.makedirs(args.shard_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    # Read the species list exactly once. It may be a process substitution or a pipe,
    # which cannot be reopened; reading it a second time silently produced an empty
    # scope and a reconciliation of nothing against nothing.
    spp = []
    wanted = None
    if args.species_list:
        with open_maybe_gzip(args.species_list, "rt") as fh:
            spp = [r for r in csv.DictReader(fh)
                   if int(r.get("isolates") or 0) >= args.min_isolates]
        wanted = {r["species"] for r in spp}

    start = time.time()
    if args.species_list:
        print(f"{len(spp)} species at or above {args.min_isolates} isolates", file=sys.stderr)
        failed = []
        skipped = []
        for n, r in enumerate(spp, 1):
            name = r["species"]
            label = name.replace(" ", "_").replace("/", "_")
            marker = os.path.join(args.shard_dir, f"{label}.DONE")
            if os.path.exists(marker):
                continue
            tx = resolve_taxid(name)
            if not tx:
                print(f"[{n}/{len(spp)}] {name}: no taxid, skipped", file=sys.stderr)
                open(marker, "w").write("no taxid\n")
                continue
            try:
                got = fetch_window(tx, "", label, args.shard_dir, args.page,
                                   args.max_depth, expected=int(r.get("runs") or 0))
            except RankArtifact as exc:
                # Recorded as done, because retrying will not help. It is a naming
                # problem in the summary table, not a transfer problem.
                print(f"[{n}/{len(spp)}] {name}: SKIPPED, {exc}", file=sys.stderr)
                open(marker, "w").write(f"skipped: rank artifact; {exc}\n")
                skipped.append(name)
                continue
            except Exception as exc:  # noqa: BLE001
                # No DONE marker, so a rerun retries this species. Pages already
                # written are kept and skipped. One bad species must not end the run.
                print(f"[{n}/{len(spp)}] {name}: FAILED, will retry on a rerun: "
                      f"{redact(exc)}", file=sys.stderr)
                failed.append(name)
                continue
            open(marker, "w").write(f"taxid {tx}, about {got} runs\n")
            print(f"[{n}/{len(spp)}] {name} (txid {tx}): ~{got} runs, "
                  f"{time.time() - start:.0f}s elapsed", file=sys.stderr)
    else:
        if not args.taxid:
            p.error("one of --taxid or --species-list is required")
        got = fetch_window(args.taxid, "", f"txid{args.taxid}",
                           args.shard_dir, args.page, args.max_depth)
        print(f"taxid {args.taxid}: ~{got} runs in {time.time() - start:.0f}s", file=sys.stderr)

    # Collect NCBI run accessions.
    ncbi = {}
    for name in sorted(os.listdir(args.shard_dir)):
        if not name.endswith(".csv"):
            continue
        with open(os.path.join(args.shard_dir, name), "rt", newline="") as fh:
            for r in csv.DictReader(fh):
                run = (r.get("Run") or "").strip()
                if run:
                    ncbi[run] = r
    print(f"NCBI distinct runs: {len(ncbi)}", file=sys.stderr)

    # Scope the ENA side to the same species. Without this, --species-list compares a
    # few species of NCBI against the whole of ENA and every other species shows up as
    # "ENA only", which is true and useless.
    if wanted is not None:
        print(f"Scoping ENA to {len(wanted)} species", file=sys.stderr)

    ena = {}
    with open_maybe_gzip(args.ena, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run = (r.get("run_accession") or "").strip()
            if not run:
                continue
            if wanted is not None:
                sp = " ".join((r.get("scientific_name") or "").split()[:2])
                if sp not in wanted:
                    continue
            ena[run] = r
    print(f"ENA distinct runs:  {len(ena)}", file=sys.stderr)

    # The NCBI side can also carry runs whose species is outside the scope, because
    # txid[Organism:exp] includes subordinate taxa whose deposited name differs.
    if wanted is not None:
        drop = {k for k, v in ncbi.items()
                if " ".join((v.get("ScientificName") or "").split()[:2]) not in wanted}
        if drop:
            print(f"  dropping {len(drop)} NCBI runs outside the species scope",
                  file=sys.stderr)
            for k in drop:
                del ncbi[k]

    ena_only = sorted(set(ena) - set(ncbi))
    ncbi_only = sorted(set(ncbi) - set(ena))

    def dump(path, accs, source, fields):
        with gzip.open(path, "wt", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(fields)
            for a in accs:
                row = source[a]
                w.writerow([row.get(f, "") for f in fields])

    stem = f"{args.taxid}" if args.taxid else "fungi_byspecies"
    dump(os.path.join(args.out_dir, f"{stem}.ena_only_runs.csv.gz"), ena_only, ena,
         ["run_accession", "sample_accession", "study_accession", "scientific_name",
          "read_count", "base_count", "center_name", "first_public"])
    dump(os.path.join(args.out_dir, f"{stem}.ncbi_only_runs.csv.gz"), ncbi_only, ncbi,
         ["Run", "BioSample", "SRAStudy", "ScientificName", "spots", "bases",
          "CenterName", "ReleaseDate", "LibraryStrategy", "LibrarySource", "Platform"])

    pct = 100 * len(ena_only) / max(len(ena), 1)
    print(f"\nENA only : {len(ena_only)} runs ({pct:.2f}% of ENA) -> unavailable_at_ncbi candidates")
    print(f"NCBI only: {len(ncbi_only)} runs")
    print(f"In both  : {len(set(ena) & set(ncbi))}")
    if args.species_list and skipped:
        print(f"\n{len(skipped)} names skipped as above species rank: {skipped}",
              file=sys.stderr)
    if args.species_list and failed:
        print(f"\n{len(failed)} species failed and have no DONE marker; "
              f"rerun to retry them: {failed[:8]}", file=sys.stderr)
    print(f"\nWrote both lists to {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
