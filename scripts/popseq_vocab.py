#!/usr/bin/env python3
"""Build the vocabularies behind the normalised country and host fields (issue #23).

  hosts      catalog/vocab/host_taxa.csv.gz: every distinct host string in the
             BioSample table, resolved to an NCBI taxid where it resolves.
             Needs the NCBI taxonomy dump (names.dmp). No network.
  countries  catalog/vocab/country_unresolved.csv: the country strings that
             pycountry plus catalog/vocab/country_aliases.csv do not resolve,
             most frequent first. This is the curator's work list.

Host resolution, in order, for each distinct host key (popseq.normalise.key):

  1. missing term                       -> missing
  2. catalog/vocab/host_overrides.csv   -> the curator's status; `resolve_as`
                                           names a taxon, resolved here through
                                           names.dmp so nobody types a taxid
  3. exact NCBI name match              -> ok (scientific name preferred; a key
                                           that matches several taxa is ambiguous)
     also tried: the text in brackets, the text without them, and the part
     before " : " ("Barley (Hordeum vulgare)", "Malus x domestica : Cox")
  4. first two words, scientific name   -> ok, method binomial
     ("Triticum aestivum cv. Chinese Spring" -> Triticum aestivum)
  5. "<Genus> sp." / "spp."             -> ok, the genus, method genus_sp
  6. the submitters' own host_taxid     -> ok, method submitter_taxid, only
     when every sample with this key that gave a taxid gave the same one
  7. otherwise                          -> unresolved

    uv run python scripts/popseq_vocab.py hosts --taxdump /path/to/taxdump
    uv run python scripts/popseq_vocab.py countries
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from popseq import normalise as nz  # noqa: E402
from popseq.query import ROOT, RUNS_TABLE  # noqa: E402

BIOSAMPLES = os.path.join(ROOT, "catalog", "kingdom", "fungi_biosamples.csv.gz")
# First non-empty wins. `lab_host` is deliberately absent: it is the organism
# used to propagate a sample in the lab, not where the isolate came from.
HOST_COLUMNS = ["host", "host_plant", "specific hosts", "host_common_name"]
NAME_CLASSES = {"scientific name", "synonym", "equivalent name", "genbank common name",
                "common name", "genbank synonym"}


def _open(path):
    return gzip.open(path, "rt", newline="") if path.endswith(".gz") else open(path, newline="")


def biosample_hosts(path=BIOSAMPLES):
    """biosample -> (raw host, submitter host_taxid)."""
    out = {}
    with _open(path) as fh:
        if not fh.readline().startswith("#"):
            fh.seek(0)
        for r in csv.DictReader(fh):
            raw = next((r[c].strip() for c in HOST_COLUMNS if (r.get(c) or "").strip()), "")
            tx = (r.get("host_taxid") or "").strip()
            out[r["biosample"]] = (raw, tx if tx.isdigit() else "")
    return out


def _read_names(dmp, want_keys, want_taxids):
    """Stream names.dmp once. Returns key -> {(taxid, class)} and taxid -> sci name."""
    hits = defaultdict(set)
    sci = {}
    with open(dmp, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\t|\n").split("\t|\t")
            if len(f) < 4 or f[3] not in NAME_CLASSES:
                continue
            tx = f[0]
            if f[3] == "scientific name":
                sci[tx] = f[1]
            k = nz.key(f[1])
            if k in want_keys:
                hits[k].add((tx, f[3]))
    return hits, sci


def _pick(matches):
    """One taxid from a set of (taxid, class), or None if ambiguous."""
    s = {t for t, c in matches if c == "scientific name"}
    if len(s) == 1:
        return next(iter(s)), "name_exact"
    a = {t for t, _ in matches}
    if len(a) == 1:
        return next(iter(a)), "name_exact"
    return None, "ambiguous"


def cmd_hosts(a):
    hosts = biosample_hosts(a.biosamples)
    n_by_key = Counter()
    example = {}
    submitter = defaultdict(Counter)
    for raw, tx in hosts.values():
        if not raw:
            continue
        for part in raw.split(";"):
            k = nz.key(part)
            if not k:
                continue
            n_by_key[k] += 1
            example.setdefault(k, part.strip())
            if tx:
                submitter[k][tx] += 1
    overrides = {}
    for r in nz._load_csv(a.overrides):
        st = r["host_status"]
        if st not in nz.HOST_STATUSES:
            sys.exit(f"{a.overrides}: bad host_status {st!r} for {r['host_key']!r}")
        overrides[nz.key(r["host_key"])] = (st, nz.key(r.get("resolve_as") or ""))

    def binomial(k):
        w = k.split()
        return " ".join(w[:2]) if len(w) >= 3 else None

    def genus_sp(k):
        m = re.fullmatch(r"([a-z][a-z-]+) (sp|spp|sp nov|species)\.?( .*)?", k)
        return m.group(1) if m else None

    def variants(k):
        """Other readings of the same string, most specific first."""
        out = re.findall(r"\(([^)]+)\)", k)
        out.append(re.sub(r"\s*\([^)]*\)", "", k))
        out.append(k.split(" : ")[0])
        return [nz.key(v) for v in out if nz.key(v) and nz.key(v) != k]

    want = set(n_by_key) | {r for _, r in overrides.values() if r}
    for k in list(n_by_key):
        for v in [k] + variants(k):
            want.add(v)
            if b := binomial(v):
                want.add(b)
            if g := genus_sp(v):
                want.add(g)
    hits, sci = _read_names(os.path.join(a.taxdump, "names.dmp"), want,
                            {t for c in submitter.values() for t in c})

    def resolve(k):
        if k in hits:
            t, m = _pick(hits[k])
            return (t, m) if t else (None, "ambiguous")
        return None, ""

    rows, status_n = [], Counter()
    for k, n in n_by_key.most_common():
        tx, method, st = "", "", "unresolved"
        if nz.is_missing(k):
            st, method = "missing", "missing_term"
        elif k in overrides:
            st, target = overrides[k]
            method = "override"
            if target:
                t, _ = resolve(target)
                if not t:
                    sys.exit(f"override {k!r}: resolve_as {target!r} does not resolve")
                tx = t
        else:
            for v in [k] + variants(k):
                t, m = resolve(v)
                if t:
                    tx, method, st = t, m if v == k else "variant", "ok"
                    break
                if m == "ambiguous" and v == k:
                    st, method = "ambiguous", "ambiguous"
            if st != "ok":
                for v in [k] + variants(k):
                    b, g = binomial(v), genus_sp(v)
                    if b and b in hits:
                        t, _ = _pick({m for m in hits[b] if m[1] == "scientific name"}
                                     or hits[b])
                        if t:
                            tx, method, st = t, "binomial", "ok"
                            break
                    if g and g in hits:
                        t, _ = _pick(hits[g])
                        if t:
                            tx, method, st = t, "genus_sp", "ok"
                            break
        if st in ("unresolved", "ambiguous") and submitter.get(k):
            c = submitter[k]
            if len(c) == 1 and next(iter(c)) in sci:
                tx, method, st = next(iter(c)), "submitter_taxid", "ok"
        status_n[st] += n
        rows.append({"host_key": k, "host_taxid": tx, "host_name": sci.get(tx, ""),
                     "host_status": st, "method": method, "n_samples": n,
                     "example_raw": example[k]})
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with gzip.open(a.out, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tot = sum(n_by_key.values())
    print(f"wrote {a.out}: {len(rows)} distinct host keys over {tot} sample values",
          file=sys.stderr)
    for s, n in status_n.most_common():
        print(f"  {s}: {n} ({100 * n / tot:.1f}%)", file=sys.stderr)


def cmd_countries(a):
    c = Counter()
    with _open(a.runs) as fh:
        for r in csv.DictReader(fh):
            v = r["country"]
            if nz.country(v)["country_status"] == "unresolved":
                c[v] += 1
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["country_raw", "runs"])
        w.writerows(c.most_common())
    print(f"wrote {a.out}: {len(c)} unresolved strings, {sum(c.values())} runs",
          file=sys.stderr)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("hosts")
    h.add_argument("--taxdump", required=True, help="directory holding names.dmp")
    h.add_argument("--biosamples", default=BIOSAMPLES)
    h.add_argument("--overrides", default=nz.HOST_OVERRIDES)
    h.add_argument("--out", default=nz.HOST_TAXA)
    c = sub.add_parser("countries")
    c.add_argument("--runs", default=RUNS_TABLE)
    c.add_argument("--out", default=os.path.join(nz.VOCAB, "country_unresolved.csv"))
    a = p.parse_args(argv)
    {"hosts": cmd_hosts, "countries": cmd_countries}[a.cmd](a)


if __name__ == "__main__":
    main()
