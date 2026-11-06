#!/usr/bin/env python3
"""Resolve any taxid to its species-rank ancestor, from the NCBI taxonomy dump.

Issue #1. The catalog was keyed on the first two words of the deposited
`scientific_name`. That key is wrong in two ways, and both corrupted real results:

  rank artifacts  A run deposited above species rank produces a row that looks like a
                  species. "Fungi" resolved to taxid 4751 and pulled the whole kingdom
                  into a job meant to survey one species.
  synonyms        One taxon under two names splits into two rows. NCBI calls 1,062 runs
                  "Candida parapsilosis" and 575 "Lodderomyces parapsilosis"; ENA calls
                  them all Lodderomyces. A name-based comparison reported 844 of them as
                  missing from NCBI when they were present under the other name.

A taxid is stable under renaming and carries a rank, so both problems disappear when the
key is the species-rank ancestor of the run's own taxid.

Builds from nodes.dmp and names.dmp. Caches the fungal subtree as a small CSV so the
300 MB dumps are read once.

Usage:
    uv run python scripts/popseq_taxonomy.py --build --root 4751
    uv run python scripts/popseq_taxonomy.py --lookup 746128 5476 4751
"""

import argparse
import csv
import gzip
import os
import sys

DUMP = "cache/taxdump"
CACHE = "catalog/taxonomy_species_map.csv.gz"

# Ranks at or below species. A strain or subspecies rolls up to its species; a genus
# or anything above does not roll down, and is reported as its own rank.
BELOW_SPECIES = {"subspecies", "varietas", "forma", "strain", "serotype", "serogroup",
                 "biotype", "genotype", "morph", "pathogroup", "forma specialis",
                 "isolate", "no rank", "clade"}


def load_nodes(path):
    parent, rank = {}, {}
    with open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t|\t")
            t = int(f[0])
            parent[t] = int(f[1])
            rank[t] = f[2]
    return parent, rank


def load_names(path, keep):
    names = {}
    with open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t|\t")
            if not f[3].startswith("scientific name"):
                continue
            t = int(f[0])
            if t in keep:
                names[t] = f[1]
    return names


def species_of(taxid, parent, rank):
    """Walk up to the species-rank ancestor. Returns (species_taxid, how)."""
    seen, cur = set(), taxid
    while cur and cur not in seen:
        seen.add(cur)
        r = rank.get(cur)
        if r == "species":
            return cur, ("self" if cur == taxid else "rolled_up")
        if r is not None and r not in BELOW_SPECIES:
            # Reached genus or above without passing a species node. This taxid IS a
            # higher rank, so there is no species to roll up to. That is the "Fungi"
            # and "Saccharomycotina" case, and it must be reported, not guessed.
            return None, f"above_species:{r}"
        nxt = parent.get(cur)
        if nxt == cur:
            break
        cur = nxt
    return None, "unresolved"


def build(root, out):
    parent, rank = load_nodes(os.path.join(DUMP, "nodes.dmp"))
    print(f"nodes: {len(parent)}", file=sys.stderr)

    # Everything under root.
    children = {}
    for t, p in parent.items():
        children.setdefault(p, []).append(t)
    sub, stack = set(), [root]
    while stack:
        t = stack.pop()
        if t in sub:
            continue
        sub.add(t)
        stack.extend(children.get(t, ()))
    print(f"taxa under {root}: {len(sub)}", file=sys.stderr)

    names = load_names(os.path.join(DUMP, "names.dmp"), sub)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    counts = {}
    with gzip.open(out, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["taxid", "name", "rank", "species_taxid", "species_name", "resolution"])
        for t in sorted(sub):
            sp, how = species_of(t, parent, rank)
            counts[how.split(":")[0]] = counts.get(how.split(":")[0], 0) + 1
            w.writerow([t, names.get(t, ""), rank.get(t, ""), sp or "",
                        names.get(sp, "") if sp else "", how])
    print(f"wrote {out}", file=sys.stderr)
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {k:<16} {v}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--build", action="store_true")
    p.add_argument("--root", type=int, default=4751)
    p.add_argument("--out", default=CACHE)
    p.add_argument("--lookup", nargs="+", type=int)
    args = p.parse_args()

    if args.build:
        build(args.root, args.out)
        return 0
    if args.lookup:
        with gzip.open(args.out, "rt") as fh:
            m = {int(r["taxid"]): r for r in csv.DictReader(fh)}
        for t in args.lookup:
            r = m.get(t)
            print(f"{t}: {r}" if r else f"{t}: not in map")
        return 0
    p.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
