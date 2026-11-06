#!/usr/bin/env python3
"""Convert the lab's Population_Genomics subprojects into PopSeq panels.

Implements spec 005 section 6. Reads the convention those projects already use
(config.txt declaring SAMPFILE, FASTQFOLDER, FASTQEXT), resolves each sample's
FileBase to real files on disk, classifies its provenance, and writes a panel in
the spec 003 form.

READ ONLY against the subprojects. It never moves, renames or deletes anything
there; those are other people's working directories. Everything it writes goes
under the PopSeq panels directory.

FileBase means four different things across these subprojects (spec 005 section 2),
and all four are handled:

    SRR1636746                                  an SRA run accession
    100D-1_S34_L001_R[12]_001.fastq.gz          a glob with a mate wildcard
    22687-1_ACTTGA_L001                         a bare prefix, no extension
    AF100.run1/lane1/x_R[12]_001.fastq.gz;...   several relative paths, one per lane

Checksums: computed only for local_only samples, which are the irreplaceable ones,
and only with --checksum. Hashing every FASTQ in these subprojects would read
terabytes to tell us something we only need for data that cannot be refetched.

Usage:
    uv run python scripts/popseq_inventory.py --list
    uv run python scripts/popseq_inventory.py --subproject B_dendrobatidis
    uv run python scripts/popseq_inventory.py --all --checksum
"""

import argparse
import csv
import hashlib
import os
import re
import sys

PG_ROOT = "/bigdata/stajichlab/shared/projects/Population_Genomics"
SRA_RE = re.compile(r"^[SED]RR\d+$")
# Mate markers, most specific first. _R1_ before _1 so that a name like
# sample_R1_001 is not split on the 1 in 001.
MATE_RE = [
    re.compile(r"^(?P<stem>.+?)[._]R(?P<mate>[12])([._].*)?$"),
    re.compile(r"^(?P<stem>.+?)[._](?P<mate>[12])([._].*)?$"),
]
FASTQ_EXT = re.compile(r"\.(fastq|fq)(\.gz|\.bz2|\.zst)?$", re.I)


def read_config(path):
    cfg = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                cfg[k.strip()] = v.strip()
    return cfg


def index_fastq(root):
    """Walk once and return [(relpath, abspath, basename)] for every read file."""
    out = []
    if not os.path.isdir(root):
        return out
    for dirpath, _, files in os.walk(root, followlinks=True):
        for f in files:
            if FASTQ_EXT.search(f):
                ap = os.path.join(dirpath, f)
                out.append((os.path.relpath(ap, root), ap, f))
    return out


def split_mate(basename):
    """Return (stem, mate) where mate is '1', '2' or None."""
    name = FASTQ_EXT.sub("", basename)
    for rx in MATE_RE:
        m = rx.match(name)
        if m:
            return m.group("stem"), m.group("mate")
    return name, None


def candidates_for(filebase, index):
    """All index entries that could belong to this FileBase.

    Handles the four forms by reducing each to a prefix test against the stable
    part of the name. A glob such as R[12] is cut at the bracket, because the
    characters after it differ between mates in some layouts and not others.
    """
    hits = []
    for piece in str(filebase).split(";"):
        piece = piece.strip()
        if not piece:
            continue
        stable = re.split(r"[\[\*\?]", piece)[0].rstrip("_.-")
        if not stable:
            continue
        base_stable = os.path.basename(stable)
        has_dir = "/" in piece
        for rel, ap, bn in index:
            if has_dir:
                if rel.startswith(stable):
                    hits.append((rel, ap, bn))
            elif bn.startswith(base_stable):
                hits.append((rel, ap, bn))
    # Stable order, no duplicates.
    seen, out = set(), []
    for h in sorted(hits):
        if h[1] not in seen:
            seen.add(h[1])
            out.append(h)
    return out


def pair_up(hits):
    """Group resolved files into (r1, r2) pairs. Unpaired files become (file, '')."""
    groups = {}
    for rel, ap, bn in hits:
        stem, mate = split_mate(bn)
        key = os.path.join(os.path.dirname(rel), stem)
        groups.setdefault(key, {})[mate or "single"] = ap
    pairs = []
    for key in sorted(groups):
        g = groups[key]
        if "1" in g:
            pairs.append((g["1"], g.get("2", "")))
        elif "single" in g:
            pairs.append((g["single"], ""))
        elif "2" in g:
            # Mate 2 with no mate 1 is a real problem, not something to paper over.
            pairs.append((g["2"], ""))
    return pairs


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def subproject_dirs(root):
    for sub in sorted(os.listdir(root)):
        p = os.path.join(root, sub)
        if not os.path.isdir(p):
            continue
        if os.path.exists(os.path.join(p, "config.txt")):
            yield sub, p
        else:
            try:
                for d in sorted(os.listdir(p)):
                    q = os.path.join(p, d)
                    if os.path.isdir(q) and os.path.exists(os.path.join(q, "config.txt")):
                        yield f"{sub}/{d}", q
            except OSError:
                pass


def process(name, cfgdir, out_root, do_checksum):
    cfg = read_config(os.path.join(cfgdir, "config.txt"))
    sheet = os.path.join(cfgdir, cfg.get("SAMPFILE", "samples.csv"))
    if not os.path.exists(sheet):
        return {"subproject": name, "status": "SAMPFILE missing", "n": 0}

    fqroot = os.path.join(cfgdir, cfg.get("FASTQFOLDER", "input"))
    index = index_fastq(fqroot)

    with open(sheet, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r and r[0].strip()]
    if not rows:
        return {"subproject": name, "status": "empty sheet", "n": 0}

    header = [h.strip().lower() for h in rows[0]]
    known = {"strain", "filebase", "runacc", "sample", "biosample", "sra"}
    has_header = any(h in known for h in header)
    data = rows[1:] if has_header else rows
    col = {h: i for i, h in enumerate(header)} if has_header else {}

    def field(r, *names):
        for nm in names:
            i = col.get(nm)
            if i is not None and i < len(r) and r[i].strip():
                return r[i].strip()
        return ""

    samples, meta, manifest = [], [], []
    counts = {}
    for r in data:
        strain = field(r, "strain", "sample") or (r[0].strip() if r else "")
        fb = field(r, "filebase") or (r[1].strip() if len(r) > 1 else "")
        acc = field(r, "runacc", "sra")
        if not acc and SRA_RE.match(fb):
            acc = fb
        if not strain:
            continue

        hits = candidates_for(fb or strain, index)
        pairs = pair_up(hits)

        if pairs and acc:
            prov = "archived_local"
        elif pairs:
            prov = "local_only"
        elif acc:
            prov = "archived"
        else:
            prov = "listed_absent"
        counts[prov] = counts.get(prov, 0) + 1

        for r1, r2 in pairs:
            samples.append({"sample": strain, "fastq_1": r1, "fastq_2": r2})
        if not pairs:
            samples.append({"sample": strain, "fastq_1": "", "fastq_2": ""})

        meta.append({
            "sample": strain, "biosample": field(r, "biosample"),
            "run_accession": acc, "provenance": prov,
            "source_subproject": name, "source_sheet": os.path.relpath(sheet, PG_ROOT),
            "filebase": fb, "n_files": sum(1 for p in pairs for f in p if f),
            "curated_notes": "",
        })

        for r1, r2 in pairs:
            for f in (r1, r2):
                if not f:
                    continue
                st = os.stat(f)
                manifest.append({
                    "local_id": f"{name}/{strain}",
                    "path": f, "bytes": st.st_size,
                    "mtime": int(st.st_mtime),
                    "sha256": sha256(f) if (do_checksum and prov == "local_only") else "",
                })

    outdir = os.path.join(out_root, "panels", name.replace("/", "__"))
    os.makedirs(outdir, exist_ok=True)

    def dump(fn, rowset, cols):
        with open(os.path.join(outdir, fn), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rowset)

    dump("samples.csv", samples, ["sample", "fastq_1", "fastq_2"])
    dump("metadata.csv", meta, ["sample", "biosample", "run_accession", "provenance",
                                "source_subproject", "source_sheet", "filebase",
                                "n_files", "curated_notes"])
    dump("local_manifest.csv", manifest, ["local_id", "path", "bytes", "mtime", "sha256"])
    return {"subproject": name, "status": "ok", "n": len(meta),
            "files": len(manifest), "fastq_indexed": len(index), **counts}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pg-root", default=PG_ROOT)
    p.add_argument("--out-root", default=".")
    p.add_argument("--subproject", action="append", help="repeatable; default is --all")
    p.add_argument("--all", action="store_true")
    p.add_argument("--list", action="store_true", help="list subprojects and stop")
    p.add_argument("--checksum", action="store_true",
                   help="sha256 for local_only files, which cannot be refetched")
    args = p.parse_args()

    found = list(subproject_dirs(args.pg_root))
    if args.list:
        for nm, d in found:
            print(f"{nm:<42} {d}")
        print(f"\n{len(found)} subprojects with a config.txt")
        return 0

    todo = found if (args.all or not args.subproject) else \
        [(nm, d) for nm, d in found if nm in set(args.subproject)]
    if not todo:
        print("No matching subproject.", file=sys.stderr)
        return 1

    results = []
    for nm, d in todo:
        res = process(nm, d, args.out_root, args.checksum)
        results.append(res)
        extra = " ".join(f"{k}={v}" for k, v in res.items()
                         if k in ("archived", "archived_local", "local_only", "listed_absent"))
        print(f"{nm:<42} {res['status']:<16} n={res['n']:<5} files={res.get('files',0):<6} {extra}")

    tot = {}
    for r in results:
        for k in ("archived", "archived_local", "local_only", "listed_absent"):
            tot[k] = tot.get(k, 0) + r.get(k, 0)
    print("\nTOTAL " + "  ".join(f"{k}={v}" for k, v in tot.items()))
    print(f"samples total {sum(r['n'] for r in results)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
