#!/usr/bin/env python3
"""Inventory the lab's existing Population_Genomics subprojects.

Reads the convention those projects already use (`config.txt` declaring SAMPFILE,
FASTQFOLDER and FASTQEXT) and reports, per subproject, how many strain rows there are,
how many name an SRA run accession, and how many resolve to a file on disk.

Read-only. It never moves, renames or deletes anything: these are other people's
working directories (spec 005 section 6).

This is the diagnostic that produced the table in spec 005. It is deliberately a
reporter, not a converter: the conversion to spec 003 sample sheets is a separate
step that should not run until the numbers here have been looked at.

Known limits, stated because the output is easy to over-read:

- Resolution is a PREFIX match against a single directory listing. It does not walk
  subdirectories, so a subproject whose FileBase carries a relative path (for example
  Afumigatus_Global) reports 0 resolved while its files are present. Spec 005 section
  2 lists the four FileBase forms that a real resolver must handle.
- A row that does not resolve is not evidence the data is lost. It is evidence this
  tool did not find it.

Usage:
    uv run python scripts/popseq_pg_inventory.py
"""

import csv, os, re
D = "/bigdata/stajichlab/shared/projects/Population_Genomics"
SRA = re.compile(r'^[SED]RR\d+')
rows = []

def cfgdirs(root):
    if os.path.exists(os.path.join(root, "config.txt")):
        yield root
    try:
        for d in sorted(os.listdir(root)):
            p = os.path.join(root, d)
            if os.path.isdir(p) and os.path.exists(os.path.join(p, "config.txt")):
                yield p
    except OSError:
        pass

for sub in sorted(os.listdir(D)):
    subp = os.path.join(D, sub)
    if not os.path.isdir(subp):
        continue
    found = False
    for cd in cfgdirs(subp):
        found = True
        cfg = {}
        for line in open(os.path.join(cd, "config.txt")):
            if "=" in line and not line.lstrip().startswith("#"):
                k, v = line.strip().split("=", 1)
                cfg[k] = v
        sf = cfg.get("SAMPFILE", "samples.csv")
        sp = os.path.join(cd, sf)
        rel = os.path.relpath(cd, D)
        if not os.path.exists(sp):
            rows.append([rel, sf, "SAMPFILE-MISSING", 0, 0, 0, "-"])
            continue
        fq = os.path.join(cd, cfg.get("FASTQFOLDER", "input"))
        lines = [r for r in csv.reader(open(sp)) if r]
        if not lines:
            continue
        hdr = lines[0]
        has_hdr = any(h.strip().lower() in ("strain","filebase","runacc","sample") for h in hdr)
        data = [r for r in (lines[1:] if has_hdr else lines) if r and r[0].strip()]
        try:
            listing = os.listdir(fq)
        except OSError:
            listing = None
        n_sra = n_res = 0
        for r in data:
            base = (r[1] if len(r) > 1 else r[0]).strip()
            if SRA.match(base):
                n_sra += 1
            if listing is None:
                continue
            for cand in base.split(";"):
                cand = cand.strip()
                if not cand:
                    continue
                pref = os.path.basename(re.split(r'[\[\*]', cand)[0]).rstrip("_.-")
                if pref and any(f.startswith(pref) for f in listing):
                    n_res += 1
                    break
        rows.append([rel, sf,
                     ("|".join(h.strip() for h in hdr)[:44] if has_hdr else "NO-HEADER"),
                     len(data), n_sra, n_res,
                     str(len(listing)) if listing is not None else "NO-DIR"])
    if not found:
        rows.append([sub, "-", "NO-config.txt", 0, 0, 0, "-"])

print(f"{'subproject':<38}{'samplefile':<15}{'header':<46}{'n':>5}{'sra':>5}{'resolv':>7}{'files':>8}")
print("-" * 124)
tot = tsra = tres = 0
for r in rows:
    print(f"{r[0]:<38}{r[1]:<15}{r[2]:<46}{r[3]:>5}{r[4]:>5}{r[5]:>7}{r[6]:>8}")
    tot += r[3]; tsra += r[4]; tres += r[5]
print("-" * 124)
print(f"TOTAL strains {tot}   FileBase matches an SRA run accession {tsra}   resolved to a local file {tres}")
