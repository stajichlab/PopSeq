#!/usr/bin/env python3
"""Import curated species labels for genome assemblies from an authority (specification 009).

    uv run python scripts/popseq_labels_import.py --source-id fusariumid_gmg_f1_2026-06-23 --format gmg_fusariumid \\
        --input genomes_labels.tsv --complexes species_complexes.tsv \\
        --authority FusariumID --authority-release 2026-06-23 --project GmG --tool "..." --method "..." \\
        --reference "..." --date 2026-11-09 --refresh "2-3 times a year, after a FusariumID release"

Writes, under --labels-dir (default catalog/labels):
    assembly_labels.csv.gz  one row per assembly and source: the authority's species, its complex, the evidence
    run_labels.csv.gz       the same labels carried to every run whose BioSample is the assembly's BioSample
    sources.csv             one row per source

Each assembly is resolved through catalog/kingdom/fungi_assemblies.csv.gz (GCA or GCF -> GCA, BioSample, taxid,
NCBI species). The rows of --source-id are replaced; rows of other sources and every curated_* column are kept,
carried over by (assembly_accession, source_id). Input may be CSV or TSV, plain, .gz or .zst. Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

LABEL_COLS = ["assembly_accession", "sample_accession", "tax_id", "species_taxid", "ncbi_species", "source_id",
              "label_species", "label_complex", "label_status", "label_agreement", "nearest_species", "distance", "differences",
              "aligned_bp", "margin", "runner_up", "loci", "curated_species", "curated_note"]
RUN_COLS = ["run_accession", "sample_accession", "assembly_accession", "source_id", "ncbi_species", "label_species",
            "label_complex", "label_status", "label_agreement", "curated_species"]
SOURCE_COLS = ["source_id", "authority", "authority_release", "project", "tool", "method", "call_rule", "reference",
               "date", "refresh", "assemblies", "runs", "notes"]
CURATED = ("curated_species", "curated_note")
ACC = re.compile(r"^(GC[AF]_\d+\.\d+)")


def gmg_fusariumid(r: dict, max_dist: float) -> dict:
    """GmG bench F1 genomes_labels.tsv (bench/F1_fusariumid/genome_resolution.py): FusariumID in-silico barcoding."""
    m = ACC.match(r["genome"])
    d = r.get("fid_dist", "")
    near = r.get("fid_call", "")
    label = near if near and d and float(d) <= max_dist else ""
    return {"accession": m.group(1) if m else r["genome"], "label_species": label, "nearest_species": near,
            "distance": d, "differences": r.get("fid_diffs", ""), "aligned_bp": r.get("fid_len", ""),
            "margin": r.get("fid_margin", ""), "runner_up": r.get("fid_other", ""), "loci": r.get("n_loci_hit", "")}


# producer format -> (row parser, the rule that turns the evidence into a label, in words)
FORMATS = {
    "gmg_fusariumid": (gmg_fusariumid, ("label when the nearest FusariumID strain over the concatenated loci is within "
                                        "distance {max_dist}; a tie between species keeps every name, ';'-joined")),
}


def open_text(path: str):
    """A text handle for a plain, .gz or .zst file."""
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), newline="")
    if path.endswith(".zst"):
        p = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True)
        return io.StringIO(p.stdout.decode(), newline="")
    return open(path, newline="")


def read_table(path: str) -> list[dict]:
    """CSV or TSV; lines starting with '#' before the header are comments."""
    with open_text(path) as fh:
        lines = [line for line in fh.read().split("\n") if line and not line.startswith("#")]
    first = lines[0] if lines else ""
    return list(csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t" if first.count("\t") > first.count(",") else ","))


def write_gz(path: str, rows: list[dict], cols: list[str]) -> None:
    tmp = path + ".tmp"
    with open(tmp, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as gz:
        fh = io.TextIOWrapper(gz, newline="")
        w = csv.DictWriter(fh, cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
        fh.flush()
        fh.detach()
    os.replace(tmp, path)


def status_of(label: str) -> str:
    if not label:
        return "no_call"
    return "tie" if ";" in label else "confident"


def agreement_of(ncbi: str, label: str) -> str:
    """How the authority's label relates to the NCBI species of the assembly."""
    if not label:
        return "no_call"
    names = label.split(";")
    if ncbi in names:
        return "same" if len(names) == 1 else "tie_includes_ncbi"
    if not ncbi or ncbi.endswith(" sp.") or " sp. " in ncbi or " cf. " in ncbi:
        return "ncbi_unnamed"
    return "differs"


def complex_of(label: str, complexes: dict[str, str]) -> str:
    """The complex all the label's names share; empty when they differ or one is not in the table."""
    if not label:
        return ""
    cs = {complexes.get(n, "") for n in label.split(";")}
    return cs.pop() if len(cs) == 1 else ""


def read_assemblies(path: str) -> tuple[dict, dict]:
    """GCA -> row, and GCF -> GCA, from catalog/kingdom/fungi_assemblies.csv.gz."""
    by_gca, gcf = {}, {}
    for r in read_table(path):
        by_gca[r["assembly_accession"]] = r
        if r.get("refseq_accession"):
            gcf[r["refseq_accession"]] = r["assembly_accession"]
    return by_gca, gcf


def runs_by_sample(path: str, samples: set[str]) -> dict[str, list[str]]:
    """BioSample -> runs, from the PopSeq run table; a multi-valued sample_accession ('S1;S2') counts for each."""
    out = defaultdict(list)
    with open_text(path) as fh:
        for r in csv.DictReader(fh):
            for s in (r.get("sample_accession") or "").split(";"):
                if s in samples:
                    out[s].append(r["run_accession"])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--source-id", required=True)
    ap.add_argument("--format", required=True, choices=sorted(FORMATS))
    ap.add_argument("--input", required=True, help="the producer table")
    ap.add_argument("--max-dist", type=float, default=0.01, help="largest distance that still gives a label")
    ap.add_argument("--complexes", default="", help="species -> complex table (columns species, complex); optional")
    ap.add_argument("--assemblies", default="catalog/kingdom/fungi_assemblies.csv.gz")
    ap.add_argument("--runs", default="catalog/runs/fungi_all_runs.csv.gz")
    ap.add_argument("--labels-dir", default="catalog/labels")
    for k in ("authority", "authority-release", "project", "tool", "date"):
        ap.add_argument(f"--{k}", required=True)
    for k in ("method", "reference", "refresh", "notes"):
        ap.add_argument(f"--{k}", default="")
    a = ap.parse_args(argv)

    parse, rule = FORMATS[a.format]
    by_gca, gcf = read_assemblies(a.assemblies)
    complexes = {r["species"]: r["complex"] for r in read_table(a.complexes)} if a.complexes else {}
    new, missing = [], 0
    for r in read_table(a.input):
        p = parse(r, a.max_dist)
        acc = gcf.get(p["accession"], p["accession"])
        asm = by_gca.get(acc)
        if asm is None:
            missing += 1
            asm = {}
        ncbi = asm.get("species_name") or asm.get("organism_name", "")
        label = p["label_species"]
        new.append({"assembly_accession": acc, "sample_accession": asm.get("sample_accession", ""),
                    "tax_id": asm.get("tax_id", ""), "species_taxid": asm.get("species_taxid", ""), "ncbi_species": ncbi,
                    "source_id": a.source_id, "label_species": label, "label_complex": complex_of(label, complexes),
                    "label_status": status_of(label), "label_agreement": agreement_of(ncbi, label),
                    **{k: p[k] for k in ("nearest_species", "distance", "differences", "aligned_bp", "margin", "runner_up", "loci")}})
    dup = [k for k, n in Counter(r["assembly_accession"] for r in new).items() if n > 1]
    if dup:
        raise SystemExit(f"{len(dup)} assemblies appear more than once in the input, e.g. {dup[:3]}")

    os.makedirs(a.labels_dir, exist_ok=True)
    lab_path = os.path.join(a.labels_dir, "assembly_labels.csv.gz")
    old = read_table(lab_path) if os.path.exists(lab_path) else []
    cur = {(r["assembly_accession"], r["source_id"]): {k: r.get(k, "") for k in CURATED} for r in old}
    rows = [r for r in old if r["source_id"] != a.source_id]
    rows += [dict(r, **cur.get((r["assembly_accession"], r["source_id"]), {k: "" for k in CURATED})) for r in new]
    rows.sort(key=lambda r: (r["assembly_accession"], r["source_id"]))
    write_gz(lab_path, rows, LABEL_COLS)

    # the run view is rebuilt from every source's rows, so it always matches assembly_labels.csv.gz
    runs = runs_by_sample(a.runs, {r["sample_accession"] for r in rows if r["sample_accession"]})
    run_rows = [dict({k: r[k] for k in RUN_COLS if k in r}, run_accession=run)
                for r in rows for run in runs.get(r["sample_accession"], [])]
    run_rows.sort(key=lambda r: (r["run_accession"], r["source_id"]))
    write_gz(os.path.join(a.labels_dir, "run_labels.csv.gz"), run_rows, RUN_COLS)

    sp_path = os.path.join(a.labels_dir, "sources.csv")
    srcs = [s for s in (read_table(sp_path) if os.path.exists(sp_path) else []) if s["source_id"] != a.source_id]
    srcs.append({"source_id": a.source_id, "authority": a.authority, "authority_release": a.authority_release,
                 "project": a.project, "tool": a.tool, "method": a.method, "call_rule": rule.format(max_dist=a.max_dist),
                 "reference": a.reference, "date": a.date, "refresh": a.refresh, "notes": a.notes})
    n_asm, n_run = Counter(r["source_id"] for r in rows), Counter(r["source_id"] for r in run_rows)
    for s in srcs:
        s["assemblies"], s["runs"] = str(n_asm.get(s["source_id"], 0)), str(n_run.get(s["source_id"], 0))
    with open(sp_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, SOURCE_COLS, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(srcs, key=lambda s: s["source_id"]))

    st, ag = Counter(r["label_status"] for r in new), Counter(r["label_agreement"] for r in new)
    print(f"{a.source_id}: {len(new)} assemblies ({dict(st)}); agreement {dict(ag)}; "
          f"{missing} not in the assembly table; {n_run.get(a.source_id, 0)} runs share a BioSample with a labelled assembly")
    return 0


if __name__ == "__main__":
    sys.exit(main())
