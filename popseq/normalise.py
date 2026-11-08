"""Normalise the free-text `country` and `host` fields (issue #23).

The archives store both as free text. *Zymoseptoria tritici* alone has 316
distinct country strings. A consumer grouping isolates by country or host must
not re-parse them, so the catalog does it once, here, and every decision is
reviewable in `catalog/vocab/`.

Every normalised field comes with a status, so "no value" has explicit,
separate states:

==============  ==============================================================
status          meaning
==============  ==============================================================
``ok``          resolved to a country (ISO 3166-1 alpha-2) or a host taxid
``absent``      the field is empty: the submitter wrote nothing
``missing``     the submitter wrote an explicit missing term ("missing",
                "not collected", "not applicable", ...). INSDC defines these;
                they are a statement, not a gap
``region``      country only: a real place above country level (an ocean, a
                continent, "West Africa"). True, but it has no ISO code
``historical``  country only: a former state with no single successor
                (Yugoslavia, USSR, Czechoslovakia)
``not_a_host``  host only: an environment, substrate or lab condition written
                in the host field ("soil", "laboratory environment", "culture")
``same_as_organism``  host only: the host field names the sequenced organism
                itself, usually a lab strain. Set when the run table is built
``ambiguous``   host only: the name matches more than one taxon
``conflicting`` a run with several BioSamples whose values disagree
``unresolved``  a value nobody has mapped yet. Add it to `catalog/vocab/`
==============  ==============================================================

The raw value is never changed. The run table keeps `country` exactly as
harvested; the normalised columns sit next to it.
"""
from __future__ import annotations

import csv
import gzip
import os
import re
import unicodedata
from functools import lru_cache
from typing import Optional

from popseq.query import ROOT

VOCAB = os.path.join(ROOT, "catalog", "vocab")
COUNTRY_ALIASES = os.path.join(VOCAB, "country_aliases.csv")
HOST_TAXA = os.path.join(VOCAB, "host_taxa.csv.gz")
# Curated. Read by scripts/popseq_vocab.py when it builds HOST_TAXA, never at
# query time: the generated table already carries every override.
HOST_OVERRIDES = os.path.join(VOCAB, "host_overrides.csv")

COUNTRY_STATUSES = ("ok", "absent", "missing", "region", "historical", "conflicting",
                    "unresolved")
HOST_STATUSES = ("ok", "absent", "missing", "not_a_host", "same_as_organism", "ambiguous",
                 "conflicting", "unresolved")

# INSDC missing-value terms, plus the variants submitters actually write.
# https://www.insdc.org/technical-specifications/missing-value-reporting/
_MISSING = {
    "missing", "not applicable", "not collected", "not provided", "restricted access",
    "unknown", "unknow", "na", "n/a", "n.a.", "none", "not available", "not determined",
    "not known", "unavailable", "null", "not recorded", "no data", "not specified",
    "unspecified", "undetermined", "not reported", "-", "--", "?", ".", "0",
    "not_applicable", "not_collected", "not_provided", "missing data",
}


def key(text: str) -> str:
    """Lower case, accents dropped, whitespace collapsed, edge punctuation stripped."""
    t = unicodedata.normalize("NFKD", text or "")
    t = "".join(c for c in t if not unicodedata.combining(c))
    t = " ".join(t.replace("_", " ").split()).lower()
    return t.strip(" .,;:'\"")


def is_missing(text: str) -> bool:
    k = key(text)
    return k in _MISSING or k.startswith("missing:") or k.startswith("not applicable:")


# --------------------------------------------------------------------------- country

def _load_csv(path: str) -> list[dict]:
    if not os.path.exists(path):
        return []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        return [r for r in csv.DictReader(fh) if not next(iter(r.values()), "").startswith("#")]


@lru_cache(maxsize=None)
def _country_index(alias_path: str = COUNTRY_ALIASES) -> dict[str, tuple[str, str, str]]:
    """key -> (iso2, name, status). Built from pycountry, then the curated aliases."""
    import pycountry

    idx: dict[str, tuple[str, str, str]] = {}
    prefix: dict[str, set] = {}
    for c in pycountry.countries:
        name = getattr(c, "common_name", None) or c.name
        val = (c.alpha_2, name, "ok")
        for n in {c.name, name, getattr(c, "official_name", None), c.alpha_3}:
            if n:
                idx[key(n)] = val
        # "Iran, Islamic Republic of" -> "iran". Only kept when it is unique, so
        # "Korea, Republic of" and "Korea, Democratic ..." do not collapse.
        if "," in c.name:
            prefix.setdefault(key(c.name.split(",")[0]), set()).add(val)
    for k, vals in prefix.items():
        if len(vals) == 1 and k not in idx:
            idx[k] = next(iter(vals))
    for r in _load_csv(alias_path):
        st = r["status"] or "ok"
        if st not in COUNTRY_STATUSES:
            raise ValueError(f"{alias_path}: bad status {st!r} for {r['alias']!r}")
        iso = r["iso2"]
        name = r["name"]
        if iso and not name:
            name = _iso_name(iso)
        idx[key(r["alias"])] = (iso, name, st)
    return idx


def _iso_name(iso2: str) -> str:
    import pycountry
    c = pycountry.countries.get(alpha_2=iso2)
    if c is None:
        raise ValueError(f"{iso2!r} is not an ISO 3166-1 alpha-2 code")
    return getattr(c, "common_name", None) or c.name


def _one_country(text: str, idx) -> tuple[str, str, str]:
    if not (text or "").strip():
        return ("", "", "absent")
    if is_missing(text):
        return ("", "", "missing")
    # INSDC form is "Country: locality". Try the part before the colon first,
    # then bracketed and comma-separated pieces ("Kerala [India]",
    # "UK (Manchester)", "Cerkvenjak, Slovenia").
    head = text.split(":", 1)[0]
    tries = [head, text]
    tries += re.findall(r"[\[(]([^\])]+)[\])]", head)
    tries.append(re.sub(r"[\[(].*?[\])]", "", head))
    parts = [p for p in re.split(r"[,/]", head) if p.strip()]
    if len(parts) > 1:
        tries += [parts[-1], parts[0]]
    for t in tries:
        hit = idx.get(key(t))
        if hit:
            return hit
    return ("", "", "unresolved")


def country(raw: str, alias_path: str = COUNTRY_ALIASES) -> dict:
    """Normalise one raw `country` value.

    A run with several BioSamples arrives with `;`-joined values. Missing parts
    are ignored when the others agree; real disagreement is `conflicting`.
    """
    idx = _country_index(alias_path)
    parts = [p for p in (raw or "").split(";")] or [""]
    res = [_one_country(p, idx) for p in parts]
    real = {r for r in res if r[2] == "ok"}
    if len(real) == 1:
        iso, name, st = next(iter(real))
    elif len(real) > 1:
        iso, name, st = "", "", "conflicting"
    else:
        # No resolved part. Report the most informative state that occurred.
        order = ["region", "historical", "unresolved", "missing", "absent"]
        seen = {r[2] for r in res}
        st = next(s for s in order if s in seen)
        hit = next(r for r in res if r[2] == st)
        iso, name = hit[0], hit[1]
    return {"country_iso": iso, "country_name": name, "country_status": st}


# --------------------------------------------------------------------------- host

@lru_cache(maxsize=None)
def _host_index(taxa_path: str = HOST_TAXA) -> dict[str, tuple[str, str, str]]:
    """key -> (taxid, name, status), from the generated vocabulary."""
    idx = {}
    for r in _load_csv(taxa_path):
        if r["host_status"] not in HOST_STATUSES:
            raise ValueError(f"{taxa_path}: bad status {r['host_status']!r}")
        idx[r["host_key"]] = (r["host_taxid"], r["host_name"], r["host_status"])
    return idx


def host(raw: str, organism_taxids: tuple = (), taxa_path: str = HOST_TAXA) -> dict:
    """Normalise one raw host value against `catalog/vocab/host_taxa.csv.gz`.

    The vocabulary is generated by `scripts/popseq_vocab.py` from the NCBI
    taxonomy, so this function needs no taxonomy dump at query time.

    `organism_taxids` are the run's own taxid and species taxid. A host equal
    to one of them is `same_as_organism`, not a host.
    """
    idx = _host_index(taxa_path)
    own = {str(t) for t in organism_taxids if t}
    res = []
    for p in (raw or "").split(";"):
        if not p.strip():
            res.append(("", "", "absent"))
        elif is_missing(p):
            res.append(("", "", "missing"))
        else:
            hit = idx.get(key(p), ("", "", "unresolved"))
            if hit[2] == "ok" and hit[0] in own:
                hit = (hit[0], hit[1], "same_as_organism")
            res.append(hit)
    real = {r for r in res if r[2] == "ok"}
    if len(real) == 1:
        t, n, st = next(iter(real))
    elif len(real) > 1:
        t, n, st = "", "", "conflicting"
    else:
        order = ["same_as_organism", "not_a_host", "ambiguous", "unresolved", "missing",
                 "absent"]
        seen = {r[2] for r in res}
        st = next(s for s in order if s in seen)
        t, n = "", ""
        if st == "same_as_organism":
            t, n = next((r[0], r[1]) for r in res if r[2] == st)
    return {"host_taxid": t, "host_name": n, "host_status": st}
