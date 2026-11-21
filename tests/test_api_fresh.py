"""catalog/api/ must match what scripts/popseq_api.py renders from the committed catalog.

Issue #21. catalog/api/ is pre-rendered and nothing else forces a regeneration after
the run table, the summaries, the QC flags, the labels or the query code change. An
external reader fetching over raw.githubusercontent.com then gets an answer from an
older snapshot than catalog/index.json advertises, with no error.

This test renders the whole directory into a temporary folder and compares it, file
by file, with the committed one. The only field ignored is `generated` in index.json,
which is the date the script ran. No network: the script reads committed files only.

When it fails, regenerate and commit:

    uv run python scripts/popseq_api.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API = os.path.join(HERE, "catalog", "api")
FIX = "uv run python scripts/popseq_api.py, then commit catalog/api/"


def _load(path: str) -> dict:
    with open(path) as fh:
        doc = json.load(fh)
    doc.pop("generated", None)
    return doc


def test_committed_api_matches_a_fresh_render(tmp_path):
    out = tmp_path / "api"
    subprocess.run([sys.executable, os.path.join(HERE, "scripts", "popseq_api.py"),
                    "--out", str(out)], check=True, cwd=HERE, capture_output=True)

    committed = sorted(f for f in os.listdir(API) if f.endswith(".json"))
    fresh = sorted(f for f in os.listdir(out) if f.endswith(".json"))
    assert committed == fresh, (
        f"catalog/api/ holds a different set of files than a fresh render. "
        f"Only committed: {sorted(set(committed) - set(fresh))}. "
        f"Only fresh: {sorted(set(fresh) - set(committed))}. Run {FIX}; delete the "
        f"files a fresh render no longer writes.")

    stale = [f for f in committed if _load(os.path.join(API, f)) != _load(out / f)]
    assert not stale, (f"{len(stale)} of {len(committed)} files in catalog/api/ are stale, "
                       f"first {stale[:5]}. Run {FIX}.")
