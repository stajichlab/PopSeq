"""An MCP server over the catalog, so an agent can query it as a tool.

Start it:

    uv run --extra mcp python -m popseq.mcp_server

Register it with Claude Code (adjust the path if the repo lives elsewhere):

    claude mcp add popseq -- uv --directory /bigdata/stajichlab/shared/projects/Fungi/PopSeq \\
        run --extra mcp python -m popseq.mcp_server

It speaks stdio and exposes the same four operations the CLI does, with the same
arguments and the same return shapes. There is no second implementation: every
tool here calls popseq.query.

Everything it reads is a committed file in this repository. It opens no socket,
needs no credentials, and writes nothing. A failed lookup comes back as an error
message, not as a guess.
"""
from __future__ import annotations

import json
from typing import Optional

# mcp 2.x renamed FastMCP to MCPServer. Pinned to >=2 in the `mcp` extra; this
# import is the 2.x name.
from mcp.server.mcpserver import MCPServer

from popseq import query as q

mcp = MCPServer(
    "popseq",
    instructions=(
        "Query the PopSeq catalog of public fungal population sequencing data. "
        "Call popseq_index first for the manifest and the caveats. Use "
        "popseq_species to find what exists, popseq_describe to judge whether a "
        "species is a usable population, and popseq_panel to select strains to "
        "profile. Counts are floors and a run count is not a strain count."
    ),
)

_SCOPE = (
    "PopSeq catalogs public whole-genome Illumina sequencing runs of fungal "
    "isolates, from an ENA snapshot of 2026-11-06 reconciled against NCBI. "
    "245,822 runs, 217,911 isolates, 6,486 species taxids."
)


def _json(obj) -> str:
    return json.dumps(obj, indent=2, default=str)


@mcp.tool()
def popseq_index() -> str:
    """Read this first. The catalog manifest: what tables exist, their grain,
    key, columns and checksums, plus the caveats that apply to every answer.
    """
    return _json(q.index())


@mcp.tool()
def popseq_species(min_isolates: int = 100, contains: Optional[str] = None,
                   limit: int = 50) -> str:
    """List the species in the catalog, most isolates first.

    Use this to find out what is available before asking for runs or a panel.
    min_isolates=100 is the threshold below which a species is usually not a
    population: 5,435 of the species in the catalog have fewer than 100 isolates.
    """
    rows = q.species(min_isolates=min_isolates, contains=contains, limit=limit)
    return _json({"scope": _SCOPE, "returned": len(rows), "species": rows})


@mcp.tool()
def popseq_describe(taxid: Optional[int] = None,
                    species_name: Optional[str] = None) -> str:
    """Everything the catalog knows about one species: run and isolate counts,
    study spread and the size of its largest study, depth percentiles, how many
    runs carry a mutant/evolved/cross/non-WGS label, and NCBI reconciliation
    status.

    Check `clonality.largest_study` against `clonality.isolates` before treating
    a species as a diverse population. Pass a taxid where you have one; a name
    that matches more than one species returns an error rather than a guess.
    """
    try:
        return _json(q.describe(taxid=taxid, name=species_name))
    except q.CatalogError as exc:
        return f"error: {exc}"


@mcp.tool()
def popseq_panel(taxid: Optional[int] = None, species_name: Optional[str] = None,
                 n: int = 100, min_depth: Optional[float] = None,
                 layout: str = "PAIRED", country: Optional[str] = None,
                 since: Optional[str] = None, include_excluded: bool = False,
                 only_exclusion: Optional[str] = None,
                 max_per_study: Optional[int] = None, seed: int = 0) -> str:
    """Select n strains of one species to profile: one run per isolate (deepest
    wins), then spread across studies by round-robin.

    This is what to call for "give me N strains of X". The result carries the run
    accessions, their ENA FASTQ paths, estimated depth, country and collection
    date, plus how many studies the panel touched.

    min_depth filters on an estimate computed from base_count and a reference
    assembly size, not from the reads; runs with no estimate are dropped by it.
    only_exclusion ("mutant", "evolution", "cross", "not_wgs") selects *only*
    that material, for a project that wants the mutant panels rather than
    avoiding them.

    Study spread is a sampling rule, not a guarantee of independent isolates.
    """
    try:
        return _json(q.panel(taxid=taxid, name=species_name, n=n,
                             max_per_study=max_per_study, seed=seed,
                             min_depth=min_depth, layout=layout, country=country,
                             since=since, include_excluded=include_excluded,
                             only_exclusion=only_exclusion))
    except q.CatalogError as exc:
        return f"error: {exc}"


@mcp.tool()
def popseq_runs(taxid: Optional[int] = None, species_name: Optional[str] = None,
                limit: int = 500, min_depth: Optional[float] = None,
                layout: str = "PAIRED", country: Optional[str] = None,
                since: Optional[str] = None, include_excluded: bool = False,
                only_exclusion: Optional[str] = None) -> str:
    """Every run of one species that passes the filters, up to limit.

    Prefer popseq_panel when the goal is a set of strains to profile; this verb
    returns many runs of the same isolate and is for inspecting what exists.
    """
    try:
        return _json(q.runs(taxid=taxid, name=species_name, limit=limit,
                            min_depth=min_depth, layout=layout, country=country,
                            since=since, include_excluded=include_excluded,
                            only_exclusion=only_exclusion))
    except q.CatalogError as exc:
        return f"error: {exc}"


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
