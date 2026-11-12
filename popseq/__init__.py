"""PopSeq: query the fungal population-sequencing catalog.

The same four operations the CLI exposes, importable:

    import popseq
    popseq.species(min_isolates=100)
    popseq.describe(taxid=1047171)
    popseq.runs(taxid=746128, min_depth=30)
    popseq.panel(taxid=746128, n=100, min_depth=30)

Everything reads committed catalog files. No network, no credentials, no
database. A call on the largest species in the catalog takes about two seconds.

The function names, their keyword arguments and the keys of what they return are
the contract, the same contract the CLI publishes. See
specifications/007-query-interface.md.
"""
from popseq.query import (  # noqa: F401
    SCHEMA_VERSION,
    CatalogError,
    EmptyResult,
    SpeciesNotFound,
    describe,
    index,
    panel,
    preload,
    preload_clear,
    resolve_species,
    runs,
    species,
    to_sample_sheet,
)

__all__ = [
    "SCHEMA_VERSION", "CatalogError", "SpeciesNotFound", "EmptyResult",
    "species", "runs", "panel", "describe", "resolve_species", "index",
    "preload", "preload_clear", "to_sample_sheet",
]
__version__ = "0.2.0"
