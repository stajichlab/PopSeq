# B1: Record whether a panel's reads are raw or trimmed

**Track:** B   **Depends on:** nothing   **Parallel-safe:** yes

**Problem.** `popseq_inventory.py` resolves `Afumigatus_Global` to `_val_1.fq.gz` files,
which are `trim_galore` output, not raw reads. That is correct behaviour, since those are
the files present, but spec 003 has no field for it.

A panel that silently mixes raw and trimmed reads will produce results that cannot be
compared across samples, and nothing in the current schema would reveal it.

**Do.** Add `read_state` to `metadata.csv` with values `raw`, `trimmed`, `unknown`.
Infer from filename conventions (`_val_`, `_trimmed`, `.trim.`) and default to `unknown`
rather than guessing `raw`.

**Acceptance.** Every panel emitted by the converter carries the field. Afumigatus_Global
reports `trimmed`.
