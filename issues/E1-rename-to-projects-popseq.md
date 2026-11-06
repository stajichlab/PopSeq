# E1: Rename to projects/PopSeq so other clades are not a fork

**Track:** E   **Depends on:** A1   **Parallel-safe:** no, serialise

Spec 006: the only fungus-specific thing in the schema is the number 4751. The repo
should become `projects/PopSeq` with a `Fungi` catalog inside it.

**This is a rename of paths other projects may already reference.** AGENTS.md section 1
and spec 005 section 6 apply: add the new path, keep the old one working, announce it,
remove the old one in a later commit.

**Acceptance.** A clade parameter drives everything. The fungal catalog still resolves
at its old path for at least one release.
