# B4: Decide what the five subprojects with no config.txt are

**Track:** B   **Depends on:** nothing   **Parallel-safe:** yes

`Population_Genomics` has 22 subprojects with a `config.txt` and five without:
`A_fumigatus`, `A_fumigatus_oldstructure`, `Afumigatus_cactus`,
`Afumigatus_deepvariant`, `H_werneckii`.

The converter cannot read them, and this needs a human to say what they are: superseded,
a different pipeline, or data that should be catalogued.

**Acceptance.** Each is either given a `config.txt`, catalogued another way, or recorded
in the catalog as deliberately excluded with a reason.
