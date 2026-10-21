# Task BD-MNPS3-L2-P1-STRUCTURE-R0

**Owner:** local DeepSeek harness (human-reviewed). **Platform:** local read-only preflight on WSL or macOS workspace. **No VASP/Japan job, no SSH required in this task.**  
**Scientific owner:** GPT B+D research notes + user approval. **Status:** `NOT_SUBMITTED`. No evidence that the local agent has received this task.  
**Full physical rationale, sources and ambiguity gates:** [MnPS3 R0 report](../docs/BD_MNPS3_REAL_MATERIAL_PREFLIGHT_20261021.md).

## Execution instruction (human can pass this file URL to their local harness)

Prepare an auditable **structure and magnetic symmetry preflight only** for the following *published* physical pair:

- Material: bilayer **MnPS₃**
- Source: Sun et al., *Nature Communications* **17**, 3101 (2026), DOI `10.1038/s41467-026-69635-2`
- First selected states: `P_1^+(L_2)` and `P_1^-(L_2)`, path 1.
- **Author connection (Eq.13):** `[E || M_z M_y]` in source's spin/space operation convention. **Do not substitute `m_z T` or L1 Eq.11.**

Perform only these steps:

1. Locate a licensable original monolayer or bilayer geometry with precise CIF/POSCAR provenance; seek the published relaxed pair and SI. If only a source monolayer is found, identify the exact source orientation operation and structural shift. Do not guess site coordinates from a schematic.
2. Validate chemical composition; in-plane lattice, vacuum and PBC; layer-site decomposition; **all four Mn sublattice spin directions** in each bilayer; exact L2 convention; mapping of site indices and spin directions from (P_1^+) to (P_1^-).
3. If original relaxed daughter structures are absent, generate **only labeled, unrelaxed** trial geometries using `extension-to-BSF` (note that this tool clones identical monolayers by default and will need an explicit top-layer orientation transform). Record exact transformation, structure hashes and any unsupported geometry convention. Never call such structures author-relaxed or DFT-stable.
4. Check nonmagnetic geometry symmetry and independently reconstruct the **SOC magnetic symmetry** for the actual moment vector. Preserve tolerances, sign conventions, and failure cases; no assumed magnetic group from spglib alone.
5. Write a **private** `structure-preflight.json` and readable `handoff.md` with citation/source URLs, frozen local Git commit SHAs, both input hashes, exact site/moment mapping, whether the partner relation passed, genuine missing data, and status `READY` or `BLOCKED`.
6. Produce only a **sanitized** handoff for later GPT physics review; licensed POTCAR and other private material stay local.

**DO NOT:** invoke VASP; submit Slurm; run Japanese HPC; claim a numerical energy, spin torque, polarization or magnetic ground state; report `RUNNING` on the strength of an old history; edit frozen paper reproduction or Wannier ranking; upload unpublished source structures/complete Wannier matrices to the public B+D repository.

**Pass criteria:** exact published-state geometry provenance and `L_2` moment map, species/layer labeling, the right Eq.13 connection, valid source/output SHA256, reproducible magnetic-symmetry checks, and no blocked structural ambiguity. If a published source file is inaccessible or conventions are ambiguous, return `BLOCKED` and its precise reason; do not invent missing input.

**Later phases are explicitly separate user decisions:** cheap paired real-material calculations on WSL, then accepted SOC Wannier spin/layer/position vertices, then frequency-resolved local spin and interlayer/exchange/SOC torque. The correct spin-response tensor component is determined only **after** a verified SOC magnetic point group, not chosen to match the synthetic toy.
