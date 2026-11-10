# BiDB layer-NMEE pilot: proposal, not a submitted campaign

2026-11-10. Sixteen proposed bilayers, eight monolayer parents. All database candidates carry Stable sliding labels, are nonmagnetic and have positive database PBE gaps. Five bilayers have historical zero-gate NMEE results; the inventory does not calculate any new NMEE.

The complete Chinese physical questions, candidate groups, protocol, advancement criteria and handoff state are in the [unified worklog](../../../docs/manuscript/worklog.md), entry “转向BiDB代表性筛选”. No job scripts or remote submissions are included.

- `inventory.json`: full stable identifiers, raw metadata, database hash, toolkit version, layer connections and spatial selection rules.
- `candidates.csv`: compact machine-readable candidate list. Formulae with doubled atoms denote the same chemical material (e.g. In2Se2 = InSe).
- `structures/`: unmodified database-exported bilayer and parent structures. Parent structures are not yet matched isolated-layer references for every bilayer.
- `inventory.py`: local read-only inventory and symmetry verification; requires ASE, spglib, NumPy and the recorded group-theory-operations-toolkit checkout. Run with --database, --toolkit and --out. It performs no DFT or remote access.

Coordinates are in the source cell. For Eip/Bz, `xx` means s_x induced by E_x B_z, `zx` means s_z induced by E_x B_z. Gate derivatives include Fz parity. These are spatially allowed components, not predictions of nonzero magnitude; no electric/magnetic input symmetrization is imposed. Connection representatives are nonunique modulo monolayer symmetry and must not be confused with symmetries of the complete bilayer. Point-orientation coset counts alone do not enumerate all registries or prove degeneracy/barriers.

GeS785 and SnS1084 recover database C2v symmetry at 1e-3 Å but not at tighter tolerances. Raw coordinates are preserved. Verify/reference this discrepancy before claiming strict tensor zeros or gate-induced spin tilting. The hypothetical tilted-spin response requires nonzero allowed coefficients; it has not been calculated.

BiDB source: https://2dhub.org/bidb/bidb.html ; cite Pakdel et al., Nature Communications 15, 932 (2024), doi:10.1038/s41467-024-45003-w. BSF/Enhanced BSF references and literature comparison are recorded in the worklog. No application performance, switching stability or universal criterion is claimed.

Current state: WAIT FOR PHYSICS/PLAN DISCUSSION, DO NOT SUBMIT. Research archive uses [skip ci]; no remote CI pass is asserted.
