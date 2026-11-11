# BD36 | BiDB PtSe2 structure-only input kit

## Provenance
Direct extraction from the user's read-only BiDB SQLite copy (11184 records, database quick_check OK). Cite Pakdel et al., *Nature Communications* 15, 932 (2024), DOI [10.1038/s41467-024-45003-w](https://doi.org/10.1038/s41467-024-45003-w). Raw SQLite database is NOT committed. See `manifest.json` and original POSCARs in `structures/`.

| BiDB ID | Geometry | Layer group | Inversion | gap_pbe (eV) | Status |
|---|---|---|---|---:|---|
| 2252 | 1T-PtSe2 monolayer, 3 atoms | 72 p-3m1 | yes | 1.165 | monolayer control |
| 2253 | parallel, unshifted homobilayer, 6 atoms | 72 p-3m1 | yes | 0.219 | **primary baseline**; `slide_stability=Stable` |
| 2255 | parallel, nonzero registry homobilayer | 72 p-3m1 | yes | 0.614 | secondary high-symmetry baseline |
| 2254 | distinct *relative-orientation* bilayer | 69 p3m1 | no | 0.384 | separate orientation family; not part of 2253 pure-sliding torus |

The two `POSCAR_PtSe2_BiDB2253_generic_...` files are **derived, rigidly shifted** by top-layer fractional displacement ±(0.08,0.19) relative to database row 2253, keeping the original layer separation. They are NOT additional database records and are NOT relaxed DFT structures. Both satisfy a site-permutation inversion test to numerical tolerance ≈1e−8 Å.

## Physics and calculation requirements

- A family of identical, inversion-symmetric parallel PtSe2 monolayers remains inversion symmetric at every rigid in-plane registry. This forbids the **layer-even** current-induced spin density `S_top^z + S_bottom^z`, but can allow the **layer-odd** `S_top^z - S_bottom^z` at general registry. C3-fixed high-symmetry registries forbid a vector `j_{x,y} -> S_-^z`.
- BiDB's `gap_pbe` values are source database properties, **not new Wannier/SOC transport results**. The 48 MB SQLite source includes atoms and scalar metadata, not `wannier90_hr.dat`, `.spn`, wavefunctions, layer-spin matrices, or an Edelstein coefficient.
- The bilayers have a gap at charge neutrality according to these data. Intraband DC Edelstein work needs clearly specified electron/hole doping and nonzero charge conductivity.
- Genuine PT spinor degeneracies make per-band diagonal spin expectation unreliable. Compute complete-degenerate-subspace projected traces or a validated gauge-invariant formulation. A bare `hr.dat` is not enough for trustworthy layer-resolved spin matrices.
- The past BD34 illustrative model encoded its registry texture by hand; this project must derive any nonzero spin response from electronic structure and real layer/atomic projectors instead.

## Safeguards
Repository `main`, pre-existing research branches and remote PBS production jobs remain untouched. User's prior research-computing permission is READ_ONLY: **no new VASP/PBS/Slurm/WSL/HPC submissions without a new, explicit bounded authorization**. Build and review the input plan first; see `DEEPSEEK_TASK.md`.
