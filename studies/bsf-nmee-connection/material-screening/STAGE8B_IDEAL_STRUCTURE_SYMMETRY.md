# Stage 8b | WSe2/MoS2 ideal AA, AB, BA structures and full-motif symmetry

**Status: idealized six-atom structural geometry gate passed; no relaxed bilayer CIF inspected, no first-principles DFT, and no material NMEE has been computed.** This is independently generated follow-up research, **not** a result or alteration of the Jia et al. PRL reproduction.

## Starting structural motif and provenance

Build two aligned 1H-MX2 monolayers, each with metal at in-plane fractional (0,0), chalcogen pair at (1/3,2/3) and vertical heights ±h. In-plane primitive vectors are a1=(a,0), a2=(-a/2,√3a/2). Separate metal planes by D=2h+g and shift only the upper layer by τ:

- AA: τ=(0,0).
- AB: τ=(1/3,2/3), top metal projection overlaps bottom chalcogen projection.
- BA: τ=(-1/3,-2/3), top chalcogen projection overlaps bottom metal projection.
- Generic test: τ=(0.13,0.20).

These names are **defined by these precise orbital registries**, not blindly inferred from the 3R stacking label. All six atomic positions (two MX2 formula units) are tested; species must map to species and the z positions must match.

| Input | WSe2 | MoS2 | Provenance |
|---|---:|---:|---|
| a (Å) | 3.28 | 3.1881 | WSe2 reported lattice; 1H-MoS2 DFT source |
| half chalcogen separation h (Å) | 1.67 | 1.560~ | WSe2 **illustrative**, MoS2 from 17.4639×0.0894 |
| interface X-X vertical gap g (Å) | 3.165 | 3.0 | WSe2 computed 3R interlayer Se gap; MoS2 **illustrative** |
| 3D periodic c (Å) | 30 | 30 | deliberately chosen computational vacuum |

References: [WSe2 monolayer lattice](https://www.sciencedirect.com/science/article/pii/S0169433221018626); [3R WSe2 interface spacing](https://doi.org/10.1039/D1NR06119D); [1H MoS2 structure](https://pubs.rsc.org/en/content/articlepdf/2021/ra/d0ra10443d); [R bilayer WSe2 AB/BA domain evidence](https://doi.org/10.1038/s41467-026-68854-x).

**These generated geometries are NOT extracted relaxed bilayer structures or database CIFs.** Changing the illustrative interlayer gap over 2.7–4.8 Å does not alter the computed point-operation counts, as long as layer planes remain separate.

## Explicit symmetry audit

Enumerate 6 in-plane rotations × 2 in-plane reflection parity × 2 z flips = 24 D6h reference operations. For each candidate operation permit a **single global in-plane shift** and require a bijection of all six atomic sites within the in-plane periodic lattice, with correct element identity and z coordinate. Then Reynolds-project a basis tensor Γ[spin=x,Zeeman=z,E=x]; spin and Zeeman are axial and field E polar, so their two determinant factors cancel and Γ transforms as an ordinary rank-three polar spatial tensor.

| Ideal structure | Number of surviving ops | Point group | C3 | mz | Inversion | Γ_xz;x |
|---|---:|---|---|---|---|---|
| WSe2 / MoS2 AA | 12 | **D3h** | yes | yes | no | **forbidden** |
| WSe2 / MoS2 AB | 6 | **C3v** | yes | no | no | allowed |
| WSe2 / MoS2 BA | 6 | **C3v** | yes | no | no | allowed |
| WSe2 / MoS2 generic sliding | 1 | **C1** | no | no | no | allowed |

The global horizontal mirror with layer interchange maps the **entire AB motif onto BA** (allowing a global in-plane origin shift) in both material constructions. Crucially, **mz is NOT a little-group symmetry of an individual AB or BA domain**. It constrains the relationship between opposite domains, not the value inside one fixed domain.

These results support a *structural symmetry* test of Stage6's C3+layer-interchanging-mirror starting-point selection rule: trace can start cubic in in-plane sliding while symmetric traceless nematic is allowed linear. The actual response magnitude, its sign beyond symmetry partner relation, gap throughout the slide and energetics remain completely untested.

## What has / has not been independently validated

- **Yes:** self-contained pure NumPy full atomic-motif operation enumeration, Reynolds tensor projection, AB↔BA mirror mapping; 3 local tests passed (6-atom stoichiometry, all four structural point groups, tolerance to illustrative gap parameters).
- **No:** no relaxed or externally downloaded bilayer CIF/POSCAR has passed the structure gate; no external \`group-theory-operations-toolkit\` cross-check yet; no DFT, phonons, total energies, polarization, gap, spinor Wannier or response run.
- **Next:** obtain real relaxed 3R-WSe2 AB/BA and AA comparator coordinates, then run the existing toolkit on the actual positions; repeat for MoS2. Carefully check species registry and actual point group for each relaxed structure. Only after an insulating full-BZ and complete-operator preflight should DFT static NMEE be planned.

Local executable deliverable: \`BSF_NMEE_STAGE8B_IDEAL_STRUCTURE_20261020.zip\` with six clearly labeled \`*_ideal_unrelaxed.cif\` scaffolds, independent \`stacking_symmetry.py\`, tests, machine-readable operation table and SHA256 checksums. These are **geometry templates**, never to be cited as actual relaxed crystallographic observations.

## Boundary with paper reproduction

The Jia et al. PRL frozen source and tests remain completely untouched. This is a symmetry gate for **potential future gapped equilibrium STATIC NMEE**, not the PRL driven-DC response or a demonstrated new effect.