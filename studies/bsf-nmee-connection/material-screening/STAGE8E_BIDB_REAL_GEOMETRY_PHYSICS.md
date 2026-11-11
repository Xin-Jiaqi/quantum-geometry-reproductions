# Stage 8e: Source-backed BiDB TMD bilayer geometries and exact symmetry gates

**Original independent BSF × NMEE research, NOT Jia et al. (PRL 2026) reproduction.** No real-material NMEE, full ionic relaxation, ferroelectric barrier or driven-DC response has been calculated.

## Source authentication, limits and numerical results

The user provided a 49,745,920-byte BiDB-format ASE SQLite file (bidb.db). Integrity PRAGMA quick_check = ok; 11,184 systems records. Local SHA256:

afe1fe7711189be5b0e189b0c101d6488a7578f4f9ba4e1bb0164534f4255984

Reference: Pakdel et al., Nat. Commun. 15, 932 (2024), DOI 10.1038/s41467-024-45003-w; official BiDB https://cmr.fysik.dtu.dk/bidb/bidb.html . Database schema and identifier conventions agree with BiDB, though no separate author-published file SHA has been compared. This 48 MB file is NOT uploaded to the repository.

Read original ASE numbers (int32), Cartesian positions and lattice (float64), and metadata in read-only SQLite mode. Require exactly 6 atoms, the correct 2:4 chemical stoichiometry, **and number_of_layers=2**. Identified 5 bilayer records each for WSe2, MoS2 and WS2: 15 structures. Excluded same-stoichiometry database single-layer records.

| Compound | ASE row | Source layer group | Independent whole-motif PG | Database slide stability | Database PBE gap (eV) | Interface vertical X–X gap (Å) | Full M_ij space dimension |
|---|---:|---|---|---|---:|---:|---:|
| **WSe2** | **2685** | p3m1 | **C3v** | Stable | 1.090 | 3.158964 | **1** |
| MoS2 | 1930 | p3m1 | C3v | Stable | 1.145 | 3.038294 | 1 |
| WS2 | 2654 | p3m1 | C3v | Stable | 1.296 | 3.086422 | 1 |

WSe2 source bilayer_uid = WSe2-1cfbe6183886-2-1_0_0_1--0.67_-0.33. Its original cell a=3.319006624 Å, c=27.717105609 Å; *upper W minus lower W* has in-plane fractional slide (1/3, 2/3) in this source basis and the bottom W is registered above the upper interface Se. Parallel 1H-TMD layers identify a **polar 3R/AB-like stacking**. Original database dipz=-0.003 (unit not inferred), slide_stability=Stable. The latter checks local frozen-layer lateral PES curvature; it does NOT establish a ferroelectric switching NEB barrier. PBE gap is original database metadata, not a new calculation or complete sliding-path band-gap proof.

**Critical qualification:** the database stores the **z-scan/frozen-monolayer geometry**, where the monolayers were individually PBE-relaxed and their interlayer spacing subsequently optimized using PBE-D3. BiDB separately compared *fully ionic-relaxed* PBE-D3 bilayers. The latter geometry files are NOT authenticated from the user-provided source. The two intralayer chalcogen thicknesses are equal to recorded precision. This is a valid initial DFT structure, NOT a claimed fully relaxed bilayer.

## Independent symmetry of ALL six source atoms

Enumerate 24 candidate hex-point operations (6 in-plane rotations × optional in-plane mirror × optional layer-normal reversal) and test all six atomic elements/Cartesian z coordinates for a one-to-one match modulo one in-plane lattice translation. Independently recovers the original database layer-group point-group classification:

| WSe2 reference | row | actual PG | source slide stability | static Gamma_xz;x |
|---|---:|---|---|---|
| Parallel high-symmetry **AA** | 2688 | **D3h**, 12 ops, horizontal m_z | **Unstable-Max** | forbidden |
| Parallel polar **3R AB-like** | **2685** | **C3v**, 6 ops | **Stable** | symmetry-allowed |
| Opposite-orientation **2H-like** | 2684 | **D3d**, 12 ops, inversion | Stable | forbidden |

**AA is an unstable sliding maximum, NOT a stable ferroic endpoint.** 2H requires relative layer reorientation and is not a simple in-plane slide of the 3R bilayer. One C3v source row was retrieved; its mirror-reflected domain partner is explicitly *derived*, NOT a second original database row. The mirror partner has reversed metal slide (2/3,1/3) and C3v symmetry, with the anticipated physical tensor relation M_AB = -M_BA in a fixed laboratory frame.

### New, stronger response test than single-component symmetry

Let M_ij = Gamma_{i,z;j}, i,j in {x,y}. Because spin and B are both axial but E is polar, their determinant factors cancel. Under an operation with in-plane R2 and normal flip fz, M transforms as

\[
M\ \mapsto\ f_z\,R_2 M R_2^\top.
\]

Reynolds-project the whole 4D vector space of 2×2 matrices using the ACTUAL atomic point-group operations:
- **C3v 3R AB or BA: dimension exactly ONE**, hence M_xx=M_yy, M_xy=M_yx=0, with unknown sign/magnitude of their common scalar.
- D3h AA and centrosymmetric D3d 2H: dimension ZERO.
- A generic displaced C1 stack would have all FOUR components allowed.

This clarifies Stages 6–7f: **nematic linear vs scalar cubic** competition is expected along **generic off-C3 slide locations**; **at an exactly C3v AB/BA endpoint the nematic MUST vanish**. This is a stringent independent validation condition for future DFT/Wannier responses. It is not a claim of nonzero material NMEE.

## Deliverable, method and status

The separate standalone Stage8e research ZIP (given to user) contains source-backed exported CIF/POSCAR for ALL 15 6-atom BiDB bilayers, one clearly named *derived mirror structure*, raw metadata+SHAs, CSV, source extractor and 8/8 passing source-backed tests. Original 48MB uploaded database is not included or published. The exported CIF files preserve six explicit atoms with P1 **data carrier** group; the independent higher point-group classification is in JSON, not silently asserted as a crystallographic CIF symmetry expansion. ASE/spglib were unavailable offline; the program runs with NumPy+SQLite and its group outcome was checked against independent BiDB original metadata. Pinned group-theory toolkit still needs an external library cross-check on these actual atoms.

Next scientific gate: **controlled actual full-bilayer relaxation and whole-BZ SOC gap**, starting from WSe2 row 2685 and its derived mirror plus AA row 2688 (reference only), then spinor Wannier H, Sxyz, physical basis overlaps/covariant velocities. Ferroelectric bistability/NEB and driven-DC NMEE are separate tasks. Published PRL source reproduction is frozen; all additions belong to the independent research study.