# v54 | Independent multiband Bloch geometry: physics-first acceptance

**Goal change:** the original Xiang et al. PRB/SM is historical evidence and a source of candidate identities, not a ground-truth acceptance standard. New physics capability is validated by independent methods, not agreement with disputed manuscript formulas. Original SM Eq.(20)/(31), Eq.(49)/(52), Fig.2(b) differences are preserved as reproduction history.

## Independently verified API

The orbital_geometry_v54.band_geometry(H, velocities) function takes any **nondegenerate multiband Hermitian H(k)** and the matched analytic velocity operators dH/dkx,dH/dky in a **fixed orthonormal orbital basis**:

- r_nm^i = i(v_i)_nm/(E_m-E_n), n≠m
- q_nm^ij = r_nm^i r_mn^j
- g_nm^ij = Re(q_nm^ij)
- Omega_nm^xy = -2 Im(q_nm^xy)
- mtilde_n^xy = (1/2) sum_m (E_n-E_m) Omega_nm^xy

Here mtilde is the **energy-weighted orbital-moment geometry convention in original SM Eq.(51)**. No charge, hbar, B, relaxation time or SI units are inserted. Per-band and per-pair metric, curvature and moment are returned with a global minimum direct band gap. Degenerate bands fail closed; no non-Abelian degenerate-subspace solver is claimed.

curvature_moment_product retains the **complete local pairwise** Omega_n*mtilde_n double sum and separates m=l from m≠l cross terms. moment_curvature_kernel calculates ONLY -mean_BZ(sum_n f'_n Omega_n mtilde_n) for a **uniform periodic momentum mesh**. This is a particular magnetic-moment×curvature *sector* of original orbital Eq.(48), **not the entire orbital MNHE or PHE**.

## Independent physics acceptance: 13 tests

1. Spectral-velocity pair formulas verified against **an independent eigenprojector finite difference**: iTr(P[Px,Py]), perpendicular-derivative metric, Im⟨∂xu|(H-E)|∂yu⟩ moment. True periodic two-band and three-band Hermitian Hamiltonians tested.
2. **Bandwise, k-dependent U(1) phases** and **constant orbital-space unitary** basis transformations leave all geometric outputs invariant. k-dependent orbital representation change is explicitly outside the API without an orbital embedding/position connection.
3. Trace curvature sum vanishes, pair exchange identities and positive metric diagonals hold; two-band upper and lower geometric orbital moments are equal; energy shifts and uniform rescalings have correct dimensional behavior.
4. Independent triple-loop m≠l expression equals vectorized product decomposition. The analytic periodic **three-band** model at mu=-1.1,T=.45 gives full -0.00079945105356, diagonal -0.00154218314981, interference +0.00074273209625. This is 92.9% of **that one complete subterm**, NOT whole MNHE; it varies with chemical potential.
5. Invalid dimensions, non-Hermitian operators, degeneracies, inconsistent eigenvectors and invalid energy units fail closed.

## Reuse

Inside papers/xiang-2026-quantum-geometric-map:

    python -m unittest tests.test_orbital_geometry_v54 -v
    python scripts/reproduce_orbital_geometry_v54.py --N 36 --mu -1.1 --temperature .45
    python scripts/reproduce_orbital_geometry_v54.py --input-npz my_model.npz --mu .1 --temperature .025

External NPZ must contain hamiltonian with shape (...,N,N) and velocities with shape (2,...,N,N). The latter must be the **operator dH/dk**, not merely band-group velocities. Loading disables pickle. Nonuniform meshes need explicit physical integration weights. Actual Wannier90 velocity and k-dependent orbital embedding must be handled before predicting observables.

## Capability states

| Capability | State | Scope |
| --- | --- | --- |
| Nondegenerate general N-band quantum metric, Berry curvature and energy-type orbital moment | **VERIFIED** | Analytic velocities in fixed k-independent orthonormal orbital basis |
| Full m=l and m≠l decomposition of one Omega×moment geometric sector | **VERIFIED** | Explicit double sum and separate three-band benchmark |
| Using the sector in a chosen semiclassical/Boltzmann current | **CONDITIONAL** | Missing response terms and preparation assumptions |
| Full untruncated gauge-covariant magnetic positional shift and total orbital MNHE/PHE | **UNRESOLVED** | Need field-coupling input, all three-band paths and independent full-response oracle |
| Degenerate non-Abelian subspaces, k-dependent orbital positions, arbitrary Wannier representation and SI transport | **UNRESOLVED** | Separate physical development and validation required |

**Next:** implement a physically defined magnetic-field perturbation/orbital-embedding API, independently derive untruncated positional shift and benchmark a real multiband response. Do not promote geometric VERIFIED status to whole transport. Original-paper scientific status PARTIAL.
