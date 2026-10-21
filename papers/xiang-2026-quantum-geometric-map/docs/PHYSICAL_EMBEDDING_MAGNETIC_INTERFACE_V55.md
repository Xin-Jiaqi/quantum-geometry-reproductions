# v55 | Orbital embedding and explicitly specified magnetic Peierls field

**Objective:** independently verified reusable physics, not forced agreement with disputed original SM expressions. The Xiang PRB original-paper science remains PARTIAL.

## Smooth k-dependent orbital basis

For a finite orthonormal orbital basis with physically supplied connection A_i, let the position operator be r_i = i d_ki + A_i. The covariant energy velocity is

D_i H = d_i H - i[A_i,H].

With smooth orbital-frame change W(k), the exact transformation law is

H' = W† H W;
A_i' = W† A_i W + i W† d_i W;
(D_i H)' = W† (D_i H) W.

The physical velocity is D_i H / ħ in wavevector inverse-length units. Functions frame_transform, covariant_velocity and connected_band_geometry implement the frame identity. Independent tests use a genuinely non-diagonal smooth orbital W(k) = exp[i f(k) K] and finite-difference the complete transformed Hamiltonian, comparing both covariant velocities and nondegenerate Berry curvature, quantum metric and energy-geometric orbital moment. **Ignoring A_i' gives wrong geometry** in a negative control.

orbital_center_connection provides point-orbital A_i = diag(r_alpha,i) only. **Wannier orbital centres alone are not generally the full projected position connection**; inter-orbital dipoles and embedding must be supplied. We do not claim arbitrary orbital position connections are physically correct, nor cover nonorthogonal bases, missing external Hilbert-space contributions, or degenerate non-Abelian subspaces.

## Real-space magnetic field from known Peierls paths

For an open finite graph of explicit site positions r_i and Hermitian hoppings t_ij, choose straight paths source site j -> destination site i:

H_ij(B) = t_ij exp[(i q/ħ) integral_(j to i) A(B)·dr]; M_B = dH/dB.

Here q is the signed particle charge; field B_z and units must be consistent. Symmetric gauge A_S = B(-y/2,x/2) and Landau gauge A_L=B(0,x) are independently implemented. Their difference is gradient Lambda=Bxy/2, so G_ii=exp[i q Lambda_i/ħ] and H_L=G H_S G†.

**Crucial magnetic gauge derivative:** because G depends on B,

M_L=G M_S G† + (d_B G) H_S G† + G H_S (d_B G†),

not only G M_S G†. In a physical 4-site square graph, naive transformation has maximum error **0.522** (chosen model units); the full derivative relation residual is **2.8e-17**. Tests also confirm exact gauge-equivalent spectra, Wilson plaquette flux exp[i q B area/ħ] relative to zero field, analytic M_B vs independent finite B differences, field reversal and origin shifts.

**Physical limits:** finite open graph only and explicitly straight Peierls paths, not a fake periodic primitive-cell d_B H(k) under finite uniform B. Spin Zeeman and atomic L·B are separate input operators. Pre-existing zero-field hopping phases also contribute to Wilson loops; isolate added magnetic flux by normalizing by the B=0 loop. No E×B response conductivity is computed.

## Acceptance: 15 independent tests

| Capability | State | Model scope |
| --- | --- | --- |
| Smooth frame-covariant velocity and geometric N-band quantities | VERIFIED | Nondegenerate finite orthonormal orbital representation with complete physically supplied A_i |
| Finite open-graph straight-path Peierls H(B), d_B H | VERIFIED | Known positions, hoppings, signed q/ħ, consistent units |
| Magnetic gauge vertex derivative for B-dependent unitary G(B) | VERIFIED | Includes all d_B G terms |
| Transport response built from these operator primitives | CONDITIONAL | Need E field preparation, current, bath and contact operators |
| Full untruncated orbital MNHE/PHE, arbitrary Wannier embedding, degenerate bands | UNRESOLVED | Independent theory and additional physical input required |

Run from paper subtree:

    python -m unittest tests.test_embedding_magnetic_v55 -v
    python scripts/reproduce_embedding_magnetic_v55.py --mode reference
    python scripts/reproduce_embedding_magnetic_v55.py --mode bloch --bloch-npz input_H_dH_A.npz
    python scripts/reproduce_embedding_magnetic_v55.py --mode peierls --graph-npz input_hoppings_positions.npz --B .32 --q-over-hbar -.87

Input Bloch NPZ keys: hamiltonian, derivatives, orbital_connection. Input finite-graph NPZ: hoppings, positions. Data loading disables pickle. See frozen v55 embedding and Peierls benchmark; unit examples are in the local handoff zip, not published in this PR.

**Next:** build independently complete finite-graph E-field and current/contact operators, verify gauge-equivalent mixed E×B response, then attempt a true orbital magnetotransport implementation. Original paper status PARTIAL.
