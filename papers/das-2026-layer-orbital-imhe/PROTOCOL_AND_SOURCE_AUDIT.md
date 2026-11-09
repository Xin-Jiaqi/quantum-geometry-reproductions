# IMHE: independent derivation and source audit

## Source and evidence boundary

Source: Das & Agarwal, arXiv:2604.20249v1, 22 April 2026, eight pages. Checked Eq. (1)–(7), Appendix A p7 and Appendix B p8 against the supplied PDF and the arXiv HTML. Supplied PDF SHA256:

`7cacdd348b3f0521d778f703aeaa6f6e9bf02ca6025341ac49976baf8800a5a1`

ArXiv HTML: https://arxiv.org/html/2604.20249v1 . No SM S1–S4 was obtained in this task. R5G's full 10×10 parameter specification is not reconstructed by guessing. Existing Table I, gate-field and gate-parity concerns remain open; they are not used as numerical acceptance criteria. All numerical values below are our independent periodic model.

Background reference for the filled-band control: Xiao, Chang & Niu, *Rev. Mod. Phys.* 82, 1959 (2010), https://arxiv.org/abs/0907.2021 .

## 1. Clearly defined input Hamiltonian

Let A=P and Q=−M₀. At each k use

$$
H(u,b)=H_0+uA+bQ+ubC.
$$

These are specified Hermitian matrices and k derivatives in a fixed orthonormal basis. A gate already present in H₀ is held fixed when differentiating with respect to u. This avoids mixing the background-gate scan with an infinitesimal perturbation.

For isolated band n and d_nm=ε_n−ε_m,

$$
\varepsilon_n^u=A_{nn},\quad\varepsilon_n^b=Q_{nn},\quad
\varepsilon_n^{ub}=C_{nn}+\sum_{m\ne n}\frac{A_{nm}Q_{mn}+Q_{nm}A_{mn}}{d_{nm}}.
$$

In particular, for C=0 the mixed energy is −2 Re Σ P_nm M_mn/d_nm. This follows from Eq. (1) as a fixed-operator Hamiltonian. It does not silently replace A7's different moment-dressing expression.

Use first-order state coefficients a_mn=A_mn/d_nm and b_mn=Q_mn/d_nm, zero for m=n. The off-diagonal mixed coefficient is

$$
c_{mn}=\frac{C_{mn}+(Ab+Qa)_{mn}-a_{mn}\varepsilon_n^b-b_{mn}\varepsilon_n^u}{d_{nm}},\quad m\ne n.
$$

Its diagonal is fixed by normalization:

$$
c_{nn}=-\operatorname{Re}\sum_m a_{mn}^*b_{mn}.
$$

This fixes the overlap with the reference state to be real. The phase-independent observables are tested under arbitrary reference-band phases. Both H·U=U·ε and U†U=1 are checked in every jet slot, including ub.

## 2. Geometry and orbital operator

Compute the fully transformed velocity V_i=U†(∂_iH)U, including derivatives of input operators. Then

$$
r^i_{nm}=\frac{iV^i_{nm}}{\varepsilon_m-\varepsilon_n},\qquad
\Omega_n=-2\operatorname{Im}\sum_{m\ne n}r^x_{nm}r^y_{mn}.
$$

The four-slot product rule retains all paths and all intermediate bands. For a scalar gap d,

$$
(d^{-1})^{ub}=2d^u d^b/d^3-d^{ub}/d^2.
$$

For the reference moment, A5 is evaluated directly for n,p, with intermediate l excluded when l=n or p. Its diagonal continuation agrees with A1. In a two-band model A5's off-diagonal sum is empty; replacing it with a physical-spin Pauli matrix would change the problem. Our nonzero-layer/orbital benchmark therefore uses three orbitals. Analytic ∂k M is checked against finite differences of the full orbital-basis M(k), including the changing eigenbasis.

## 3. Transport: partial versus complete product rule

With f=f(ε), the Berry-current mixed coefficient has five terms:

$$
\partial_u\partial_b(\Omega f)=
\Omega^{ub}f+\Omega\varepsilon^{ub}f'
+(\Omega^u\varepsilon^b+\Omega^b\varepsilon^u)f'
+\Omega\varepsilon^u\varepsilon^b f''.
$$

The first two are Eq. (6)'s structure. All five reproduce Appendix B A16's structure. Return each separately. The antisymmetric coefficient is [[0,−H],[H,0]] with H equal to their integrated sum.

For the specified drift functional −∫v_a v_b f', the full coefficient is

$$
-\int\!\Big[
(v_a^{ub}v_b+v_a v_b^{ub}+v_a^u v_b^b+v_a^b v_b^u)f'
+\{(v_a^u v_b+v_a v_b^u)\varepsilon^b
+(v_a^b v_b+v_a v_b^b)\varepsilon^u
+v_a v_b\varepsilon^{ub}\}f''
+v_a v_b\varepsilon^u\varepsilon^b f'''
\Big].
$$

Eq. (7)'s selected sector omits the single-field cross products; both are returned. An f''' term is required in the full drift expansion. This function is symmetric in a,b and is per τ. Magnetic phase-space measure/velocity modifications and Lorentz-force terms are outside this specifically defined functional.

## 4. Independent numerical checks

The finite-field oracle diagonalizes H at four corners ±h_u, ±h_b, calculates the undiscretized-in-k spectral velocities/curvature at each corner, integrates the entire occupation-weighted response, then divides its alternating sum by 4h_u h_b. It does not use the eigenstate perturbation jets.

A second curvature oracle uses finite-k spectral projectors:

$$
\Omega_n=i\operatorname{Tr}P_n[\partial_xP_n,\partial_yP_n].
$$

The supplied periodic three-band model has H(k,η)*=H(−k,−η), P*=P, and M(k,η)*=−M(−k,−η). Numerical tests therefore check Ω^ub even, ε^ub odd, integrated intrinsic even under η reversal and drift odd. Gate parity is not inferred from time parity and is not hard-coded.

`results/validation.json` records actual field-step, mesh and projector residuals. The recipe emits `mu_scan.csv` when `--outdir` is specified. The frozen numerical values are checked in addition to independent oracles, not used instead of them.

## 5. Executable discrepancies / no forced source agreement

### A4 inverse-gap sign

With Eq. (1) and constant P,

$$
\partial_u[-iV_{nm}/d_{nm}]
=-i[V^u_{nm}d_{nm}-V_{nm}(P_{nn}-P_{mm})]/d_{nm}^2.
$$

A4 prints a plus between the two numerator terms. For H=k_xσ_x+k_yσ_y+(m+u)σ_z at k=0, r_x=i/(2(m+u)) for the lower-to-upper pair. Its derivative is −i/(2m²); A4 gives +i/(2m²). The test uses m=1.2. This is a local inconsistency under the explicitly stated convention, not an empirical fitting correction.

### Second-order states in an A8-like velocity expansion

The exact mixed matrix element contains U_ub†VU + U†VU_ub as well as U_u†VU_b+U_b†VU_u and explicit operator derivatives. The first-order ket in Eq. (2) cannot by itself generate the former terms. A diagonal P counterexample has U_u=0 at the expansion point but nonzero U_ub and mixed interband velocity. A8's first-order-cross construction alone misses it. The normalized full-Hamiltonian derivative is the accepted branch here.

### Filled-band DC control

For a smooth periodic two-dimensional Bloch Hamiltonian, a completely filled isolated occupied subspace has an integer total Chern number. While the gap stays open, ∂u∂b of its full-BZ Berry integral is zero. A local Ω^ub can nevertheless be nonzero. The model's zero-T in-gap response converges to zero; no offset is added to reproduce Fig. 2(d)/4's claimed nonquantized plateau.

A finite valley patch, a non-equilibrium occupation or a different magnetic/boundary protocol requires its own analysis. This test does not establish the origin of the paper's plateau; it shows that its mere presence cannot be an acceptance target for the present full-BZ static effective-band calculation. The earlier intuition “there is a Fermi sea, so its integrated mixed Hall coefficient must be finite” is insufficient.

## 6. Scope exclusions

No R5G curve-fitting, no presumed gate-even extrinsic law, no inferred SOC or spin, no physical ℰ/B SI calibration, no ordinary Lorentz Hall, no all-scattering-order theory, no arbitrary k-dependent Wannier gauge, and no unresolved-degeneracy band formula. None of these is silently marked complete. A supplied C,dC can represent a separately derived mixed contact, but its microscopic validity remains the caller's responsibility.
