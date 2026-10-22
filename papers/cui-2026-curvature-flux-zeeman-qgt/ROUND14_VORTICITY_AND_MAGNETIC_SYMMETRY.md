# Round 14 | Spin vorticity and 2D magnetic point groups

**Reference:** Cui et al., *Curvature-Flux Representation of Local Dirac-Node Topology from Non-Hermitian Zeeman Quantum Geometry*, arXiv:2604.09725v3. Main-text Eqs. (20)–(28) and Table I. Original PDF is NOT bundled.

**Status:** Spin-vorticity Eq.(27) **VERIFIED for nondegenerate fixed-basis synthetic models**; 2D MPG candidate tensor spaces **CONDITIONAL (assumed Jahn transformation rule)**; full-paper **PARTIAL**. These are deliberately distinct scientific acceptance gates.

## Physical identity

Let \`O_a^c = {dH/dk_a, sigma_c}/4\`, where \`sigma_c\` is a genuine spin Pauli operator and not an orbital pseudospin. In a k-independent orthonormal basis with fixed physical-spin matrix:

\`T_nm^{ab,c} = O_{a,nm}^c v_{b,mn}/(E_n-E_m)^2\`.

The anomalous curvature vector is \`Omega_nm,z^A = -Re[T_xy-T_yx]\`. Energy weighting gives

\`barOmega_n,z^A = sum_(m!=n) (E_n-E_m) Omega_nm,z^A = [d_kx <O_y>_n - d_ky <O_x>_n]/2\`.

This is exactly source Eq.(27) in the specified nondegenerate setting. It has been verified through THREE independent routes:
1. all band-pair spectral Zeeman-spin matrix elements;
2. direct eigenstate diagonalization and central finite differences of the spin-current expectation;
3. the full operator derivative with explicit Hessian spin-vertex from the prior Appendix C module.

The antisymmetric mixed Hessian cancels because \`d_a d_b H = d_b d_a H\` and \`sigma_c\` is momentum-independent. The **symmetric Hessian term does not cancel**; consequently the reduced geometry-only RTA expression obtained from Eq.(28)/(133) is NOT an exact complete response for a generic periodic tight-binding model.

Local 4-band periodic example \`k=(0.31,-0.22)\`, lowest-band \`barOmega_z^A=0.006122411735850528\`, spectral vs finite difference max error across bands \`4.79e-12\`, explicit Hessian curl residue exactly zero, maximal symmetric spin-vertex \`0.1766604381\`.

A smooth periodic band has \`int_BZ curl<O> =0\`, as in a torus Stokes identity. All-band unweighted separate-band full BZ midpoint averages on this model are ~\`10^-18\`, **even though** local vorticity is nonzero. This is NOT proof of a nonzero macroscopic spin current; occupation/symmetry/scattering matter.

## 2D magnetic point groups: necessary constraints, not a dynamic Kubo theorem

The file \`data/toolkit_2d_generators.json\` contains **31 generator-only** 2D magnetic point groups derived from the user's BSD-3-Clause \`group-theory-operations-toolkit\` pinned at commit \`551a0daf8d4ba5a5f24b2811694bfe438616de48\`. Proper 2D group embedding uses a VERTICAL mirror \`m_x\` and \`C_n\` about z. The original toolkit's standalone mirror convention differs; it is explicitly re-embedded. All 31 groups independently pass orthogonality, closure and group order tests.

Candidate transformation of a current (in-plane polar) versus magnetic field (in-plane axial) static **real coefficient** is \`X -> eta^theta det(R3) R2 X R2^T\`, with \`eta=-1\` for source's \`ae\` sectors and \`eta=+1\` for \`e\`. Both the Reynolds group projector and a separate stacked \`(D(g)-I)X=0\` nullspace determine the same invariant dimension for **31 × 4 = 124** candidate subspaces.

Caution: at finite frequency, antiunitary symmetry acts on the complex response and frequency analytic structure, and needs a separate microscopic Kubo/Liouville oracle. **These 124 projector matches do not certify the paper's full Table I as a dynamical-response theorem.**

For reproducible provenance the printed Table I is transcribed without silent changes; six CONDITIONAL differences with the specified static rule occur for \`2'\`, \`m'\`, \`m'm2'\` in source-labeled \`gN\` and \`OmegaA\`. They are documented for follow-up and not asserted as errors in the physical paper. Crucially, a symmetry-allowed diagonal symmetric metric does NOT imply an allowed symmetric **transverse xy** response (e.g. 3, 4, 6 and related black-white mirror groups).

Microscopic positive control: the planar massless Dirac \`H=v_x k_x sigma_x + v_y k_y sigma_y\` obeys \`C2z*T\` represented by \`sigma_x K\` at fixed k; its in-plane mixed QGT real part and out-of-plane anomalous curvature vanish, consistent with its local tensor constraints. This positive control does **not** constitute 31 direct Hamiltonian calculations.

## Scientific boundaries

- Verified: Eq.(27) and its spin/orbital operator assumptions, independent numerical derivatives, reduced-vs-complete Hessian distinction.
- Computationally accepted with assumptions: 31 generator snapshots, 124 real static candidate Neumann subspaces, group closure, direct nullspace, selected Bloch spinor check.
- Not established: complete complex finite-frequency magnetic-group Table I, full Eq.(28)/(134) for generic lattice models, non-Abelian Kramers multiplets, SOC spin conservation or torque dipole, scattering vertices, actual Wannier same-gauge physical-spin matrices, SI spin currents and real material response.
- No private AutoWann files, source PDF or licensed DFT material are included.

## Reproduce

~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/demo_vorticity_mpg.py
~~~

**80 independent tests** = 56 existing + 24 new (spin vorticity and MPG audit). Results pinned in \`results/frozen/round14_vorticity_mpg_reference.json\`.
