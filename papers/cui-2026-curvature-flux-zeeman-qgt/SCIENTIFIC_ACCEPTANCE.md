# Scientific acceptance · public geometry checkpoint (2026-10-22)

**Scientific status: MODEL/METHOD VALIDATED (bounded); PAPER PARTIAL.**

## Source and scope

- Source: Cui, Xiang, Xu, Wang, *Curvature-Flux Representation of Local Dirac-Node Topology from Non-Hermitian Zeeman Quantum Geometry*, arXiv:2604.09725v3. Main text and Appendix A–C were reviewed in the offline research archive; the publisher PDF is **not** redistributed.
- This public checkpoint transplants the independently verified two-band geometry kernel and 13 matching tests from the local Round 11 research archive (the archive itself remains outside this GitHub PR).
- Adds a separately tested nondegenerate N-band Zeeman QGT / position-connection primitive with eight tests. Inputs are explicitly supplied matched physical operators, **not reconstructed from Wannier90 files**.
- Reproduces Eq. (2)–(15) model subset, **not** a full transport calculation. The local broader 185-test archive does not become an accepted remote paper calculation by citation alone.

## Verifiable gates

| Capability | Independent checks | Status |
|---|---|---|
| Two-band full complex mixed Zeeman QGT and four sectors | eigenvectors, spectral projectors, Pauli analytic identity; signs/band parities | PASS, nondegenerate |
| Local anomalous curvature flux and winding | direct QGT loop integral, Hodge winding, mass regularization, independent Gauss surface integral | PASS, selected Dirac model |
| Berry-phase counterpart | independent Wilson spinor overlaps | PASS, selected loops |
| Four-band physical spin | explicit spin operators, spectral vs projector, basis covariance and wrong-gauge negative control | PASS, model only |
| Non-diagonal k-space orbital embedding | covariant Hamiltonian derivative under a smooth orbital gauge | PASS, synthetic control |
| Degenerate bands and missing same-gauge operators | explicit exceptions | PASS, fail closed |
| Real-material Zeeman response | no same-gauge spin and full position operators | BLOCKED |
| Finite-frequency conductivity, KME, Table I, Appendix C, SI response | not imported in this PR | PARTIAL / NOT ACCEPTED |

## Model conventions

- Energies and momentum are natural units, no lattice length or SI conductivity is asserted.
- The spin input is a triplet of **dimensionless physical Pauli** matrices. Do not substitute pseudospin, hbar*Pauli/2 or an untransformed Bloch-gauge SPN file.
- Derivative `velocity` is the **covariant derivative of H with respect to k**, not group velocity (divide by hbar separately if required).
- The nondegenerate resolvent is undefined at band touchings and does not cover Kramers/PT doublets; no hidden broadening is applied.
- For finite mass and finite contour radius R, the anomalous flux is -pi R^2/(R^2+m^2) in the isotropic convention, not an integer-valued invariant.

## Reproduction

~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
~~~

`21` tests pass locally. GitHub Actions runs the same isolated tests on Python 3.11 and 3.12; CI success must be inspected separately and cannot stand for full paper acceptance. Other paper directories, root indexes and historical frozen results are not modified.


## v2 limited response capability

A standalone physical-spin local **interband-only** causal-response method was added in [CAUSAL_RESPONSE_ACCEPTANCE.md](CAUSAL_RESPONSE_ACCEPTANCE.md): full complex band Kubo versus orbital-basis Liouville; exact same-regulator interband KME reciprocity; O(omega_source²) Eq.(16) Zeeman four-sector expansion; frozen synthetic BZ and local benchmarks. This does **not** verify complete Eq.(19), experimental Onsager, material magnetoelectric transport, metallic intraband Drude or SI units. 18 new tests, 39 overall.


## v3 · Appendix C conventional spin-current method

Verified for the stated nondegenerate fixed-basis spin × orbital **synthetic model only**: Eq.(127)/(128) intrinsic conventional spin Kubo (independent H+E r finite-field state oracle), derivative identity including the necessary band-diagonal Hessian spin-vertex, and Fermi-surface vs Fermi-sea RTA equivalence. 17 Appendix C tests; 56 total. The geometry-only Eq.(133) is NOT general on a periodic lattice. Real spin vs orbital pseudospin, SOC nonconservation, unknown same-gauge Wannier spin, degeneracies, collisions and SI response remain explicit limits. See [APPENDIX_C_ACCEPTANCE.md](APPENDIX_C_ACCEPTANCE.md).


## v4 · main-text spin-vorticity and conditional magnetic symmetry

Eq.(27), for a smooth nondegenerate fixed-basis Hamiltonian with *physical spin* operator independent of k, is VERIFIED through spectral pair summation, direct finite-k eigenstate curl and exact Hessian-aware derivative. The explicit Hessian cancels in the antisymmetric curl, not in the full symmetric lattice RTA. Local vorticity does not establish macroscopic transport.

The 31 2D magnetic-point-group generators and 124 static real Jahn-character tensor candidate subspaces are CHECKED via independent group closure, Reynolds projector and direct nullspace. Their use for a complex finite-frequency antiunitary conductivity remains **CONDITIONAL**; original Table I is not declared fully resolved. See [ROUND14_VORTICITY_AND_MAGNETIC_SYMMETRY.md](ROUND14_VORTICITY_AND_MAGNETIC_SYMMETRY.md). 24 new tests; 80 total on this paper package.


## v5 · three microscopic antiunitary examples, not all magnetic groups

The 1', C2zT, and m_xT true-spin toy Hamiltonians obey exact
antiunitary k/velocity/spin, band-resolved Zeeman QGT, full complex-frequency
GHE/KME interband Kubo and independent N² Liouville constraints; the same
result survives separate physical spin × orbital coupling. Exact reversal
takes z to -z (advanced at eta>0). Opposite-real-frequency **retarded**
substitution demonstrably fails. A 16² BZ exercise includes symmetry-forbidden
zeros and permitted responses, not converged material numbers.
20 additional tests / 100 total. The other magnetic groups, SI observable
KME, dissipative collisions and actual Wannier spin operators remain
unverified. [Scope](ROUND15_ANTIUNITARY_MICROSCOPIC.md).


## v6: all-31 local generator-to-microscopic bridge

31/31 2D MPG projective real-spin generator checks plus generic-point
4-band Kubo χ_JB, α_ME symmetry covariance, 62/62 even/odd in-plane
candidate-rank comparisons, seven independent orbital-basis Liouville
cross-checks and five selected square BZ acceptance models. Crucial:
the catalog real-space R is distinct from **Q=(-1)^anti R** for Bloch
momentum; the pure TR gray group must invert k. Clean interband analytic
branch under antiunitary is full -z, NOT fixed positive eta at -omega.

18 new tests / 118 total. Full Table I, source Eq.(19) measurable KME,
non-Abelian degenerate bands, hexagonal crystal BZ, physical scattering,
real Wannier spin/position and SI response remain PARTIAL or blocked.
[Scope](ROUND16_ALL31_MICROSCOPIC.md).


## v7 · Round 17: non-Abelian clusters and first hexagonal BZ

The new cluster module is **VERIFIED for the declared clean interband cluster-flat protocol**: exact Kramers doublets, global basis covariance, arbitrary within-cluster gauge rotations, projector QGT trace/eigenvalue invariance, independent orbital Liouville, and TRIM even/odd frequency checks. Near-degenerate clusters are flattened to their mean energy; do not generalize this to intraband metallic response.

The new hexagonal module is **VERIFIED for synthetic periodic models and symmetry checks**: triangular reciprocal periodicity, Cartesian velocity conversion, actual first Wigner-Seitz hexagon (area \(2\sqrt3\pi^2\)), fully D6-symmetric quadrature, local complex-frequency Kubo covariance and integrated even/odd tensor projection for selected groups. BZ order-3 output is not convergence evidence.

**137/137 tests pass** in Python 3.11 and 3.12. Real Wannier spin/position, collision/contact terms, non-Abelian metallic KME, converged observable, source Table I and SI response remain **UNRESOLVED/PARTIAL**. See ROUND17_NONABELIAN_HEX_BZ.md.


## Round20: slow-convergence mechanism diagnostics (bounded)

Pair-resolved synthetic clean-cluster response is separated exactly into cross-mu and same-side energy-transition tensors and checked pointwise against existing Kubo for both complex frequency branches. Overlapping near-pole/small-gap subsets track denominator and Fermi factors but are not additive channels. Weighted magnitude, signed net and cancellation metrics are computed *before* any magnetic tensor projection. In the group-6 hexagonal order-12 model, cancellation ratio in allowed even xx is ~6921, whereas group 6m'm' allowed even xy is ~1.44. Strong near-pole conditions may carry Fermi differences ~1e-6. These are explanatory model diagnostics, not finite-order physical-convergence claims. See ROUND20_HOTSPOTS_OCCUPATIONS.md.


Round20 CI-confirmed model subterm: 155/155 tests pass on Python 3.11 and 3.12 (Actions run 39746294291), and full audit JSON is frozen. Sixfold order8 even xx k-integrand absolute-to-net cancellation ratio ~8265; anti-mirror 6m'm' allowed xy ratio ~1.49. Pointwise near-pole occupation differences ~1e-6; fixed-k thermal scan isolates same-side occupation changes. No BZ convergence or real-material KME claim.
