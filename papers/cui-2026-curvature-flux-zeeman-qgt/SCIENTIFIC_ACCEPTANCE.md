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


## v9 · Round20: k-space cancellation and occupied-transition diagnostics

VERIFIED for the bounded clean cluster-interband synthetic model: independent analytic even/odd band-pair formula, agreement with the previous causal Kubo at single k and throughout both BZ domains, weighted absolute-local integrals, selected-component cancellation, sampled gaps and diagnostic occupancy/resonance distances. The nonphysical response-truncation threshold error and symmetry-component floating tie were caught and repaired. In group 6 the specified yy even component cancels by ~98.6% at hexagonal order5; group 6m'm' xy cancels by ~21.1%. These are **illustrative finite-grid** results, not a universal symmetry theorem nor a converged material susceptibility. 159/159 tests PASS in Python 3.11/3.12; full paper PARTIAL. [Acceptance and limits](ROUND20_HOTSPOTS_AND_CANCELLATION.md).


## Round22: empirical grid disagreement is NOT certified integration accuracy

The mathematical conditional inequality \`|δI|/|I|<=eps/r\`, given
\`∫|δF|<=eps ∫|F|\`, is verified, as are four-shift midpoint/primitive
Gauss/hexagonal Gauss comparisons and an exact Fourier aliasing
negative control for the synthetic clean cluster-interband model.
The measured finite-grid discrepancies (~98.87% of hex5 for group6 even yy,
~52.83% for group6m'm' even xy) are empirical disagreements, NOT
strict errors against the unknown continuum value. Both Python
3.11 and 3.12 passed **178/178** tests (Actions 39750100446).
Numerical BZ convergence remains **UNVERIFIED** and the original paper
**PARTIAL**. [Complete bounded protocol](ROUND22_ERROR_BUDGET_AND_ALIASING.md).


## Round25 bounded k-dependent Bloch-spin Ward

Static finite-T gauge-covariant Zeeman-current response, including \(-\kappa Tr[f(H) D_k S]\) explicit spin-current contact, is proposed for a synthetic periodic k-dependent basis/vertex and checked by independent finite-B and finite-k differences. Exact Kramers degeneracy is handled by existing Round24 Fréchet method. No finite-frequency metallic KME, true Wannier operator embedding or source-paper completion. [Protocol](ROUND25_COVARIANT_SPIN_EMBEDDING.md). Initial 213/213 Python 3.11/3.12 CI PASS (Actions 39877965441); adding one frozen audit regression for final 214-test gate.
