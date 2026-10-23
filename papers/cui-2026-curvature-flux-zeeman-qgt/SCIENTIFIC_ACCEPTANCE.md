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

## Round26: Wannier-operator input gate (synthetic, restricted)

The new [same-gauge operator contract](ROUND26_SAME_GAUGE_WANNIER_OPERATORS.md) checks explicit physical Pauli spin, H(R), position connection A(R), orbital-center consistency, Fourier derivatives, basis provenance declarations and independent Wilson/unitary QGT controls. It is a **synthetic method verification**, NOT a genuine Wannier .spn converter or authenticated physical-spin matrix for a real material. A falsified basis ID can pass metadata checks but yield a demonstrably wrong QGT: this is an explicit **physical negative control**, not a guarantee of automatic gauge discovery. No new original-paper figure, Eq.(19) SI material KME or finite-q conservation is accepted. **PAPER PARTIAL**. The prior Round21/Round25 accepted modules and audit remain unchanged.

## Round27: finite-q longitudinal U(1) charge-Ward prerequisite

**Model/method VALIDATED only after the remote 3.11/3.12 isolated gate. Original paper remains PARTIAL.** The [Round27 physical report](ROUND27_FINITE_Q_GAUGE_WARD.md) tests finite-(q,omega) charge-gauge invariance in a finite 2D periodic SOC lattice: matching longitudinal Peierls and scalar source, explicit bond-current contact, local continuity, spectral vs orbital-basis Liouville, finite-amplitude gauge and negative omissions. The physical *pure-gauge* current/density vanish but the omitted-contact current is nonzero. No finite-q Zeeman KME, scattering conservation, real Wannier material, SI normalization or source Eq.(19) becomes accepted. Round21 and Round26 remain intact.

## Round28: site-charge conserving, gauge-covariant phenomenological relaxation

See [Round28 method](ROUND28_LOCAL_CONSERVING_COLLISION.md). A local-projector collision \({\cal Q}=1-\sum_iP_i(\cdot)P_i\) with a **gauge-covariant affine reference** is designed to preserve finite-(q,omega) local charge continuity and electromagnetic pure-gauge invariance at gamma>0 on a synthetic 3×3 spinful SOC torus. The separate naive-global-RTA and omitted-reference negative controls explicitly distinguish charge conservation from gauge covariance. Physical Pauli-spin Zeeman response and two independent linear solves are included. **This is not a microscopic conserving bath, a thermodynamic material KME, or original Eq.(19) in SI. PAPER PARTIAL.** Acceptance of numerical outputs depends on actual dual Python CI logs; previous 243 tests remain intact.
 
### Round28 actual first-gate numerical evidence

Python 3.11 and 3.12 both passed 256 isolated tests and the Round28 audit. A new frozen number-by-number test and mandatory audit reference must also pass the final 257-test CI. The sitewise collision source is zero, the local scalar-response Ward residual is 7.71e-15, and an independent matrix Liouville solve agrees at 2.59e-15. The two **separate** negative controls give 0.01171157 local-charge violation with global RTA and 0.02747677 fake pure-gauge current with a site-charge-conserving but nongauge-covariant collision. The original paper remains PARTIAL; CP-bath, microscopic self-energy/vertex and material KME are explicitly unverified.

## Round29 | Finite GKSL thermal balance versus local continuity

[Round29](ROUND29_THERMAL_LINDBLAD_WARD.md) distinguishes **local site-charge conserving** CP dephasing from **canonical one-particle pairwise detailed-balanced** CP energy transitions. This is NOT a universal no-go theorem: the local-dephasing toy does not preserve finite-beta Gibbs, and the synthetic Gibbs-preserving jump bath has an explicit local charge source, while globally preserving trace. The correct sitewise Ward must include that source. Bath gauge covariance requires **joint** unitary transformation of the jump operators, with a finite-difference-verified bath-contact/vertex; frozen jumps are an explicit negative control. CI acceptance requires actual GitHub 3.11/3.12 tests and audit. No microscopic material bath, fermionic metal or Eq.(19) SI, and whole Cui paper remains **PARTIAL**.

### Round29 bounded acceptance (remote observed, not entire-paper)

The canonical one-particle thermal GKSL and local dephasing contrast is independently audited with the KMS rate master equation, matrix Liouvillian, RK4 finite-time dynamics, and the explicit local source versus current balance. Frozen synthetic numbers are `results/frozen/round29_thermal_lindblad_ward_audit.json`; final status requires dual-interpreter CI on the frozen regression. The locally dephasing bath has exactly zero site-charge bath source but does not preserve finite-β Gibbs, while the thermalizing jump bath preserves Gibbs but exchanges site population and requires an **explicit site bath source**. The bath gauge vertex is independently finite-differenced and is not optional. No fermionic many-electron KME, real material, or Ward-Takahashi microscopic self-energy+vertex is accepted. **PAPER PARTIAL.**

## Round30: interacting fermionic thermal bath, magnetic collision contact

[Round30 physics and nonclaims](ROUND30_INTERACTING_FOCK_THERMAL_BATH.md) builds explicit four-mode **fermionic CAR** (16 many-body Fock states), spinful SOC with Hubbard interactions, a true grand-canonical Fock Gibbs density and independently validated single-particle Fermi-Dirac zero-interaction limit. A prescribed number-conserving KMS-sector GKSL bath preserves total N and stationary interacting grand Gibbs but **cannot fix mu**; local site density has a nonzero **bath source**, so a closed-electron sitewise Ward is invalid. Missing bath gauge and equilibrium Zeeman derivative contacts produce explicit failures. This is not microscopic electron-phonon scattering, a conserving metal response or source Eq.(19) in SI. Actual Python 3.11/3.12 CI and frozen-data regression are required for bounded acceptance. **Full Cui paper PARTIAL.**

### Round30 actual numerics and exclusion contract

First remote CI passed 291/291 tests in Python 3.11 and 3.12 with independent Round30 audit. Final frozen regression raises expected count to **292**, to be accepted only after exact final-HEAD CI. Interacting fermionic Fock U,V change orbitals by 0.075984; 54 thermal jumps preserve total N but create local site bath source 0.271191 for the selected two-electron state. Detailed balance residual 1.73e-18 and Gibbs stationarity 2.36e-17. Omitting the magnetic bath derivative creates false response 0.005563 despite full derivative residual 7.33e-11. Frozen raw metrics are in `results/frozen/round30_interacting_fock_bath_audit.json`. Fermionic Pauli and grand KMS are verified for this finite model, NOT microscopically derived electron-phonon vertices or original SI KME; **PAPER PARTIAL**.

## Round31 | Explicit bosonic correlation and driven finite-z conserving-bath prerequisites

[Round31 physical protocol](ROUND31_BOSONIC_DAVIES_DYNAMIC_WARD.md) **derives within a declared harmonic-bath model** KMS absorption/emission from thermal oscillator correlation with specified Ohmic continuum \(J(\omega)\), and checks a finite oscillator independently of the interacting Fock jump construction. Strict Born/Markov/secular assumptions and a separate collective zero-frequency channel are explicit. In the four-mode interacting fermion Fock toy, exact grand Gibbs and total N conservation coexist with nonzero local bath site charge source. Co-rotated bath jumps, explicit magnetic bath derivative and **retarded finite-z** electromagnetic gauge state/current contact enforce Ward constraints, with omitted source and contacts as negative controls. **No actual ab-initio boson spectrum, Lindblad Lamb shift, interacting metallic Ward-Takahashi self-energy/vertex, real-material SI KME or source Eq.(19) accepted. Whole Cui paper PARTIAL.** Dual-version CI and frozen data must be checked before acceptance.

### Round31 first real dual-interpreter numerical evidence

The corrected initial synthetic source model passed **308/308** in GitHub-hosted Python 3.11 and 3.12, all prior audits and Round31 audit. The subsequent frozen regression is one additional test; final acceptance requires **309/309** plus all audits on both interpreters and a post-merge main check. Correct model has finite-z exact gauge-state error ~1.07e-16 and contact-complete physical current ~5.89e-17, local bath charge source 0.406784, and finite-z source-inclusive magnetic charge Ward ~1.5e-17. Missing U(1) bath vertex gives spurious response 0.0107735, and missing magnetic bath vertex 0.0150318. Oscillator Fourier correlation normalization was independently corrected after an initial fourfold error, with regression retaining strict original thresholds. Spectral density is specified **Ohmic model**, not a real material; no source SI KME. **PAPER PARTIAL.**

## Round32: real-input eligibility is not real-material reproduction

[Round32 input contract](ROUND32_EPW_MATERIAL_SCATTERING_GATE.md) introduces **a fail-closed normalized export reader, not a universal EPW raw binary parser**. It requires explicit matched H/g/physical spin/position metadata, exactly closed k+q and paired inverse phonons, eV/Å units, phonon-active mask, and original source SHA256 across electronic SCF, DFPT phonons, complex g, Wannier H, *physical* spin, position and exact exporter mapping. **Hashes and gauge IDs cannot prove ab-initio origin or true common spin gauge**: even a fully checksummed deliberately fake bundle is returned only as \`CANDIDATE_BYTES_VALIDATED_PHYSICAL_AUTHENTICITY_UNVERIFIED\`. The synthetic fixture remains synthetic even under a real material name. Golden-rule finite-Gaussian Pauli blocking is not exact on-shell detailed balance, a normalized momentum integral, or Ward-consistent conductivity. Actual source data remain absent; original Cui PAPER PARTIAL. Round32 bounded CI acceptance requires both remote Python versions and frozen negative-control audit.

### Round32 numerical evidence and strict distinction from material verification

First corrected remote CI completed **326/326 in Python 3.11 and 3.12**, plus all historical audits and the Round32 source gate. One added mandatory frozen regression yields expected final count **327**. Frozen numerical **synthetic** max rate is **4.182e12 s^-1** and matched k-dependent H/g gauge error **0.00244 s^-1**, while incorrectly rotating H alone changes rates by **2.716e13 s^-1**; these are **not** real material transport lifetimes or KME. Complete SHA256 evidence for explicitly fake raw files returns PHYSICAL_AUTHENTICITY_UNVERIFIED and never `real_material_verified=true`. The present repository has no complete authentic EPW complex g and same-gauge physical spin/position source bundle. Entire Cui PAPER PARTIAL.

## Round33: original Wannier90 H/r file decoding, not physical-spin certification

[Round33](ROUND33_NATIVE_W90_HR_R_SPIN_GATE.md) adds explicit native dense \`hr.dat\`/\`r.dat\` parsers using the **official Wannier90 raw text formats** and crucial shared Wigner–Seitz degeneracies (r lacks its own weights), hermiticity and independent H/A Bloch Fourier checks. \`write_ndegen_applied\` is handled without double division, but **unexpanded wsvec** and unspecified orbital-shift conventions are rejected. Test files are generated in proper Wannier90 syntax, **not authentic materials**. Parsing H and projected position r does NOT authenticate source physical spin, because raw \`spn\` is Bloch-gauge and needs actual disentanglement/Wannier rotations. No real EPW g, SI KME, material data or source Eq.(19) is accepted. Return \`PARSED_H_R_SPIN_NOT_AUTHENTICATED\`, \`kme_ready=False\`, full Cui PAPER PARTIAL. Dual-version GitHub CI and frozen numeric data required for accepted parser scope.

### Round33 initial numerical acceptance and exact scientific stop

Initial remote Python 3.11 and 3.12 each completed **344/344** and every Round17–33 audit. Mandatory frozen reference adds one final test, expected 345, and requires exact-final-SHA CI. Under the generated native Wannier90 fixture, correct shared hr/r Wigner–Seitz weights reconstruct H with \`1.55e-17 eV\` analytic Fourier residual, while omitting hr degeneracy yields \`0.1059858 eV\` error; unsupported unexpanded wsvec is rejected. The \`spn\` gauge transform is deliberately not implemented, and **no actual first-principles Wannier source has been authenticated**. **Cui paper PARTIAL**.

## Round34: native physical Pauli Bloch-to-Wannier spin, not yet source authentication

[Native Round34 protocol](ROUND34_NATIVE_SPN_SOURCE_CHAIN.md) parses official formatted \`spn.fmt\` Bloch Pauli (Hermitian upper triangle, xyz fastest), native **Fortran column-major** \`_u.mat\`, \`_u_dis.mat\`, optional independent \`_v.mat\`, native eigenenergies and Round33 Wigner–Seitz-corrected \`hr/r\`. It checks \(V=UdisU\), true \(S^W=V^\dagger\sigma V\), independent \(H_\mathrm{eig}^W=H_\mathrm{hr}^W\) at matching k, physical spin Pauli bound and the BLOCH phase covariance. Stale .spn while changing Bloch gauge produces a finite wrong spin; row-major matrix reading is a separate negative control. Local \`collect_round34_local_sources.py\` and [one-shot DeepSeek handoff](ROUND34_LOCAL_DEEPSEEK_HARNESS.md) discover source file roles **read-only** and never claim authentication or upload raw data. Generated source-format fixtures prove math, not authentic DFT. Actual first-principles spin provenance and source Cui KME still absent; full paper **PARTIAL**. Bound acceptance requires exact final-head dual-interpreter CI and frozen data.
