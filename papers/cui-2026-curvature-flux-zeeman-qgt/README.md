# Cui et al. (2026): non-Hermitian Zeeman QGT and local Dirac flux

> **Codex research handoff (reviewed 2026-10-23):** [Start work immediately](CODEX_START_HERE.md) · [Independent full-reproduction and research-idea review](FULL_REVIEW_2026-10-23.md). Distinguishes independently green main (371/371 dual Python), user-reported local Fe demo, and currently unmerged scientifically contested [Round36 PR #107](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/107). **Whole Cui theoretical paper remains PARTIAL.** No material EPW archive is required by the source paper.


**Scope: VERIFIED NONDEGENERATE MODEL-METHOD SUBSET; FULL PAPER PARTIAL.**

Cui, Xiang, Xu and Wang, *Curvature-Flux Representation of Local Dirac-Node Topology from Non-Hermitian Zeeman Quantum Geometry*, arXiv:2604.09725v3 (September 2026). Independently implemented; not affiliated with the authors.

This intentionally small public checkpoint includes the authors' **geometry / local topology** through Eqs. (2)–(15), reproduced from three independent two-band formulations and a four-band physical-spin interband extension. `core.py` computes both chiralities and finite-mass flux; `multiband.py` computes the four Zeeman sectors for any nondegenerate Hermitian Bloch Hamiltonian with same-gauge physical-spin operators. The test suite compares eigenstate and spectral-projector oracles, full complex QGT, Wilson phase, winding, finite-mass Gauss flux, k-dependent basis covariant derivatives and failure at degeneracy.

**Not included / not yet certified in this GitHub checkpoint:** the paper's full finite-frequency conductivity Eq. (16)–(18), reciprocal KME Eq. (19), Table I's disputed magnetic point-group classification, Appendix C spin currents, experimental SI response, real Wannier material data and complete continuum/collision effects. The separate 11-round local research archive retains broader numerical/model validations, with 185 local tests; those results are NOT promoted to remote full-paper acceptance by this publication.

~~~bash
python -m pip install -e .
python -m unittest discover -s tests -v
~~~

For user-supplied Wannier data, `wannier90_hr.dat` and centers alone are insufficient: `spin` must be the three *physical* Pauli matrices projected into **the same gauge** as H, and `velocity` must use a declared position connection when the basis is k-dependent. A `.spn` in the original Bloch gauge is NOT directly interchangeable. The public API refuses missing operators and unresolved degeneracies; it does not infer Pauli spin matrices from the Wannier orbital count.

No publisher PDF, personal archive, restricted VASP files, raw AutoWann matrices or private user notes are distributed. CI PASS validates the named methods only. Future paper response branches require separate protocol-specific verification and versioned numerical baselines.


## v2: finite-frequency causal *interband* response (2026-10-22)

The separately verified [causal_response.py](src/flux_lab/causal_response.py) computes magnetic-field-induced current, electric-field-induced spin magnetization in a declared model normalization, causal Liouville/Kubo and static-subtracted interband reciprocity. It also independently evaluates the **four** gN, gA, ΩN and ΩA contributions to the source Eq.(16) to second order. Source Fourier omega and Liouville drive frequency have opposite signs: **omega_source=-omega_drive**. Full scope, results and open physics are in [CAUSAL_RESPONSE_ACCEPTANCE.md](CAUSAL_RESPONSE_ACCEPTANCE.md).

~~~bash
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/demo_causal.py --kmesh 48
~~~

**Not full Eq.(19) KME or a real-material transport package.** The v2 suite covers 21 geometry and 18 new response tests; actual CI status must be checked separately. Intraband/collision/magnetization contacts, real Wannier spin and position, and SI conversion stay open.


## v3: Appendix C conventional physical-spin current (2026-10-22)

A deliberately bounded [spin_current.py](src/flux_lab/spin_current.py) implementation adds clean intrinsic Eq.(127)/(128), the complete physical-spin velocity derivative and Fermi surface / Fermi sea identity. It requires the explicit Hessian of the Hamiltonian for general periodic models, rather than silently discarding the Eq.(132) spin-current vertex as in the special linear-H limit of Eq.(133).

Details and actual numerical values: [APPENDIX_C_ACCEPTANCE.md](APPENDIX_C_ACCEPTANCE.md), with a synthetic fixed-basis four-band benchmark and frozen 64×64 BZ results. **17 new independent tests, 56 total.** The current operator is conventional and generally NOT conserved under spin-orbit torque. This does not establish an experimental spin Hall conductivity, real-material position/spin operators, or the full Eq.(134) in arbitrary systems.


## v4: spin-vorticity Eq.(27) and magnetic-group screening (2026-10-22)

[spin_vorticity.py](src/flux_lab/spin_vorticity.py) validates the energy-weighted anomalous spin curvature against **half the curl of the band spin-current expectation** by spectral sums, independent eigenstate finite differences and the exact Hessian-aware derivative. The mixed Hessian cancels in the curl but not in the symmetric RTA response. Full-zone integrals of this local curl can still cancel in a periodic band.

[magnetic_symmetry.py](src/flux_lab/magnetic_symmetry.py) supplies 31 2D magnetic point groups and **124 conditional Jahn-character real-coefficient tensor projections** independently checked by Reynolds averaging and nullspace. This is NOT a claim that all 31 finite-frequency antiunitary Kubo responses or original Table I are fully validated. See [ROUND14_VORTICITY_AND_MAGNETIC_SYMMETRY.md](ROUND14_VORTICITY_AND_MAGNETIC_SYMMETRY.md) for conventions, six recorded conditional source-table differences, attribution and exclusions.

~~~bash
python -m unittest discover -s tests -v
python scripts/demo_vorticity_mpg.py
~~~

**80 total tests** (56 preceding + 24 new). Original full-paper reproduction remains **PARTIAL**.


## v5: microscopic antiunitary response selection (2026-10-22)

[antiunitary_response.py](src/flux_lab/antiunitary_response.py) supplies exact H/velocity/true-spin transforms, Zeeman QGT pair identities, and nondegenerate frequency-domain interband Kubo/Liouville for TR (1'), C2zT (2'), and mxT (m'). It includes four-band spin x orbital coupling controls. Crucially the complex-frequency identity connects z to **-z**, not -omega+i eta at finite causal broadening. A clean paired BZ verifies forbidden and permitted odd/even channels. [Scientific details](ROUND15_ANTIUNITARY_MICROSCOPIC.md). **100 tests** (20 new); no full 31-group dynamical proof or real material observable.


## v6: full 31-group local microscopic bridge and progress dashboard

[Group-to-Kubo bridge](src/flux_lab/group_microscopic.py) connects the
generator-based 31 magnetic point groups to a 4-band orbital x physical-spin
Bloch model, explicit projective spinor representations, proper real-space
versus momentum-space antiunitary maps, full complex-frequency local Kubo,
and 62 even/odd tensor rank cross-checks with the existing Jahn projector.

[Round16 scientific acceptance](ROUND16_ALL31_MICROSCOPIC.md) | [Three-track
progress ledger](REPRODUCTION_PROGRESS.md). The 31 checks are local generic-k
synthetic models, NOT a 31-material or complete finite-frequency bulk-KME
demonstration. 18 new unit tests bring the limited paper package to 118.


## v7: non-Abelian clusters and first hexagonal BZ (Round 17)

- nonabelian_response.py: projector-resolved Zeeman-QGT blocks, cluster-pair interband Kubo and independent Liouville for exact/near degeneracies.
- hexagonal_bz.py: triangular reciprocal periodic model, Cartesian velocities, first Wigner-Seitz hexagon, D6-symmetrized quadrature.
- Round 17 acceptance: ROUND17_NONABELIAN_HEX_BZ.md; progress ledger: REPRODUCTION_PROGRESS.md; audit script: scripts/audit_round17.py.

19 new tests bring this isolated package to **137 tests**. Python 3.11 and 3.12 both pass after correcting a mirror-asymmetric quadrature rule. Near-degenerate cluster flattening is explicit; the implementation is not a conserving metallic or SI transport theory. Full-paper status remains **PARTIAL**.


## v8: reciprocal primitive cell cross-check (Round19)

[primitive_bz.py](src/flux_lab/primitive_bz.py) integrates the same periodic model on a reciprocal parallelogram using independent Gauss-Legendre and shifted midpoint schemes, without projection of the result. Compare full complex even/odd tensors to the D6-symmetrized hexagonal Wigner-Seitz calculation. [Scope](ROUND19_PRIMITIVE_BZ_CROSSCHECK.md) and [audit](scripts/audit_round19.py). Eight new tests bring the isolated suite to **145/145 PASS on Python 3.11 and 3.12** (Actions 39742758842), with reproducible audit results in results/frozen/round19_primitive_bz_audit.json. Numerical BZ convergence remains **UNRESOLVED**.


## v9: Round20 k-resolved gap and hotspot diagnostics

[bz_hotspots.py](src/flux_lab/bz_hotspots.py) independently evaluates the even/odd complex interband pair response and diagnoses band gaps, occupation contrasts, weighted k hotspots, and the cancellation of local contributions in the BZ integral. **14 new tests; 159 total PASS** in Python 3.11/3.12. [Scientific scope](ROUND20_HOTSPOTS_AND_CANCELLATION.md), [frozen audit](results/frozen/round20_khotspot_cancellation.json), [audit script](scripts/audit_round20.py). These model-method tests do **not** show the net BZ response is numerically converged.


## Round 21: Fermi occupation-sector decomposition

The accepted Round20 bz_hotspots.py stays the authoritative k-space
hotspot/cancellation tool. The separate
[occupation_partition.py](src/flux_lab/occupation_partition.py) resolves
clean directed interband Kubo terms into disjoint **cross-mu** and
**same-side** energy pairs, checks their exact reconstruction of the
complex-frequency response, and tracks tiny thermal occupations even
near small-gap resonances. The integration uses existing hexagonal and
primitive BZ nodes, not a new symmetry-projected quadrature. See
[Round21 scientific acceptance](ROUND21_OCCUPATION_PARTITION.md) and
[reproducible sector audit](scripts/audit_occupation_partition.py).
Eight new tests; actual CI **167/167 PASS** on Python 3.11/3.12 (Actions run 39748027692). Frozen data: results/frozen/round21_fermi_sector_audit.json. BZ convergence remains UNRESOLVED.


## Round22: cancellation-conditioned numerical accuracy and aliasing controls

A new [BZ error-budget method](src/flux_lab/bz_error_budget.py) quantifies the
**conditional** amplification of local weighted-L1 integrand error by strong
response cancellation and compares fixed complex response components across
independent primitive-cell midpoint/Gauss and first-hexagon BZ sampling.
The empirical multi-grid spread is explicitly **not a strict quadrature error
bound**, and Fourier aliasing has a negative control.
[Scientific acceptance](ROUND22_ERROR_BUDGET_AND_ALIASING.md) |
[frozen data](results/frozen/round22_bz_error_budget_audit.json).
**178/178 PASS** on both Python versions, Round21 occupation
partition retained; complete-paper status **PARTIAL**.


## Round 23: equilibrium static Zeeman-current Ward control

[Equilibrium Ward kernel](src/flux_lab/equilibrium_ward.py) resolves the
difference between the earlier clean interband Zeeman current response
and the **complete finite-temperature equilibrium static response**:
the diagonal Fermi-occupation derivative is required. In a fixed Bloch
basis with k-independent physical spin, their sum is
minus kappa times the momentum derivative of the equilibrium spin density.
The continuum BZ integral vanishes by periodicity (no uniform equilibrium
bulk charge current). Independent finite-B and finite-k differences test
the identity. See [Round23 acceptance](ROUND23_EQUILIBRIUM_STATIC_WARD.md)
and [audit](scripts/audit_equilibrium_ward.py).

This does not certify dynamic metallic KME, scattering, gauge contacts,
degenerate bands or source-paper completeness. **188/188 tests passed** on Python 3.11 and 3.12 (Actions run 39756688525). Frozen executed audit: [Round23 JSON](results/frozen/round23_equilibrium_ward_audit.json).


## Round24: exact Kramers-degenerate equilibrium static Ward identity

[degenerate_static_ward.py](src/flux_lab/degenerate_static_ward.py)
extends Round23's static Zeeman–current Ward identity to **exact
degeneracies** without choosing an arbitrary single-band gauge. The
Fermi-matrix divided difference takes the regular limit f'(eps) on a
degenerate cluster, including n != m entries *inside* the cluster.
A separately contracted spectral-projector oracle, finite uniform
Zeeman field, finite momentum gradient, and within-Kramers U(2) basis
controls verify this equilibrium-static relation. Unlike the earlier
nondegenerate formula, this method does not reject a Kramers doublet.

[Scientific scope](ROUND24_DEGENERATE_STATIC_WARD.md) ·
[executable audit](scripts/audit_degenerate_static_ward.py).
15 new tests (including frozen-audit regression) target **203 paper-specific tests** on the isolated Cui
package; require actual Python 3.11/3.12 CI before promoting to VERIFIED.
Not dynamic metallic KME; original full paper **PARTIAL**.


## Round25: k-dependent Bloch gauge and effective-spin static Ward

[covariant_spin_ward.py](src/flux_lab/covariant_spin_ward.py) extends the accepted Round24 Kramers static Fréchet identity to a smooth k-dependent Bloch basis and a deliberately synthetic effective k-dependent spin vertex. It tests the connection commutator and explicit Zeeman-current contact against independent finite-B current and finite-k spin-density differences. [Scientific scope](ROUND25_COVARIANT_SPIN_EMBEDDING.md). 11 new tests, target 214 total including frozen CI regression; first 213 succeeded in both Python versions. Equilibrium static only, not dynamic/material KME.

## Round26: same-Wannier-gauge physical operators

The [in-memory operator contract](src/flux_lab/wannier_operator_contract.py) checks declared same-Wannier-basis H(R), *physical Pauli* spin S(R), and full projected position A(R) including orbital centers; consistent eV/Å units, R/-R hermiticity and reciprocal Fourier derivatives. Its synthetic SOC benchmark uses independent Wilson finite differences and projector Zeeman QGT, plus deliberately wrong spin basis / omitted connection negative controls. See [bounded physical acceptance](ROUND26_SAME_GAUGE_WANNIER_OPERATORS.md) and [frozen synthetic audit](results/frozen/round26_same_gauge_operator_audit.json). The interface does not parse .spn or certify any material's spin operators, SI KME or finite-q transport. **Full original paper PARTIAL.**

## Round27: finite-q longitudinal electromagnetic Ward (synthetic model)

[Finite-q U(1) gauge-Ward kernel](src/flux_lab/finite_q_gauge_ward.py) proves that a time-dependent **pure gauge** built from matched scalar potential and Peierls hopping produces zero induced *physical* density and bond current when the current-operator contact is retained. This is checked using independent spectral/orbital Liouville, exact finite-amplitude gauge, explicit per-site continuity and negative controls. See [scope and residuals](ROUND27_FINITE_Q_GAUGE_WARD.md), [audit](scripts/audit_finite_q_gauge_ward.py) and [frozen data](results/frozen/round27_finite_q_gauge_ward.json). **No conserving bath, genuine metal, real Wannier source or source Eq.(19) SI KME is verified. Full paper PARTIAL.**

## Round28: local-density-conserving collision with U(1) gauge reference

[Linearized collision model](src/flux_lab/conserving_collision.py) adds a finite-torus, finite-q local site-charge conserving relaxation \(-\gamma{\cal Q}[\delta\rho-i[G,\rho_0]]\). Its two intentionally separate gates are **local continuity** and **Peierls-gauge covariance**, with naive RTA and omitted covariant reference as negative controls. An independent dense orbital Liouville solve checks the spectral fixed-point solution, including physical Pauli Zeeman driving. [Scope and derivation](ROUND28_LOCAL_CONSERVING_COLLISION.md) | [audit](scripts/audit_conserving_collision.py). This is an *affine linearized model* with an explicit external-gauge reference, **not** a microscopically thermalizing CP Lindblad bath or a conserving metal KME in SI. **Full paper PARTIAL.**

## Round29: CP thermal balance versus local charge Ward, and a bath gauge vertex

[Two finite Lindblad bath kernels](src/flux_lab/thermal_lindblad_ward.py) contrast local-site charge-preserving CP dephasing (not Gibbs stationary) with deliberately constructed canonical detailed-balance eigenstate jumps (Gibbs stationary but with nonzero local *bath charge source*). The thermal bath's local U(1) gauge covariance requires *rotating all jump operators*, and a directly finite-differenced bath-response contact cancels the pure-gauge state response. See [derivation and hard exclusions](ROUND29_THERMAL_LINDBLAD_WARD.md) and [audit](scripts/audit_thermal_lindblad_ward.py). Synthetic two-site **one-particle** canonical model only, not finite-frequency metallic KME or a derived electron-phonon bath. **Full original paper PARTIAL.**

## Round30: interacting fermionic Fock-sector thermal bath

[Exact CAR Fock operators and many-electron GKSL](src/flux_lab/fermionic_fock_bath.py) now test **four physical-spin fermion modes (16 Fock configurations)**, Hubbard/density interaction, Fermi-Dirac noninteracting control and Pauli blocking, exact 0–4 particle sectors and number-preserving grand-canonical KMS jumps. It distinguishes **total charge** conservation from a genuine local **bath source** and independently differentiates electromagnetic and magnetic bath contacts; omitting those terms fails. [Derivation, negative controls, scope](ROUND30_INTERACTING_FOCK_THERMAL_BATH.md) and [executable audit](scripts/audit_interacting_fock_bath.py). The chosen small-system rates are *not* a microscopic material electron-phonon kernel or Eq.(19) SI KME. **Whole Cui paper PARTIAL.**

## Round31: harmonic-boson KMS spectrum, interacting Fock bath, finite-frequency Ward

[Weak-coupling bosonic model](src/flux_lab/bosonic_davies_vertex.py) explicitly introduces electron density coupling to a thermal harmonic oscillator **continuum** with declared Ohmic \(J(\omega)\). Bose emission/absorption weights are checked by an independent truncated oscillator matrix trace, then used to construct number-preserving interacting fermionic GKSL jumps instead of Round30's arbitrary exponential rates. A zero-frequency collective jump, local bath-charge source, **retarded finite-\(z\)** electromagnetic pure-gauge charge/current Ward with Peierls contact, and a field-differentiated **magnetic** bath vertex are audited. [Full physical protocol and hard limits](ROUND31_BOSONIC_DAVIES_DYNAMIC_WARD.md) · [audit](scripts/audit_bosonic_davies_dynamic_ward.py). This is Born-Markov-secular **toy-model-derived** scattering, not ab-initio electron-phonon matrix elements, Ward-Takahashi metal self-energy/vertex or Cui Eq.(19) SI. **Whole paper PARTIAL.**

## Round32: first-principles scattering source and gauge readiness gate (no real input yet)

[EPW normalized export contract](src/flux_lab/epw_material_contract.py) specifies **complex** g(k,q,phonon,m_final,n_initial) in eV, exact k+q and +/-q phase-paired mesh, phonon masks, H, **physical** Pauli spin and projected position connection in a declared shared Wannier gauge. It verifies NPZ and seven original-source evidence hashes without pretending SHA256 proves physics authenticity. A synthetic-only finite-width Pauli-blocked golden-rule diagnostic checks fully matched **k-dependent H/g gauge** versus a wrong-unmatched-g negative control. A fake fully hashed source bundle deliberately **cannot** earn a real-material claim. [Scope, original EPW sources, hard blockers](ROUND32_EPW_MATERIAL_SCATTERING_GATE.md) | [executable audit](scripts/audit_epw_material_scattering_contract.py). **No real EPW files, authentic spin provenance, conserving metal collision vertex, SI response, or full Cui-paper reconstruction; PAPER PARTIAL.**

## Round33: native Wannier90 H/r read-only parser and Wigner–Seitz weighting gate

[Native \`hr.dat\`/\`r.dat\` dense ASCII parser](src/flux_lab/wannier90_raw_gate.py) handles official one-based matrix rows, **hr-only Wigner–Seitz degeneracy header**, shared r-file weights, two-dimensional \(\pm R\) Hermiticity, Fourier H/A, and already expanded \`write_ndegen_applied\`. Rejects unsupported **unexpanded wsvec** corrections, unmatched source blocks, bad weights and center double-counting. Synthetic source-format fixtures prove that dropping degeneracy changes H(k), while correct decoding reproduces independent complex analytic Fourier oracles. A parsed \`spn\` Bloch→Wannier **physical** spin projection is NOT available: \`physical_spin_available=False\` and \`kme_ready=False\`. [Round33 source scope](ROUND33_NATIVE_W90_HR_R_SPIN_GATE.md) · [audit](scripts/audit_native_wannier90_hr_r.py). No genuine source files received, full Cui **PAPER PARTIAL**.

## Round34: native formatted physical spin, disentangled Wannier gauges and one-shot local handoff

[Native input parser](src/flux_lab/native_spin_chain.py) reads Wannier90 formatted \`.spn.fmt\`, \`_u.mat\`, \`_u_dis.mat\`, \`_v.mat\`, \`.eig\`, constructs physical \(S^W=V^\dagger S^\mathrm{Bloch}V\) with \(V=U^\mathrm{dis}U\), and independently checks Bloch-eigenenergy H against Round33 native Wigner–Seitz-corrected \`hr.dat\`. Independent band-phase matched/unmatched spin and column-major U negative controls prevent common false results. Real spinor sources absent: \`physical_source_provenance_authenticated=false\`, and the original Cui **PAPER PARTIAL**. [Scientific acceptance](ROUND34_NATIVE_SPN_SOURCE_CHAIN.md) · [local one-shot DeepSeek harness instruction](ROUND34_LOCAL_DEEPSEEK_HARNESS.md) · [read-only local source receipt](scripts/collect_round34_local_sources.py). **No automatic raw data upload or large DFT execution.**
