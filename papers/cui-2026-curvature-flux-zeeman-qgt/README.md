# Cui et al. (2026): non-Hermitian Zeeman QGT and local Dirac flux

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
