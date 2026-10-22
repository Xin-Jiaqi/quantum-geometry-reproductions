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
