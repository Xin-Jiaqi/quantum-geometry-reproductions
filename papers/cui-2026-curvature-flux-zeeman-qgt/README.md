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
