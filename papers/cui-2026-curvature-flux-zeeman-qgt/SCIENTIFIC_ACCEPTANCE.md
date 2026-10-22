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
