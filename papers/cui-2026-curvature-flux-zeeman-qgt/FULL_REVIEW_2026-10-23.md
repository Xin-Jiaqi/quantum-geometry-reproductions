# Cui arXiv:2604.09725v3 | Independent end-to-end research review
**Reviewed:** 2026-10-23. **Source paper:** Cui, Xiang, Xu & Wang, *Curvature-Flux Representation of Local Dirac-Node Topology from Non-Hermitian Zeeman Quantum Geometry*, arXiv:2604.09725v3 (revised 2026-09-17). **Correct target:** a theoretical, 2D Dirac / non-Hermitian Zeeman-QGT / GHE-KME / 2D magnetic-point-group paper. The printed Table I contains *symmetry classes*, **not** a table of computed materials. Never make real-material EPW/DFPT source collection a paper-reproduction acceptance prerequisite.

**Start the next local Codex session at [CODEX_START_HERE.md](CODEX_START_HERE.md)**. That file is the executable handoff, this file is the scientific ledger.

## Verified / not verified | a strict evidence hierarchy

**Green, independently checked on GitHub main:** SHA **ac9b0cb9d5dce19baf55e6fe95d00428bab5ee53**, merged [CI run 39916655754](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39916655754), Python 3.11 **371/371 PASS** and 3.12 **371/371 PASS**, all Round17–34 audits. This validates the **bounded test protocols**, not the complete original paper.

**Local-only report, not independently accessible from GitHub:** Round35 DeepSeek harness says official Wannier90 v3.1.0 Fe spinor postprocessor example (28 Bloch bands, 18 Wannier, 8 k points) passed `Udis U`/spin/Hr consistency; approximate HR/eig matrix error **3.1×10^-6 eV**, local **371 tests**. Local receipt and original files are git-ignored and uninspected by this reviewer. Classification: **USER_REPORTED_LOCAL_FE_POSTPROCESS**, not authenticated from GitHub, and not a Cui source-paper reproduction.

**Open, **NOT** accepted:** [PR #107](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/107) at SHA **c97dee0b0a39747b17451b53f70652dfbdf2dc7e**, [CI run 39942431072](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39942431072). Python 3.12 **383/383 PASS**; Python 3.11 **FAILED** on import (unparsable f-string at `src/flux_lab/zeeman_v3.py:741`), **372 tests and one import error**; no 383/383 dual-pass claim is permissible. Scientific review and exact fixes: [PR #107 review comment](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/107#issuecomment-6364242247). DO NOT MERGE until corrected and retested.

## 1. What was established, from the beginning

| Phase | Accepted bounded result (source-aligned versus optional) | Never infer |
|---|---|---|
| Early geometry, Round8–13 | Conventional QGT, two-band Zeeman-QGT split, analytic Dirac mass/flux, Berry phase and winding; eigenstate versus projector/analytic comparisons, clean **interband** causal Kubo versus orbital Liouville and reciprocal channel tests; Appendix C conventional spin current including the necessary Hessian vertex. See **README.md**, **CAUSAL_RESPONSE_ACCEPTANCE.md**, **APPENDIX_C_ACCEPTANCE.md**. | Full metallic Eq.(19) SI KME; all-band contacts/collisions; real materials. |
| R14 | Nondegenerate **physical-spin** vorticity derivative and curl Eq.(27) with spectral, central-difference, Hessian-aware independent routes; 31 groups and **124 conditional static** tensor spaces. See **ROUND14_VORTICITY_AND_MAGNETIC_SYMMETRY.md**. | Static Jahn projectors as a complete **finite-frequency antiunitary** theorem. |
| R15–16 | Correct spinful antiunitary transformations; complete complex-frequency `z -> -z` relation, not an incorrect `-omega+i eta` retarded substitution; 31/31 group generators at **generic k** in synthetic four-band local Kubo and selected full-BZ examples. See **ROUND15_ANTIUNITARY_MICROSCOPIC.md**, **ROUND16_ALL31_MICROSCOPIC.md**. **CRITICAL earlier resolved rule:** for spatial R and antiunitary theta, **Q=(-1)^theta R**, with physical axial spin action **L=(-1)^theta det(R) diag(R,1)**. | A full material-based or all-group BZ classification; every allowed tensor being realized by every model. |
| R17–22 | Covariant degenerate projector/cluster synthetic response, six-triangle mirror-symmetrized hexagon BZ quadrature, independent primitive-cell comparisons, hot-spot/cancellation budgets, explicit occupied/unoccupied Kubo partition; frozen acceptance and negative controls. | Real material or genuinely converged bulk KME. |
| R23–27 | Equilibrium static no-bulk-charge-current Ward and occupation derivative, exact-degenerate and k-dependent spin-gauge cases, matched H/S/position operator interface, Peierls finite-(q,ω) charge Ward and current-contact negative control. | Static interband-only band sum equals a complete equilibrium observable; all corrections captured in a real metal. |
| R28–31 | Local charge-preserving affine collision; thermal Lindblad / interacting Fock and bosonic bath KMS models; conservation, detailed balance and finite-frequency gauge/magnetic vertices *in explicitly synthetic finite models*. | Microscopic real-metal phonon/impurity self-energy with conserving vertex; paper Eq.(19) derived from these baths. |
| R32–34 | Normalized complex EPW g source **eligibility** checks, official native text Wannier90 hr/r degeneracy parsing; formatted Bloch physical Pauli spn + Udis U projection and eig/HR oracle; careful no-fake-provenance negative controls. | Authentic first-principles material source simply because schemas/hashes/fixtures pass; this is **optional material-extension infrastructure**. |
| R35 local | **User-reported**, official Fe spinor postprocessing files, valid-looking HR eig match and U/spn chain in independent local harness. | QE-SCF to spn full origin chain, material Cui data, CI-verified tests from private local artifacts. |
| R36 open #107 | **Candidate** v3-specific independent `dhat` geometry, massless Hodge/flux/winding, GHE/KME four-frequency-channel tests, v3 claimed sign mismatch and Table I microscopic oracles. | Treat 3.12-only success as dual-CI, or claim manuscript Table I error from a model not honoring its antiunitary generator. |

**High-confidence conceptual result:** the central physical route is
\[
 T^{Z,ab}_{nm}(\mathbf k)=r^a_{nm}(\mathbf k)\,\sigma^b_{mn}(\mathbf k),\qquad
 \{\text{normal metric},\text{anomalous metric},\text{normal curvature},\text{anomalous curvature}\},
\]
and in a 2D Dirac model the **anomalous curvature radial flux** supplies a Hodge-dual view of the **local winding** (\(\pi_1\)); the four sectors feed distinct symmetry/frequency components of the gyrotropic conductivity and reciprocal KME. See v3 arXiv abstract and the repository's early accepted oracles. **This is a model/topology claim, not a real-material DFT claim.**

## 2. Essential scientific corrections before accepting Round36

### (A) Antiunitary mirror: the proposed m′ paper contradiction is currently INVALID

The unmerged PR's `_probe_models()[('Mx', 1)]` uses
\[
 d=(k_x,\ k_x k_y,\ m+k_x^2+k_y^2).
\]
But \(m'_x=m_x\mathcal T\), with \(R_{m_x}=\operatorname{diag}(-1,+1)\), acts on momentum as **\(Q=-R=\operatorname{diag}(+1,-1)\)** and on physical spin as **\(\Lambda=\operatorname{diag}(-1,+1,+1)\)**. Required invariance:
\[
 d(k_x,-k_y)=(-d_x(k_x,k_y),\ d_y(k_x,k_y),\ d_z(k_x,k_y)).
\]
The candidate fails visibly (at \((0.37,-0.52)\), the d-vector mismatch has first two components \(+0.74,+0.3848\), norm ~0.834). It therefore **cannot** establish the reported Table-I `m′ g^N` source error. A correct positive example is \(d=(k_y,k_x,m+k_x^2+k_y^2)\), which obeys that invariance; its relevant integrated \(g^N\) was numerically ~0 in a separate local check, but **one example is not a theorem**. Recompute before stating any source contradiction.

The baseline **R15–16** already implemented the correct \(Q=(-1)^\theta R\), and its antiunitary response complex \(-z\) branch. Reuse this instead of reinventing an inconsistent antiunitary representation.

### (B) C3 symmetry does NOT force all two-band in-plane d to be radial

The PR claims universal `d_x+i d_y=z f(|z|²)` and universal zero integrated normal curvature. **False.** A valid spinful \(C_3\)-invariant counterexample:
\[
d=(-k_y,k_x,m),\qquad d_+=i(k_x+i k_y).
\]
It is **tangential** and invariant under simultaneous C3 k/spin rotation. By source's normal-sector analytic Eq.(4),
\[
 \Omega^N_z = \frac{n(r^2+2m^2)}{2(r^2+m^2)^{3/2}},
\]
nonzero at finite \(r\) and with nonzero radial integral (for \(m\ne0\)). General trigonal term \(\bar z^2\) is also C3 covariant, so the single radial test in PR #107 cannot establish a universal theorem. This invalidates the reported blanket no-realization claim for \(\Omega^N\); other channel assertions require separate verification.

### (C) Manuscript O(ω) sign difference remains a FLAG, not a proven correction

PR #107 records the main-text Eqs. labeled `exp` / `exp1` O(ω) blocks opposite to the derivation at Appendix B `supp:eq:sigma`; the O(ω²) terms agree per its local v3 source transcription. **The v3 LaTeX input is git-ignored**, so the exact equation text was not independently verified from those private bytes by this reviewer. Before labeling a source error, cross-check conjugation, Fourier sign `omega_source=-omega_drive`, order of pair indices, retarded \(\pm i\eta\), static contact and Onsager definitions using independent Kubo and Liouville routes.

### (D) CI / governance must be fixed

Current PR #107 has a Python 3.11 syntax error. Replace the f-string key by `name + ("'" if theta else "")`. Its **Round36 frozen JSON is not executed as a standalone workflow audit step**. A CI test count is necessary but not sufficient: build corrected explicit **antiunitary invariance, tangential C3, independent v3 sign** negative controls and freeze only their measured outputs. Preserve **main** baseline 371/371 and untouched earlier frozen data. Do not promote the paper beyond **PARTIAL** without closing actual v3 equations/table checks.

## 3. Three worthwhile research ideas, prioritized

**P0 — A clean finite-complex-frequency MPG theorem for the response, not sectors by fiat.** Derive/verify how \(\chi_{JB}(Qk,z)\) and \(\alpha_{ME}(Qk,z)\) transform under all unitary/antiunitary generators, including the advanced \(-z\) partner; then project the **proper transverse symmetric/antisymmetric observable**, not arbitrary rank-2 tensor entries. Compare each Table-I row to explicit model responses and independent Liouville. This may resolve several outstanding symmetry claims with actual publishable theoretical content.

**P1 — Stability and scope of the curvature–flux local invariant.** With correct Dirac cones, vary velocities/chirality/mass/warping, verify the Hodge relation and contour/IR dependence, and separate topological contour charge from nonuniversal pointwise curvature. Check lattice regularization and node annihilation; do not infer an invariant from a single continuum cutoff.

**P2 — One sign-convention ledger from microscopic Kubo to reciprocal KME.** Explicitly lay out \(H_B\), \(H_E\), Fourier signs, \(r_{nm}\), energies, contact/static subtraction, and Onsager normalization. Derive all \(O(\omega)\) and \(O(\omega^2)\) channels twice (spectral and Liouville) before concluding the paper appendix has a sign error. Distinguish intrinsic clean **interband** response from equilibrium bulk static or dissipative metallic KME.

## 4. Interpretation discipline

- **GREEN** = remotely verified bounded method with frozen numbers on exact commit, not whole-paper acceptance.
- **AMBER** = unmerged or local-only results, convention-sensitive interpretation, representative-not-universal model.
- **RED/BLOCK** = fails physics symmetry covariance or supported CI; must be fixed before merging.
- Preserve source-derived statements, our calculations, inferred possibilities and independent criticisms as **four different categories**. Do not write a speculative erratum to authors.
- This branch changes **documentation only**; no source, frozen results, old PRs, solvers or registry status is overwritten.

**Immediate operations:** [CODEX_START_HERE.md](CODEX_START_HERE.md) · [Round36 PR #107](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/107) · [review notes](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/107#issuecomment-6364242247) · [Round36 original issue #105](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/issues/105).
