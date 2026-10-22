# Cui Flux — three-track progress checkpoint

**Date: 2026-10-22. Scientific full-paper status: PARTIAL.**
This is a gate-based ledger, not a percent-complete estimate.
“Passed tests” means the listed model protocol is numerically validated,
not that a real experimental observable or a disputed source claim is resolved.

| Workstream | Accepted gate | Still missing |
| --- | --- | --- |
| **1. Paper reproduction** | Eqs.(2)–(15) two-/multi-band Zeeman QGT, local Dirac flux and Wilson loop; Eq.(16) restricted clean interband causal response, fourfold sector expansion; local GHE/KME reciprocity; Appendix C conventional spin current, Hessian and Eq.(27) spin vorticity; 31-group local antiunitary Kubo covariance | Entire source Fig set and Table I six convention-sensitive entries not certified; full Eq.(19) material-normalized KME, genuine many-band degeneracies, contacts, collisions, SI response |
| **2. Physics capabilities** | Reusable nondegenerate Zeeman QGT; spectral Kubo plus independent Liouville; conventional spin-current derivative; finite-frequency antiunitary selection; static Jahn rank projector and all 31 local synthetic microscopic constraints | Conserving metallic electric drive in general Bloch/Wannier gauge, non-Abelian bands, full spin-torque and material-operator interfaces |
| **3. GitHub integration** | Merged PRs #47/#49/#52/#55/#57 plus physics-registry PRs #48/#50/#53/#56/#59; published source+tests+limited scientific-acceptance docs; Round16 PR #60 merged with 18 additional tests and acceptance documentation | Round16 PR #60 **merged**, double-version CI **118/118 PASS**; remaining offline Round9–11 research protocols deliberately not upstreamed as production |

## Verifiable scope counters

- 31/31 **MPG generator/real-spin projective representation + generic-k four-band local Kubo identities** (new).
- 62/62 **even/odd in-plane Jahn candidate tensor-space rank comparisons** (new).
- 7 representative groups independently tested via 16x16 orbital Liouville in Round16.
- 5 representative square-compatible groups with 8x8 BZ selection checks; C3/C6 hexagonal BZ **not validated**.
- Local tests: 100 previous + 18 new = 118; GitHub Actions run 39722106861: Python 3.11 and 3.12 jobs both success; actual logs each show Ran 118 tests.
- Full-paper acceptance still **PARTIAL**, irrespective of these counters.

## Real-material data gate

Historical AutoWann position-block candidates include HgI2, AlAgS2 and Hg3O3,
but no same-Wannier-gauge **physical spin** operator has been verified.
The Japan system search remains incompletely verified. Do not submit new VASP
runs until the nondegenerate physics interface is reviewed; a purposeful
simple-SOC control run is the next option if verified spin matrices remain absent.

## Next independent milestones

1. Demonstrate BZ convergence of symmetry-allowed response components across quadrature orders on the synthetic hexagonal model.
2. Extend cluster response to general non-Abelian/intraband metallic and conserving-collision regimes; the current cluster-flat interband protocol is not that theory.
3. Integrate conserved finite-q longitudinal electric and gauge-covariant position-operator protocols from the earlier local Round9–11 research.
4. Validate same-gauge Wannier Hamiltonian, physical spin and position matrices for a material; only then attempt macroscopic SI KME claims.

No global completion percentage is stated because test count, physics
generality, and true-paper observables do not share a defensible denominator.


## Round 17 — non-Abelian clusters + hexagonal first BZ

| Workstream | Round 17 update |
|---|---|
| **1. Paper reproduction** | Added exact Kramers-cluster clean interband response and first Wigner-Seitz hexagon integration in periodic triangular coordinates. Still PARTIAL: full source Table I, converged material response, intraband/contacts/collisions and SI KME. |
| **2. Physics capabilities** | Added projector Zeeman-QGT block, cluster-pair Kubo, independent cluster Liouville, D6-symmetric hexagonal quadrature, and Cartesian-velocity conversion. Near-degenerate groups are flattened to their mean energy. |
| **3. GitHub integration** | PR #63 is **merged**. Two modules, 19 tests, a reproducible audit script, scientific acceptance and the CI audit step are on main. Python 3.11 and 3.12 both report **137/137 PASS**, and the audit script ran successfully (Actions run 39730445788). |

### Round 17 measurable gates

- 9 non-Abelian tests: Kramers doublets, global/unitary subspace gauge, projector QGT invariants, Kubo/Liouville equality and TRIM frequency parity.
- 10 hexagonal BZ tests: lattice automorphism, reciprocal periodicity, exact hexagon area, D6 weighted-node symmetry, local Kubo covariance and selected integrated tensor projectors.
- First CI caught a genuine quadrature defect: C6 rotation symmetry alone did not enforce mirror symmetry. The fixed rule averages with a mirror image; forbidden integrated tensor components now pass.
- Package count: **137 tests** (118 prior + 19 new); full paper remains **PARTIAL**.

### Updated next milestones

1. Increase hexagonal BZ quadrature order and demonstrate numerical convergence of allowed response components (still synthetic).
2. Extend cluster response to the general non-Abelian/intraband metallic and conserving-collision regime; the present cluster-flat interband protocol is not that theory.
3. Validate same-gauge Wannier Hamiltonian, physical spin and position operators for a material.
4. Resolve the six conditional source Table I differences and complete material-normalized Eq.(19)/SI checks.


## Round 17 final CI and physics audit (2026-10-22)

- GitHub PR #63 merged; main commit 2083a063bda0f030dd357b553c86fd6265ca97c7.
- Actions run 39730445788: Python 3.11 and 3.12 both passed 137 tests and ran scripts/audit_round17.py.
- Non-Abelian doublets: cluster dimensions [2, 2]; spectral cluster Kubo vs independent orbital Liouville maximum residual \(1.92\times10^{-16}\); maximum trace magnitude of the selected Zeeman-QGT block 0.09184 (synthetic natural units).
- Hexagonal local covariance maximum across H, Cartesian velocity, Kubo and ME tensors: \(4.46\times10^{-16}\) for the four tested groups \(6,6mm,6m'm',3m'\).
- First Wigner-Seitz hexagon area \(2\sqrt3\pi^2=34.1893125466\); order-4 rule uses 192 weighted nodes, weight sum 1.0000000000000002.
- Order-3 integrated even/odd projector residuals: group \(6\), \(2.17\times10^{-18}\) and \(4.34\times10^{-19}\); group \(6m'm'\), \(2.43\times10^{-18}\) and \(5.96\times10^{-20}\).
- These are synthetic model values and symmetry checks, not a converged material prediction. The full paper remains **PARTIAL**.


## Round19: independent primitive-cell BZ cross-check

| Workstream | Current gate |
| --- | --- |
| Paper reproduction | Same clean cluster-interband 4-band model integrated over hexagonal Wigner-Seitz and reciprocal primitive fundamental domains; exact continuum equality, finite-order convergence not assumed. Whole-paper PARTIAL. |
| Physics capabilities | Independent reciprocal primitive-cell Gauss and shifted midpoint quadratures, full complex even/odd tensor comparison without projecting away numerical errors. |
| GitHub | Round19 PR #67 **merged into main**; independent primitive method, 8 tests, full audit and frozen reference data. Python 3.11/3.12 **145/145 PASS** with successful final audit in Actions run 39743169120. |

Round18 demonstrated order20-to-order24 relative differences ~6% for group 6 and ~3.9% for group 6m'm' in even response; projector residuals close to machine precision do NOT establish convergence. Next priority: primitive/hex domain numerical agreement, spectral-gap and temperature sensitivity. No additional VASP work at this gate.


### Round19 independent audit result (2026-10-22)

- Frozen raw comparison: results/frozen/round19_primitive_bz_audit.json.
- Group 6: hex12 even max 4.85529e-7, primitive Gauss20 5.14225e-7,
  primitive midpoint32 5.47477e-7, but max full-tensor
  Gauss20-minus-hex12 is 2.78444e-7 (57% relative to the comparator).
- The unprojected primitive grid has symmetry residuals of 1.25e-7
  (Gauss20 group6) and 1.34e-8 (midpoint32 group6), versus ~1e-18
  in the symmetry-enforced hex quadrature. This diagnoses numerical
  integration, not fundamental symmetry breaking.
- CI only certifies the *independent integration method and fixed tests*.
  The BZ tensor remains numerically UNCONVERGED and full paper PARTIAL.
- Next: k-resolved resonant/near-degenerate/occupation diagnostics and
  a higher-order independent method comparison, before material DFT.


### Round19 GitHub handoff finalized (2026-10-22)

- Source PR #67 merged; main squash commit `34ef72a469e61d7dc8efcccc737cb86d745e0834`.
- Final head before squash `bd82b5a6103b729f7d8242df1d870d0943182efe`; Actions 39743169120 successfully ran **145 tests and Round19 independent primitive-cell audit** in both Python 3.11 and Python 3.12.
- Central Physics Registry explicitly labels the primitive-cell method VERIFIED but the numerical BZ tensor convergence UNRESOLVED. Full paper remains PARTIAL.


## Round20 — microscopic origin of slow convergence

| Workstream | Current scope |
| --- | --- |
| Paper reproduction | Pair-resolved clean interband response, temperature/occupation and k-space cancellation diagnostics of the paper-inspired synthetic model. Original paper still PARTIAL; no converged material KME. |
| Physics capabilities | New bounded occupation partition, resonance/gap subset analysis, k-hotspot and absolute-vs-net cancellation metrics, with independent reconstruction against existing Kubo. |
| GitHub | Round20 branch adds 10 new regression tests (155 target), an audit script, and scientific findings. Count actual PASS only after dual-version CI. |

Key preliminary independent observation: group 6 hexagon-order12 net even xx ~4.85e-7 but weighted absolute kernel ~3.36e-3, cancellation factor ~6921. In this finite sample, cross-mu absolute integrand dominates (~3.36e-3) yet integrated cross-mu net ~3.77e-9, while same-side thermal transitions give ~4.82e-7 net. Group 6m'm' allowed xy has only ~1.44 cancellation ratio. These are not converged material quantities.
