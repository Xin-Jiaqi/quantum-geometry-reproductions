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


## Round 20: spectral k-hotspot and cancellation diagnosis (2026-10-22)

| Workstream | Accepted method / remaining gate |
| --- | --- |
| **1. Paper reproduction** | Exact algebraic interband even/odd pair decomposition and finite-grid BZ cancellation studied independently. This explains potential sensitivity of tiny net values to discretization; full source Table I, Eq.(19) material normalization and full paper still PARTIAL. |
| **2. Physics capabilities** | Added k-resolved gap, active-occupation detuning, weighted hotspot concentration, complex integral absolute magnitude, cancellation ratio, and fixed-observable cross-domain comparisons. Diagnostic thresholds never truncate physical terms. 14 new tests, 159 total. |
| **3. GitHub and CI** | PR #69 **merged** (main squash commit b60956d18596a1b99bcb7cd885ac6273c2b20444); both Python 3.11 and 3.12 **159/159 PASS** with full Round20 audit (GitHub Actions run 39746316981). Scientific audit and frozen data are on main; central Physics Registry updated via companion documentation PR. |

**Important:** Group 6 yy even cancellation 98.575% (hex order5), compared with group 6m'm' xy 21.12%; neither physical net has been demonstrated fully converged. Group 6 active-transition detuning min ~4.44 while min sampled intercluster gap ~0.011; 'active' uses a threshold only in diagnostic reports. Next test high-order adaptive sampling focused on cancellation error bars and true global gap extrema, before material computation.


## Round21 — thermal band-pair partition supplementary acceptance

| Workstream | Current status |
|---|---|
| Paper reproduction | Existing clean interband Kubo response further decomposed by exact cross-mu versus same-side transitions; whole original paper remains PARTIAL, BZ convergence UNVERIFIED. |
| Physics capability | New occupation_partition.py partitions every directed intercluster term while preserving arbitrarily small occupation differences; fixed-H thermal scans and independent prior Kubo are the validation gates. Prior bz_hotspots.py remains authoritative for spatial hotspots. |
| GitHub | Round21 PR #74: 8 new tests (159 prior → **167 total**); Python 3.11 and 3.12 both **PASS**, plus reproducible audit in Actions run 39748027692. Frozen synthetic JSON in results/frozen/round21_fermi_sector_audit.json. |

The model group-6 hexagon-order8 allowed xx has cancellation factor ~8265, with cross-mu net -2.007e-7 and same-side net 6.025e-7, adding to total 4.018e-7. Thermal factor scans at fixed Hamiltonian give strong same-side response variation, but none of these are converged real-material observables. Next physics gate: analytic interband cancellation sum rule and consistent BZ convergence tests.


### Round21 GitHub final status

- PR #74 merged (squash commit b17f043972425804ceb680fd5a5839b67808e392).
- Final double-version GitHub Actions 39748495850: Python 3.11 and 3.12 each ran **167 tests PASS**, including the reproducible occupation-sector audit.
- Physics Registry lists this as a VERIFIED bounded model method, NOT a converged BZ response or real material KME. The previously merged Round20 bz_hotspots code remains the separate authoritative spatial hotspot tool.


## Round22 final checkpoint — cancellation-conditioned BZ error budget

| Workstream | Verified today | Open gate |
| --- | --- | --- |
| **1. Source-paper reproduction** | Conditional net-error amplification identity for a fixed complex response component and model-level multi-method discrepancy validation | True absolute/relative BZ integration error, original material-specific Eq.(19), full Table I and SI KME; entire Cui paper still **PARTIAL** |
| **2. Reusable physics** | Bounded \`bz_error_budget.py\` method: four shifted midpoint meshes at coarse and fine resolution, primitive Gauss, independent hexagon Gauss, complex drift, aliasing negative control | Rigorous/adaptive integration bounds, real Wannier spin and position operators, conserving collision model |
| **3. GitHub** | **Round22 PR #75 merged** at \`db0e10d9a5921027f0b95c16c9d1ff4a3e58d9bf\`. Prior Round21 occupation-sector PR #74 remains intact. Final Actions run **39750100446**, Python 3.11 and 3.12 both **178/178 PASS**, including all paper-specific audits | Numerical convergence certification deliberately absent |

### Quantified finite-mesh evidence

- Previous Round20 group 6 even yy *finite-grid* retention: 0.01424847. Under the unverified weighted-L1 perturbation premise, obtaining ≤1% relative net error requires ≤0.014248% local weighted-L1 error.
- Group 6m'm' even xy retention: 0.788792, giving a corresponding ≤0.788792% conditional local accuracy target.
- Four shifted 9×9 midpoint grids for group 6 have maximum complex spread 2.37929e-7; cross-method discrepancy relative to finite hex5 is ~98.87%. Group 6m'm' xy corresponding shift spread 4.43700e-7, relative disagreement ~52.83%.
- These numbers are **empirical finite-grid disagreements**, neither rigorous bounds on the unknown exact BZ integral nor probabilistic confidence intervals.
- The Round21 occupation-sector tool and Round22 error-budget tool are complementary: thermal sector decomposition does not in itself establish quadrature convergence.
- Next physical milestone: controlled band/occupation-sector-resolved adaptive quadrature, treating cancellation before asserting small nonzero responses. No new VASP run is required at this point.


## Round23: equilibrium no-current Ward (under dual-version CI review)

| Stream | New bounded gate |
| --- | --- |
| Paper | The original source is still PARTIAL. Independent equilibrium Ward provides a necessary static-limit negative control, not a new claim of source Eq.(19) reproduction. |
| Physics | Explicit intraband Fermi-derivative + interband static Kubo = minus kappa momentum derivative of equilibrium spin density in fixed-orbital k-independent-spin basis. The BZ integral vanishes for a complete periodic equilibrium crystal. |
| GitHub | New equilibrium_ward.py, 10 tests, analytic scope note and finite-B/finite-k audit in an isolated PR; accepted count requires real remote CI results. |

A discrete periodic spin-density gradient telescopes to zero at every grid
resolution, while an independently sampled spectral derivative may have
nonzero finite-grid integration error. Neither calculation is a finite-
frequency metallic transport response. True Wannier embedding and
collision/contact terms remain outstanding.


### Round23 actual equilibrium static Ward acceptance

- Dual-version Actions 39756688525: Python 3.11/3.12 both **188/188 PASS**, including finite-B/finite-k independent Ward audit.
- Fixed-k maximum discrepancies: analytic full static vs finite B = **7.88e-11**, vs finite-k spin-density gradient = **1.25e-10**.
- Primitive midpoint n=4,8,12 maximum sampled full static BZ response residuals: **1.1573e-2**, **7.0627e-5**, **8.6852e-7**; periodic discrete spin-density derivative telescopes below **8e-18**.
- Model method VERIFIED in fixed-basis, spin-k-independent clean static protocol; whole original paper **PARTIAL**, dynamic metallic KME unresolved. Frozen output in results/frozen/round23_equilibrium_ward_audit.json.


## Round24: degenerate finite-temperature equilibrium Ward

| Stream | Method gate |
| --- | --- |
| Paper | Exact Kramers-degenerate **equilibrium static** no-current Ward as an independent negative control; whole paper PARTIAL. |
| Physics | Fréchet Fermi divided difference, spectral projector / within-doublet gauge invariance, local field and momentum derivative oracles, smooth near-degenerate limit. |
| GitHub | Round24 PR #83 **merged**: 15 new tests, frozen numerical reference, static audit and scope report. GitHub Actions 39780389586: **203/203 PASS** independently on Python 3.11 and 3.12, including the executable audit. |

The earlier Round23 static formula deliberately rejected degeneracies;
Round24 will extend the equilibrium method without calling any
finite-frequency metallic response complete. The original 188 verified
paper-specific tests remain unmodified.


### Round24 actual independent acceptance (2026-10-22)

- Stage PR #83 merged to main at commit 5d8c7814f9b5834131ccf10d62f282182e340ef6.
- Fixed 4-band spin × orbital Kramers toy model at k=(0,0), mu=0.15, T=0.22, kappa=0.7 (natural units): two exact doublets [2,2].
- Full static field response and independent spectral projector contraction residual <= 7.02081e-18; independent finite B difference <= 1.17030e-11; independent finite k spin-density gradient <= 5.47129e-11.
- Max intracluster static response component 0.00733779; max full response 0.139812. f'[E] explicitly includes n != m within an exact Kramers cluster; a diagonal-only population derivative is basis-dependent and fails its negative-control test.
- Initial CI failed five tests due an incorrect U†VU / U†SU Einstein-contraction index ordering; the implementation was corrected, and the final code passed **203/203** on both Python 3.11/3.12 with a frozen-data regression. The error is part of the visible audit history, not erased.
- Repository file: results/frozen/round24_degenerate_static_ward_audit.json. The analytic whole-BZ equilibrium current vanishes on a periodic torus; this is NOT a demonstration of dissipative finite-frequency KME or an experimental material response.
- Full author paper **PARTIAL**. Remaining: dynamic metallic charge conservation, electromagnetic contacts, Wannier physical spin+position in the same gauge, source Eq.(19), real material SI and Table I reconciliation.


## Round25: Bloch-gauge covariant static spin Ward (in independent CI)

| Main workstream | New restricted gate |
| --- | --- |
| Paper reproduction | Derived gauge-covariant **equilibrium static** Zeeman current Ward as an independent physical negative control. Full original Cui paper remains PARTIAL; this does not reproduce a new source Eq.(19) transport measurement. |
| Physics capability | Bloch connection A=W†dW, covariant physical velocity/spin-operator derivatives and explicit Zeeman-current contact for k-dependent synthetic effective-spin vertex; Kramers Fréchet response reused unchanged. |
| GitHub | Round25 PR: 10 new tests and scripted audit. Count 203 existing + 10 target = 213 only after both 3.11 and 3.12 remote CI pass. |

Next required step is a realistic same-Wannier-gauge H, physical spin and position-operator interface; current results are neither SI nor dissipative KME.
