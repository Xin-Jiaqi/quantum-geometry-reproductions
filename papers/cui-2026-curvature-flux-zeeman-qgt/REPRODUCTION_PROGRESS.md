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


## Round25: Bloch-gauge covariant static spin Ward (merged)

| Main workstream | New restricted gate |
| --- | --- |
| Paper reproduction | Derived gauge-covariant **equilibrium static** Zeeman current Ward as an independent physical negative control. Full original Cui paper remains PARTIAL; this does not reproduce a new source Eq.(19) transport measurement. |
| Physics capability | Bloch connection A=W†dW, covariant physical velocity/spin-operator derivatives and explicit Zeeman-current contact for k-dependent synthetic effective-spin vertex; Kramers Fréchet response reused unchanged. |
| GitHub | Round25 source PR #87 **merged** (squash 2bce0fb8f15b288cd864189fe55d19cf794a3cf2). Eleven new tests including frozen numerical regression: Python 3.11/3.12 **214/214 PASS** and reproducible audit on Actions run 39878438728. |

Next required step is a realistic same-Wannier-gauge H, physical spin and position-operator interface; current results are neither SI nor dissipative KME.


### Round25 physical acceptance and residuals (2026-10-23)

- Synthetic k-dependent Bloch gauge: full 2×3 static susceptibility basis covariance error ≤2.4e-16 at sampled k, including exact Kramers doublets.
- At k=(0.36,-0.29), effective projected-spin variation 0.16, T=0.22, mu=0.15, kappa=0.7: complete Ward vs independent finite-B 2.65e-11; vs finite-k spin-density gradient 4.16e-11.
- Explicit Zeeman current-vertex contact max = 0.0421665; population-only omission has error 0.0421665. The raw k derivative of the rotated Hamiltonian differs from the physical covariant velocity by 0.425737 in the chosen toy gauge.
- With purely k-dependent gauge but k-independent physical spin, the **covariant** spin derivative/contact is zero. All values are synthetic natural-unit model checks, NOT real material or finite-frequency KME.
- Frozen raw audit: results/frozen/round25_covariant_spin_ward_audit.json. Full source-paper reproduction **PARTIAL**.

## Round26 | same-Wannier-gauge physical-operator prerequisite (2026-10-23)

| Track | New bounded milestone | Still unresolved |
|---|---|---|
| **1. Paper reproduction** | A four-band SOC *synthetic* same-gauge operator input gate precedes existing Zeeman-QGT calculations; existing scientific baselines retained. **Full original Cui paper PARTIAL.** | Source Eq.(19) SI KME/material figures, original Table I and actual same-gauge material spin/position. |
| **2. Physics capability** | Validated Fourier H/S/A, orbital embedding and units; independent Wilson finite-difference covariant dH/dS, projector QGT covariance; R/-R, physical spin torque and bad-gauge/missing-connection controls. | Declared matching gauge ID is **not** independent spin-provenance authentication; true material operator ingestion, finite-q EM Ward, collisions, metallic transport. |
| **3. GitHub/CI** | Isolated Round26 13 new tests, executable audit, frozen synthetic JSON, scoped Physics Registry entry and Python 3.11/3.12 Cui CI audit step. CI/merge status must be recorded after checking actual jobs and logs. | No Round21 repeat, no rewriting the existing 214-case baseline, no inference of Cui tests from 1105 repository tests. |

At fractional q=(.183,-.271), Wilson dH residual 1.51e-10 eV·Å; spin gradient 5.83e-16 Å; QGT covariance residual 2.09e-15; **wrong spin gauge** QGT discrepancy 0.148; **missing position connection** vertex error 0.263 eV·Å. All **synthetic**, not macroscopic KME.

## Round27 | Finite-(q,omega) longitudinal charge-gauge Ward (2026-10-23)

| Progress line | New physically bounded result | Remaining |
|---|---|---|
| **1. Cui source-paper reproduction** | A necessary finite-q charge-gauge negative control for clean microscopic response on a *3×3 synthetic 2D SOC torus*. It does not reproduce source KME or revise original Table I. **Whole paper PARTIAL.** | Magnetic Zeeman KME in a conserving, contact-complete metallic continuum; materials, Eq.(19) SI and full figures. |
| **2. Physics capability** | Exact scalar/vector pure-gauge decomposition, Peierls current contact, sitewise continuity, spectral Kubo vs independent N² Liouville, nonlinear Peierls/time evolution, strong contact-omission counterexample. | Ward-consistent collision vertex/self-energy, transverse finite-q EM response, thermodynamic response, real Wannier spin and position. |
| **3. GitHub / CI** | Isolated Round27 source, 16 new tests, frozen numerical audit, scoped Physics Registry entry, new workflow step. **CI/merge acceptance determined from actual remote logs only.** | Unrelated open PR #65 untouched; old Round21/26 frozen code intact. |

At q=(1,1) on the 3×3 torus, z_drive=.42+.09i: spectral vs Liouville 3.73e-15, charge 2.22e-15, contact-complete bond current 3.43e-15. Dropping the bond-current contact yields a spurious current **1.34659** and finite-q local continuity violation **3.59468**. These are natural-unit synthetic benchmarks, not SI observables.

## Round28 | Charge-conserving linearized collision and Zeeman-response gauge guard

| Progress line | Actual bounded new capability | Still unresolved |
| --- | --- | --- |
| **1. Original Cui paper** | Clean finite-q Ward upgraded to a *synthetic linearized* finite-gamma local continuity/gauge-covariance check, plus magnetic Pauli-Zeeman drive. **Whole original paper PARTIAL.** | Microscopic collision self-energy+Ward vertex, metallic transport/Drude, real material same-gauge spin and position, source Eq.(19) SI, Table I. |
| **2. Physics capability** | Site-local projector Q conserves every charge density; external Peierls gauge covariant reference; independent orbital Liouville versus spectral fixed-point; naive RTA violates local continuity, omitted reference violates gauge covariance despite local conservation. | No autonomous CP thermal bath, detailed balance, real impurity/phonon kinetics, thermodynamic BZ convergence. |
| **3. GitHub/CI** | Isolated Round28 tests, executable independent collision audit, frozen reference after dual CI, Physical Capability Registry and a new workflow step. PR/CI status must be verified on the exact tested SHA. | Round21/26/27 and parallel PR #65 untouched; paper never auto-upgraded by test count. |

**Acceptance gate:** genuine remote Python 3.11/3.12 entire Cui suite, all Round17–28 scripts and fixed numerical negative-control benchmark. Scientific scope and derivation: [Round28](ROUND28_LOCAL_CONSERVING_COLLISION.md).
 
### Round28 measured first dual-interpreter acceptance

Initial isolated Cui remote CI: Python 3.11/3.12 **256/256 PASS**, independent audit successful. Reported gamma=0.12 benchmark: scalar density 0.1041668, scalar current 0.0389254, independent Liouville residual 2.59e-15, pure-gauge current 1.15e-15, wrong global RTA charge-Ward violation 0.01171157, and missing gauge-reference spurious current 0.02747677. The finalized frozen numerical regression brings the expected suite to **257**; final status must be inspected in GitHub Actions independently before merging. Full Cui source **PARTIAL**.

## Round29 | Canonical thermal GKSL versus local charge conservation (2026-10-23)

| Track | Independently testable bounded milestone | Remaining |
|---|---|---|
| **1. Cui paper reproduction** | Compare two explicit CP toy baths and expose which Ward/source terms each supports; not a published Eq.(19) real-material or experimental response. **Full paper PARTIAL.** | Microscopic electron/phonon/disorder self-energy and conserving vertex, fermionic metallic occupations, SI response, original Table I and materials. |
| **2. Physics capability** | Local GKSL dephasing conserves onsite bath charge but is not finite-T Gibbs preserving; canonical KMS pair-transition GKSL preserves Gibbs but creates a local bath source. Joint jump/H gauge covariance and explicit bath-vertex finite difference distinguish correct and stale jump bases. | Only two-site one-particle canonical synthetic baths; no thermodynamic limit, actual material scattering or current vertex reconstruction. |
| **3. GitHub** | Round29 isolated code/tests/audit, scoped Physics Registry, explicit negative controls and full pre-existing Cui CI matrix. Final test count and acceptance determined from exact remote logs. | Round21,26–28 and parallel numerical BZ PR #65 left unchanged. |

Derivation: [Round29 scope](ROUND29_THERMAL_LINDBLAD_WARD.md). Full source Cui paper **PARTIAL**.

### Round29 first frozen-audit checkpoint

The first remote 273-case test and audit run exposed the finite two-site bath tradeoff without changing original Cui physics: KMS pair residual 3.47e-18; Gibbs thermal-jump stationarity 3.69e-17; local CP dephasing Gibbs failure 0.023248; thermal bath site-charge source 0.074817; gauge-stale jumps false covariance 0.010845; omission of bath vertex false pure-gauge residual 0.017868. Newly frozen data and one numerical regression increase the final requested package count to **274**. Final merge requires real dual-version CI on the exact SHA. **Whole source paper remains PARTIAL.**

## Round30 | Interacting many-fermion thermal sectors and collision Ward (2026-10-23)

| Stream | New physically bounded milestone | Still missing |
|---|---|---|
| **1. Original Cui reproduction** | Exact four-mode interacting fermionic Fock test with Pauli blocking and fixed-N energy-pair KMS bath; an upstream real Fermi many-body prerequisite, **not** original KME. **Whole paper PARTIAL.** | Material-specific finite-q dissipative KME, source Eq.(19) SI, Table I and published figures. |
| **2. Reusable physics** | Exact CAR + 16-state grand Gibbs, number-sector jumps, Pauli/no-double-occupancy, interacting vs independent one-body Fermi limit, total-number Ward vs explicit local bath source, gauge and **magnetic** bath vertices with independent finite differences. | Microscopically derived phonon/disorder self-energy and Ward ladder, fixed-mu reservoir, real Wannier spin source, thermodynamic metallic response. |
| **3. GitHub** | Isolated Round30 module/tests/audit/frozen regressions, scoped Registry, dual Python 3.11/3.12 workflow and all Round17–29 historical steps. Final counts/status require exact CI logs. | Parallel open PR #65 and prior Round21/26–29 untouched; no full-paper promotion. |

**Bounded method**: [Round30 physics](ROUND30_INTERACTING_FOCK_THERMAL_BATH.md). Number-preserving KMS jumps do **not** select a unique chemical potential; the synthetic Fock ensemble is not yet a bulk metal.

### Round30 remote first-run baseline and final frozen gate

Dual Python 3.11/3.12 **291/291 PASS** in first CI run 39899779828, including all past audits and the new many-fermion test. Adding an independent frozen numerical regression gives the final expected Cui count **292**; final PASS/merge claims require fresh exact-SHA remote evidence. Model diagnostic: 54 Fock number-preserving KMS jumps, interaction-dependent occupations (difference 0.075984), grand KMS 1.73e-18, Gibbs stationarity 2.36e-17, total N source ~7e-18 but local bath site source 0.271191, and missing Zeeman bath contact error 0.005563. Explicit finite interacting Fock only, **Cui source PARTIAL**.

## Round31 | Harmonic-boson spectrum to interacting Fock collision (2026-10-23)

| Track | Bounded new verifiable advance | Not accepted |
|---|---|---|
| **1. Cui source** | First *explicit model bath correlation* for energy-dependent fermion transition rates; dynamic finite-z pure electromagnetic gauge Ward and magnetic field collision derivative on interacting synthetic Fock. **Whole paper PARTIAL.** | Source-specific metal KME Eq.(19), SI units/material figures and Table I; real phonon Wannier coupling. |
| **2. Physics capability** | Finite oscillator explicit exact correlation/absorption/emission KMS versus continuous Ohmic J; Born-Markov-secular number-preserving jumps including zero-frequency collective bath, grand-Gibbs, finite-z charge gauge and physical current contact, bath source and Zeeman bath-vertex negative omissions. | Spectral J is declared, not ab-initio; no actual phonon eigenmodes, Lamb shift, nonsecular clusters, unique chemical potential or translationally conserving electronic metal. |
| **3. GitHub & audit** | Isolated source, added tests, executable independent Round31 audit, frozen model output after first remote run, Registry and 3.11/3.12 Cui workflow audit. | Previous Round17–30 frozen outputs never rewritten; parallel BZ convergence PR #65 untouched. Final counts/CI/merge require logs on exact SHA. |

[Complete Round31 derivation and scoped claims](ROUND31_BOSONIC_DAVIES_DYNAMIC_WARD.md).

### Round31 remote baseline (first clean dual CI)

After fixing an incorrect oscillator field alias that multiplied the independent correlation test by four, the first corrected GitHub-hosted [Round31 run 39904654968](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39904654968) passed **308/308** on both Python 3.11 and 3.12 and all past audit steps. Frozen first-run metrics show oscillator thermal KMS ~2.2e-16, Fock grand KMS 3.5e-18, local bath charge source 0.406784, finite-z full pure-gauge current 5.9e-17, missing dynamic gauge bath vertex 0.010774, source-complete finite-z magnetic Ward 1.5e-17 and omitted magnetic collision vertex 0.015032. One additional strict frozen-data regression yields expected **309** final tests, to be checked on exact final HEAD. No promotion of source Cui **PARTIAL**.

## Round32 | Material electron-phonon input provenance and gauge gates (2026-10-23)

| Track | Bounded work implemented | Unresolved |
|---|---|---|
| **1. Original Cui reproduction** | No actual material EPW/DFPT/Wannier spin or position source files currently present; source status **PARTIAL**. Formalize what must be provided rather than fabricate source tensors. | Full original real-material spectra, g, metal KME Eq.(19) SI, figures and Table I. |
| **2. Reusable physics** | Explicit normalized external export (not raw EPW binary reader): H/g/spin/connection and phonon energies/phase/units, k+q closure, +/-q, SHA256 evidence, source-stage labels. Synthetic Pauli-blocked finite-width golden-rule rates have independent matched H/g gauge positive oracle and unmatched H gauge negative control. | Source file mappings/electron/phonon joint gauge cannot be authenticated from IDs and checksums; position connection needs k-gauge derivative; Gaussian broadening not exact detailed balance/vertex; no acoustic/polar divergence or BZ convergence. |
| **3. GitHub** | New isolated Round32 source/tests/audit, synthetic-only frozen benchmark after first remote run, scoped Registry and dual Python workflow. Exact merge and test count require actual remote logs. | All Round17–31 preserved, open BZ convergence PR #65 unaffected. |

[Round32 scientific and source-integrity boundary](ROUND32_EPW_MATERIAL_SCATTERING_GATE.md).
