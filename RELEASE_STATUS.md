# Quantum geometry releases (2026-10-20)

| Paper | Remote module | Software gate | Scientific status |
| --- | --- | --- | --- |
| Jia 2026 PRL NMEE | `papers/jia-2026-nonlin-magnetoelectric-edelstein` | `scripts/verify_release.py` | PARTIAL |
| Xiang 2025 PRL IGMC | `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current` | `scripts/verify_round9.py` | PARTIAL |
| Xiang, Jin & Wang 2025 PRL spin quantum geometry | `papers/xiang-2025-spin-quantum-geometry` | 70 tests incl. V32 SM Eq.(2) four sources, V31 Laurent and V30 finite-eta/PT; Python 3.11/3.12 CI | PARTIAL; curated paper model plus verified method subsets; physical metallic DC and published absolute factors OPEN |
| Xiang 2026 PRB magnetotransport | `papers/xiang-2026-quantum-geometric-map` | `scripts/verify_release.py` | PARTIAL; GitHub **curated core + migrated nondegenerate spin/orbital method subset** |

For the PRB, the complete 30-round original **offline** v30 checkpoint is separate from the independently readable/callable core published in this repository. It includes additional model branches, scans, plots, and scientific logs; those are **not** claimed to have been uploaded here. The core implements selected independently checked methods and leaves absent advanced interfaces explicitly unimplemented.

Open physics:
- Spin QGT PRL 135: Fig.1(a) ~2× axis, SM Eq.55/56 and SM Eq.69/70 algebra, finite-eta DC interpretation; selected finite-eta/PT/Laurent methods are public, but original V14 four-source audit and full V17–V19 non-Abelian differential derivations remain offline.
- NMEE: static-versus-driven density-matrix response and CuMnAs PT-degenerate DC.
- IGMC: source Fig.2(d) normalization and response convention.
- PRB: Eq.(11) vs Fig.2(b) ~2π SI amplitude, physical steady-state PHE collision protocol, SM Eq.(52)/(56) internal algebra, multiband orbital response and Rashba QMQ attribution.

Do not infer a full-paper science PASS from the software CI result.

## 2026-10-20 advanced method migration

The PRB module additionally contains reusable four-band spin-MNHE (SM 16–19 vs 29), orbital PHE Eq.56, orbital MNHE 2D SBC and 3D Eq.49/52 diagnostics. These are callable via `python scripts/reproduce_advanced.py --case ...` and covered by independently scoped GitHub CI. The disagreement branches remain OPEN and are never re-labelled as full-paper scientific acceptance. The remaining v30 offline exact finite-system Kubo continuity/contact benchmark, limit-order scans, full Kubo–Středa terms and historical raw figures remain outside this remote subset.


## 2026-10-20 SCBA/Kubo/Ward migration (v34)

Migrated the original v11/v17/v18/v19 kernels and 48 regression tests: matched SCBA self-energy/ladder, RA control, four-channel symmetric spectral Kubo, and finite-(q,omega) resolvent Ward. Latest method scope: `papers/xiang-2026-quantum-geometric-map/docs/SCBA_KUBO_WARD_MIGRATION_V34.md`. Full Kubo–Středa, real-frequency physical continuity, and original paper's driven-DC PHE remain open.


## 2026-10-20 v35 exact Kubo / SCBA regulator migration

The original v20 exact finite-disorder real-frequency longitudinal Kubo continuity/contact tests and independent SCBA regulator-order scan are now in the PRB module, with 18 original tests and JSON entry points. Their PASS status is a mathematical/method validation, not acceptance of the original paper's transverse steady-state PHE. See `docs/EXACT_KUBO_LIMITS_MIGRATION_V35.md`.

## 2026-10-20 v36 second-order Liouville migration

Migrated unchanged original v21–v23 physics modules and 66 historical tests, plus frozen JSON numbers, source-file SHA256 manifest and three callable JSON cases. Finite-q scalar electric and homogeneous length-gauge electric remain distinct protocols. Original continuum steady-state DC PHE and complete Kubo–Středa remain OPEN.

## 2026-10-20 v37 continuum/SI/Fig. 2 migration

Migrated the original v24–v27 continuous tilted Dirac, published Eq. 11/SM Eq. 20, independent direct SI Liouville and approximate published Fig. 2(b) digitization. All **62 original tests** are in the PRB module, with exact-source Git blob hashes and four JSON cases. The ~2pi Eq.11-versus-published-ordinate discrepancy and full physical steady-state PHE remain OPEN.

## 2026-10-20 Xiang 2025 spin-QGT Eq.(13) migration

The Xiang/Jin/Wang PRL 135, 146303 paper-model module now includes the **printed reduced Fig. 1(a) intrinsic Eq. (13)** in addition to the two-valley NDSC and supplemental Dirac models. It freezes the Fermi-sea and Fermi-surface branches at μ=0,12,20,40 meV and verifies Pauli-algebra kernels against Bloch-eigenbasis matrix elements. GitHub Actions Python 3.11/3.12 pass 20 local-scoped tests. The source Fig. 1(a) factor ~2, omitted four-band terms and metallic driven-DC response remain **OPEN**. The offline PT/multiband methods have *not* been migrated in this update.


## 2026-10-20 v38 clean antisymmetric sea benchmark

A new `streda_clean_sea_v38` method computes the *bare clean* two-band antisymmetric Green-function Fermi-sea integral and independently checks it with direct real-energy quadrature, Chern link phases and Berry-curvature integration. Twelve tests are included. The metallic rational pole boundary is shown separately and is NOT claimed as a complete disorder-dressed Středa surface term. Original paper physical DC EB and generic Kubo–Středa remain PARTIAL/OPEN.


## 2026-10-20: v39 energy-dependent SCBA model audit

Published a named QWZ noncrossing scalar-disorder SCBA method with self-energy energy derivative, same-disorder RR density Bethe–Salpeter vertex, finite-difference proof, 14 method/frozen tests, an exact finite-disorder torus spectral control and a JSON entry point. \(\partial_E G^R=-G^R(I-\partial_E\Sigma^R)G^R\) is essential for the defined sea kernel; the disorder-dressed transverse current vertex, full Středa I sector and the original paper's nonlinear EB response remain OPEN.

## 2026-10-20 Xiang 2025 spin-QGT V30 methods migration

PR #12 merged: the Xiang/Jin/Wang 2025 PRL paper module now contains independent N-band (N=2–6 regression) **finite-eta** Liouville band/orbital solvers, exact synthetic PT rank-2 projector checks, U(2) gauge covariance, a 6-band three-doublet test, and a V16 offline staggered-spin baseline (2.0232340477418247 at the declared synthetic-model parameters). Paper-core+methods scope: **32 tests, GitHub Actions Python 3.11 and 3.12 PASS**; see `papers/xiang-2025-spin-quantum-geometry/docs/METHODS_ACCEPTANCE_V30.md`. This is a curated/reimplemented subset, **not** a full upload of the V14–V19 offline derivations or a validated intrinsic tau^0 driven-DC solver. The original Fig. 1(a) ~2x absolute discrepancy, metallic DC protocol and the supplemental algebra issues remain open. No stacking, intercalation or real-material DFT claims.

## 2026-10-20: v40 current-vertex Hall method benchmark

In the PRB quantum-geometric-map module, the noncrossing four-Pauli-basis current ladder is now paired with the v39 SCBA self-energy, across RA/AR/RR/AA channels. The result includes the self-energy-matched RA density Ward, periodic-momentum RR current Ward, exact finite-torus bare RA control, and the **vanishing antisymmetric 4-channel symmetric-spectral** diagnostic. The *isolated* antisymmetric RA bubble is not identified as a physical Hall conductivity. Full disordered Kubo-Bastin/Středa plus the original second-order EB DC remain OPEN.

## 2026-10-20: v40b independent RR spatial Ward

Added independent uniform-vector-potential finite-difference tests of the RR SCBA current ladder in both x and y, checking \(\partial_A\Sigma^R=C^{RR}\) and \(\partial_A G^R=G^R\Gamma^{RR}G^R\). Four new tests were added. This is separate from the v40 RA current response and does not constitute physical Hall transport. Full disordered Kubo–Středa and original two-field EB DC remain OPEN.

## 2026-10-20 Xiang 2025 V31 Laurent

Migrated the **adiabatic finite-eta** Laurent decomposition A/eta^2+B/eta+C using isolated-energy-cluster projectors, including exact PT three-doublet six-band benchmarks and independent dH/dk population-gradient identity. Original offline V19 nonzero extra term is preserved and frozen. A controlled counterfactual outer regulator ratio shows the constant C depends on the specified adiabatic protocol, so **C is not claimed to be a general intrinsic driven-DC tau^0 coefficient**. PR #13 merged after 46/46 tests succeeded on Python 3.11 and 3.12. SM Eq.(2) four-source migration and full physical DC limit are still open.

## 2026-10-20: v41 exact finite-torus linear Hall benchmark

Merged PR #14 (head `70b3d8d`, squash commit `b3d2821`) after the **291/291 PRB software tests PASS** at Actions run `39435808942`, including eight new exact-Hall gates. The independent same-Hamiltonian full-spectrum **linear** Hall Kubo response on finite disordered QWZ tori is checked against a separate occupied-subspace flux-twist Chern number, finite-frequency retarded correlator, model-specific zero mixed diamagnetic contact and current finite differences. Frozen representative: `L=4,m=-1,W=.12,seed=1947,ntheta=7,mu=0`, fixed-flux Hall `-1.287619933`, twist-averaged Hall `-1.000000019`, FHS Chern `-1` in units `e²/h`. This is **not** separate physical Středa surface/sea, thermodynamic disorder-averaged conserving SCBA or the original EB PHE. Scientific status is **PARTIAL**.

## 2026-10-21 Xiang 2025 V32: original SM Eq.(2) four-source audit

PR #16 merged. An independent nondegenerate 2–5-band explicit evaluation of all *four* SM Eq.(2) source groups agrees with fixed-orbital Liouville (max absolute error < 8e-9 in the frozen reference), with a further separate V30 solver crosscheck. The source Eq.(2) eta^0 result differs from printed SM Eq.(6) by the occupation-gradient term -i*r_mn*(f'_n-f'_m)/(epsilon_m-epsilon_n)^2; this is recorded as a **source-algebra discrepancy under the specified eta/2eta prescription**, not a universal claim of incorrect physical DC response. 70/70 tests and both Python 3.11/3.12 GitHub CI pass. Previously unresolved Fig.1(a) normalization, SM Eq.55/56, Eq.69/70 and metallic collision/limit-order physics remain OPEN. No DFT, stack or intercalation simulation claimed.


## 2026-10-21: v42 exact finite-disorder Streda linear Hall

PR #18 merged as squash commit `c39aed71` after remote PRB Actions run `39530869244` confirmed **301/301 PASS**, including ten new v42 tests. For exactly the v41 frozen-disorder QWZ finite torus, standard linear Smrčka–Středa I(a), Hall-odd I(b)=0, and II are independently evaluated and their **sum** converges to v41's full Lehmann Hall in the external regulator `eta->0` limit. The II eigenpair primitive agrees with a separate direct real-energy dense Green matrix integral and the v38 clean Bloch result; I(a) matches the independent inverse-matrix current trace. The surface/sea split has overlap, hence the separate terms are not independently observable intrinsic/extrinsic mechanisms.

For `L=4,m=-1,W=.12,seed=1947,theta=(.11,.33),mu=1.1`, `eta=.003`: I(a) `-0.115816488`, II `-0.191177158`, sum `-0.306993647`, v41 exact `-0.307005888` e²/h. Near finite sample Fermi eigenlevels, `eta` regularization is highly sensitive. Six exact-disorder seeds are frozen as a **finite sampling diagnostic**, not infinite-size or SCBA acceptance.

**Scientific status PARTIAL.** Conserving energy-dependent SCBA/Středa I+II with matched current/density vertices, contact and thermodynamic dc order of limits remain OPEN, as do original PRB nonlinear EB response, material-specific calculations, and other historical discrepancies. See PRB paper-local `docs/EXACT_STREDA_V42.md`.

## 2026-10-21: v43 matched SCBA *surface I(a)* versus exact disorder ensemble

PR #20 was merged as `cd262de2` after Actions run `39533046760` confirmed **312/312 PRB tests PASS**, including 11 new tests. The v43 module reuses v39 energy-point SCBA and v40 one-leg RA/AR current BSE with identical W, and independently compares to v42's dense real-space exact disorder-averaged **surface I(a)**. It tests left/right dressed-vertex equivalence, separately iterated BSE, RA density resolvent Ward, and seed—not twist—uncertainty. At W=.12, μ=1.1, η=.12: bare -1.27807, single RA ladder -1.43032, exact L4 24-seed 7x7 twist mean -1.50028±.05502 in e²/h. L3–L5 scans of 64–96 seeds retain a visible discrepancy. **Scientific status PARTIAL:** this is not a complete conserving averaged Středa I+II, transverse DC, or original nonlinear mixed EB PHE. Next scientific gate: independently derive and validate disorder-dressed sea II/contact/field-density responses on the same protocol before any full-Hall claim.

## 2026-10-21: v44 provisional SCBA I(a)+II versus exact disordered Hall

PR #22 merged as squash `5348caed` after GitHub Actions run `39535057861` confirmed **322/322 PASS**, including ten v44 tests. Added a same-W energy-dependent RR SCBA density-vertex `dG/dE` sea integrator (420 Gauss nodes, E_lower=-12) and a **provisional** combination with the v43 one-leg RA surface I(a). Independent numerical Sigma/G finite differences, exact clean two-band pole primitive and finite-disorder v42 full Středa I+II sample averaging verify the *declared model approximation*. Frozen m=-1,mu=1.1,eta=.12,SCBA 35x35,exact L4/24 independent seeds/7x7 twists: W=.03 candidate -0.376123538 vs exact -0.389852385±.021380310; W=.12 candidate -0.456059226 vs exact -0.483689664±.027787510 e²/h. **Not conserving-field theory acceptance**: cross diagrams, transverse Ward/contact, thermodynamic metallic DC and the published second-order EB steady-state PHE remain OPEN. Overall science **PARTIAL**.

## 2026-10-21: v45 original symmetric PHE physics realignment

PR #23 merged as `9a696fac` after full PRB Actions run `39537421781` confirmed **332/332 PASS**. The independent direct main Eq8 angular integral, eigenvector-checked main Eq10 Zeeman curvatures, full-plane fixed-tilt Bx->ky translation, and separately cancelling velocity/population derivatives now have callable Python, 10 tests and frozen JSON. At beta=1/8, mu>0, Eq8-normalized coefficient **+0.00250641590669**; main Eq11 geometry is **−0.00250641590669** under the explicitly printed band/field conventions. Fully equilibrated static-B Dirac Drude is field-invariant by a momentum shift. This exact model identity does **not** prove paper Eq11 is physically wrong: its Liouville magnetic-field preparation, selected τ-scaling and contact/collision terms must be derived *as the same protocol*. Prior QWZ v41–v44 SCBA/Středa antisymmetric linear Hall remains a distinct methodological extension, **not** the published symmetric mixed EB spin PHE. The paper itself uses phenomenological RTA and leaves microscopic bath theory to future work. Fig2(b) magnitude discrepancy (~2pi) still OPEN. Overall **PARTIAL**; next original-paper gate is the complete same-protocol two-field RTA response, not another unrelated Hall benchmark.

## 2026-10-21: v46 original mixed PHE preparation-dependent RTA

PR #24 merged as squash commit `4847a319` after PRB Actions run `39539072615` confirmed **343/343 PASS** (11 new tests). The original symmetric mixed EB PHE of Xiang PRB is now tested on the **same fixed-tilt Dirac Hamiltonian** in two deliberately different scalar RTA protocols. An equilibrium-following magnetic bath target `f(H_B)` produces complete mixed derivative at most numerical residual (~4.57e−10) from source preparation (+0.076275298) and B commutator (−0.076275298) under beta=.125,mu=.3,T=.055,eta=.065. This exactly follows the `B_x -> k_y` shift Ward identity for the infinite-UV, fixed-tilt model. A different relaxation target fixed at zero-field `f(H0)` gives B-then-E (+0.110050304), E-then-B (−0.076275298), combined +0.033775006; as eta tends to zero its tau-linear coefficient has the same sign and approximate magnitude as the original formal Eq.(8) at small beta, but it is **not identical** and the beta-dependent deviation grows strongly (ratio .764 at beta=.6). No fitted prefactors or silent sign changes.

**Acceptance is only the two named RTA model protocols, not published original PHE.** Neither one has been identified with the authors' full Supplemental Material collision prescription. The v45 printed Eq8/11 sign branch, v37 Fig2(b) ~2pi normalization, and physical preparation/limits remain OPEN. The earlier v41–v44 QWZ charge Hall study is distinct method development, not a prerequisite for the authors' phenomenological original spin PHE. Whole-paper scientific status **PARTIAL**. Next priority: read/verify the original SM collision/field-frequency prescription and reproduce it *unchanged* as a third protocol.
