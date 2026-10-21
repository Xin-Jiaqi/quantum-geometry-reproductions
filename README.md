> **v54 PHYSICS-FIRST CAPABILITY RELEASE:** The project now has a separate [physics capability registry](physics/README.md), explicitly distinct from article-specific reproduction history. A reusable **nondegenerate N-band** Bloch geometry core accepts H(k) and matched dH/dk, and independently calculates per-pair/per-band metric, Berry curvature and energy-geometric orbital moment. **VERIFIED** within fixed orthonormal k-independent orbital basis by an independent finite-difference projector oracle, band U(1) rephasing, constant orbital-basis unitary invariance, and true two-/three-band Hamiltonians. A complete m=l/m≠l decomposition of **one named curvature×orbital-moment subterm** is likewise verified, independently reproducing v53 full -0.00079945105, diagonal -0.00154218315, and cross +0.00074273210 for a selected model. **This is NOT full magnetic positional shift, orbital MNHE or magnetic response; those remain UNRESOLVED.** [v54 API and scope](papers/xiang-2026-quantum-geometric-map/docs/INDEPENDENT_MULTIBAND_ORBITAL_GEOMETRY_V54.md) · [PR #34 merged](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/34) · [438/438 remote tests PASS](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39574187471). Original paper scientific audit stays PARTIAL.

> **PRB v49, 2026-10-21: Authentic original Supplemental Material now verified.** The user supplied the 16-page SM and 8-page published PRB (SHA256 source provenance pinned in the isolated PRB subtree). The **actual source** regulates one electric-field adiabatic pole as tau (SM Eq15/20). Independent direct eigenbasis calculations find a **literal overall sign discrepancy between SM Eq20 and printed Eq31** after applying their printed velocity/Berry relation, and an **opposite spin helicity** between main Eq10 and SM Eq62. At beta=1/8,mu=.3,T=.03, on main Eq10 helicity direct Eq20 gives −0.002506188247 while printed Eq31 gives +0.002506188247 in normalized physical-B xy; Eq32 with the SM helicity agrees with signed Eq11. The previous v47 numerical eta-scaled symmetric coefficient approaches the **Eq20 + main helicity** branch. No source formula has been silently corrected. These findings are **local source-level inconsistencies, not proof the full physical phenomenon is invalid**. Complete original-paper status **PARTIAL**. [v49 audit](papers/xiang-2026-quantum-geometric-map/docs/ORIGINAL_SM_SOURCE_AUDIT_V49.md) · [PR #28](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/28) · [376/376 PRB regression PASS](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39556012282).

> **SCIENTIFIC ERRATUM, 2026-10-21, v47:** An incorrect sign in the v46 electric-y equilibrium source invalidated its originally frozen mixed-EB results despite green CI. Old v46 frozen numbers and XY-only Eq8 comparison are **RETRACTED**. Corrected source, BOTH XY/YX response channels and independent finite-B Liouvillian now recover the signed *main Eq11* symmetric spin-PHE geometry to ~5×10⁻⁵ relative across beta=.025–.6 in a declared frozen-H0 RTA model. The author Supplemental Material bath prescription, literal Eq8 signed conventions and Fig2(b) absolute scale remain OPEN. See [v47 technical correction](papers/xiang-2026-quantum-geometric-map/docs/ORIGINAL_PHE_SYMMETRIC_SOURCE_CORRECTION_V47.md). CI [353/353](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39543408276), [PR #26](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/26). Overall scientific status PARTIAL.

# Quantum Geometry Reproductions

Independent, **paper-by-paper** reproductions of quantum geometry and magnetotransport. Each paper has an isolated Python package, response protocol, tests, and scientific acceptance status. **CI PASS is not the same as full scientific reproduction.** These are not official author implementations.

## Paper index

| Paper | Available independent capability | Scientific status |
| --- | --- | --- |
| [Jia et al., PRL **136**, 226302 (2026): NMEE](papers/jia-2026-nonlin-magnetoelectric-edelstein/) | Kane–Mele/CuMnAs, spin response, printed vs full-static expressions, non-Abelian diagnostics, Wannier input preflight | **PARTIAL:** equilibrium vs driven DC remains open |
| [Xiang et al., PRL **134**, 116301 (2025): IGMC](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/) | Rashba/Dirac, Zeeman geometry, finite-frequency response, 122 magnetic point groups | **PARTIAL:** Fig. 2(d) amplitude and some conventions remain open |
| [Xiang, Jin & Wang, PRL **135**, 146303 (2025): Spin Transport Revealed by Spin Quantum Geometry](papers/xiang-2025-spin-quantum-geometry/) | Two-valley Dirac NDSC/reduced intrinsic Eq.13, independent SM integrals, finite-eta N-band Liouville and PT projector/U(2) method subset | **PARTIAL:** Fig. 1(a) absolute factor, metallic DC protocol and SM algebra remain open; initial curated core only |
| [Xiang et al., PRB **113**, L201406 (2026): Quantum Geometric Map of Magnetotransport](papers/xiang-2026-quantum-geometric-map/) | Dirac/Rashba, Liouville, spin/orbital responses, SCBA + ladder, four-channel spectral Kubo, finite-(q, omega) Ward, **exact finite-disorder linear Hall + twist Chern (v41)** | **PARTIAL:** full Kubo-Streda, physical dc continuity and real-material multiband acceptance remain open; this is a curated subset of the v30 archive |

The two Xiang papers are **different physics**: 2025 IGMC concerns magnetic-field-driven gyrotropic AC current; 2026 PRB concerns spin/orbital magnetotransport under electric and magnetic perturbations. Do not merge them simply because they share authors and quantum geometry terminology.

## Layout

```text
papers/
  jia-2026-nonlin-magnetoelectric-edelstein/         # NMEE, existing frozen module
  xiang-2025-spin-quantum-geometry/                    # Xiang PRL 135 (2025), curated Dirac paper core
  xiang-2025-intrinsic-gyrotropic-magnetic-current/   # IGMC, existing frozen module
  xiang-2026-quantum-geometric-map/                   # PRB magnetotransport, independent
studies/                                               # Cross-paper research
.github/workflows/                                     # Path-scoped CI for each module
```

## Run

Use Python 3.11; create a separate virtual environment per paper.

```bash
cd papers/jia-2026-nonlin-magnetoelectric-edelstein
python -m pip install -e .
python scripts/verify_release.py
```

```bash
cd papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
python -m pip install -e .
python scripts/verify_round9.py
```

```bash
cd papers/xiang-2026-quantum-geometric-map
python -m pip install -e .
python scripts/verify_release.py
```


```bash
cd papers/xiang-2025-spin-quantum-geometry
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/reproduce.py
```

See each module's README and [PRB physical integration contract](papers/xiang-2026-quantum-geometric-map/docs/INTEGRATION.md). The published PRB core uses its own import namespace `qgmap` and never modifies NMEE or IGMC frozen benchmarks.

## Reuse without overstating physics

- Choose the physical protocol first: intrinsic vs extrinsic, static equilibrium vs driven DC vs AC, magnetic-field preparation, temperature, relaxation, and degeneracy.
- Preserve source equations, units, signs, and convergence evidence. An independently generated model figure does not equal original publication data.
- Never silently replace the printed formula with an independently preferred correction. Label branches and keep disputed factors visible.
- For a real Wannier-material response, `wannier90_hr.dat` alone is **not sufficient**. Same-gauge projected positions/Berry connection, velocities, spin/magnetic moment operators, derivatives where necessary, occupations, and validated collision assumptions are required.
- Publisher PDF/SM, unlicensed figures, private research archives and unpublished author communications are not mirrored here.

## Workflows and licensing

See [.github/workflows](.github/workflows). Workflows are scoped to their own module. The existing [IGMC module's BSD-3-Clause license](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/LICENSE) does not automatically license the entire repository. No root license is asserted.

See [RELEASE_STATUS.md](RELEASE_STATUS.md) for the cross-paper release status and outstanding scientific gates.


The PRB module now also includes an independently tested SCBA/ladder, four-channel spectral Kubo and finite-q, finite-ω resolvent-Ward method subset; see [PRB integration notes](papers/xiang-2026-quantum-geometric-map/docs/SCBA_KUBO_WARD_MIGRATION_V34.md). These do not constitute full Kubo–Středa or scientific acceptance of the original dc PHE.

The PRB module now also exposes the independently benchmarked v20 full-eigenstate longitudinal Kubo continuity/contact checks and SCBA eta limit-order diagnostics. These are *not* a replacement for a full electromagnetic Kubo–Středa calculation.

The PRB module also exposes v21–v23 second-order EB response methods (finite-q exact mixed Kubo, dynamic two-color orderings and homogeneous length-gauge Liouville with the SM Eq.20 pole). All are independently verified synthetic-model subsets, not full steady-state PHE acceptance.

The PRB reproduction also contains the original v24–v27 continuum Dirac, direct-SI Liouville and Fig. 2(b) source-curve diagnostics. The unresolved factor of approximately \(2\pi\) between Eq. (11) and the plotted ordinate is preserved as an OPEN issue, not silently corrected; see [v37 method notes](papers/xiang-2026-quantum-geometric-map/docs/CONTINUUM_SI_FIG2_MIGRATION_V37.md).

The PRB module now includes an independent **clean antisymmetric Hall Fermi-sea** benchmark against integer Chern numbers; it is separate from v18's symmetric dissipative Kubo and from the original paper's EB response. See [v38 scientific scope](papers/xiang-2026-quantum-geometric-map/docs/V38_CLEAN_STREDA_AUDIT.md).

The PRB module also includes the v39 energy-dependent SCBA derivative and RR density Ward benchmark, with separate finite-torus exact disorder comparisons. This is not a full disordered Hall conductivity or original mixed-field EB response. [v39 method limits](papers/xiang-2026-quantum-geometric-map/docs/SCBA_ENERGY_DERIVATIVE_V39.md).

The PRB module now provides a v40 same-disorder RA/AR/RR/AA **current-vertex** calculation, plus a Hall-zero check proving that the antisymmetric RA correlator alone is **not** a physical Hall conductivity. [Scientific scope](papers/xiang-2026-quantum-geometric-map/docs/CURRENT_VERTEX_V40.md).

The PRB v41 module adds the **full exact finite-torus linear Hall Lehmann response** and an independent disorder-compatible occupied-subspace **twist Chern** reference. The PRB CI tested 291 cases at the PR head. It remains **PARTIAL**: neither a conserving disorder-averaged Středa I/II decomposition nor the original nonlinear EB PHE is complete. [v41 physical protocol](papers/xiang-2026-quantum-geometric-map/docs/EXACT_HALL_V41.md).

The PRB v42 module adds the **exact finite-sample linear Hall Smrčka–Středa I(a)+I(b)+II decomposition** in the **same** v41 disordered QWZ torus and current convention. Independent inverse-matrix current traces, direct energy integration, v38 clean-Bloch sea, and the full Lehmann limit were verified; 10 new tests brought remote PRB GitHub Actions to **301/301 PASS**, run [39530869244](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39530869244). The standard surface/sea contributions cancel strongly and are *not* uniquely physical mechanisms. **Still PARTIAL:** conserving disorder-averaged Středa with matched self-energy/vertices, physical metallic DC, and the published second-order EB response remain open. See [v42 protocol](papers/xiang-2026-quantum-geometric-map/docs/EXACT_STREDA_V42.md).

The PRB **v43** study now checks the **matched noncrossing SCBA Středa I(a) *surface* current vertex only**, versus v42's independent frozen-disorder, twist-averaged finite-torus surface correlator using **the same W and bare eta**. At W=.12, mu=1.1, eta=.12, L4/24 seeds/7x7 twists: bare SCBA -1.27807, one-ladder -1.43032, exact -1.50028±.05502 e²/h, **a mismatch remains**. GitHub Actions [39533046760](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39533046760) passed **312/312** PRB tests (11 new). [v43 scope](papers/xiang-2026-quantum-geometric-map/docs/SCBA_SURFACE_V43.md). **Do not promote to complete conserving disorder-averaged Středa Hall or published EB response; overall PARTIAL.**

The PRB **v44** method now calculates a *provisional, model-only* SCBA Středa I(a)+II sum, with the energy-dependent RR density vertex solved anew at every integration energy. Independent finite-difference energy derivatives and a W=0 analytic two-band Green function verify it; full v42 finite-disorder I+II samples provide a separate benchmark. At W=.03,.12 the candidates are -0.376124 and -0.456059 versus exact -0.389852±.021380 and -0.483690±.027788 e²/h. **The result is NOT a demonstrated conserving disorder-averaged Hall conductivity or the original PRB second-order EB response. Overall scientific state PARTIAL.** [v44 protocol](papers/xiang-2026-quantum-geometric-map/docs/SCBA_SEA_V44.md). PRB CI [39535057861](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39535057861) **322/322 PASS** (10 new model tests).

## 2026-10-21 v45: Restore original spin-PHE response target

The Xiang 2026 PRB central response is **symmetric** mixed-field spin PHE (Eq.8/11, τ-linear) in a **phenomenological RTA** model. The extensive v41-v44 QWZ **antisymmetric linear charge Hall** SCBA/Středa work is useful but an **independent extension**, not an original-paper reproduction prerequisite. Independent v45 source and eigenbasis checks reproduce the main Eq8 signed geometry coefficient +0.00250641590669 versus literal main Eq11 −0.00250641590669 at beta=1/8 and mu>0; this is a labelled convention question, not a proven correction. The fixed-tilt surface Dirac with uniform Zeeman Bx obeys an exact ky-translation identity, so a **fully equilibrated static-B constant-tau Drude** response over all momentum is invariant, with nonzero velocity/population derivative pieces cancelling. This **does not establish** that the source's two-color/adiabatic Liouville PHE is wrong: a same-protocol full EB collision/preparation comparison is still OPEN. The Fig2(b) approx. 2pi factor remains unresolved. [v45 original-PHE protocol](papers/xiang-2026-quantum-geometric-map/docs/ORIGINAL_PHE_PROTOCOL_V45.md). GitHub Actions [39537421781](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39537421781): **332/332 PASS** (10 new). Scientific state **PARTIAL**.

## 2026-10-21 v46: Same tilted Dirac Hamiltonian under two different EB RTA preparations

The v46 PRB PHE work deliberately addresses the **symmetric, tau-linear mixed EB response** of original Eqs. (8), (10), (11), unlike the separate v41–v44 antisymmetric charge-Hall methodology. For the exact fixed-tilt surface Hamiltonian `H_B(k)=H_0(k_x,k_y+b)`, a stationary RTA whose bath targets the magnetically **equilibrated** `f(H_B)` produces large but cancelling first-order B contributions to the E response: source preparation +0.076275298 and B-on-E commutator −0.076275298, full <1e−9 under the declared model parameters. A mathematically distinct RTA whose collision target remains **frozen** at `f(H0)` instead yields nonzero B-then-E +0.110050304 and E-then-B −0.076275298, sum +0.033775006. Its low-eta tau coefficient is **close but not exactly equal** to source formal Eq. (8), with larger differences at greater tilt. Neither collision rule is yet verified as the authors' Supplemental Material prescription; the original paper PHE and Fig.2(b) amplitude issues therefore **remain OPEN**. [v46 independent physical protocols](papers/xiang-2026-quantum-geometric-map/docs/ORIGINAL_PHE_MIXED_RTA_V46.md). PR #24 merged after [343/343 PRB CI PASS](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39539072615). Scientific status **PARTIAL**.

## 2026-10-21 v48: symmetric original spin PHE robustness and original-SM access boundary

PRB [PR #27](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/27) added 10 independent model-physics tests, raising confirmed [Actions 39551301028](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39551301028) to **363/363 PASS**. The corrected v47 **symmetric** tensor, using a specifically **frozen-H0 scalar RTA collision target**, approaches signed main Eq.(11) without fitted normalization: η=.016→.001 gives ratio 1.003215→1.000012885; fitting three smallest η² points extrapolates to **1.0000002435**. Mu=±.15,±.3 odd parity, finite-T tanh smoothing, UV=1–4 robustness, and vanishing η-scaled antisymmetric remainder are independently gated.

**Source boundary:** the original source Supplemental Material could not be retrieved from accessible APS/arXiv routes, despite direct main-text verification. The paper-specific RTA denominators and bath/field-preparation protocol **are NOT source-verified**; no equivalence to the authors' exact SM is asserted. The main Eq.(8) signed branch and published Fig.2(b) ~2π SI scale remain OPEN. See [v48 detailed scientific note](papers/xiang-2026-quantum-geometric-map/docs/ORIGINAL_PHE_LIMIT_SOURCE_GATE_V48.md) and [frozen summary](papers/xiang-2026-quantum-geometric-map/results/frozen/v48_original_phe_limit_summary.json). Overall original-paper science **PARTIAL**; previous v46 sign-error numerical outputs remain RETRACTED. 

## 2026-10-21 PRB v50: Direct authentic SM finite-two-frequency Liouville source check

[PR #29](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/29) has merged after full PRB Actions [run 39558836890](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39558836890) confirmed **388/388 PASS** (12 new source-level tests). We independently solved **the actual original SM Eq.(4)/(6)/(7)/(9)** with **independent omega_a/omega_b and eta_a/eta_b**, *both E/B perturbation orderings*, complete fixed-spin-basis matrix k derivatives, and **no invented uniform bath**. Extracting the electric eta_a DC pole **after summing** both orderings converges to **SM Eq.(20) direct current**, but not the *opposite sign* printed in **SM Eq.(31)**. At beta=.125,mu=.3,T=.03,100x128 grid, eta_a=eta_b=.0005 gives eta_a full response **−0.00250621056**, versus source Eq20 **−0.00250618825** and Eq31 **+0.00250618825** (per tau*gmuB, main source helicity). The two individual orderings have substantial ratio-dependent contributions; holding omega_a nonzero eliminates the selected DC 1/eta_a pole. Source Eq20/31 signed mapping, main/SM helicity, the published Fig2(b) SI magnitude and steady-state collision physics remain OPEN. **Whole paper scientific state PARTIAL.** Full details: [v50 primary source audit](papers/xiang-2026-quantum-geometric-map/docs/AUTHENTIC_SM_TWO_COLOR_LIOUVILLE_V50.md). Source PDFs were user-uploaded and remain private; only hashes and equations are tracked.

## 2026-10-21 PRB v51: Authenticated SM Eq20/31 basis-covariant sign audit

[PR #30](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/30) merged after complete PRB [CI 39561999512](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39561999512) **400/400 PASS**, 12 new source-level tests. Independent **generic random Hermitian two-band** algebra proves literal SM Eq20 yields +tau sum Z df per E*Bbar while printed Eq31 yields -tau sum Z df, using the author's own velocity/Berry-connection relation. The main Eq10 vs SM Eq62 helicities are related by U=sigma_z, but **the physical magnetic source sigma_x and both current operators must also transform**; consistent rotation preserves both coefficients and cannot explain their *relative* sign. Original SM Eq33 completeness (sum Z=d_k <sigma_x>) independently matches eigenstate spin finite differences at ~1e-10; the momentum-integrated xy/yx identity also holds. No author sign has been silently changed. Original signed source inconsistency, Fig2(b) SI magnitude and orbital multiband/three-band physics remain **OPEN**; **scientific overall PARTIAL**. [Source audit](papers/xiang-2026-quantum-geometric-map/docs/AUTHENTIC_SM_SIGN_COVARIANCE_V51.md). Original publisher PDFs private; SHA256 tracked only.


## 2026-10-21 PRB v52: printed Fig2(b) and quoted ~43 microV share the same 2pi SI discrepancy

[PR #31](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/31) squash merged after [full PRB Actions CI 39565928035](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39565928035) **413/413 PASS** (13 new science gates). The **original published 8-page PRB page4** plots Fig2(b) T=0 plateau at roughly **+1.75** under a printed **e²/hbar** ordinate. Literal main Eq11 at beta=1/8,tau=.5ps,g=10,Bx=1T instead predicts **+11.0208361** with those units. The **separate independent publisher paragraph** estimates Hall voltage **~43 microV** for Ey=10^4 V/m,width=100 um,R~10^3 ohms, whereas literal Eq11 SI yields **268.262 microV**. Evaluating the *alternative* e²/**h** conductance quantum instead of e²/**hbar**, WITHOUT fitted factors, yields **+1.7540205** and **42.6952 microV**, matching **both** publisher approximations. Original Fig2(b) 20/40 K raster-digitized traces agree with a separately stated no-fit tanh finite-T model at this scale within ~0.02 ordinate units (source Eq11 itself is T=0). This is an **unresolved source normalization/unit-convention inconsistency**, NOT verified author plotting intent or an author erratum; the SM Eq20/31 source sign and orbital physics remain OPEN. [v52 source SI report](papers/xiang-2026-quantum-geometric-map/docs/ORIGINAL_FIG2B_ABSOLUTE_SI_V52.md). Main+SM original PDFs remain private and SHA256-pinned; scientific overall PARTIAL.


## 2026-10-21 PRB v53: original orbital quantum metric / three-band source verification

[PR #33](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/33) squash merged after **425/425 PASS** in [Actions run 39569348968](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39569348968) (12 new SM physics tests). The authentic original SM pp11–12 contains an explicit **factor 2** in Eq49's velocity–quantum-metric term, but **printed Eq52** uses factor **1**. A genuine 3D two-band periodic local quantum metric and direct Eq48(F49) BZ integration select **2**: at N=48³ source Eq48 −1.701587e−5, printed Eq52 −2.380431e−5 and factor-2 −1.687477e−5 (converging residual); missing velocity +6.92954e−6. **Eq52 explicitly drops** m≠l three-band orbital-moment/Berry interference: a real periodic three-band Hamiltonian finds omitted +0.000742732, 92.9% **of this subterm only** at mu=−1.1 (5.6% at mu=1.5); this is NOT the fraction of whole MNHE. Positive control: **both printed Eq56 orbital PHE forms agree** under BZ integration by parts, with N^-2 convergence. This is a bounded original-source algebra/approximation audit, **not an author erratum**. Spin SM Eq20/31, Fig2 SI 2pi and full untruncated orbital transport remain OPEN; overall science **PARTIAL**. [Source v53 report](papers/xiang-2026-quantum-geometric-map/docs/ORIGINAL_ORBITAL_SM_EQ49_52_56_V53.md) · [frozen JSON](papers/xiang-2026-quantum-geometric-map/results/frozen/v53_original_orbital_source_audit.json). Original PDFs remain private.

## 2026-10-21 v55 physics-first: physical orbital embedding and finite Peierls magnetic vertex

[PR #35](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/35) merged after [Actions run 39584072202](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39584072202) verified **453/453 PASS**, 15 new original-source-independent physical tests. **VERIFIED, bounded:** nondegenerate band quantum geometry with a physically specified complete finite-orbital position connection A_i and covariant energy velocity D_iH=d_iH-i[A_i,H] under non-diagonal smooth k-dependent orbital W(k); and finite **open** real-space graph Peierls H(B), dH/dB with supplied site positions/straight paths, symmetric and Landau magnetic gauge spectra and Wilson flux checks. **Important physical contact:** B-dependent gauge G requires derivative terms in magnetic vertex transformation, not merely G M G†; naive residual 0.522 versus complete ~2.8e-17 in the model. Real Wannier inter-orbital dipoles, non-Abelian degeneracies, primitive translational B field response and complete orbital E×B transport remain **UNRESOLVED**; no fabricated experimental observable. See [v55 scientific details](papers/xiang-2026-quantum-geometric-map/docs/PHYSICAL_EMBEDDING_MAGNETIC_INTERFACE_V55.md), [physics capability registry](physics/README.md), and [frozen reference](papers/xiang-2026-quantum-geometric-map/results/frozen/v55_embedding_magnetic_reference.json). Whole original PRB science **PARTIAL**.
