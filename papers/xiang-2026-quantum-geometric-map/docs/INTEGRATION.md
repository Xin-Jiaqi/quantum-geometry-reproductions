> **v49 direct source update:** Authentic published PRB main and 16-page Supplemental Material now directly inspected and SHA256-pinned. Earlier v47–v48 caveats saying the original SM is unavailable are historical and superseded. Do not conflate SM's selected electric (1/eta_a) pole turned into tau with the uniform finite eta all-sector frozen-H0 RTA. An independently calculated literal Eq20 current and printed Eq31/32 curvature current have opposite signs; main Eq10 and SM Eq62 flip the spin helicity. See [source-level v49 audit](docs/ORIGINAL_SM_SOURCE_AUDIT_V49.md). Overall scientific status PARTIAL.

## v47 critical erratum: v46 source sign and PHE symmetry

Before using v46, read [v47 corrected physics](docs/ORIGINAL_PHE_SYMMETRIC_SOURCE_CORRECTION_V47.md). The analytic y-source in v46 had the sign of d_y(kx/r) wrong, so former v46 frozen outputs and near-Eq8 XY-only comparison are **RETRACTED**. Corrected source is finite-difference checked. `qgmap.dirac_source_tensor_v47.tensor` evaluates both XY/YX orders and their **symmetric** part, compares eta*symmetric (NOT XY alone) to the signed original main Eq.(11) and separates an antisymmetric remainder. `independent_mixed_B_derivative` recalculates the finite-b stationary density through an independently diagonalized Liouville superoperator for both electric directions. The successful five-beta numerical Eq.(11) comparison is bounded to the **researcher-defined frozen-H0 collision target**. Exact original Supplemental Material collision input and Fig.2(b) absolute units remain unverified; entire paper PARTIAL.

# Callable physical interfaces and fail-closed material input contract

Import from `qgmap`, **not** from the other two independent modules. This module deliberately avoids a speculative universal NMEE/IGMC/magnetotransport solver.

| Request | Available API | Units, boundary |
| --- | --- | --- |
| Dirac `H(k)`, bands | `qgmap.dirac.hamiltonian`, `energies` | energy eV, k Å^-1, `vF_eVA` explicitly supplied |
| Zeeman geometry | `local_zbc_analytic`, `local_zbc_eigensolver` | no k=0 nondegenerate formula |
| Printed Eq11/SM63 | `eq11_kernel`, `eq63_literal_integral` | separate, sign-disagreeing protocols |
| Thermal Dirac | `continuum_temperature_kernel` | ideal unbounded clean Dirac; not a full transport calculation |
| Rashba OHE | `rashba.zero_temperature`, `direct_fermi_ring_sum` | e=ħ=1, constant τ, T=0, outside band-edge/crossing |
| Classical OHE | `orbital.ohe_parabolic_2d` | trivial one-band model, not multiband quantum geometry |
| Two-color density matrix | `liouville.two_color_response` | finite matrices, explicit `V_E,V_B,η_E,η_B,ω_E,ω_B` |
| Magnetic groups | `symmetry.table_i_acceptance` | 11 explicit generator settings, not all 122 MPG |

## Required for material transfer

Do NOT infer a quantitative magnetic response from `wannier90_hr.dat` alone. Stage-1 read-only contracts for other workflows can help validate structural dimensions, SOC and checksum completeness, but physical inputs also include:

- lattice and reciprocal units, dimensionality and occupancy;
- same-gauge Hamiltonian H(k), position/Berry connection `A_a(k)`, velocity `D_a H`;
- projected spin and/or orbital moment `S_i(k)`, `m_i(k)`, and covariant derivatives needed by the target response;
- specified field-frequency/collision/thermalization protocol and near-degenerate subspace treatment.

**Without these, real-material magnetotransport is BLOCKED.** No mock HR-only real-material result is claimed. The protocol distinction between `pre-equilibrated-B` and `adiabatic-B` is not an implementation detail; it can change observables.

## Cross-paper reuse

NMEE, IGMC and PRB have different perturbations and observables. Reuse common linear-algebra, integration, tensor symmetry and provenance helpers only after agreeing on observables, units, gauge and limiting procedure. Never copy a response coefficient across papers because the API name sounds similar.

## Newly migrated advanced physical methods (v31–v33)

| Calculable method | Python function | Independent physics gate |
| --- | --- | --- |
| Four-band intrinsic spin MNHE, SM (16–19) | `qgmap.spin_multiband_v5.four_band_audit(N=23)` | sum of four density-matrix currents vs independently evaluated SM Eq.(29) |
| 2D orbital PHE, SM Eq.(56) | `qgmap.orbital_phe_v4.orbital_phe_eq56(N=61)` | literal first/second lines and independent finite-field Boltzmann derivative |
| Strict 2D orbital MNHE SBC | `qgmap.orbital_mnhe_v5.orbital_mnhe_2d(N=61)` | orbital-moment Fermi-surface term vs squared Berry curvature |
| 3D two-band orbital MNHE, SM Eq.(52) | `qgmap.orbital_mnhe_v5.orbital_mnhe_3d(N=17)` | positional-shift Eq.(48)/(49) versus literal Eq.(52), factor-two issue retained |

Each function returns a `protocol` label and normalized model-unit results. The synthetic spin and orbital protocols must **not** be combined into a physical SI coefficient without separately defining the perturbation, magnetic moment units, collision model, dimensionality and band/occupancy convention.

Unified JSON commands (run from this paper module after `pip install -e .`):

```bash
python scripts/reproduce_advanced.py --case spin-four-band
python scripts/reproduce_advanced.py --case orbital-phe
python scripts/reproduce_advanced.py --case orbital-mnhe-2d
python scripts/reproduce_advanced.py --case orbital-mnhe-3d --output results/local-orbital3d.json
```

### Disorder, Kubo and Ward methods (v34)

| Method | API | What is tested |
| --- | --- | --- |
| SCBA + ladder | `qgmap.microscopic_ward_v11.scba`, `ladder_vertex` | same scalar disorder variance in self-energy and Pauli-basis Bethe-Salpeter vertex |
| RA control | `qgmap.warped_disorder_v17.ra_observable` | RA bubble and exact momentum-shift control in an independent periodic model |
| Spectral Kubo | `qgmap.spectral_kubo_v18.spectral_bubble` | four (RA,AR,RR,AA) channels, each with its own noncrossing ladder |
| Finite-q, finite-ω Ward | `qgmap.finite_qomega_ward_v19.finite_ward` | matrix resolvent identity with separately solved energies and finite-difference current insertion |

Run `python scripts/reproduce_advanced.py --case scba-ra`, `--case kubo-spectral`, or `--case finite-qomega-ward`. See [v34 method scope and limits](docs/SCBA_KUBO_WARD_MIGRATION_V34.md).

A real Wannier input pipeline is still **BLOCKED** until consistent magnetic/position operators and response limits are provided. The finite-eta Ward identity is a resolvent identity, not a proof of physical dc charge conservation.


## Exact finite-system charge continuity and regulator limits (v35)

- `qgmap.real_frequency_continuity_v20.kinematic_ward`: **exact** finite-disorder longitudinal Kubo Ward identity, including the equal-time density-current contact term. The q=0 dynamical density response and equilibrium compressibility are intentionally separate.
- `qgmap.limit_order_v20.scan_limit`: **SCBA** finite-regulator limit, contrasting nonzero probe ((q,\omega)) against static uniform origin. This is not the physical charge Ward identity at finite eta.

```bash
python scripts/reproduce_advanced.py --case exact-kubo-v20
python scripts/reproduce_advanced.py --case limit-order-v20
```

See [v35 scientific gate and numbers](docs/EXACT_KUBO_LIMITS_MIGRATION_V35.md). These methods do not upgrade the paper's static PHE or full Kubo–Středa status.

## v36: independently callable second-order EB protocols

| Workflow | API | Protocol |
| --- | --- | --- |
| Finite-q E / static B | `qgmap.electromagnetic_eb_v21`; `--case eb-finite-q-v21` | disordered finite system, scalar-electric and Zeeman, contact terms and thermal equilibrium density derivative |
| Finite-q two-color Liouville | `qgmap.two_color_exact_v22`; `--case eb-two-color-v22` | dynamic E/B, two separate action orderings, independent Sylvester and time-domain validation |
| Uniform E length gauge | `qgmap.homogeneous_length_gauge_v23`; `--case eb-uniform-v23` | periodic synthetic two-band E/B, FFT covariant derivative and SM Eq.20 pole |

The finite-q scalar-electric models and homogeneous-E length-gauge model must **not** be merged as if they used the same experimental electric-field protocol. All results retain model units; neither static RTA `tau` nor (B) in tesla is automatically inferred. See [v36](docs/SCIENCE_MIGRATION_V36.md).

## v37: continuum tilted Dirac, physical SI and published-figure checks

| Purpose | Independent API | Science boundary |
| --- | --- | --- |
| Original main Eq. (10) Dirac two-color kernel | `qgmap.continuum_ti_v24`, `--case continuum-dirac-v24` | \(\bar B_x\) units, finite UV cutoff and adiabatic regulators, NOT transport DC |
| Eq. (11), Fig. 2(b), voltage dimensional audit | `qgmap.si_fig2_v25`, `--case fig2-si-v25` | \(\hbar\) vs \(h\) hypothesis separately named |
| SI current from full two-band Liouville | `qgmap.physical_si_liouville_v26`, `--case si-liouville-v26` | [A/m], [S], tesla and J units, still adiabatic response |
| Fig. 2(b) fitted curves | `qgmap.fig2_source_audit_v27`, `--case fig2-curve-v27` | Approximate publisher raster pixels, source raw data not available |

See [v37 provenance and physical limitations](docs/CONTINUUM_SI_FIG2_MIGRATION_V37.md). No original PDF or figure is included. Do not silently replace Eq. (11) with its \(1/(2\pi)\) diagnostic branch.


### v38 reusable clean Hall Fermi-sea calculation

`qgmap.streda_clean_sea_v38.bulk_sea(n,mass,mu,eta)` runs a clean nondegenerate QWZ Bloch model and returns `sea_sigma_over_e2h` and the independent sharp occupied-band Berry-Chern value. `fhs_chern(H)` supplies a separate velocity-free topological check; `direct_energy_sea` checks the matrix-energy integral against the analytic primitive. `metallic_pole_audit` **only diagnoses** the sea integral's rational pole boundary and does not compute a complete physical Středa Fermi-surface term. See [full science and formulas](docs/V38_CLEAN_STREDA_AUDIT.md). The original PRB mixed \(EB\) PHE and real Wannier are not in this routine.

## v39: QWZ SCBA with energy-dependent self-energy and exact finite disorder

- `qgmap.scba_streda_v39.derivative_audit`: SCBA self-energy in the independent QWZ model, `Gamma_0^{RR}=I-dSigma^R/dE`, full four-Pauli-basis RR density-ladder solve, and finite-difference derivatives of both `Sigma` and `G`.
- `qgmap.scba_streda_v39.sea_energy_scan`: independent real-energy quadrature of the **bare-current** antisymmetric sea kernel using either the consistent self-energy derivative or a deliberately naive constant-`Sigma` derivative, with both branches named.
- `qgmap.scba_streda_v39.exact_finite_disorder_sea`: exactly diagonalized finite torus with sampled site disorder and complete single-particle eigenstates, **not** SCBA.

Use `python scripts/reproduce_advanced.py --case scba-energy-sea-v39` or the detailed `scripts/reproduce_scba_streda_v39.py`; both return `PARTIAL` science status. See [v39 standalone protocol](docs/SCBA_ENERGY_DERIVATIVE_V39.md). General metallic Hall, current vertices, Středa I and nonlinear EB remain OPEN.


## v40 matched noncrossing current vertex

The exact four-Pauli-basis solver `qgmap.current_vertex_v40.pauli_vertex` returns the current vertex and its Bethe–Salpeter residual for a named RA, AR, RR or AA branch with the SAME disorder variance as SCBA. `current_audit` also checks the RA density Ward identity and the RR periodic-momentum Ward limit, and enforces that the four-channel symmetric spectral Hall antisymmetry vanishes. `exact_disorder_ra_torus` provides an independent full-eigenspectrum **bare** RA benchmark. Use `--case scba-current-vertex-v40`; see [method limitations](docs/CURRENT_VERTEX_V40.md). Neither the isolated RA Hall-odd correlator nor the dressed symmetric spectral 4-channel sum is the complete physical disorder Hall conductivity.

## v41 exact disordered linear Hall reference

`qgmap.exact_hall_v41.torus`, `hall_lehmann`, `hall_retarded_crosscheck`, `twist_chern` and `audit` accept a declared finite-size QWZ torus with **frozen site disorder**, Peierls derivatives and zero-T occupancy. `scripts/reproduce_exact_hall_v41.py` returns a protocol-labelled JSON object with `scientific_status=PARTIAL`. For named gapped finite tori, independently computed occupied-subspace FHS Chern and flux-averaged exact Hall agree. There is no claim of separate Středa I/II, SCBA conserving vertex completion, material-ready operators or original EB DC reproduction. Details: [v41](docs/EXACT_HALL_V41.md).

## v42: Exact finite-sample charge Hall Středa split

Callable: `qgmap.exact_streda_v42.spectral_streda(H,Jx,Jy,mu,eta)` for surface I(a) and sea II plus full total; `direct_surface_trace`, `direct_sea_quadrature` are intentionally independent checks; `streda_audit`, `ensemble_audit` provide JSON-ready finite-sample audits. The external positive `eta` is a Green regulator, not a conserving collision operator. The method reuses `qgmap.exact_hall_v41.torus`; its exact Hall sum agrees with the independent v41 Lehmann as `eta -> 0` at fixed finite gapped samples. Do **not** use I or II alone as physical conductivity or treat this as a disorder-averaged SCBA solver or nonlinear EB PHE. See [v42](docs/EXACT_STREDA_V42.md).

## v43 same-disorder Hall I(a) (not a total Hall coefficient)

`qgmap.scba_surface_v43.scba_surface(mesh,mass,W,mu,eta)` returns the matched v39/v40 noncrossing RA **surface** I(a) with one current ladder and checks left/right reciprocity and RA resolvent Ward. `exact_ensemble_surface(L,mass,W,mu,eta,seeds,ntheta)` calculates the **same surface** using independent v42 exact dense Green traces and sample-level twist integration. `comparison` returns seed SEM and remaining model discrepancy; only independent disorder seeds contribute statistical degrees of freedom. The v43 module deliberately has **no total Středa I+II API**. Conserving averaged sea II, crossed diagrams, physical metallic DC and original nonlinear EB PHE are OPEN. See [v43 protocol](docs/SCBA_SURFACE_V43.md).

## v44: Provisional SCBA Hall sea II and candidate I(a)+II

`qgmap.scba_sea_v44.sea_kernel` solves the same-W RR density Bethe–Salpeter derivative insertion for a chosen energy, with optional independent centered finite differences of Sigma and G. `sea_scan` integrates the bare-spatial-current Středa-II *candidate* over real energy, solving a fresh self-energy at each node; default frozen regulator eta=.12, 420 Gauss nodes, lower=-12. `candidate_total` combines that sea with **one** v43 RA current ladder on I(a) and explicitly labels the result **provisional**, with the complete conserving electromagnetic response still OPEN. `exact_ensemble_total` compares to the v42 full I+II finite-sample model with sample-level uncertainties; `clean_bloch_pair` avoids the W=0 finite-torus degeneracy and independently checks the exact analytic clean pole primitive. The public JSON always uses `scientific_status=PARTIAL`. See [v44](docs/SCBA_SEA_V44.md).

## v45 Original RTA symmetric spin PHE vs fixed-tilt equilibrium B field

`qgmap.dirac_phe_protocol_v45` exposes `hamiltonian`, `zeeman_curvature`, `angular_eq8_geometry`, `printed_eq11_geometry`, `equilibrium_drude_xy`, `drude_delta_derivative_terms`, `shift_identity` and `audit`. It independently audits original Eq.(8)/(11) and the exact pure ky shift of the fixed-tilt Dirac model under Bx, with separate fully equilibrated static-B Drude populations. The original theoretical scope is phenomenological RTA and **symmetric EB PHE**; v41-v44 conserving SCBA-Streda **antisymmetric linear Hall** is an independent extension. Different magnetic preparation/UV/collision protocols must never be silently merged. See [v45](docs/ORIGINAL_PHE_PROTOCOL_V45.md). Fig2(b) ~2pi scaling and stationary mixed-field response remain OPEN.

## v46 original symmetric mixed EB under two selected scalar RTA targets

`dirac_mixed_rta_v46.Spec`, `current`, `decompose`, `unprepared_two_field`, `audit` build the same H_B and complete two-band source matrix under two distinct collision targets. `decompose` checks the source-preparation and magnetic insertion against an independent full B finite difference, obeying the Bx-to-ky translation identity in the B-prepared equilibrium target. `unprepared_two_field` freezes the collision target at f(H0), includes B-then-E and E-then-B orderings and returns a generally nonzero response. It is a protocol **comparison**, not an implementation claimed to match the author's unspecified full SM collision operator; the unequal beta dependence of the tau-linear coefficient proves these are not interchangeable. See [v46](docs/ORIGINAL_PHE_MIXED_RTA_V46.md). Whole paper PARTIAL.

## v48: Controlled original Eq11 PHE model limits (NOT source-SM certified)

`qgmap.phe_limit_audit_v48` provides `eta_limit`, `chemical_potential_temperature`, `cutoff_sweep` and `full_audit`. These repeatedly call **the corrected v47 symmetric (xy+yx)/2** frozen-H0 RTA response, never v46 retracted XY-only values. The eta scan fits 3/4/5 smallest η² points and checks symmetric Eq11 ratio intercept and vanishing η-rescaled antisymmetric remainder. The mu/T protocol verifies odd electron-hole response and tanh smoothing; the cutoff protocol checks a finite range of continuum UV values. Output explicitly labels author original SM collision source UNAVAILABLE, and all results remain `scientific_status=PARTIAL`. [v48 details](docs/ORIGINAL_PHE_LIMIT_SOURCE_GATE_V48.md). Do not infer the source-specific relaxation, measured Hall voltage or Fig2 normalization.

## v50 independently solved original SM two-colour Liouville

`qgmap.sm_two_color_liouville_v50.Tone` defines omega_a,omega_b,eta_a,eta_b and k derivative; `two_color` inverts L(z)=z+i[H0,.] in the physical-spin basis for both E-before-B and B-before-E and returns complex physical B-current per complex Fourier amplitude, with no bath collision. `pole_only` isolates the original electric-field intraband 1/eta_a residue; `run_limit` scans eta_a while preserving eta_b. First-order eigenbasis source Eq6/7, different regulator paths, original SM Eq20 sign and printed SM Eq31 opposite sign, finite omega electric pole disappearance and complex-conjugation are independently tested. The original Eq9 **now** supports literal Eq20; the source relative sign and Fig2 remain unresolved. See [v50 note](docs/AUTHENTIC_SM_TWO_COLOR_LIOUVILLE_V50.md). All science PARTIAL.

## v51: Primary source sign and basis covariance audit

`qgmap.sm_sign_covariance_v51.source_sign_identity` verifies the SM Eq20/printed Eq31 opposite signs on arbitrary nondegenerate Hermitian two-band data, independently of Dirac topology, and handles independent band U(1) phases. `integrated_table` returns all four options (main/SM helicity, ±sigma_x magnetic source) in **per E Bbar** units; these show exact physical-spin representation covariance when H and the magnetic source rotate together, with the relative Eq20/31 sign unchanged. `finite_difference_spin_ward` validates the separate SM Eq33 identity via direct spin expectation derivatives, and the XY/YX equality is tested after integration (not falsely pointwise). An explicit test cross-checks the original v50 two-frequency Liouville DC sign. Read [v51](docs/AUTHENTIC_SM_SIGN_COVARIANCE_V51.md) and frozen `results/frozen/v51_sm_sign_covariance_summary.json`. Overall paper PARTIAL, original Eq20→Eq31 source-level sign OPEN.

## v52 original published Fig2(b) and Hall voltage SI audit

`fig2_si_v52.Fig2Config`, `eq11_f`, `absolute_fig2_audit`, `finite_temperature_shape_audit` give a source-pinned **literal printed Eq11/e²-hbar SI** prediction versus an explicit **alternative e²/h** unit convention, without fitting any 2pi or signs. At published beta=.125,tau=.5ps,g=10,Bx=1T, literal normalized Fig2b amplitude is 11.0208361 and Hall voltage 268.262 microV; alternative yields 1.7540205 and 42.6952 microV, agreeing with publisher T0 plateau about 1.75 and prose about 43 microV. CODATA joule and eV routes are independently tested; 8 source-page pixel readings of T20/40 K match a **declared model-only** tanh finite-T continuation within ~.02 axis units; source Eq11 was at T=0. [v52 source SI report](docs/ORIGINAL_FIG2B_ABSOLUTE_SI_V52.md). Figure-generation code and the true author normalization convention are unavailable; original signed SM Eq20/31 source and orbital multiband remain OPEN. Whole-paper PARTIAL.

## v53: authentic orbital SM quantum metric, three-band curvature and Eq56

`qgmap.orbital_sm_v53.metric_eq49_to_eq52` constructs a nondegenerate **physical 3D two-band** metric and compares direct source Eq48(F49) with printed Eq52 and the algebraically factor-2 velocity counterpart. The difference is controlled by actual band-sum velocity and BZ finite-difference convergence; flat d0 recovers identical forms. `orbital_moment_three_band` diagonalizes a genuine 2D periodic **three-band** Bloch Hamiltonian, differentiates it analytically for velocities and compares the **full m,l interband product** in Eq48 using Eq50–51 against the **m=l subterm retained in Eq52**. Crucially the omitted fraction refers **only to this curvature-orbital-moment subterm**, not the total MNHE; gauge phases and Fermi energy dependence are tested. `orbital_phe_eq56` confirms the **two printed Eq56 forms** agree after periodic BZ integration by parts (second-order finite-difference convergence). Exact original source PDF SHA pinned, other papers untouched. [v53 evidence](docs/ORIGINAL_ORBITAL_SM_EQ49_52_56_V53.md), frozen JSON. Overall paper **PARTIAL**, untruncated orbital transport still OPEN.
