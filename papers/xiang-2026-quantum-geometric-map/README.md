> **v49 ORIGINAL PRIMARY SOURCE AVAILABLE:** Both user-provided PDFs have now been SHA256-verified: published PRB (8 pages) and authentic Supplemental Material (16 pages). Prior v45–v48 statements that the SM was unavailable are superseded. The original SM Eq15 inserts phenomenological tau at a selected electric-field adiabatic pole, while Eq20 and printed Eq31 have a literal **overall sign mismatch** under the velocity-r identity quoted in Eq31. Main Eq10 and SM Eq62 use opposite spin helicities for fixed physical sigma_x. Independent both-direction eigenbasis calculations keep all branches, no signs fitted. See [v49 source-certified scientific audit](docs/ORIGINAL_SM_SOURCE_AUDIT_V49.md). The article is still **PARTIAL**, with the source-level sign bridge OPEN.

> **v47 SCIENTIFIC CORRECTION:** The original v46 electric-y equilibrium source contained a wrong Pauli derivative sign. The v46 frozen outputs and claim that unsymmetrized xy approximately matches Eq.(8) are **RETRACTED**. Code and tests are corrected. The properly symmetrized frozen-target response `(xy+yx)/2` matches signed main Eq.(11) across beta=.025–.6 at eta=.002, with relative residual about 5×10^-5; the author's Supplemental Material relaxation operator is still **not independently checked**. See [v47 correction](docs/ORIGINAL_PHE_SYMMETRIC_SOURCE_CORRECTION_V47.md). Whole paper scientific status **PARTIAL**.

# Xiang et al. (2026): Quantum Geometric Map of Magnetotransport

**Reference:** L. Xiang, J. Jia, F. Xu and J. Wang, *Physical Review B* **113**, L201406 (2026). [DOI: 10.1103/kd12-4ycj](https://doi.org/10.1103/kd12-4ycj).

**Scientific status: PARTIAL.** Independent reproducible calculations; not an official author implementation.

## GitHub scope versus full offline research

This repository module publishes **directly importable, tested core capabilities** from the v30 research program, not the entire 421-file v30 checkpoint. The complete offline archive (including historical 30-round numerical audits, detailed response branches, and figures) remains separate. The CI here certifies only the particular code checked in, not every historical research result.

Supported:

- Tilted Dirac band spectrum, Hamiltonian conventions (main Eq. 10 vs SM Eq. 62), Zeeman Berry curvature via analytic and gauge-invariant matrix elements, Eq. (11) and literal SM Eq. (63), thermal-step and finite-field Drude diagnostics.
- Rashba orbital Hall response Eq. (9) at `T=0`, with distinct one-ring/two-ring occupancy sectors and independent Fermi-ring calculation.
- Classical orbital Hall effect (spinless parabolic band, finite temperature).
- Finite-dimensional **two-color Liouville resolvent** with explicit `E→B` and `B→E` orderings, independent frequencies, broadenings, and a **separately labelled** equilibrium population correction. This is a method module, not the complete published driven-DC formula.
- SM Table I 11 representative magnetic-group tensor spaces using Neumann operators.
- SCBA self-energy and a noncrossing Bethe-Salpeter ladder using the same scalar-disorder variance.
- Symmetric four-channel spectral Kubo-Greenwood response with separate RA, AR, RR, and AA ladders.
- Finite-(q, omega) resolvent Ward identity with the exact finite-difference longitudinal insertion.
- SCBA self-energy and noncrossing Bethe-Salpeter ladder in an independent periodic two-band regulator.
- Symmetric four-channel spectral Kubo-Greenwood response (`RA+AR-RR-AA`) with a distinct ladder in each channel.
- Finite-((q,\omega)) SCBA resolvent Ward identity using the exact finite-difference longitudinal insertion `[H(k+q/2)-H(k-q/2)]/q`.

## Start

```bash
python -m pip install -e .
python scripts/verify_release.py
```

```python
from qgmap.dirac import Parameters, eq11_kernel, local_zbc_eigensolver
from qgmap.rashba import zero_temperature, direct_fermi_ring_sum
from qgmap.liouville import two_color_response
print(eq11_kernel(1/8))
print(zero_temperature(1.0, 0.4, -0.04))
```

The source and benchmarks use explicit model units; they are not automatically a SI transport or materials solver.

## Scientific limits

- **Fig. 2(b):** temperature dependence verified in the original research; original plotted plateau is about `2π` below the literal SI Eq.(11) result. Do not fix by arbitrary scaling.
- **Steady-state PHE:** adiabatic Liouville pole and equilibrium/static full response differ under some preparation/relaxation protocols.
- **Orbital MNHE/PHE:** SM Eq.(52)/(56) discrepancies remain open. General multiband quantum geometry, omitted three-band terms, SCBA/Středa full transport, and first-principles calculations are **not** established by this published core.
- **Rashba OHE:** model calculation is available, QMQ attribution as a whole-paper statement remains open.
- Degeneracies and response limits must not be treated as covered by nondegenerate expressions.

See [scientific acceptance](docs/SCIENTIFIC_ACCEPTANCE.md), [integration contract](docs/INTEGRATION.md), and [provenance](docs/PROVENANCE.md). Do not redistribute publisher PDF, original figure images, or private logs.

## Migrated advanced-method subset (v31)

The public core now includes independently callable spin-MNHE, orbital PHE/MNHE, SCBA/ladder, four-channel spectral Kubo, and finite-(q, omega) resolvent-Ward method subsets. See the [v31 spin](docs/SPIN_MULTIBAND_MIGRATION_V31.md), [v32 orbital](docs/ORBITAL_MIGRATION_V32.md), and [v34 SCBA/Kubo/Ward](docs/SCBA_KUBO_WARD_MIGRATION_V34.md) protocols. These are model-level methods, not a general degeneracy-robust real-material solver.

### Orbital advanced-method subset (v32)

The migrated `qgmap.orbital_phe_v4` and `qgmap.orbital_mnhe_v5` modules independently check 2D/3D model responses and preserve the literal SM Eq. (52)/(56) and finite-field comparison branches. [Methods and limitations](docs/ORBITAL_MIGRATION_V32.md). They **do not** implement the paper's omitted general multiband three-band orbital terms.


### Disorder/Kubo/Ward subset (v34)

The original v11/v17/v18/v19 SCBA, ladder, four-channel spectral Kubo, and finite-((q,\omega)) resolvent-Ward modules are now migrated with 48 original regression tests. See [the v34 protocol and boundaries](docs/SCBA_KUBO_WARD_MIGRATION_V34.md). This remains a method subset: it is not full Kubo–Středa, physical finite-(eta) continuity, or the paper's driven-DC PHE.


### Exact Kubo and limit-order subset (v35)

The original v20 finite-disorder full-Lehmann Kubo continuity/contact-term benchmark and the separate SCBA external-regulator limit-order scan are now migrated with their 18 original tests. Two new JSON cases are `exact-kubo-v20` and `limit-order-v20`. See [v35 methods, reference numbers and physical limits](docs/EXACT_KUBO_LIMITS_MIGRATION_V35.md). Neither constitutes a complete transverse or mixed-field DC response.

### Two-color and length-gauge EB migration (v36)

Three more independent original research modules are published with 66 historical tests: `electromagnetic_eb_v21` (finite-q scalar E and static-B Liouville/Kubo derivative), `two_color_exact_v22` (both dynamic E/B switch orderings), and `homogeneous_length_gauge_v23` (periodic uniform E, gauge-covariant derivative and SM Eq.20 pole). Use `python scripts/reproduce_advanced.py --case eb-finite-q-v21`, `--case eb-two-color-v22`, and `--case eb-uniform-v23`. See [v36 migration scope](docs/SCIENCE_MIGRATION_V36.md) and [frozen provenance](results/frozen/v36_source_provenance.json). The paper's physical steady-state DC response remains OPEN.

### Continuum Dirac and physical SI completion (v37)

The original v24–v27 continuum tilted Dirac, restored-SI Liouville, Eq. (11)/Fig. 2(b) dimensional audit and approximate publication-curve digitization are now migrated with **62 original model tests**, numeric snapshots and exact Git blob SHA provenance. The main Eq. (10) helicity, Supplemental Material helicity and \(\bar B_x=-g\mu_B B_x\) sign conventions remain explicitly distinct. No publisher PDF or original figure is mirrored.

```bash
python scripts/reproduce_advanced.py --case continuum-dirac-v24
python scripts/reproduce_advanced.py --case fig2-si-v25
python scripts/reproduce_advanced.py --case si-liouville-v26
python scripts/reproduce_advanced.py --case fig2-curve-v27
```

The literal SI formula gives an ordinate near **11.02**, while the published figure's approximate digitized plateau is **1.75**. The \(1/(2\pi)\) branch fits the ordinate but remains a **diagnostic hypothesis, not a proven correction or confirmed source of the difference**. Source curve coordinates are approximate PDF pixel readings, not author raw data. For full details see [v37 migration and verification](docs/CONTINUUM_SI_FIG2_MIGRATION_V37.md).

## v38: clean antisymmetric Fermi-sea Green-function benchmark

In a *separate, intentionally topological* QWZ two-band model, `qgmap.streda_clean_sea_v38` evaluates the antisymmetric retarded–advanced Fermi-sea trace via both an analytic two-pole primitive and independent direct real-energy quadrature. A gauge-invariant link Chern invariant and clean Berry-curvature integration verify the insulating limit. Unlike v18's **symmetric dissipative** spectral kernel, this is an **antisymmetric Hall** test. At finite doping, the module explicitly exposes the rational boundary-pole contribution and avoids identifying it with the physically complete disorder-dressed Středa surface term. [Derivation and 12 test gates](docs/V38_CLEAN_STREDA_AUDIT.md).

```bash
python scripts/reproduce_streda_v38.py --mesh 61 --mass -1 --mu 0 --eta .02
python scripts/reproduce_streda_v38.py --mesh 101 --mass -1 --mu 1.1 --eta .02
python scripts/reproduce_advanced.py --case clean-streda-sea-v38
```

All resulting conductances are dimensionless `sigma/(e²/h)`, not an SI prediction of the PRB steady-state \(EB\) PHE. Full finite-energy disorder self-energy, conserving vertices and Kubo–Středa remain OPEN.

### v39: energy-dependent SCBA Green function and RR density Ward derivative

The new `qgmap.scba_streda_v39` method solves `Sigma(E)=W<GR(E)>` together with the **same-`W` RR density Bethe–Salpeter vertex** `Gamma_0=I-dSigma/dE`, verified by finite differences in a QWZ band model. It integrates the explicitly defined antisymmetric bare-current Green-function sea kernel both with and without the selfenergy derivative, and separately checks an **exact finite static-disorder Lehmann sum**. These are intentionally different algorithms, not assumed to produce identical full Hall conductivities.

```bash
python scripts/reproduce_scba_streda_v39.py --mesh 21 --W .12 --eta .12
python scripts/reproduce_advanced.py --case scba-energy-sea-v39
```

14 new tests plus frozen JSON checks enforce the method boundary. The v39 module is **not** a complete disordered Kubo–Středa conductivity, an electromagnetic transverse Ward identity, or the original paper's nonlinear \(EB\) DC PHE. See [v39 numerical audit](docs/SCBA_ENERGY_DERIVATIVE_V39.md).

### v40: Current-vertex Bethe–Salpeter and Hall zero controls

`qgmap.current_vertex_v40` solves the four Pauli-basis noncrossing current vertices separately in RA, AR, RR and AA channels, using **the same scalar-disorder variance** as `scba_streda_v39`. In a doped QWZ model the RA antisymmetric correlator changes when the ladder is included, but the **four-channel symmetric spectral combination has zero antisymmetric Hall part**. The RR current correction converges to zero on a full periodic BZ (momentum Ward identity), and the RA density ladder matches the SCBA self-energy difference.

```bash
python scripts/reproduce_current_vertex_v40.py --mesh 21 --energy 1.1 --W .12 --eta .12
python scripts/reproduce_advanced.py --case scba-current-vertex-v40
```

The independent exact disordered finite torus is a **bare RA control**, not an ensemble-converged SCBA transport coefficient. [v40 model proof and benchmarks](docs/CURRENT_VERTEX_V40.md). Full Hall Kubo–Bastin/Středa surface and sea, transverse gauge contact terms and the PRB's second-order physical EB PHE remain **OPEN**.

### Independent finite-twist RR Ward control (v40b)

The `qgmap.rr_spatial_ward_v40` module independently computes `Sigma^R(A_x, A_y)` at small positive and negative **uniform vector potential** and verifies that its finite difference matches the RR current Bethe–Salpeter correction. The same SCBA propagators are also checked against `dG^R/dA=G^R Gamma^{RR} G^R`. This **independent twist test** is different from the RA transport ladder and the v39 energy/density Ward; both momentum directions and the zero-disorder limit are tested.

```bash
python scripts/reproduce_rr_ward_v40.py --mesh 21 --W .12 --energy 1.1 --eta .12
```

Details and limitations: [v40b RR spatial Ward](docs/RR_SPATIAL_WARD_TWIST_V40.md). Passing the two spatial derivative tests **does not** establish a complete disordered Hall/Středa conductivity or the physical mixed-EB DC protocol.

### v41: Exact finite-torus Hall reference (not full Kubo–Středa)

`qgmap.exact_hall_v41` computes the **full exact finite-system single-particle linear Hall Lehmann response** for a fixed clean or scalar-disordered QWZ torus. Its twist average is independently compared with the FHS occupied-subspace Chern number, including large-gauge seams; a separate retarded correlator and the QWZ-specific zero mixed diamagnetic contact are tested. Neither the finite-sample exact Hall nor v40's RA ladder is a complete disorder-averaged Středa surface/sea calculation, and the original second-order `EB` PHE is still **OPEN**.

```bash
python -m unittest tests.test_exact_hall_v41 -v
python scripts/reproduce_exact_hall_v41.py --L 4 --mass -1 --W .12 --seed 1947 --ntheta 7
```

See [v41 exact Hall protocol](docs/EXACT_HALL_V41.md).

### v42: Exact finite-disorder *linear* Hall Středa I+II

`qgmap.exact_streda_v42` evaluates the **standard Smrčka–Středa** antisymmetric I(a)+I(b)+II terms for the **same frozen scalar-disorder QWZ torus** and current operators as v41. It checks a direct matrix-trace I(a), a separate real-energy dense Green-function integral II, the analytic pair primitive, the v38 clean Bloch sea and convergence of the **complete sum** to the v41 independent Lehmann DC Hall. The Hall-odd I(b) is zero for these Hermitian charge-current operators. The individual surface and sea terms strongly cancel at finite `eta` and are not uniquely separable measurable mechanisms.

```bash
python -m unittest tests.test_exact_streda_v42 -v
python scripts/reproduce_exact_streda_v42.py --L 4 --mass -1 --W .12 --mu 1.1
```

[Full v42 derivation and reproducibility](docs/EXACT_STREDA_V42.md). **The self-consistent disorder-averaged conserving Středa I/II, physical thermodynamic DC, and Xiang original mixed second-order EB PHE are still OPEN; whole-paper state PARTIAL.**

### v43: Matched SCBA Hall I(a) versus exact sample averaging

The **linear charge Hall Středa I(a) RA surface channel only** can now be computed using the same scalar disorder variance, Green regulator, and QWZ currents on both sides: `qgmap.scba_surface_v43.scba_surface` reuses v39 SCBA and v40 **single-leg** current ladder, while `exact_ensemble_surface` independently averages the v42 **dense real-space** surface trace over disorder seeds and twists. `comparison` returns bare/dressed/exact values and **seed-based standard errors**. At W=.12, L=4, 24 seeds and 7x7 twists (μ=1.1, η=.12), bare SCBA is -1.27807, single ladder -1.43032, exact finite mean -1.50028±.05502 e²/h: the residual is explicitly **not closed**. Eleven new physics gates include independent iterated BSE and left/right-vertex identity. See [v43 scope](docs/SCBA_SURFACE_V43.md).

```bash
python -m unittest tests.test_scba_surface_v43 -v
python scripts/reproduce_scba_surface_v43.py --W .12 --seeds 24
```

**This is NOT a full conserving disorder-averaged Středa I+II, thermodynamic DC Hall, or original nonlinear EB PHE. Status remains PARTIAL.**

### v44: Energy-dependent RR density Ward and provisional noncrossing Hall sea

The standalone `qgmap.scba_sea_v44` integrates a **candidate Středa II sea** from self-consistent Sigma(E) and the same-W RR density vertex at every real-energy node, then reports the **provisional** sum with v43's one-leg RA I(a) surface. It independently validates `dG/dE` and `dSigma/dE` by finite differences, the W=0 Bloch two-pole analytic formula, cutoff/quadrature convergence, and v42 exact finite-disorder **full I+II**. At W=.03,.12 with mu=1.1, eta=.12 the candidate totals are respectively -0.376124 and -0.456059 e²/h, versus exact 24-seed L4 results -0.389852±.021380 and -0.483690±.027788. A numerical match within finite seed SEM **is not a conserving transport proof**.

```bash
python -m unittest tests.test_scba_sea_v44 -v
python scripts/reproduce_scba_sea_v44.py --W .12 --mesh 35 --L 4 --seeds 24 --ntheta 7 --lower -12 --quadrature 420
```

See [v44 scientific limits](docs/SCBA_SEA_V44.md). Remains `PARTIAL`: transverse electromagnetic Ward/contact, crossed disorder sectors, thermodynamic metallic dc order of limits and original paper nonlinear mixed EB are OPEN.

### v45: Original symmetric spin PHE (distinct from QWZ linear Hall)

The target Eq.(8)/Eq.(11) is **symmetric** in xy/yx, linear in tau, and formulated by the authors in **phenomenological RTA**. The v41-v44 QWZ *antisymmetric* linear Hall/Streda calculations are a useful **separate method extension**, not acceptance criteria for the original PHE. The `dirac_phe_protocol_v45` module independently verifies the printed Zeeman curvatures, main Eq.(8) angular integral and Eq.(11) signed coefficient; at beta=1/8, mu>0 these are +0.002506415907 and -0.002506415907, respectively, under the explicitly named normalized branches. The same fixed-tilt Dirac Hamiltonian with spin Bx is exactly related to zero B by **ky translation**, so the *fully equilibrated, full-plane static-B Drude protocol* has zero linear Bx response after velocity and occupation derivatives cancel. This **does not automatically refute** the published different Liouville/adiabatic RTA branch. Field preparation/collision and Fig2(b) normalization remain OPEN.

```bash
python -m unittest tests.test_dirac_phe_protocol_v45 -v
python scripts/reproduce_dirac_phe_protocol_v45.py --beta .125 --mu .3 --temperature .055 --zeeman-shift .1
```

[Scientific scope](docs/ORIGINAL_PHE_PROTOCOL_V45.md), whole-paper PARTIAL.

### v46: HISTORICAL, numerics superseded by v47 (same Hamiltonian, distinct magnetic preparation)

`qgmap.dirac_mixed_rta_v46` independently calculates the **full E-linear, B-derivative density matrix** of the same fixed-tilt Dirac Hamiltonian under two explicitly different scalar RTA collision targets. If the bath targets the B-equilibrated `f(H_B)`, the complete field response obeys the exact ky-translation identity and source-preparation (+0.07627530) cancels the magnetic commutator (-0.07627530) for beta=1/8, eta=.065. If the collision target remains **frozen at f(H0)**, both B-then-E and E-then-B orderings instead yield total +0.03377501 under the same model conditions. Its low-eta tau-linear coefficient is **near but not equal to** the formal original Eq.(8), and deviates further at larger tilt. Neither selected collision target is yet identified with the full author Supplemental Material RTA. See [v46 detailed method](docs/ORIGINAL_PHE_MIXED_RTA_V46.md).

```bash
python -m unittest tests.test_dirac_mixed_rta_v46 -v
python scripts/reproduce_dirac_mixed_rta_v46.py --eta .065 --beta .125 --nr 96 --nphi 144 --h .00003
```

Original physical PHE, Fig2 normalization and source-specific RTA/bath remain OPEN; scientific status PARTIAL.

### v47: Physically correct electric source and symmetric planar Hall tensor

A centered finite difference of the complete Fermi matrix f(H_B) caught the wrong sign of d_y(kx/r) in the old v46 source. Both xy and yx are now computed, including both EB orderings, and independently reproduced by an entire **finite-B Liouvillian** solve. At eta=.065, beta=1/8, frozen-target xy is **−0.042199786819**, yx **−0.038976453956**, symmetric part **−0.040588120388**; old v46 xy **+0.033775006344** is invalid. At eta=.002, temperature=.03, the signed eta*symmetric responses agree with literal main Eq.(11) to about 5e-5 relative at five tilt values from beta=.025 to .6. Earlier v46 frozen JSON is marked RETRACTED and retained only for archive. See [v47](docs/ORIGINAL_PHE_SYMMETRIC_SOURCE_CORRECTION_V47.md).

```bash
python -m unittest tests.test_dirac_mixed_rta_v46 tests.test_dirac_source_tensor_v47 -v
python scripts/reproduce_dirac_source_tensor_v47.py --beta .125 --eta .002 --temperature .03
```

This is acceptance of an **independently constructed frozen-target model RTA**, **not** confirmation of the author's original Supplemental Material collision prescription or the absolute published Fig.2(b) SI scale. Overall scientific status PARTIAL.

### v48: Model PHE Eq.(11) survives eta, mu, temperature and UV tests

For the corrected v47 *symmetric* spin PHE with **frozen-H0 scalar RTA**, `qgmap.phe_limit_audit_v48` sweeps η (.016→.001), electron/hole doping, two temperatures and multiple ultraviolet cutoffs. At beta=1/8, mu=.3, T=.03 the ratio `eta*PHE / signed main Eq11` decreases from 1.003215 to 1.000012885; a fit against eta² with the three smallest η points extrapolates to **1.0000002435** without any arbitrary factor or sign fit. The η-scaled antisymmetric remainder also tends to zero approximately as η². Across |mu|=.15,.3 and T=.02,.055 the response is odd under electron/hole reversal and follows the model's tanh finite-T smoothing; UV 1–4 has ratio spread below 3e-8. Ten model tests check all limits.

**Critical source gate:** The author's Supplemental Material has **not** been retrieved via the available publisher/arXiv endpoints. v47–v48 reproduce the **main Eq.(11) model geometry**, *not an independently verified original SM relaxation operator*. The literal Eq.(8) sign branch and Fig.2(b) absolute ~2π ordinate issue remain OPEN. See [v48 provenance and limit audit](docs/ORIGINAL_PHE_LIMIT_SOURCE_GATE_V48.md).

```bash
python -m unittest tests.test_phe_limit_audit_v48 -v
python scripts/reproduce_phe_limit_v48.py --mode all
```

Scientific state remains **PARTIAL**.

### v49 authentic original Supplemental Material and main-article audit

For the **actual user-provided primary PDFs**, SM Eq.(4)–(20), Eq.(31)–(33) and Eq.(62)–(63) have been read directly and linked to the published Eq.(10)–(11). `sm_literal_spin_phe_v49` evaluates the literal Eq20 tau-density-current and the Eq31/32 Zeeman-curvature expressions **independently** on both printed Hamiltonians with exactly the same physical sigma_x. At beta=.125,mu=.3,T=.03, Eq20/main gives -0.002506188247 versus Eq31/main +0.002506188247; Eq32/main gives +0.002506188497 whereas Eq32/SM gives -0.002506188497, agreeing with the signed Eq11 (-0.002506188335 finite-T model continuation). The earlier v47 eta-scaled frozen-H0 response approaches the **literal Eq20/main** coefficient, not Eq32/main, under fixed spin/B conventions. Separate adiabatic eta_a/eta_b regulators and one selected tau pole are explicitly distinguished from the v47 all-sector finite-eta bath model. No material experimental correction is claimed. Source PDF raw bytes are **not** committed; hashes and source lines are frozen.

Run `python -m unittest tests.test_sm_literal_spin_phe_v49 -v` and `python scripts/reproduce_sm_source_v49.py --beta .125 --mu .3 --temperature .03`. The original-publisher Eq20/31 and helicity sign issues remain **OPEN**, whole-paper scientific acceptance **PARTIAL**.

### v50: Original Supplemental Eq.(4)/(9) full two-frequency Liouville, independently solved

`qgmap.sm_two_color_liouville_v50` now computes the **actual source SM** first/second-order response using **separate** omega_a, omega_b, eta_a, eta_b and **both** field orderings in a fixed spin basis, with no imagined uniform finite collision operator. First-order Eq.(6)/(7) are verified as eigenbasis fractions at finite frequencies; full second-order Eq.(9) is solved with entire matrix k derivatives. Taking the electric dc **eta_a** pole only *after summing both orderings*, at beta=.125,mu=.3,T=.03,100x128, eta_a=eta_b=.0005 gives **-0.00250621056**. The same-basis authentic SM Eq.(20) direct current gives **-0.00250618825**; **printed SM Eq.(31) gives +0.00250618825**. Two field-ordering terms are individually large and strongly regulator-path dependent. Keeping omega_a nonzero removes the selected 1/eta_a dc divergence. Hence the original Liouville Eq9 independently supports the Eq20 sign, while source Eq20/31 sign disagreement remains OPEN. No published Fig2 normalization adjustment is applied. [v50 scientific protocol](docs/AUTHENTIC_SM_TWO_COLOR_LIOUVILLE_V50.md).

```
python -m unittest tests.test_sm_two_color_liouville_v50 -v
python scripts/reproduce_sm_two_color_v50.py --eta-a .002 --eta-b .001 --nr 100 --nphi 128
```

The original main+SM source SHA256 are retained but raw PDFs are private and not published. Whole paper scientific status **PARTIAL**.

### v51: Original SM Eq20/31 signed-source discrepancy survives complete spin-basis covariance

An independent **generic Hermitian two-band** matrix check (120 random cases, arbitrary U(1) band phases) shows that direct current from source SM Eq.(20), using the very velocity–Berry-connection identity printed below Eq.(31), yields **plus tau sum Z d f per E Bbar**, whereas printed Eq.(31) yields **minus**. This is an explicitly bounded local source formula inconsistency. The opposite spin helicities in main Eq.(10) and SM Eq.(62) are basis-related by **U=sigma_z**, but that same transformation must flip **the physical sigma_x magnetic operator** and the current velocity, leaving all correctly paired results invariant. Flipping the Hamiltonian alone changes both signs but **cannot** cure their relative disagreement. Separately, original SM Eq.(33) `sum Z_nm = d_k <sigma_x>_n` passes independent finite-difference spin-texture and integrated xy/yx checks (~1e-10 at generic momentum). See [v51 primary-source signed-covariance report](docs/AUTHENTIC_SM_SIGN_COVARIANCE_V51.md).

```bash
python -m unittest tests.test_sm_sign_covariance_v51 -v
python scripts/reproduce_sm_sign_covariance_v51.py --beta .125 --mu .3
```

No original PDF is redistributed, source SHA256 hashes remain pinned. This does not declare any physical experiment invalid. Whole-paper science **PARTIAL**, original signed Eq20/31 source **OPEN**.
