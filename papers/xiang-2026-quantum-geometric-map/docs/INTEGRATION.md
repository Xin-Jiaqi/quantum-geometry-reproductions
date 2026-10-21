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
