# SOURCE MAP · Xiang et al., PRL 134, 116301 (2025)

Date: 2026-10-19. Status: source inventory checked against **all 7 pages of article and all 11 pages of SM**, including article's two composite figures and SM tables. This file maps evidence, not validation.

## Source provenance

- Xiang, Jia, Xu, Qiao, Wang, *Intrinsic Gyrotropic Magnetic Current from Zeeman Quantum Geometry*, PRL **134**, 116301 (2025), DOI [10.1103/PhysRevLett.134.116301](https://doi.org/10.1103/PhysRevLett.134.116301).
- Article locally reviewed: 7-page PDF, SHA256 `6b61f8a1ddd37123053d0de744795f986587d837e2d9706cc5640580d94080c2`.
- SM locally reviewed: 11-page PDF, SHA256 `9efa730976083b892e201659509e75cdc640e88a5065abd190524a99bc82b073`.
- Neither PDF is redistributed in this reproduction package. Consult the user's local originals.
- Public source code baseline: [`Xin-Jiaqi/quantum-geometry-reproductions`](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions), `main` commit `78a34431bd4a165484ac3f714786b98120132f81`, prior module `papers/jia-2026-nonlin-magnetoelectric-edelstein/`.

## Main text, equation by equation

| Eq. | Content | Independent numerical/theory check |
|---|---|---|
| (1)-(2) | Momentum quantum distance; conventional band-resolved QGT `g=rr` | gauge phases; metric/curvature Hermiticity |
| (3)-(6) | Combined spin-rotation/momentum quantum distance; Zeeman QGT `z=rσ`, Zeeman quantum metric Q and curvature Z | mixed parameter Berry geometry; gauge invariance; P/T |
| (7) | AC Zeeman-induced linear density matrix, infinitesimal `η`, equilibrium occupations | direct finite-ω Liouville/Kubo susceptibility |
| (8) | Current computed from interband velocity and `ρ^(1)` | independent frequency-domain band summation |
| (9)-(10) | Nonresonant split into `J_conduction` and polarization derivative `∂t P` | frequency parity, Fermi surface versus sea |
| (11) | `σ^C_ab = ∑_nm∫ f_n Z_nm^ab` | Fermi surface identity Eq. (13), Kubo ω→0 |
| (12) | `σ^D_ab = ∑_nm∫ f_n 2ω Q_nm^ab / ε_mn`; neglect O(ω²) in P | low-ω Kubo slope, insulating nonzero response |
| (13) | `σ^C_ab= -∑_n∫ f'(ε_n) v_n^a σ_nn^b` | Fermi-surface integral vs Eq. (11), periodic integration by parts |
| (14) | Magnetic point-group tensor constraint, `η_T=+1` unitary / `-1` antiunitary **as printed for C/D** | audit against P/T table, TR-invariant Rashba and SM Jahn symbols **(open issue)** |
| (15) | 2D continuum Rashba spin Hamiltonian | analytic eigenvalues, symmetry, spin geometry |
| (16) | Zero-T closed form `σ^C_xy(μ)` for Rashba | independent occupied-interval and Fermi-surface integrations; correct e/ℏ units |
| (17) | Square-lattice regularized massive Dirac spin Hamiltonian | band spectrum, period BZ, symmetries, two responses |

## Main-text Figure acceptance inventory

| Panel | Observable and input | What exactly must be reproduced |
|---|---|---|
| Fig. 1(a) | Rashba eigenenergies vs `k_x`, `k_y=0`, chemical potential dashed | two shifted parabolas; minimum at `-λ_R²m/(2ℏ²)` |
| Fig. 1(b) | Conduction Hall response `σ^C_xy(μ)` | Eq. (16): zero below band bottom; square-root window; constant above crossing |
| Fig. 1(c) | `Z^{-+}_{xy}(k_x,k_y)` | two-lobe even/anisotropic texture; correct sign and singular `k→0` behavior |
| Fig. 1(d) | `Z^{-+}_{xx}(k_x,k_y)` | alternating four-lobe texture, angular parity and sign |
| Fig. 2(a) | Dirac-lattice bands along displayed `k_x` cut | gap `2β=0.5 eV` at Γ and high-k structure |
| Fig. 2(b) | `σ^C_xy(μ)` and inset `Z` distribution | *exactly zero in the zero-T band gap*, peaks near band edges; tensor antisymmetry |
| Fig. 2(c) | `2Q^{-+}_{xx}/ε_{+-}` over k, inset bare `Q^{-+}_{xx}` | band-gap-weighted central Γ hotspot **separately** from unweighted metric |
| Fig. 2(d) | `σ^D_xx(μ)` at finite driving frequency | in-gap response, linear low-ω scaling, resonant caveat |

Main Fig. 1 parameters (caption): `λ_R=0.6 eV Å`, `m=1.1m_e`.

Main Fig. 2 parameters (caption): `β=0.25 eV`, `a=1 Å`, `v_F=0.825 eV Å`, `M=0.1 eV Å²`, `ω=10^12 Hz` as printed. Note: reported `Hz` versus angular `ω` must be disambiguated before absolute SI-value acceptance. Article plots carry prefactor `-g μ_B σ` with the axes' chosen physical units; do not confuse this plotted observable with bare `σ`.

**Digitized numeric data:** not provided in PDFs. Any reconstructed point values from figures must be labeled independently digitized and given extraction errors. First acceptance checks analytic shapes and signs, not pixel-perfect overlay.

## Supplemental Material (complete section coverage)

| Section / SM equations | Actual content | Acceptance |
|---|---|---|
| I, (1)-(5) | expansions of both quantum distances and cross `z` term | independent expansion / invariant Q/Z |
| II, (6)-(11) | Bloch-basis Liouville, E and B, first-order `ρ` | recover main Eq. (7), carefully retain diagonal field/occupation term if E is switched on |
| III, (12)-(34) | semiclassical wavepacket, interband `M_mn`, positional shift, velocity, both IGMC currents | independent wavepacket derivation of low-ω C and D |
| IV, (35)-(47) | lesser Green's function Kubo derivation | recover main Eq. (8) and full-ω benchmark |
| V, Tables I, II, III | allowed magnetic point groups, separate conduction/displacement tensor forms and classification | independently generate selected unitary/antiunitary subgroup constraints, later exhaustive coverage |
| VI, no numbered main figure | electromagnetic Faraday-induced electric field, field geometry, frequency/gap discussion | quantify applicability of purely Zeeman driving; not a further numerical plot |

**SM numerical figure inventory: none** (SM has three tables, derivation equations and experimental discussion). No SM DFT datasets or material-dependent calculations are supplied. Two continuum/tight-binding example Hamiltonians in the main text are the entire article's explicit computational benchmarks.

## Scope and reuse

Reuse from prior NMEE: independently checked matrix/spectrum infrastructure, numerical tests, SHA256 provenance, separate named physics protocols, offline `verify_release`, CI/handoff design. **Do not reuse** the NMEE's static EB susceptibility, CuMnAs/PT degenerate backend or its response normalization as an IGMC implementation. `wannier90_hr.dat` alone is insufficient: projected Pauli-spin operator in a consistent Bloch/Wannier gauge is required.
