# Scientific acceptance (public core, 2026-10-20)

**v49 primary-source update:** The original 8-page published PRB and authentic 16-page Supplemental Material were directly supplied by the user and SHA256-verified. Statements from v45–v48 that SM was inaccessible are superseded. The original mixed-field τ pole (SM Eq15/20) is now source-verified, but a literal SM Eq20/31 sign difference and opposite main Eq10/SM Eq62 helicities are OPEN.

**Overall: PARTIAL.** The original paper works with phenomenological RTA; conserving SCBA-Streda is a separate extension, not a prerequisite for the symmetric PHE. Software regression PASS verifies only the published code and named numerical protocols, not every claim in Xiang et al., PRB 113 L201406 (2026).

| Branch | Public capability | Scientific gate |
| --- | --- | --- |
| Spin Zeeman geometry | Direct eigenvectors versus analytic two-band response, helicity and U(1) gauge tests | PASS for stated nondegenerate model |
| Four-band spin MNHE | SM Eqs. (16–19) density-matrix currents compared with Eq. (29) | PASS synthetic nondegenerate model; physical driven DC and degeneracy OPEN |
| Tilted Dirac Fig. 2 | Spectrum, local Zeeman geometry, formal Eq. (11), literal SM Eq. (63) branch | PASS for stated model definitions |
| Fig. 2(b) physical ordinate | Literal SI Eq. (11) remains about (2\pi) above the published plateau in independent audits | OPEN |
| Two-color Liouville | Sequential resolvents, both orderings and separately defined equilibrium population derivative | PASS method subset; no general steady-state claim |
| Steady-state static PHE | Collision operator, bath preparation, Ward/contact terms and order of limits | OPEN |
| Orbital OHE | Classical one-band and analytic T=0 Rashba Fermi-ring benchmarks | PASS for named models |
| Orbital PHE (2D Bz) | Literal SM Eq. (56) compared with finite-field semiclassical derivative | PASS model method; printed-line sign difference OPEN |
| Orbital MNHE (strict 2D and 3D) | SBC, Eq. (49) positional shift and literal Eq. (52) comparison | PASS named model subset; factor-two difference OPEN |
| SCBA + ladder | Matched scalar-disorder self-energy and noncrossing Bethe-Salpeter vertex | PASS independent periodic model; not full PHE |
| Four-channel spectral Kubo | (RA+AR-RR-AA), with channel-specific ladders | PASS symmetric finite-eta spectral sector; not full Kubo–Středa |
| Finite-(q,omega) Ward | Matrix resolvent identity with exact finite-difference longitudinal insertion | PASS finite-eta resolvent identity; physical continuity remains OPEN |
| Exact real-frequency Kubo charge continuity | `real_frequency_continuity_v20`: full finite-system Lehmann sum and lattice contact terms | PASS finite-disordered longitudinal benchmark; mixed-EB transverse PHE OPEN |
| SCBA regulator order | `limit_order_v20`: eta leakage at fixed nonzero probe and at uniform static origin | PASS numerical model limit-order diagnostic; physical dc Ward OPEN |
| Finite-q exact mixed EB | `electromagnetic_eb_v21`: Kubo/Liouville static-B derivative, density-current contact and pure gauge | PASS independent finite-q test; uniform DC PHE OPEN |
| Finite-q two-color | `two_color_exact_v22`: both E/B orderings; resolvent, Sylvester and time-domain | PASS two-color finite-system method; steady-state OPEN |
| Homogeneous length gauge | `homogeneous_length_gauge_v23`: FFT covariant derivative, equilibrium preparation and SM Eq.20 residue | PASS periodic two-band method; continuum Fig.2/physical DC OPEN |
| Continuous tilted Dirac SM Eq. (20) | `continuum_ti_v24`: independent projector vs eigensolver; analytic angular kernel and two-color limits | PASS model-adiabatic source formalism; transport DC OPEN |
| SI Liouville + Eq. (11) | `si_fig2_v25` and `physical_si_liouville_v26`: joules, magnetic moment, electric current and conductance | PASS standard-SI literal coefficient; published Fig. 2(b) normalization OPEN |
| Published Fig. 2(b) digital fit | `fig2_source_audit_v27`: approximate 20/40 K published-raster curve digitization and T fits | PASS approximate curve-shape check; author raw data unavailable |
| Clean antisymmetric Kubo–Středa sea | `streda_clean_sea_v38`: exact two-pole energy primitive, independent matrix energy integral and FHS Chern | PASS only clean two-band gapped benchmark; metallic/disordered full Hall OPEN |
| Energy-dependent SCBA RR Ward + sea | `scba_streda_v39`: independent self-energy finite difference versus RR density ladder; clean W=0 and exact finite-disorder controls | PASS for model derivative and defined bare-current sea; complete disordered Středa Hall and EB OPEN |
| Same-disorder RA/AR/RR/AA current ladder | `current_vertex_v40`: four-Pauli-basis BSE, finite-eta RA density Ward, RR momentum Ward, spectral Hall-zero, finite-torus bare RA | PASS QWZ model method; physical disordered Středa and EB DC OPEN |
| Independent RR spatial Ward by vector-potential twist | `rr_spatial_ward_v40`: finite differences of separately converged \(\Sigma^R(A_x,A_y)\) and \(G^R(A_x,A_y)\) versus RR current ladder | PASS same-model two-direction twist method; RA Hall, full Středa and EB DC OPEN |
| Exact finite-disorder *linear* Hall Kubo (v41) | Full finite-torus occupied/empty Lehmann current response; independent retarded current-current Kubo; twist-averaged Hall versus occupied-subspace FHS Chern; model-specific zero mixed diamagnetic contact | PASS for named gapped QWZ sample and finite-size method; disordered Streda I/II, conserving disorder average and nonlinear EB OPEN |
| Exact finite-disorder linear Středa I(a)+I(b)+II (v42) | Same v41 frozen-sample Peierls QWZ torus; independent dense surface Green trace, independent direct sea energy quadrature, all-pair energy primitive, clean v38 Bloch check, finite-seed ensemble | PASS for explicitly named finite nondegenerate QWZ samples as external eta->0, against independent v41 Hall; conserving disorder-averaged SCBA Středa and nonlinear EB OPEN |
| Same-disorder SCBA noncrossing RA Hall I(a) vs exact finite-sample ensemble (v43) | v39 self-energy + v40 one-leg current ladder, independent iterative BSE, RA density resolvent Ward, v42 dense Green exact 24-seed/7x7-twist disorder ensemble, finite-L sample scaling | PASS *method and explicitly reported approximate comparison*; nonzero residual remains, matched conserving sea II and full Hall/EB OPEN |
| SCBA RR density-vertex energy-dependent sea II + RA one-ladder surface (v44) | W-matched selfenergy at each energy, RR density derivative verified by independent finite differences, independent W=0 Bloch-pole primitive, fixed-seed v42 full I+II ensemble benchmark | PASS **model-level provisional SCBA I+II candidate** with visible remaining difference; complete conserving electromagnetic Ward/contact, thermodynamic disorder and original nonlinear EB all OPEN |
| Original symmetric spin PHE Eq.(8)/(11) vs equilibrated Bx Dirac shift (v45) | Eigenvector Zeeman curvature, independent angular Eq.(8); exact uniform Bx-to-ky translation, static-B Drude velocity/occupation cancellation and finite UV control | PASS named original model identities; **same-preparation complete RTA PHE still OPEN**, Eq.(8)/Eq.(11) printed sign difference remains |
| v46 historical model RTA (source derivative retracted by v47) | The original k_y electric source contained an incorrect sign in d_y(k_x/r); previous v46 frozen numbers and its single-xy Eq8 match are INVALID | **RETRACTED numerics**; repaired implementation and all new physical gates in v47; original author SM still OPEN |
| v47 repaired original symmetric PHE tensor vs main Eq.(11) | Independent finite-difference of full f(H_B) electric source; corrected sign; both xy and yx field orderings; independent direct finite-B Liouville superoperator; eta->0 numeric beta scan | PASS **named frozen-H0 model**: symmetric eta*response matches signed Eq.(11) at beta=.025-.6 to about 5e-5 relative; author SM collision operator and Fig2 SI amplitude OPEN |
| v48 original PHE Eq. (11) physical limit robustness and source-provenance gate | Corrected v47 symmetric tensor across eta, eta-squared extrapolation, mu electron/hole reversal, thermal tanh crossover, four UV cutoffs, full Fermi source finite differences | PASS the **named frozen-H0 model**: eta² extrapolated ratio 1.00000024; source-author Supplemental Material relaxation inaccessible, hence **author-method identity OPEN** |
| v49 authentic source PDF/SM provenance | Independently inspect original 8-page PRB and 16-page SM, freeze both SHA256; read SM Eq4–20,31–33,62–63 rather than infer from historic notes | PASS original source now available and directly checked |
| v49 literal SM Eq20 → Eq31/32 sign | Independent two-band eigenbasis current from actual SM Eq20 versus printed Eq31/32, both XY/YX, U(1) gauge, multiple tilts and raw signed physical B | **OPEN**: literal Eq20 contraction and Eq31 are equal magnitude/opposite sign; must resolve the manuscript's index/velocity conventions |
| v49 main Eq10 versus SM Eq62 helicity | Direct comparison of two opposite spin helicities; sigma_z rotation also changes Zeeman sigma_x; source Eq32/SM Eq62 agrees with Eq11 whereas main Eq10 gives negative | **OPEN** physical spin/B sign convention; no silent correction |
| v50 authentic full two-colour source SM Eq4–9 | Two independently regulated complex frequencies, fixed spin-basis Liouvillian, independently verified SM Eq6/7, both EB orderings, finite-omega conjugation, eta_a/eta_b limit path | PASS original source-method protocol for named 2-band model; physical bath/experiment OPEN |
| v50 independent original Eq9 -> Eq20/31 source sign | Extract eta_a times complete current along multiple regulator paths, compare to source literal Eq20 and printed Eq31 *without inserting either in source solver* | **Eq20 sign supported; Eq31 opposite sign remains an OPEN source-algebra inconsistency** |
| v51 SM Eq20/31 intrinsic signed-algebra and spin-rotation covariance | Random generic Hermitian two-band Eq20 contraction vs printed Eq31; physical U=sigma_z transformation of H, v, spin source; independent Eq33 spin derivative and integrated xy/yx; v50 source Liouville cross-check | **PASS named algebra/covariance gates**; Eq20→printed Eq31 remains **OPEN local source sign discrepancy**; no full physical PHE acceptance |
| v52 published Fig2b + printed 43 microV absolute SI | Source-pinned original 8-page PRB Fig2b ordinate in e²/hbar; independent Eq11 J/s and eV/s units, exact e²/hbar vs e²/h, Hall voltage Ey width R, 20/40 K pixel traces, parameter sweeps | **PASS bounded source-level 2pi discrepancy**: literal Eq11 Fig2b 11.0208 and 268.26 microV; alternate e²/h scale 1.754 and 42.70 microV matches published approximately 1.75/43. Actual author plotting convention **OPEN**, no inserted correction |
| v53 orbital SM Eq49→52 metric velocity coefficient | Source Eq49 explicit 2; genuine periodic 3D two-band quantum metric, source Eq48 direct F derivative vs printed Eq52 and 2-variant, BZ mesh convergence and zero-tilt control | **OPEN printed factor 1**: source Eq48 selects factor 2; local source-algebra discrepancy, not an author-certified erratum |
| v53 orbital SM Eq52 explicitly omitted three-band processes | Genuine periodic three-band Hermitian Bloch velocity/Berry geometry, full curvature-orbital-moment double sum vs retained m=l, U(1) gauge and Fermi-energy checks | **PASS exact decomposition; physical approximation OPEN**: m≠l is 92.9% of selected full subterm in one model, NOT total MNHE |
| v53 orbital PHE Eq56 printed line equality | Actual smooth periodic two-band Bloch moment/Berry curvature, both printed line integrands, independently integrated BZ and N^-2 convergence | **PASS named Eq56 integration-by-parts equivalence**; entire orbital transport remains PARTIAL |
| v54 independent reusable N-band Bloch geometry | Per-band/pair Berry curvature, quantum metric and energy-type orbital moment in nondegenerate fixed orthonormal k-independent orbital basis; independently cross-check analytic-velocity spectral results with projector finite differences, U1 phases, constant orbital unitaries, two/three-band models and degeneracy rejection | **VERIFIED named primitive only**, no full E×B response |
| v54 all m=l / m≠l decomposition of Omega × orbital moment geometric sector | Exact pairwise product, explicit independent band-triple loops; genuine periodic three-band example reproduces full -0.00079945105, retained -0.00154218315, crossed +0.00074273210 | **VERIFIED complete geometric subterm**, **NOT total orbital MNHE** |
| v54 actual magnetic-field positional shift, non-Abelian/Wannier embedding and full orbital MNHE | Need physically specified B-coupling, orbital positions, all omitted three-band diagrams, contact/bath terms and independent full-response oracle | **UNRESOLVED**; prior v53 paper-level PARTIAL remains |
| General orbital QMQ/BCQ | General omitted three-band pieces and complete Kubo–Středa physical DC response | OPEN |
| Real Wannier materials | Requires consistent projected magnetic/position operators, Berry connection and validated response protocol | BLOCKED |

This matrix applies only to the code in this GitHub repository. The original v30 research archive includes additional offline programs and results that are not represented by remote CI. In particular, a finite-eta resolvent Ward identity is not a proof of physical charge conservation at finite eta, and a symmetric spectral Kubo sector is not a complete Kubo–Středa calculation.

## v55 physical capabilities, independent of source paper

| Module | Acceptance | Physical gates |
| --- | --- | --- |
| k-dependent orbital basis/position connection and D_iH=d_iH-i[A_i,H] | **VERIFIED** for complete supplied projected A_i and nondegenerate finite bands | Non-diagonal W(k), independent transformed-H k derivative and geometry invariance; wrong no-A control |
| Finite open graph Peierls H(B), d_BH | **VERIFIED** for explicit site coordinates, straight paths and signed q/hbar | Magnetic gauge spectra, Wilson flux, B finite differences, field reversal/origin shift |
| B-dependent magnetic gauge transformation | **VERIFIED** | Full d_B G derivative: residual 2.8e-17 versus naive omission 0.522 |
| Full Wannier embedding or total orbital E×B magnetotransport | **UNRESOLVED** | Position dipoles, current/contact terms, collision physics not yet certified |

15 new independent physics tests; original-paper scientific status **PARTIAL**.


## v56: complete finite-open-model mixed EB versus conditional transport physics

| Independent acceptance gate | Result |
| --- | --- |
| Finite Peierls H(B), real-space P and current operator J=i[H,P], full f(H_B) preparation | **VERIFIED named finite-model operators** |
| Complete analytic d_B sigma, including current contact, Fréchet equilibrium matrix derivative and Liouvillian term | **VERIFIED numerically** versus independently evaluated 5-point full finite-B differences (error around 1e-13) |
| Magnetic gauge covariance, Onsager reciprocity, ac reality, q/B inversion, origin change and zero static equilibrium current | **VERIFIED in specified finite open graph** |
| Mixed ac EB response with scalar gamma relaxation toward field-equilibrated f(H_B) | **CONDITIONAL bath-dependent physics** |
| Actual intrinsic infinite-system orbital MNHE and conserving scatterer/current contact theory | **UNRESOLVED** |

New independent tests 19/19 local PASS; GitHub CI required before merge. Xiang original source paper scientific state remains PARTIAL.


## v57 physics-first independent four-terminal open-device acceptance

| Physical capability | Acceptance | Independent validation |
| --- | --- | --- |
| Local wide-band coherent two-/four-lead Caroli scattering matrix and lead magnetic derivative | **VERIFIED specified finite-model implementation** | Independent Caroli trace and local-site formula, retarded/advanced optical theorem, same T_ab in two gauges, analytic d_B T against five-point finite B |
| Floating voltage probes (zero I_top/I_bottom), Hall ratio and magnetic slope | **VERIFIED lead current-conserving solution** | Independent complete five-point B derivative of *all* contacts and voltages, Onsager T_ab(B)=T_ba(-B), total/current-probe conservation |
| Two-probe B-even longitudinal conductance; symmetric four-probe B-odd Hall voltage | **VERIFIED model symmetries** | Explicit B→-B and B=0 tests, contact/temperature/coordinate dependence |
| Four-terminal device voltage for particular finite graph and wideband leads | **CONDITIONAL physically modeled coherent device** | Lead Γ, temperature, geometry materially change response |
| Intrinsic thermodynamic bulk MNHE/PHE, scattering/microscopic leads and bulk size limit | **UNRESOLVED** | Sparse diamond ring enlargement is *not* controlled thermodynamic scaling |

20 new independent physical tests PASS locally; full GitHub CI required before accepting merged module. Original-paper source status remains PARTIAL.


## v58: fixed-site-density planar family and verified adaptive energy integration

| Gate | Status | Independent evidence |
| --- | --- | --- |
| Fixed a,t,onsite,mu and qBa²/hbar across real Nx×Ny rectangular devices, centered multi-site wideband leads | **VERIFIED model construction** | True site/hopping checks, disjoint symmetric contacts, per-site vs fixed-reservoir-total Γ |
| Multi-site Caroli T and d_BT, floating voltage Hall ratio and true Rxy=VH/I | **VERIFIED mathematical algorithm** | Independent Caroli trace; v57 single-site limit; Onsager/current conservation and magnetic-gauge invariance |
| Complete magnetic derivative of Rxy and Hall voltage | **VERIFIED within finite model** | Independent five-point full B/transmission and re-solved floating probes; matched rotated rectangles |
| Narrow coherent resonances and energy integration | **VERIFIED convergence gate** | Deliberate 80/120/180-point Gauss mismatch; fail-closed adaptive split at H eigenvalues, separately tightened integration, adaptive two magnetic gauges |
| Size/contact-dependent spinless coherent Peierls Hall voltage | **CONDITIONAL device observable** | Nonmonotonic finite-size response and large changes with lead coverage/Γ normalization |
| True bulk orbital MNHE/PHE or fixed carrier-density thermodynamic limit | **UNRESOLVED** | This fixes lattice-site density, NOT electron filling, and has no microscopic inelastic/dephasing/lead bulk oracle |

25 new independent local tests. Source Xiang original manuscript remains **PARTIAL**; v58 is a physics tool, not source-equation correction.


## v59 exact finite-band coherent leads + independent periodic bulk charge Kubo

| Scope | Status | Science gate |
| --- | --- | --- |
| 1D semi-infinite chain retarded g, Σ and Γ(E) | VERIFIED | Dyson equation, causal ImΣ, band edges and high energy asymptote |
| Full four-probe Caroli magnetic derivative, including optional dBΣ and dBΓ | VERIFIED named mathematical algorithm; CONDITIONAL physical device | Explicit Caroli trace and optical theorem, both gauges, independent complete five-point magnetic field derivatives, v57 wideband limit |
| Periodic bulk charge AHE QWZ | VERIFIED bounded charge Hall benchmark | Separate spectral Berry versus Fukui-Wilson links, 0,+1,-1,0 Chern, gapless phase refused |
| Genuine orbital bulk mixed E/B, realistic matrix Wannier electrodes and conserving microscopic bath | UNRESOLVED | Original-paper science remains PARTIAL |

21 new independent scientific tests. Spinless Onsager verified only for B-independent lead onsite; a B-odd scalar lead shift may violate it.


## v60 independent physics acceptance

| Component | State | Independent physical evidence |
| --- | --- | --- |
| General n-orbital semi-infinite block lead retarded g(z), finite positive eta | **VERIFIED** | Matrix Dyson residual, spectral causality, 160-layer independent sparse inversion and complex one-orbital analytic root |
| Full matrix/nonlocal contact C, Sigma, Gamma and magnetic derivative d_B C | **VERIFIED specified finite-model algorithm** | Gauge transform with nontrivial G(B), negative omitted-contact control, full independent finite B Hall derivative |
| Same QWZ Fourier H for periodic Bloch Chern and finite-open local projector Chern marker | **VERIFIED topology only** | Exact hopping Fourier reconstruction, center marker tends to periodic Chern, vanishing finite total marker with edge compensation |
| Specific coherent two-orbital finite Hall voltage | **CONDITIONAL** | Depends on eta, electrode band structure and contact C |
| Intrinsic bulk orbital EB and same-H **transport** finite-to-bulk limit | **UNRESOLVED** | Topology agreement is not an electrical conductivity equality |

21 new tests for independently derived physics. Original Xiang paper scientific audit remains **PARTIAL**.


## v61: same-H QWZ coherent four-terminal charge Hall, NOT original orbital mixed E/B

| Component | Science status | Independent gates |
| --- | --- | --- |
| QWZ Bloch ↔ finite sample ↔ matching semi-infinite leads | **VERIFIED** for the named exact QWZ Fourier Hamiltonian | Analytic Fourier reconstruction, Hermitian finite matrix, correctly oriented 2-orbital transverse strips and matrix H01; lead coupling matches actual cut bonds |
| Zero-T four-probe charge Hall Landauer-Büttiker implementation | **VERIFIED mathematical scattering problem; CONDITIONAL finite device** | Independent S†S scattering oracle and probability blocks versus Caroli; positive Γ, floating I_T=I_B=0, total conservation, orbital gauge rephasing; wrong unchanged contacts negative control |
| Single-channel charge Hall plateau vs periodic Chern | **CONDITIONAL clean-gapped model bulk-edge consistency** | m=-1/+1 yield σxy=+1/-1 q²/h, Rxy≈+1/-1 h/q², correct longitudinal current; 3/5/7 size, η/μ/contact sensitivity |
| Trivial insulating m=+3 | **VERIFIED negative control** | Chern=0, tiny η-dependent current; Rxy deliberately undefined when |I|<10^-7 rather than fabricated very large resistance |
| Same microscopic H full Kubo conductivity / four-terminal resistance equality | **UNRESOLVED** outside clean chiral plateau | These are distinct observables; equality requires tensor inversion and carefully ordered contact, finite-size and transport limits |
| Intrinsic multiband orbital mixed E/B of Xiang et al. | **UNRESOLVED** | No substitution of charge Hall topology or finite Hall voltage for orbital magnetoelectric transport |

28 new scientific tests, independently exercised locally. Reference [v61 physical report](SAME_QWZ_FOUR_TERMINAL_BULK_EDGE_V61.md). Source paper status remains **PARTIAL**.


## v62: periodic finite-T charge Kubo vs real-space NEGF edge and contact geometry

| Response/algorithm | Science status | Independent evidence |
| --- | --- | --- |
| Periodic QWZ intrinsic charge σxy(μ,T) | **VERIFIED specified 2-band band intrinsic charge** | Independent analytic Pauli d-vector curvature and eigenbasis velocity Kubo; signed opposite Chern, metallic μ=1.5 and thermal suppression, invalid gapless phases rejected |
| Finite QWZ real-width lead and local spectral/bond currents | **VERIFIED specified finite-system identities; CONDITIONAL physical device** | v61 full-width limit, actual truncated strip H00/H01, A=i(Gr−Ga)=Σ GrΓGa, current commutator, zero interior divergence, local orbital-gauge covariance and contact position controls |
| Four-probe voltage plateau vs spatial current | **CONDITIONAL** | At L7,w3,m=-1: Rxy≈0.99999783, but central x-section current≈0.0228959 vs injected 0.9975553, due to locally invasive floating probes |
| Nontrivial metallic Kubo→Hall resistivity inversion | **NEGATIVE/UNRESOLVED** | μ=1.5: σxy=-0.07741266 and finite Rxy≈0.01001052, so Rxy ≠ 1/σxy; nonzero σxx and contact/bath needed |
| Full original orbital mixed EB | **UNRESOLVED**, manuscript audit remains **PARTIAL** | None of these charge-Hall source-independent tests replaces complete original orbital E/B Ward/contact method |

30 new independent local tests, plus v60/v61 49 existing tests. [v62 physics report](QWZ_BULK_KUBO_SPATIAL_CHARGE_V62.md).


## v63: non-invasive voltage-probe coupling limit (charge transport only)

| Question | Status | Evidence / boundary |
| --- | --- | --- |
| Top/bottom coupling scaling | **VERIFIED** | (C_{T,B}	olambda C), fixed lead (g^r), therefore (Sigma,Gammaproptolambda^2); left/right drive leads are held fixed |
| Local probe-current invasiveness | **VERIFIED for this finite QWZ device** | (L=7,w_d=w_p=3): (Delta J_{cut}) falls from 0.974659 at (lambda=1) to (7.58\times10^{-8}) at (10^{-4}), with interior continuity and zero net floating-probe current |
| Weak-probe (R_{xy}) | **CONDITIONAL** | At fixed geometry (R_{xy}	o1.00215746,h/q^2); do not claim exact quantization, because it differs from the strong-probe result and is size/position dependent |
| Probe width | **CONDITIONAL geometry result** | With drive width fixed at 3, wider probes perturb the local current more strongly at unit coupling; weak-coupling current-profile perturbation vanishes |
| (eta) and size limit | **UNRESOLVED generally** | (eta=10^{-3}ldots10^{-6}) varies (R_{xy}) by ~(9.2\times10^{-6}) for this device; sizes 5–11 remain measurably size-dependent and do not establish a thermodynamic limit |
| Kubo–NEGF same-observable limit | **UNRESOLVED** | No general proof for metallic (sigma_{xx}), disorder, bath, Ward/current-vertex or ordered ((lambda,eta,L)) limits |
| Xiang orbital mixed (E/B) | **UNRESOLVED**, paper **PARTIAL** | The probe-limit study is charge-only and does not complete Wannier orbital B perturbation or conserving mixed Kubo–Středa response |

32 new v63 tests; local v60–v63 standalone **111/111 PASS**. [v63 report](QWZ_NONINVASIVE_PROBE_LIMIT_V63.md).
