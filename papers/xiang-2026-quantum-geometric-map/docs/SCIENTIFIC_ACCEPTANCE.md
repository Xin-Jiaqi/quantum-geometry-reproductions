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


## v64: independent two-terminal virtual floating probes, NOT orbital mixed EB

| Capability | Status | Independent science gate |
| --- | --- | --- |
| Virtual lambda→0+ voltage probe | **VERIFIED mathematical finite-QWZ operator-level limit** | Two-drive-only G0; probe weights w_pL/R from independent semi-infinite QWZ strip Gammas; no top/bottom Sigma inserted; compare v63 lambda² approach |
| Local current disturbance and Hall finite-lambda deviations | **VERIFIED specified coherent finite model; CONDITIONAL measurement** | Leading lambda² scaling, exactly fixed driving contacts; distinct weak positive coupling vs undefined lambda=0 voltage |
| Probe geometry, finite size and eta | **CONDITIONAL** | L=5,7,9,11 still drifts; offsets shift Rxy despite invariant drive T; width spectral weights alter even when selected Hall ratios nearly cancel |
| Trivial insulator and metallic fillings | **NEGATIVE controls** | Under-threshold I yields Rxy=None; metallic virtual Hall ≠ reciprocal intrinsic bulk Hall |
| Thermodynamic bulk electrical transport and Xiang orbital E/B | **UNRESOLVED** | No controlled L→∞/eta→0 or contact-fraction family, nor orbital position operators, full conserving mixed-field response |

26 new independent local v64 tests, 137 total v60–v64 standalone. [Details](QWZ_VIRTUAL_PROBE_ORDERED_LIMITS_V64.md). Original-paper source audit remains **PARTIAL**.


## v65 | Static projected-orbital E/B Maxwell response, NOT original transport

| Mathematical / physical scope | Gate | Status |
|---|---|---|
| Supplied non-diagonal finite Wannier-projected P = q diag(position)+D | Independent Peierls H(B), full P_B and d_B P_B, degenerate Fermi matrix Fréchet; B-dependent symmetric/Landau and noncommuting orbital representations | **VERIFIED named finite static model** |
| Complete fixed-mu equilibrium mixed grand-potential derivative | 3 independent fully rediagonalized finite-difference Maxwell routes; omitting position magnetic contact fails gauge covariance | **VERIFIED equilibrium thermodynamics** |
| Fixed-N Helmholtz response and particle conservation | mu_B Ward, numerical chemical potential refit at every perturbed E/B, 3 independent canonical derivatives, P+cI origin negative control | **VERIFIED static canonical sector** |
| Microscopic orbital Peierls coupling / position dipoles | Finite four-orbital graph, supplied D, no real material / no complete field-dependent Wannier embedding | **CONDITIONAL physical approximation** |
| Original Xiang nonlinear orbital E*B transport, conserving Kubo-Streda, bath and thermodynamic limit | NOT calculated by equilibrium Maxwell derivative; SM Eq49/52 and sign/3-band issues remain | **UNRESOLVED; source paper PARTIAL** |

40 new independently checked v65 physics tests; stand-alone v60-v65 baseline 177 tests. [Scientific details](STATIC_ORBITAL_EB_MAXWELL_V65.md).

## v66: projected orbital mixed E/B at finite frequency, conditional gamma bath

| Mathematical/physical sector | Independent science evidence | Status |
| --- | --- | --- |
| Full projected positions Pa(B), Pb(B) and explicit electric current contact -i[Pb,Pa] | Supplied noncommuting orbital dipoles, fixed B Peierls finite graph, electric perturbation H-E Pb; omissions fail static gamma→0 equilibrium current | **VERIFIED named finite projection** |
| Complete analytic magnetic derivative of AC response | Seven named terms: magnetic current H, current dPa, drive dPb, equilibrium Fermi, magnetic Liouvillian, electric-current contact dB, electric-current contact rhoB | **VERIFIED symbolic/matrix finite-model identities** |
| Magnetic/orbital frame covariance | Symmetric, Landau and noncommuting U(B)=exp(iBK) give identical complete complex response; 5-point full B differences error ~1e-12; independent interband eigenbasis response | **VERIFIED specified AC model** |
| Independent finite E+B dc superoperator current | Dense nonlinear finite-E Liouvillian solved directly, 4-corner mixed E/B finite differences, 3 frames agree with analytic omega=0 gamma model to <2e-9; spinless TR reciprocity | **VERIFIED same phenomenological model, NOT conserving microscopic oracle** |
| Original full Xiang orbital E/B transport | No conserving bath, complete Kubo-Streda, infinite-system limit, realistic magnetic Wannier embedding or source-equation resolution | **UNRESOLVED; paper audit PARTIAL** |

40 new v66 science gates plus one numerical frozen-regeneration test (total standalone v60-v66 219). See [scientific v66 report](DYNAMIC_ORBITAL_EB_SEVEN_CONTACT_V66.md). The frozen GitHub JSON is acceptance-compared *numerically* across supported Python/JSON serializers; no unverified byte-identical claim.


## v67: physically independent elastic reservoir Ward benchmark, NOT original orbital Hall

| Target | Evidence | Acceptance |
| --- | --- | --- |
| Genuine microscopic semi-infinite, noninteracting tight-binding leads | Retarded surface Dyson root, energy-dependent exact selfenergies, 2-channel Caroli and independent S†S=I | **VERIFIED specified conserving elastic scattering** |
| Magnetic current/contact Ward and gauge covariance | Both Sigma_B/Gamma_B, two-channel d_B(S†S)=0, d_B optical spectral identity, d_B T from scattering amplitudes vs Caroli, uncontacted local Q(B) continuity | **VERIFIED finite-model scattering and magnetic connections** |
| Electrostatic E gate, projected orbital B, reservoir voltage | Full d_B H, d_B P_b, 8 mixed differentiated Caroli terms and independent 16-corner d_B d_E gate finite differences, 3 frames | **VERIFIED named gate/bias model, CONDITIONAL physical embedding** |
| Wrong contact transformation or omitted d_B coupling | Deliberate gauge-breaking negative controls, origin shift corrected only by moving lead onsite and chemical potential | **NEGATIVE tests passed** |
| Original homogeneous electric-field-driven Xiang orbital mixed E/B | Landauer gate derivative is a **DIFFERENT observable**; no same-observable Keldysh/SCBA or field-complete Wannier | **UNRESOLVED; original source PARTIAL** |

v67 has 45 new physical + 1 frozen tests, standalone v60–v67 **265/265 local PASS**. [Detailed v67 report](CONSERVING_RESERVOIR_WARD_V67.md).

## v68: same finite device polarization-current operator, microscopic elastic-lead Keldysh AC E/B

| Named claim | Independent gate | Status |
| --- | --- | --- |
| Finite-frequency, device-local E, real energy-dependent semi-infinite reservoir lesser Keldysh | Physical lead band surface branches, positive/complete spectral weight, 512-to-640 resonance-resolved energy integration, exact two-time lesser Dyson identity and projected J(E)=i[H-E Pb,Pa] | **VERIFIED bounded noninteracting finite-model identity** |
| Same *device current operator* as v66, different bath | Both models use supplied Pa,Pb and explicit electric current vertex -i[Pb,Pa]; v68 bath from semi-infinite noninteracting leads, v66 scalar gamma; their numerical susceptibilities differ | **VERIFIED comparison; NOT bulk equivalent** |
| Microscopic lead charge Ward at finite omega | Construct lead collision from independently evaluated two-energy Sigma-r/a and Sigma-less, not from residual; full matrix continuity residual ~9e-11, uncontacted site ~1e-11; global lead injection balances device charge accumulation | **VERIFIED current/charge continuity, specified model** |
| Magnetic orbital and lead contact completeness | Full dB H,Pa,Pb,Sigma,Gamma and both energies; three B-dependent frames agree; analytic vs independent five-point full B numerical derivative ~1e-11 | **VERIFIED named finite operator derivative** |
| Strict DC alternative finite-E/B oracle | Recompute full equilibrium lead lesser density with H-E Pb and J(E) at each gate, compare 3-frame finite E derivatives and four-corner mixed E/B to omega=0 Keldysh | **VERIFIED same physical elastic-lead bath, independent full-field computation** |
| Full original Xiang bulk intrinsic orbital mixed E/B | Device-only field, boundary response, finite arbitrary off-diagonal projected dipoles, no real-material magnetic Wannier, no bulk thermodynamic normalization or original SM Eq49/52 three-band closure | **UNRESOLVED; source paper PARTIAL** |

40 v68 independent physics tests + 1 frozen-data reprovenance gate, total standalone v60–v68 306 PASS. See [v68 physics report](FINITE_FREQUENCY_SAME_DEVICE_KELDYSH_EB_V68.md).


## v69: finite-lead Lehmann independent oracle

| Gate | Evidence | Status |
| --- | --- | --- |
| Full device+two-lead Hermitian Hamiltonian | Direct eigendecomposition; no surface self-energy/Keldysh integral in the response | **VERIFIED finite model** |
| Same v68 device polarization-current observable | Identical embedded J0 and noncommuting electric contact JE=-i[Pb,Pa]; device-local field only | **VERIFIED observable matching at finite-device protocol level** |
| Finite closed-system Ward | (omega+i eta)rhoE-[H,rhoE]+[Pb,rho0]=0, residual 7.8e-17; induced device charge exactly balanced by explicit lead charge | **VERIFIED finite Hamiltonian** |
| Static equilibrium response | Lehmann Frechet including diagonal -f'Pnn agrees with full finite-E rediagonalization to 4.2e-13 | **VERIFIED isothermal finite system** |
| Ordered finite lead/broadening limit | (L,eta)=(250,.02),(350,.02),(500,.01),(1000,.005) discrepancy to v68 falls 9.27e-4 -> 2.12e-4 | **IMPROVING, NOT YET CONVERGED** |
| Original Xiang bulk intrinsic orbital mixed E/B | Source equation sign/factor/three-band issues, real Wannier magnetic embedding and physical bulk orbital-current correspondence remain open | **UNRESOLVED; source PARTIAL** |

45 v69 physics gates plus one quick frozen-data regeneration gate. The high-resolution L=1000 point is independently regenerated in the main v69 test class. See [v69 report](FINITE_LEAD_LEHMANN_ORACLE_V69.md).


## v70: multi-frequency ordered finite-lead/continuum benchmark

| Gate | Acceptance | Scientific state |
| --- | --- | --- |
| Same finite-device observable at ω=.20/.29/.40 | Independent full-Hermitian eigenstate Lehmann versus two-energy semi-infinite-reservoir Keldysh, shared electric-current contact and device-only projected current | **VERIFIED specified model** |
| Joint length and η path | L250→350 (η=.02)→500 (η=.01)→1000 (η=.005); all three production errors smaller than low-L errors | **IMPROVING; NONZERO DISCREPANCY** |
| Nonmonotonic finite-spectrum negative control | At ω=.20, error rises 7.77e-4→7.93e-4 for L250→350 | **VERIFIED no pointwise monotonic guarantee** |
| Far-end boundary-cap sensitivity | ω=.29, L500 η=.01, cap 0→±1.1 changes complex response ≈2.18e-5 | **VERIFIED finite-size artifact, not universal independence** |
| Bulk orbital current / original Xiang PRB | Current defined as device-only i[H−EP_b,P_a], no general bulk q→0 mapping, full Wannier B embedding or closure of printed source sign/factor/three-band sectors | **UNRESOLVED; source PARTIAL** |

35 new v70 numerical checks. See [multi-frequency report](MULTIFREQUENCY_ORDERED_LIMITS_V70.md).

## v71: current-observable identity and literal source branch

| Scientific target | Evidence | Status |
| --- | --- | --- |
| Device current vs full device-supported polarization derivative | Exact Jfull = Jdevice + i[Hcontact,Pdevice], nonzero contact; independently decomposed optical Lehmann susceptibility and Ward | **VERIFIED finite bounded operator** |
| Finite equilibrium isothermal E derivative | Jdevice +0.002996376694233, lead-contact −0.002996376694233, Jfull=0; independent finite-E/B rediagonalization and 3 B-dependent frames | **VERIFIED finite-equilibrium contact cancellation; NOT biased DC transport** |
| Printed SM Eq20/31 spin PHE signs | Direct interband velocity contraction equals inferred +tau sum Z df, while literal printed Eq31 has opposite sign, independent Hermitian/U1/band relabel controls | **VERIFIED literal source-index discrepancy; author intent OPEN** |
| Main Eq10 vs SM Eq62 spin helicity | H transforms with sigma_z only if physical Zeeman sigma_x also changes sign | **VERIFIED conditional basis covariance; original printed branch OPEN** |
| PHE/MNHE capability inputs | Separate spin Zeeman versus orbital Wannier moment/full positional shift; missing inputs rejected, presence alone never establishes reproduced physics | **VERIFIED input inventory only, source response UNRESOLVED** |
| Original paper bulk orbital tensor and Fig2 SI | Finite device is not bulk; Eq49/52 factor issue, m!=ell three-band omissions and full magnetic Wannier embedding still missing | **UNRESOLVED; original paper PARTIAL** |

42 v71 named physical checks plus one frozen numerical regeneration gate; exact PR-head and postmerge CI require independent verification. [v71 details](SOURCE_OBSERVABLE_BOUNDARY_AND_WANNIER_GATE_V71.md).

## v72: periodic finite-q spin Zeeman charge-Ward method, NOT source dc transport

| Gate | Method | Status |
| --- | --- | --- |
| Periodic finite-q scalar source and k-dependent frame | Full density U(k+q)†U(k), longitudinal J=[H(k+q)-H(k)]/q; commensurate BZ global occupation cancellation | **VERIFIED finite-q lattice model** |
| Exact longitudinal charge and magnetic derivative Ward | (omega+i eta)chi_density-q chi_Jy=0 and dB identity (residual O(1e-18)) | **VERIFIED defined finite-q current** |
| Mixed spin Zeeman transverse response | H(k,B)=H0(k)+B sigma_x; 4x4 orbital Sylvester independently matches eigenstate sum, analytic Fermi-Fréchet+Liouvillian dB matches 5-point full B difference | **VERIFIED declared finite-q model** |
| Gauge/vertex negative controls | Replacing U+†U- by I violates Ward ~8.24e-4; endpoint-averaged vy instead of divided-difference fails ~1.66e-5 | **REJECT incorrect finite-q Ward methods** |
| Original source spin PHE and orbital MNHE | Selected eta_a adiabatic pole, both E/B orderings, q/omega/eta to dc, magnetic moment and intercell Wannier position, transport-magnetization subtraction and SM disputed equations not implemented | **UNRESOLVED; original source PARTIAL** |

37 independently motivated physics checks plus one frozen numeric provenance gate. [Full report](BULK_FINITE_Q_SPIN_ZEEMAN_V72.md).
