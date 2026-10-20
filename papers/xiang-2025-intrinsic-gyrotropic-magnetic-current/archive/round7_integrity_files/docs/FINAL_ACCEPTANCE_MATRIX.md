# Complete status matrix after seven research rounds

**Overall:** `PARTIAL` / validated model-method subset; **never** 100% reproduction of published panels. Source: Xiang et al., PRL 134, 116301 (2025), all 7 main and 11 SM pages read. No author source figure data or code received.

Legend: `PASS` independently verified named physical/model *subset*; `PARTIAL` incomplete absolute curve or secondary method; `BLOCKED` missing external provenance; `NOT STARTED` not implemented.

| Target | Current gate | Evidence | Missing gate |
|---|---|---|---|
| Main Eqs. (1)–(6) geometry | PASS model algebra | spectral vs independent d-vector Pauli traces; 12 component checks | handling band degeneracy beyond nondegenerate two-band cases |
| Main Eqs. (7)–(8) Liouville/Kubo | PASS model | independent orbital-basis Liouville; spectral Kubo finite-ω; 12 TR complex cases | material dissipative collisions not claimed |
| Main Eq. (9)–(13) C/D split | PASS nonresonant 2-band protocol | Wilson Berry-phase, d-vector, Liouville, Fermi line Eq.13, 1 THz frequency | author Fig.2(d) plotted amplitude remains separate question |
| Main Eq. (14) magnetic antiunitary | PARTIAL | all 122 MPG channel-resolved Neumann/character + complex-ω TR | paper-printed single antiunitary factor inconsistent with naive real C channel; original author's frequency-response convention unconfirmed |
| Fig.1(a) bands | PASS model | analytic Rashba two branches, spectral tests | exact layout only |
| Fig.1(b) conduction | PASS model | Eq.16 analytic vs numerical + independent SI prefactor | author original curve data unavailable |
| Fig.1(c–d) k-resolved Z | PARTIAL | analytic angular texture and symmetry | singular k→0 error/uniform plot masking not exhaustively validated |
| Fig.2(a) lattice Dirac bands | PASS model | exact Γ gap, vector spectrum | author plot digitization unavailable |
| Fig.2(b) metallic C_xy | PASS model / PARTIAL author curve | full-domain boundary partition + Eq13 delta Fermi line agree ≤5.35e−10 Å⁻¹ across sample µ; sample peaks | entire original plot digitization/source unavailable |
| Fig.2(c) weighted and bare Q | PARTIAL | Pauli d-vector, Γ-weighted hotspot and raw Q | no quantitative original k-map comparison |
| Fig.2(d) insulating and metallic D_xx | PASS Eq12 method / BLOCKED original amplitude | 4 independent formulas/numerics, in-gap 0.01775389 and µ=±1 −0.06417845 A/(mT), converged metal | author graph amplitude visually about ×2; cannot close until input/code clarified |
| SM §I distance expansion | PARTIAL | definitions and gauge identities inspected | fully independent Fubini–Study subtraction proof not a separate release gate |
| SM §II Liouville | PASS selected magnetic-response branch | direct first-order Liouville, full Kubo | mixed electric perturbation not in this paper's final IGMC numerical targets |
| SM §III semiclassical | PARTIAL | Wilson wavepacket-position proxy validates D | full original wavepacket Lagrangian derivation not independently rederived line-for-line; SM Eq34 normalization issue |
| SM §IV Green function | PARTIAL | Kubo/finite-ω response agreement | independent lesser Green's contour method not coded |
| SM §V Tables I–III | PASS 243/244 source tensor spaces; PARTIAL total | 122 MPG generator closure, 244 SVD/character, full source transcription | single SM Table II `−42′m′` missing equality |
| SM §VI induced electromagnetic E | PARTIAL | source frequency and ideal sheet geometry checked | realistic finite sample and Faraday-field quantitative model not developed |
| Real material and Wannier pipeline | BLOCKED | prior NMEE package contract exists; not an IGMC real-material calculation | same-gauge projected spin/position operators and DFT data absent |
| Reproducible code | PASS offline | parametrized reproducibility, frozen numeric baselines + hashes, test suite | GitHub remote folder and CI not yet published |
| Publisher source plot digitization/data | BLOCKED | PDF visual estimates explicitly not taken as ground truth | author-supplied raw data, plotting scripts or transparent normalization |

**Current scientific decision:** keep the physically consistent main-text branch as the published-formula implementation; keep SM Eq. (34) literal ×2 as a named *invalidated for physical acceptance* diagnostic; never tune g, omega, cell length or numerical broadening to force plot agreement. User must authorize author correspondence or remote repository push separately.
