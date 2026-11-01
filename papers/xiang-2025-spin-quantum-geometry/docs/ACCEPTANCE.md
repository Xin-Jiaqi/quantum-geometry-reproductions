## Current core-access decision — 2026-11-01

Model-level spin-current and spin-polarization capabilities are callable via
`methods.response`, with both linear and quadratic finite-frequency fields,
all intermediate bands, the full physical velocity vertex and an explicit
position connection. 76 tests pass on Python 3.11/3.12. The independent
four-source SM(2) audit covers 2/4/6 bands. At eta=.02 the nested k derivatives
need a converged step: .004 gives <2e-7 absolute source residual; taking an
arbitrarily tiny step amplifies eigenproblem roundoff and is not improvement.

`results/CORE_RESPONSE.json` is a nonzero four-band BZ benchmark. Its 25² to
41² mesh difference is 1.091e-8 and .003 to .006 derivative-step difference
is 1.430e-11. It is full finite-frequency response, not tau^0 extraction.
A shifted-inversion-symmetric initial toy gave zero and was rejected as a
positive transport benchmark; an independent sin(kx) hopping breaks that
symmetry in the saved example.

Source recheck: SM p12 Eq.(55) integrates Delta² over [0,2pi] but prints
Delta² instead of 2pi*Delta² in its last line. Existing printed/corrected
branches remain separate. SM p6 explicitly drops terms from Eq.(26) for
main Eq.(13), and p9 drops three-band processes for main Eq.(14). The core
access decision does not require those approximations to become universal.
The physical metallic tau^0 preparation and original absolute axes remain
open; neither prevents use within the validated named protocols.

---

# Scientific acceptance, 2026-10-20 (curated V28 subset)

| Item | Gate | Meaning |
| --- | --- | --- |
| Main Fig. 1(a), Eq. 9 NDSC | PARTIAL | Two valleys, local spin-QGT and shape checked; paper ordinate remains approximately twice independent SI calculation. |
| Main Fig. 1(a), Eq. 13 intrinsic | PASS reduced mathematical kernel / PARTIAL physics | V29 frozen two-band sea+surface, independent eigenbasis matrix-element check, cutoff/mesh convergence. **Not full tau^0 driven DC**; factor ~2 remains open. |
| Main Fig. 1(b) | PASS for printed closed-form curves only | SM Eq. 67/main text/SM Eq. 70 agree with original PDF vectors offline; printed Drude sign and gap formula disagree across source sections. |
| Main Fig. 1(b) independent response | OPEN | SM Eq. 70 literal integral and source density-matrix response disagree with published closed form; physical DC protocol remains unresolved. |
| SM Fig. 1(a), (c), (d) | PASS for analytic kernels; PARTIAL published maps | Band and metric/curvature kernels computed; colorbar maps lack absolute calibration. |
| SM Fig. 1(b) | PASS published expression; OPEN integral | Printed Eq. 56 disagrees with full polar integral of Eq. 55. |
| General finite-eta method | OFFLINE validated subset | Not included here. Not a validated intrinsic tau^0 tool. |
| Real material / stacking | OUT OF SCOPE / DEFERRED | Not represented as completed paper calculations. |

**Open source algebra:** SM Eq. (2) to Eq. (6) metallic finite-eta expansion; Eq. (43)–(45) covariant spin-matrix derivatives; SM Eq. (69)–(70) closed integral; SM Eq. (55) angular integration. Keep these as questions rather than automatic error corrections.

**No-fit rule:** source-printed formulas and independently derived expressions remain separately callable with explicit names. CI only proves the specific frozen model tests. Full V27 archive remains locally frozen and is not part of this GitHub upload.

## V29 extension

V29 migrates the frozen intrinsic Eq. (13) two-band reduced model. The eigenstate matrix oracle validates local Q and g_yy terms independently; the Fermi-pocket and Fermi-shell integrations agree with offline V27 frozen reference at μ=0,12,20,40 meV. This does not resolve the Fig. 1(a) ordinate normalization, finite-eta/DC interpretation, or omitted four-band terms. Multiband and PT solvers remain unported.

## V31 isolated-cluster Laurent audit

PASS for 46 scoped unit tests (32 old + 14 V31), with 3/4/5 nondegenerate bands, exact PT four/six bands, an independent dH population-gradient derivative, true three-cluster intermediate term, U(6) fixed-orbital covariance, finite-eta O(eta) remainder and a converged one-dimensional BZ diagnostic. The adiabatic finite constant is **NOT** established as a universal intrinsic metallic driven-DC tau^0. No source equations were silently corrected. See LAURENT_ACCEPTANCE_V31.md.


## V32 — SM Eq.(2) four-source algebra

PASS: 70 scoped tests on nondegenerate source Eq.(2), independent finite-eta orbital recursion (N=2..5), gauge choices, printed Eq.(6) algebra and frozen original V14 values. Max matrix difference < 8e-9 for the four original band-versus-orbital comparisons. Printed Eq.(6) has an additional source-equation occupation-gradient difference: deltaC=-i*r_mn*(f'_n-f'_m)/(epsilon_m-epsilon_n)^2. **Scientific status stays PARTIAL; no metallic intrinsic physical DC gate is closed.** See SM_EQ2_ACCEPTANCE_V32.md.
