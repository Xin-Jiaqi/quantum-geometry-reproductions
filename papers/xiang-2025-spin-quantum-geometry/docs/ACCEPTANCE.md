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
