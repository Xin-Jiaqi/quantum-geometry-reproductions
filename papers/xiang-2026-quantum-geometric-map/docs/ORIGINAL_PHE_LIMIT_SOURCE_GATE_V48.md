> **Source status superseded by v49:** The user has now supplied the original 16-page SM and 8-page published main PDF, both directly inspected and SHA256-pinned. Earlier paragraphs claiming SM is unavailable accurately describe v48 only, NOT current access. See [v49 primary-source audit](ORIGINAL_SM_SOURCE_AUDIT_V49.md) for source Eq15/20/31 sign and Eq10/62 helicity discrepancies.

# v48 | Original symmetric spin PHE: robust RTA limit and source-provenance audit

**Paper:** Xiang et al., *Quantum geometric map of magnetotransport*, PRB **113**, L201406 (2026), DOI [10.1103/kd12-4ycj](https://doi.org/10.1103/kd12-4ycj). Entire-paper scientific acceptance remains **PARTIAL**.

## Main paper actually accessed vs missing source

**Directly checked:** the authors' public seven-page [arXiv:2510.02661v2](https://arxiv.org/pdf/2510.02661). Main Eqs. (2)–(6) distinguish symmetric PHE from antisymmetric MNHE/OHE, Eqs. (7)–(8) specify Zeeman geometry, Eq. (10) fixes the tilted Dirac cone, and Eq. (11) gives its zero-temperature step-like spin PHE. The paragraph before Eqs. (7)–(8) uses an iterated Liouville density matrix. The authors explicitly call the relaxation-time approximation **phenomenological** in their final discussion.

**Not retrieved:** the original Supplemental Material. The APS Supplemental Material pages and the arXiv source archive were inaccessible from this runtime, and there is no mounted original SM. Archived Eq.-by-Eq reconstructions in the GitHub paper module are *not* substitutes for source inspection. It is therefore **incorrect to claim** a source-certified RTA denominator, collision insertion, or exact field-preparation equivalence in v48.

The independently constructed model remains the v47 **frozen zero-field collision target** \`f(H0)\` with a sign-corrected electric source, *both* E/B orderings, and the **symmetric** PHE projection \`(xy+yx)/2\`, not the old retracted v46 unsymmetrized XY result.

## Stable original Eq. (11) geometric coefficient

For model \`e=hbar=vF=1\`, fixed \`beta=.125, mu=.3, T=.03\`, 140 radial ×160 angular points, UV=2.5, compare \`eta*sigma_symmetric\` to the signed finite-T low-energy continuation
\[
F_{11}(\beta,\mu,T)=
-\frac{\beta}{4\pi[1+\sqrt{1-\beta^2}]^2}
\tanh\left(\frac{\mu}{2T}\right).
\]
The paper's *main Eq. (11)* is explicitly presented at zero temperature; the finite-temperature tanh continuation in this test is an independently defined model check, **not** an SM-derived assertion.

| eta | eta*sigmaS / F11 | eta*sigmaA / eta² |
| ---: | ---: | ---: |
| 0.016 | 1.0032150358 | −0.01600496 |
| 0.008 | 1.0008104803 | −0.01606647 |
| 0.004 | 1.0002038834 | −0.01614689 |
| 0.002 | 1.0000512677 | −0.01622829 |
| 0.001 | 1.0000128853 | −0.01633109 |

For small eta the relative symmetric discrepancy and the eta-scaled antisymmetric remainder both shrink approximately as **eta²**. Fits of the ratio against eta² on the three/four/five smallest eta points give intercepts **1.0000002435**, **1.0000007042**, **1.0000025998** respectively. **No fitting of 2π, overall sign, or additive response.**

## Extra physical gates

- **Electron/hole reversal:** \`mu=±.15,±.3\` at \`T=.02,.055\`, eta=.002, gives opposite symmetric response to within numerical quadrature. The Eq11-tanh ratios range from **1.00004638** (|mu|=.3,T=.02) to **1.00068813** (|mu|=.15,T=.055); accuracy decreases as the Fermi window approaches the Dirac point.
- **Thermal smoothing:** the |mu|=.15 response is smaller at higher T; the tanh continuation is genuinely checked rather than copied into the output.
- **UV cutoff:** cutoffs \`1,1.5,2.5,4\` have ratios \`1.0000512796,1.0000512618,1.0000512677,1.0000512889\`, spread below **3×10^-8**. This is a continuum model control, *not* a lattice completion.
- The complete Fermi matrix source is still independently finite-differenced; the erroneous v46 analytic derivative and its archived beta scan remain **RETRACTED**.

## Reproduce

From this paper's directory:

~~~bash
python -m unittest tests.test_phe_limit_audit_v48 -v
python scripts/reproduce_phe_limit_v48.py --mode all --output results/frozen/v48_phe_source_limit.json
~~~

Ten new physics tests validate eta² scaling/extrapolation, antisymmetric residue, mu odd parity, thermal smoothing, UV stability, independent source differences, correct xy/yx symmetry, input rejection and explicit source noncertification.

**Still OPEN:** author-specific SM collision/relaxation and both perturbation orderings; literal main Eq. (8) versus main Eq. (11) signed-branch audit; Fig. 2(b) published absolute SI amplitude (~2π discrepancy); remaining multiband orbital response discrepancies. The separate v41–v44 QWZ SCBA/Středa **antisymmetric charge Hall** remains an extension, not a prerequisite for the original spin PHE. Next gate is *access the actual original SM and implement its formula unchanged*. Do not contact authors without user request.
