# v52 | Direct publisher Fig.2(b) / Eq.(11) / claimed 43 microV: absolute SI audit

**Source**: Xiang et al., *Quantum geometric map of magnetotransport*, *Physical Review B* **113**, L201406 (2026), genuine user-provided published 8-page PDF, page **L201406-4**; original user-provided 16-page SM. Main PDF SHA256: b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730. SM SHA256: 9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4. Original copyrighted PDFs are NOT uploaded. **Entire-paper scientific acceptance PARTIAL.**

## 1. Exactly what the publisher prints

Published page 4 Fig. 2(b) has vertical axis

\[
-10^4\,\sigma^{(1)}_{xy,x}\,g\mu_B B_x\quad(e^2/\hbar)
\]

with a **positive T=0 plateau approximately 1.75** at positive chemical potential. The figure caption states beta=1/8, Bx=1 T, tau=0.5 ps, g=10. The T=20 K green and T=40 K blue traces exhibit thermally smoothed steps. No figure-generation data or code was supplied by the authors.

Main Eq. (11), after its displayed restoration of e and hbar, gives the signed coefficient for positive mu:

\[
\sigma_{xy,x}^{(1)}
=\frac{\tau e^2}{\hbar^2}\,
F(\beta),\qquad
F(\beta)=\frac{\beta^2+2\sqrt{1-\beta^2}-2}{4\pi\beta^3}
=-\frac{\beta}{4\pi[1+\sqrt{1-\beta^2}]^2}.
\]

The physical spin field contributes energy g mu_B Bx. The authors additionally state, **on the same printed page**, Ey=10^4 V/m, 100 micrometer lateral size, 10^3 ohm Hall resistance and claimed voltage about **43 microV**. The dimensional Hall-bar estimate is |delta_sigma_xy| Ey width R, with delta_sigma a 2D sheet conductance (Siemens). This is an estimate, not a directly measured voltage.

## 2. Direct SI calculation, no free prefactors

Use CODATA physical constants as returned by scipy.constants, with proper joule and second units, g=10, tau=0.5 ps, Bx=1 T:

- F(1/8) = **-0.00250641590669**
- g mu_B Bx = **5.788381798e-4 eV**
- tau*g*mu_B*Bx/hbar = **0.43970500209**
- e^2/h = **3.87404586493e-5 S**, e^2/hbar = **2.43413480579e-4 S**
- Ratio (e^2/hbar)/(e^2/h) is **2pi EXACTLY**.

Literal Eq11:

\[
\frac{\Delta\sigma}{e^2/\hbar}
=\frac{\tau g\mu_B B_x}{\hbar}F(\beta)
=-0.00110208361150,
\quad |\Delta\sigma|=2.682620078\times 10^{-7}{\rm S}.
\]

Therefore the Fig2(b) ordinate **should be +11.0208361**, if both the printed Eq11 and the printed axis e^2/hbar are interpreted literally.

Since Ey*width = (10^4 V/m)(100 micro m) = 1 V, the source Hall-bar estimate would be **268.2620 microV**, not approximately 43 microV, under this literal equation-to-axis normalization.

Now evaluate a strictly **separate** alternative conductance-unit convention: preserve the same beta, tau, field and dimensionless F, but multiply by e^2/h instead of e^2/hbar before expressing the result in the figure's printed e^2/hbar units.

\[
\Delta\sigma_{\rm alt}/(e^2/\hbar)
=[\tau g\mu_B B_x/\hbar]F(\beta)/(2\pi).
\]

This gives **ordinate +1.754020545**, sheet conductance **4.269522458e-8 S** and **Hall voltage 42.6952 microV**.

| Independent reference | Fig.2(b) T=0 ordinate | Hall voltage |
| --- | ---: | ---: |
| Printed Eq11, literal e^2/hbar | +11.020836 | 268.262 microV |
| Separate alternative e^2/h conductance normalization | +1.754021 | 42.695 microV |
| **Published** plot readout and paragraph | approximately +1.75 | approximately 43 microV |

**Two independent publisher claims (printed figure amplitude AND claimed 43 microV) coincide with the e^2/h alternative, smaller by 1/(2pi) than literal Eq11/e^2/hbar.** This is a high-confidence **original-source absolute-normalization tension**. It is NOT yet proof which code or convention the authors intended, and no silent division by 2pi is accepted as a manuscript correction. The literal vs approximately read publisher figure ratio is 6.298 rather than exactly 2pi because the graph has only ~two-digit readout; the exact comparison of independent conductance quanta is 2pi.

## 3. Independent original finite-temperature curve digitization

The publisher PDF page4 was rasterized at 8x and cropped in PDF point coordinates (140,55,285,215), with NO OCR and no image recreation. Coordinate calibration: zero chemical potential x=718.5; mu ±40 meV ticks x=379,1059; ordinate zero y=418.5; each ordinate unit 170.5 pixels. Digitized source green T=20 K and blue T=40 K at mu=±5 and ±10 meV. This is approximate image sampling (ordinate precision about ±0.03), not author data.

An independently formulated **finite-temperature model continuation** uses the predicted e^2/h alternative amplitude 1.754020545 * tanh[mu/(2k_BT)] with physical CODATA kB. **Eq11 itself is given for T=0; this finite-T function is not claimed to appear as an author-supplied analytic source.** No additive offset, temperature rescaling or amplitude fit was made.

| T | mu (meV) | Original figure digitized | Model no-fit prediction |
| ---: | ---: | ---: | ---: |
| 20 K | -10 | -1.739 | -1.743 |
| 20 K | -5 | -1.551 | -1.571 |
| 20 K | +5 | +1.569 | +1.571 |
| 20 K | +10 | +1.739 | +1.743 |
| 40 K | -10 | -1.563 | -1.571 |
| 40 K | -5 | -1.070 | -1.088 |
| 40 K | +5 | +1.094 | +1.088 |
| 40 K | +10 | +1.575 | +1.571 |

Maximum pixel-digitized ordinate minus continuation **0.01994**, within the declared ±0.03 figure-reading resolution. This corroborates that the published finite-temperature traces also have the smaller absolute-amplitude scale; it cannot certify the author's finite-T integration algorithm.

## 4. Independent physics and numerical gates

Thirteen new regressions: Eq11 original vs rationalized geometry; SI e^2/hbar vs e^2/h EXACT 2pi; independently encoded CODATA Bohr magneton in J/T and eV/T; literal Eq11 ordinate; alternative original figure amplitude and 43 microV using Ey width R; beta, tau, g and geometry parameter scaling; source digitized green/blue thermal branches; symmetry of the finite-T model; explicit rejection of invalid physics inputs; and strict source provenance plus PARTIAL scientific acceptance. Some independently tabulated CODATA eV/T constants have fewer published digits than their J/T counterpart, leading to ~1e-8 harmless representation differences: this is bounded in tests, not a 2pi source.

## 5. Unresolved science; no fabricated source repair

**OPEN**: whether authors meant e^2/h or e^2/hbar on the printed y axis; how the original Eq20/31 opposite signed-current branches (v49-v51) were selected in generating Fig2b; original figure-generation code is unavailable; multiband orbital response, physical steady-state field-preparation/bath/contact terms. We do not contact authors or claim their measured voltage is incorrect. The source's quoted 43 microV is an *estimate*, and source Eq11/SM signs must remain separately tracked.

**Whole original paper remains PARTIAL.** Proposed v53: source orbital multiband Eq49/52/56 and 3-band terms, prioritizing checks that can distinguish a genuine physical contribution from an index/sign or contact convention. Do not create another unrelated Hall benchmark.

### Reproduce from isolated PRB subproject

    python -m unittest tests.test_fig2_si_v52 -v
    python scripts/reproduce_fig2_si_v52.py --beta .125 --tau-ps .5 --g 10 --bx-tesla 1 --ey-vm 10000 --width-um 100 --resistance-ohm 1000

You can pass --main-pdf and --sm-pdf to verify original private byte hashes. The original PDF image is not redistributed. Full-precision model values and original axis pixel sample coordinates reside in results/frozen/v52_original_Fig2b_SI_audit.json.
