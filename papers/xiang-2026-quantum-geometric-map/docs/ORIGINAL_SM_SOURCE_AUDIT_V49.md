# v49 | Primary-source scientific audit: published PRB and original 16-page SM

**Scientific state: PARTIAL.** Both authentic user-supplied PDFs have now been directly checked; v45–v48 descriptions of the Supplemental Material as inaccessible were historical, now superseded. **The original copyrighted PDF files are not uploaded to this public repository.**

| Source | Pages | SHA256 |
| --- | ---: | --- |
| Published PRB 113 L201406 main PDF | 8 | b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730 |
| Supplemental Material to Quantum geometric map of magnetotransport | 16 | 9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4 |

## 1. What the source *actually* does with relaxation

**SM pp.3–4 Eqs.(4)–(10)** use two separately regulated time-dependent external fields with frequencies omega_a and omega_b, positive adiabatic infinitesimals eta_a and eta_b, eta_Sigma=eta_a+eta_b, and magnetic coupling Bbar=-g*mu_B*B. Both E/B perturbation orderings are included, with Bloch-band covariant derivatives and off-diagonal Berry connections.

**SM p.5 Eq.(15)** introduces the phenomenological tau by regularizing a **selected electric-field 1/eta_a pole**, **not** by applying a uniform finite relaxation eta to every matrix-element or frequency channel. In particular, **SM Eq.(20)** for m != n reads:

\[
\rho^{(EB;5)}_{mn}=-\tau\frac{\sigma^c_{mn}\partial_b f_{nm}}{\varepsilon_{mn}}\,E^b\bar B^c,\qquad f_{nm}=f_n-f_m.
\]

**SM p.8 Eqs.(31)–(33)** contract the density matrix with velocity, define the Zeeman curvature, and symmetrize the resulting PHE. The independently constructed **frozen-H0 finite-eta all-sector RTA** used in v47/v48 is **not source-identical** to this adiabatic pole prescription, even though their leading tau coefficients can agree.

## 2. Literal SM Eq.(20) -> Eq.(31) algebra has an opposite overall sign

The source prints \(v^a_{nm}=i\varepsilon_{nm}r^a_{nm}\), where \(\varepsilon_{nm}=-\varepsilon_{mn}\), and defines \(Z^{ac}_{nm}=i(r^a_{nm}\sigma^c_{mn}-r^a_{mn}\sigma^c_{nm})\). Thus contracting **the literal Eq.(20)** gives:

\[
\begin{aligned}
\frac{j_a^{(EB;5)}}{E^b\bar B^c}
&=-\tau\sum_{m\ne n}\frac{\sigma^c_{mn}v^a_{nm}}{\varepsilon_{mn}}\partial_b f_{nm}\\
&=+i\tau\sum_{m\ne n}\sigma^c_{mn}r^a_{nm}(\partial_b f_n-\partial_b f_m)\\
&=+\tau\sum_{m\ne n}Z^{ac}_{nm}\partial_b f_n.
\end{aligned}
\]

**But the printed SM Eq.(31) says MINUS tau times the same sum.** The mismatch appears already in the first line of Eq.(31), where minus i*sigma*r appears instead of the plus sign implied by Eq.(20) and the printed velocity identity. Since Bbar=-gmuB*B multiplies both with the same minus, changing to physical B does not remove their relative sign.

This is a **literal source-algebra inconsistency**, not automatic proof that a measured phenomenon or the complete original paper is wrong. Both intermediate routes are reported; no unproved repair is applied.

## 3. Main and SM print opposite spin helicity

**PRB main p.4 Eq.(10)** prints \(H_{\rm main}=t_x k_xI+v_F(k_x\sigma_y-k_y\sigma_x)\).

**SM p.13 Eq.(62)** prints \(H_{\rm SM}=t_x k_xI+v_F(k_y\sigma_x-k_x\sigma_y)\).

For the same physical spin Zeeman sigma_x held fixed, these have opposite helicity and opposite signed local spin PHE. Indeed \(\sigma_z H_{\rm main}\sigma_z=H_{\rm SM}\), but **the same spin rotation also sends sigma_x to MINUS sigma_x**. A basis equivalence only applies when all spin and magnetic perturbation operators are transformed consistently; it is invalid to silently keep them both unchanged.

Independently diagonalized two-band eigenstate integrations, beta=1/8, mu=.3, T=.03, 100 radial by 128 angular Gauss mesh, UV=2.5, units e=hbar=vF=tau=gmuB=1:

| Fixed spin/B convention | Direct Eq20 physical-B XY | Printed Eq31 physical-B XY | Printed Eq32 symmetric | Main Eq11 signed geometry (model tanh factor) |
| --- | ---: | ---: | ---: | ---: |
| **Main Eq10 helicity** | −0.002506188247 | +0.002506188247 | **+0.002506188497** | −0.002506188335 |
| **SM Eq62 helicity** | +0.002506188247 | −0.002506188247 | **−0.002506188497** | −0.002506188335 |

The **printed Eq32 evaluated with SM Eq62** agrees with signed main Eq11/SM Eq63. With the **main Eq10** model it gives the opposite sign. At eta=.002, the earlier v47 full symmetric frozen-H0 response gives eta*sigmaS=**−0.002506316935**, close to **literal Eq20+main Eq10** (−0.002506188497) and opposite **Eq32+main Eq10** (+0.002506188497). This identifies the source/response branch involved in the former numerical agreement.

The tanh(mu/2T) used in this numeric comparison is a *model-derived finite-temperature continuation*; the authors' Eq11/Eq63 is stated at **T=0**.

## 4. Acceptance and remaining source-level physics

**PASS:** both original PDFs independently read and cryptographically pinned; correct SM tau-pole structure; direct Eq20 and printed Eq31/32 independent eigenbasis calculations with both xy/yx channels, both helicities, electron/hole reversal, gauge phase invariance and multiple tilts; direct v47 small-eta limit relative to literal Eq20; original Eq63 analytic expression; Eq15 exact partial-fraction identity. **13 new physical tests.**

**OPEN:** a fully consistent signed mapping of **all** source equations SM Eq20/31/62 and main Eq10/11; finite-frequency source Liouville including eta_a,eta_b, both field orderings and all contact terms; previously identified published Fig2(b) absolute SI ordinate scale about 2*pi different from the literal Eq11 conversion; original orbital Eq49/52/56 and omitted three-band processes; real-material predictions. Software CI PASS is never whole-paper scientific PASS.

Reproduce from the PRB paper directory by running:

    python -m unittest tests.test_sm_literal_spin_phe_v49 -v
    python scripts/reproduce_sm_source_v49.py --beta .125 --mu .3 --temperature .03

To independently verify the exact source PDF bytes, pass both optional main-pdf and sm-pdf flags with the original uploaded PDF paths. **v50:** evaluate the original two-field complex-pole Liouville directly with the source index conventions, then distinguish which sign conflict is purely convention and which is irreducible. Whole paper remains **PARTIAL**.
