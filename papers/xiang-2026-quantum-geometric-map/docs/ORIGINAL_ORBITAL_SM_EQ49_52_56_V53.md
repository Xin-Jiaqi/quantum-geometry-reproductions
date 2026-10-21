# v53 | Authentic orbital SM Eq. (49)/(52)/(56): source-level physics audit

**Original source:** User-provided, SHA256-pinned 16-page Supplemental Material of Xiang et al., *Quantum geometric map of magnetotransport*, PRB **113**, L201406 (2026), especially pp. **11–12** (printed Eqs. 48, 49, 50, 51, 52, 56). Original publisher PDFs not redistributed. **Entire-paper science remains PARTIAL.**

## 1. The source Eq49→Eq52 metric velocity term has a coefficient discrepancy

Original **SM Eq. (49)** writes

\[
F_n^{bc}=-\sum_m\epsilon_{bij}
\left[\boxed{2}\frac{(v_n^i+v_m^i)g_{nm}^{jc}}{\varepsilon_{nm}}+\partial_i g_{nm}^{jc}\right].
\]

Original **SM Eq. (48)** gives the magnetic positional-shift sector \(\sum_n\int_{\rm BZ}(\partial_aF_n^{bc}-\partial_bF_n^{ac})f_n\).
For a smooth, periodic BZ, integration by parts makes its velocity–local-quantum-metric term

\[
\boxed{2}\sum_{nm}\int_{\rm BZ}
\frac{[v_n^a\epsilon_{bij}-v_n^b\epsilon_{aij}]
(v_n^i+v_m^i)g_{nm}^{jc}}{\varepsilon_{nm}}f'_n,
\]

plus the printed minus-metric-Hessian sector. Yet the corresponding velocity-metric term in **printed SM Eq. (52)** has **coefficient 1**, not 2. Curvature–magnetic-moment terms are an independently identifiable sector and cannot generically repair this relative metric coefficient.

**Independent physical, non-invented local metric control:** use a nondegenerate *genuine* smooth 3D 2-band periodic Bloch Hamiltonian \(H(k)=d_0(k)I+\mathbf d(k)\cdot\boldsymbol\sigma\), with physically valid \(g_{nm}^{ij}=\frac14\partial_i\hat d\cdot\partial_j\hat d\) and analytically differentiated dispersions. Compute actual F^{yz},F^{xz} from printed Eq49; finite-difference (∂x F^{yz}−∂y F^{xz}) f. Compare printed Eq52 metric terms versus coefficient-2 form after integrated BZ finite differences. Component (a,b,c)=(x,y,z), mu=2.55,T=.34 model units, BZ mean \((2\pi)^{-3}\int d^3k\):

| Mesh | Direct Eq48 using Eq49 | Printed Eq52 metric | Coefficient-2 variant | Missing velocity term |
|---:|---:|---:|---:|---:|
| 20³ | −1.53668e−5 | −2.15018e−5 | −1.45723e−5 | +6.92954e−6 |
| 28³ | −1.63178e−5 | −2.28365e−5 | −1.59069e−5 | +6.92954e−6 |
| 36³ | −1.67315e−5 | −2.34111e−5 | −1.64815e−5 | +6.92954e−6 |
| 48³ | −1.70159e−5 | −2.38043e−5 | −1.68748e−5 | +6.92954e−6 |

The coefficient-2 residual contracts toward zero with mesh size while the literal Eq52 residual approaches a nonzero limit. **Independent negative control**: set d0(k)=0, so \(v_n+v_m=0\), and both formulas coincide at floating-point precision; this isolates precisely the missing velocity factor. This is a **bounded source-algebra inconsistency**, NOT an author-certified erratum or an experimental invalidation. Eq49 itself expressly drops a class of three-band terms; the two-band control has no distinct third intermediate band.

## 2. Source Eq52 explicitly drops the m≠l interband Berry–orbital-moment term

Original **SM Eq. (50)–(51)** imply the untruncated moment–Berry sector of **SM Eq. (48)** contains a double interband sum

\[
-\frac14\sum_n\int_{\rm BZ} f'_n
\left(\sum_m\Omega^{ab}_{nm}\right)
\left(\sum_l\varepsilon_{nl}\epsilon_{ijc}\Omega^{ij}_{nl}\right).
\]

For (a,b,c)=(x,y,z), \(\epsilon_{ijz}\Omega_{nl}^{ij}=2\Omega_{nl}^{xy}\). Printed SM Eq52 retains only the **m=l** diagonal of this double sum, with the **m≠l** terms explicitly dropped by the authors. This is a declared approximation, **not** an algebra mistake. We quantify its potential importance with a genuine periodic **3-band Hermitian** tight-binding Hamiltonian, analytical Bloch velocity operators and separately diagonalized interband Berry connections, including an independent band-U(1) gauge phase test.

For mu=−1.1,T=.45, 56² mesh (all arbitrary model energy units; BZ average):

| Berry curvature × orbital magnetic moment subterm | Result |
|---|---:|
| Full m,l double sum | −0.00079945105 |
| Retained m=l term of printed Eq52 | −0.00154218315 |
| **Omitted m≠l interband interference** | +0.00074273210 |

The magnitude of the omitted process is about **92.9% of this complete subterm** in this example; at mu=1.5 the ratio is only **5.6%**. This is **not 92.9% of the total orbital MNHE** and is not a universal material estimate. It does not include other dropped processes in the original Eq49.

## 3. SM Eq56: the two printed orbital-PHE expressions pass

For (a,b,c)=(x,y,z), the common Berry term is \(-v_xv_y\Omega_z f'_n\). The first printed line contains

\[
-\partial_x m_z v_y f'_n-v_x\partial_y(m_z f'_n),
\]

and the second contains

\[
2m_z\partial_{xy}^2 f_n-v_xv_y m_z f''_n.
\]

On a smooth periodic physical two-band Bloch model with orbital moment \(m_z=\frac12\varepsilon_{nm}\Omega_{nm}^{xy}\), independent numerical BZ integration supports their equivalence by integration by parts. Absolute first-minus-second residuals at grids 24²,32²,48²,64² are respectively **1.960e−5,1.104e−5,4.913e−6,2.765e−6**, approaching zero approximately as \(N^{-2}\). **PASS bounded integrated identity**; we do not falsely flag Eq56 as inconsistent.

## Engineering and acceptance

New model **12 physics tests**: authentic g_nm constructed from eigenprojector d, source F49 vs F48, Eq52 coefficient-2 convergence, flat-d0 no-velocity negative control, 3-band full/m=l/interference separated, mesh and U(1) gauge convergence, Fermi-window sensitivity, original Eq56's two independent forms and integration-by-parts rate, fail-closed inputs and SHA-pinned source. Reproduce:

```bash
python -m unittest tests.test_orbital_sm_v53 -v
python scripts/reproduce_orbital_sm_v53.py --output results/frozen/v53_original_orbital_source_audit.json
```

**Open:** reconstruct the untruncated orbital positional shift including Eq49's *other* omitted three-band processes, verify full orbital MNHE (not just subterm), original SM spin Eq20/31 sign, original Fig2 e²/h versus e²/ħ SI ambiguity, and fully physical bath/contact or material response. **Overall paper PARTIAL**. v54 should work on the complete gauge-covariant orbital magnetic positional shift, followed by realistic multiband checks. No source signs/coefficient silently altered; no author contact.
