# BD31 | Conditional local equivariant-index rule for layer-resolved current-induced spin response

**Status:** exact symmetry derivation conditional on the assumed \(D_{3h}\)-type registry point and smooth response field; not a new mathematical theorem by itself. Direct LEE companion-paper overlap remains unresolved.

## 1. Define the field and the operation laws

Let the relative in-plane layer displacement near a high-symmetry registry be
\[
z=u_x+i u_y,\qquad \omega=e^{2\pi i/3}.
\]
For the response of the layer-even/odd out-of-plane spin density to in-plane current, define the complex response vector
\[
F_\pm(z,\bar z)=\chi_{\pm,x}(u)+i\chi_{\pm,y}(u),\quad
S_\pm^z=S_{\rm top}^z\pm S_{\rm bottom}^z.
\]
This is a **current-to-spin-density response field over registry space**, not a spin-current field, Berry curvature, or quantum metric.

Assume the high-symmetry reference structure has \(C_{3z}\), a vertical mirror \(m_x\), and a horizontal mirror \(m_z\) that exchanges layers. Since \(S_z\) is an axial component, it is even under \(m_z\) but odd under \(m_x\); the extra layer-odd factor makes \(S_-^z\) odd under layer exchange. The relative displacement also changes sign under layer exchange. Thus the equivariance conditions are
\[
F_\pm(\omega z,\omega^*\bar z)=\omega F_\pm(z,\bar z),
\]
\[
F_\pm(-\bar z,z)=\overline{F_\pm(z,\bar z)},
\]
\[
F_-(-z,-\bar z)=-F_-(z,\bar z),\qquad
F_+(-z,-\bar z)=+F_+(z,\bar z).
\]
The mirror relation is written for \(m_x:z\mapsto-\bar z\). These laws need to be re-derived if the material has a different layer-exchange Seitz operation, if the response is defined under fixed electric field instead of fixed current, or if magnetic/antiunitary symmetries are present.

## 2. Allowed local Taylor terms

Expand a smooth field as \(F(z,\bar z)=\sum_{mn}c_{mn}z^m\bar z^n\). Threefold covariance requires
\[
m-n\equiv1\pmod 3.
\]
Layer-exchange parity requires odd total degree for \(F_-\) and even total degree for \(F_+\).

### Layer-odd output

The lowest permitted term is degree one, \(z\). The vertical mirror forces its coefficient to be purely imaginary:
\[
F_-(z,\bar z)=i\gamma z+O(|z|^3),\quad \gamma\in\mathbb R.
\]
If \(\gamma\ne0\), the response has an isolated zero at \(z=0\) with local winding \(+1\).

### Layer-even output

No linear term is allowed. The lowest permitted term is degree two, \(\bar z^2\); the vertical mirror makes its coefficient real:
\[
F_+(z,\bar z)=\beta\bar z^2+O(|z|^4),\quad \beta\in\mathbb R.
\]
If \(\beta\ne0\), the zero at \(z=0\) has local winding \(-2\).

**Critical logical qualifier:** symmetry allows these coefficients but does not force \(\gamma\) or \(\beta\) to be nonzero. If a leading coefficient vanishes accidentally, higher-order terms control the local index. The index is defined only on a loop where \(F\ne0\). A finite-radius winding can differ from the local index if extra zeros lie inside the loop or cross it as parameters vary.

## 3. A general stabilizer test for global forced zeros

For any registry \(\mathbf u_*\), compute its actual stabilizer \(H_{\mathbf u_*}\) on the periodic registry torus, including fractional translations/Seitz parts. For every \(g\in H_{\mathbf u_*}\), derive the response representation \(D_\pm(g)\) on the current-response vector, including axial-spin sign and layer parity. Then
\[
F_\pm(\mathbf u_*)=D_\pm(g)F_\pm(\mathbf u_*)\quad\forall g\in H_{\mathbf u_*}.
\]
Therefore the allowed value lies in the common invariant subspace
\[
\mathrm{Fix}\,D_\pm(H_{\mathbf u_*})
=\bigcap_{g\in H_{\mathbf u_*}}\ker[D_\pm(g)-I].
\]
If this common fixed subspace is \(\{0\}\), the response zero is symmetry-forced. If it is one-dimensional, symmetry only constrains the direction, not the amplitude. This criterion is generalizable across layer groups and is a more precise statement than inferring every zero from the high-symmetry label alone.

The current Pd2Te2 toy numerically reproduces:
- layer-odd zero at G/A/B/C/M/N;
- layer-even zero at G/M/N but nonzero at A/B/C;
- conditional local index +1/−2 near G for the tested parameters.

The six-point pattern is a **model result**; the stabilizer derivation should be implemented independently for all six points before claiming a general classification.

## 4. What might be publishable, and what is not

Not novel by itself:
- Neumann's principle / second-rank axial-tensor classification;
- generic layer Edelstein response;
- a linear versus quadratic Taylor expansion derived from point-group covariance;
- ordinary winding numbers or equivariant-map index theory.

Potential S2-level contribution, **only if not already in the submitted sliding-REE companion**:
- a systematic registry-torus stabilizer algorithm for operator-resolved response fields;
- a theorem linking layer-exchange parity to forced zero sets and conditional local indices;
- proof that these constraints differ for \(S_+\) and \(S_-\) in the same bilayer, with a material-grounded, gauge-covariant response calculation.

The published LEE article already covers general bilayer LEE tensor criteria over all 80 layer groups. It cites a same-group work labelled submitted: *Control of the Rashba-Edelstein effect via sliding ferroelectricity in bilayer systems*. Its full text must be obtained before treating continuous-sliding response as an independent novelty. See [BD30 scope and robustness audit](BD30_LEE_SCOPE_AND_MODEL_ROBUSTNESS_AUDIT_20261022.md).

## 5. Next gate

1. Derive stabilizers and response representations for G/A/B/C/M/N from explicit fractional coordinates and Seitz operations; do not infer them from names.
2. Search for / obtain the submitted sliding-REE manuscript and compare exact propositions.
3. Replace band-diagonal Boltzmann expectations by a gauge-covariant formulation around degenerate bands; fit a spinful Pd-d + Te-p Hamiltonian before interpreting the field physically.
4. Do not submit DFT/HPC while prior art and the response formulation remain unresolved.
