# Stage 7f: irreducible tensor competition between orbital NMEE mechanisms

**Scope:** synthetic 8-band sliding Kane–Mele bilayer, 0 K, gapped *equilibrium/static* finite-Zeeman mixed spin–momentum curvature. Neither a claim about the paper's unsettled driven-DC NMEE, nor proof of sliding ferroelectric energetics, DFT, or a material prediction.

## Physical discovery

The previous exact occupied-projector derivative provided the weak interlayer tunneling coefficients for one response component. Now the same **unfitted projector derivatives** calculate all four \(M_{ij}=\Gamma_{iz;j}\) for \(i,j=x,y\), for three orbital-channel combinations:

\[
M(\mathbf d;s,c)=s^2\mathsf A(\mathbf d)+c^2\mathsf B(\mathbf d)+sc\mathsf C(\mathbf d)+O(t_\perp^4).
\]

Each is resolved into scalar \(m_X=\operatorname{Tr}\mathsf X/2\), symmetric traceless \(N_X=(\mathsf X+\mathsf X^T)/2-m_X I\), and antisymmetric sector. The **positive cross-sublattice-only \(\mathsf B\)** is dominated (away from extremely small sliding) by a rapidly growing *isotropic scalar*, while the **negative same/cross interference \(\mathsf C\)** is dominated by a *linearly induced nematic anisotropy*. It is therefore a competition between **different tensor irreducible representations**, not merely two arbitrary signed numbers.

**Qualification:** \(\mathsf B\) still has a symmetry-allowed but tiny *linear nematic*, and \(\mathsf C\) has a cubic scalar; do not claim either sector identically vanishes.

## Symmetry theorem checked independently for each channel

For this specifically declared A/B orbital-channel partition, both hopping perturbations respect the requisite transformed \(C_3\) rotation and layer-exchanging \(m_z\) reference symmetries. The independently recomputed matrices satisfy, separately for \(X=A,B,C\),

\[
\mathsf X(R_3\mathbf d)=R_3\mathsf X(\mathbf d)R_3^T,\qquad
\mathsf X(-\mathbf d)=-\mathsf X(\mathbf d).
\]

In Cartesian complex displacement \(z=d_x+id_y\),

\[
\operatorname{Tr}\mathsf X = a_X\operatorname{Re} z^3+b_X\operatorname{Im}z^3+O(|z|^5),
\quad n_X=(N_X)_{xx}+i(N_X)_{xy}=\gamma_X\bar z+O(|z|^3).
\]

Thus the **trace begins at cubic order or higher** and the symmetric traceless sector **may begin linearly**. Identical cubic selection holds for the antisymmetric in-plane scalar. Under a 60-degree rotation of the real displacement, \(\operatorname{Tr}\mathsf X(R_{60}\mathbf d)=-\operatorname{Tr}\mathsf X(\mathbf d)\), separately for all three channel coefficients. Direct recomputation after physically rotating stacking by 120 and 60 degrees and reversing its sign, at u=.1, gives covariance residual around \(10^{-15}\) in the original program units. This is **not** a theorem for an arbitrary orbital splitting that is not preserved by the spatial symmetries.

## Exact projector numerics: two K/K' disks

Radius .75 nm^-1; radial Gauss × angular quadrature 18×48, with u=.1,.4 refined to 26×72. Coefficients in the original Gamma/eV² units.

| u | coefficient | scalar trace/2 | nematic Nxx |
|--:|:--:|--:|--:|
| .02 | A | +8.63e-10 | +1.718e-6 |
| .02 | B | +2.4855e-5 | +1.3416e-5 |
| .02 | C | -8.357e-7 | -6.0983e-3 |
| .10 | A | +1.0678e-7 | +8.4225e-6 |
| .10 | B | **+0.00307607** | +6.5746e-5 |
| .10 | C | -0.00010343 | **-0.02988470** |
| .40 | A | +5.8328e-6 | +2.4417e-5 |
| .40 | B | **+0.16804370** | +0.00018913 |
| .40 | C | -0.00564872 | **-0.08590333** |

Small-slide asymptotics along d=u(.23*a1+.41*a2): \(\operatorname{Tr}B/2\simeq3.11u^3\), \(N^C_{xx}\simeq-0.305u\), \(N^B_{xx}\simeq+6.7\times10^{-4}u\), \(N^A_{xx}\simeq+8.6\times10^{-5}u\). The scalar B term increasingly overcomes the negative C nematic at larger sliding. Actual finite hopping and full BZ need their own convergence verification, as in Stages 7c–7e.

For equal infinitesimal hopping \(s=c=t\), \(M_{xx}/t^2\) is -0.026838 at u=.1 and +0.076711 at u=.4. **Directional counterfactual:** at u=.4 and the *same magnitude* of d and s=c, independently rotating the displacement gives (angle, \(M_{xx}/t^2\)) = (0°, +0.07667), (60°, -0.21692), (120°, +0.19341), (180°, -0.07667). Trace changes sign under 60°, while the nematic transforms nontrivially; \(M_{xx}\) need not simply reverse under 60°. There is no demonstration that these stacking geometries are accessible stable ferroic partners.

## Closed source-allocation audit: a complementary, not uniquely physical partition

The exact source projector curvature \(i\mathrm{Tr}\{P[P_\lambda,P_k]\}\) expanded in hopping sources s,c and B has 27 placements of those three sources over the three projector factors. Grouping by whether s,c share a factor and which factor receives B gives five contributions to C. At u=.4 (12×36 valley quadrature), contributions are [-.0524645, -.0049929, -.0146662, -.0154644, -.0039371], totaling -0.0915251. A negative 57.3% contribution comes from all three sources inserted into the same factor; all five grouped terms are negative here and also at u=.1. This is an **algebraic grouping under a fixed curvature representation**, not a unique observable decomposition into individual virtual transitions or a general sign theorem.

Stage7e showed that B arises mainly from *Bloch state dressing read out through the intralayer velocity*, while C is dominated by *explicit cross-sublattice interlayer velocity interference*. Stage7f clarifies that these mechanisms mainly occupy the **scalar cubic** and **nematic linear** tensor sectors, respectively. All remote-band virtual processes remain included in the eight-band projector.

## Boundaries and next decisive test

This establishes a **symmetry-governed mechanism in a validated toy model**. It does not yield a universal sign of B or C, a response-preserving 2-/4-band downfolding, a first-principles material, an insulating BZ error theorem, or equivalence with the unresolved driven-DC PRL coefficient. The most useful next question is whether a minimal *valley spin–orbital* matrix-element model with remote-band dressed vertices can reproduce both the positive cubic trace and negative linear nematic, including their relative phase and angular dependence.

Repository files: bsf_nmee/irrep_channel_tensor.py, tests/test_stage7f_irrep_tensor.py, scripts/run_stage7f_irreps.py. An **independent standalone** archive contains the full 2×2 matrices, all five u samples, 27 source placements, plots, code and their SHA256 manifest. Study branch only; frozen paper reproduction is untouched.
