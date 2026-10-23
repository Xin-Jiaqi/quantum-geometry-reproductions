# Stage 7d: valley chirality, sliding velocity vertex, and channel-ratio sign test

**Protocol:** synthetic gapped eight-band sliding Kane–Mele, T=0 equilibrium/static finite-Zeeman NMEE. NOT the unsettled driven-DC NMEE, DFT, or proven ferroelectric.

## Exact microscopic sliding/velocity relation

With D=R+r_b-r_a+d and the original hopping T_ab(k,d)=sum_R t exp(-D²/xi²) exp(i k·D), differentiation gives

\[
\boxed{\partial_{d_j} T_{ab}=ik_j T_{ab}+(2i/\xi^2)\partial_{k_j}T_{ab}.}
\]

This is specific to **Gaussian hopping and the chosen Bloch embedding**. Independent checks: maximum identity residual ~4.2e-17 eV/nm and finite slide derivative error ~4.8e-10 eV/nm.

At u=0, t_perp=5 meV, K/K' have T_AB≈T_BA≈0 while |T_AA|=3.611744 meV and |v_AB,x|=|v_AB,y|=0.389169 meV nm. For valley index nu=+1 (first K), -1 (second K'):

\[
v_{AB,y}^\nu=i\nu v_{AB,x}^\nu,\qquad
v_{BA,y}^\nu=-i\nu v_{BA,x}^\nu.
\]

Thus for small q,d,

\[
T_{AB}(K_\nu+q,d)\simeq
v_{AB,x}^\nu [q_x+i\nu q_y+(2i/\xi^2)(d_x+i\nu d_y)].
\]

A cross-orbital hopping vanishing at K still has a **nonzero velocity vertex**. This local expansion is not valid for arbitrary large slides.

## Exact Rashba-helicity sign constraint

For U_z=I_(layer,sublattice)⊗sigma_z, U_z H(alpha_R) U_z†=H(-alpha_R), while S_x is odd and S_z even; velocity transforms covariantly. Therefore

\[
\boxed{\Gamma_{xz;x}(-\alpha_R)=-\Gamma_{xz;x}(+\alpha_R),\qquad
\Gamma_{xz;x}(0)=0.}
\]

At u=.1/.4, t=5 meV and +20 meV Rashba the two-valley disk response is -6.887029e-7/+2.166877e-6; reversing Rashba exactly reverses the signs (numerical residual <1e-15), while alpha_R=0 gives ~1e-18. Hence the microscopic signs of the hopping coefficients are NOT fixed by C3+m_z. This operation is a spin-helicity counterfactual at fixed laboratory Sx, not a real material switching claim.

## Testable channel-ratio criterion

At fixed A/B basis and analytic weak hopping, write H=H0+s V_same+c V_cross with c=rho*s. Layer-sign parity gives

\[
\Gamma/s^2=A(u)+C(u)\rho+B(u)\rho^2+O(s^2).
\]

For A,B>0, C<0, C²>4AB this polynomial is negative strictly between its two positive roots and positive outside. It does NOT predict those coefficient signs from symmetry alone. At u=.4, two K/K' disks (radius .75 nm^-1, 26x72):

| s (meV) | A | B | C | positive upper root |
|---:|---:|---:|---:|---:|
| 2.5 | 3.02906e-5 | .167600 | -.091756 | .54714 |
| 1.25 | 3.02602e-5 | .168079 | -.091606 | .54468 |
| .625 | 3.02527e-5 | .168200 | -.091568 | .54407 |

An **independent finite-s scan** at s=2.5 meV finds Gamma/s² negative for rho=.53 (-.000583) and positive for rho=.54 (+.000331), root ~.5364. As s decreases the root moves towards ~.544. A separate adaptive full-BZ sign check (u=.4, s=2.5 meV, valley disks plus 151² exterior) yields Gamma=-1.916e-8 at rho=.50 and +2.086e-8 at rho=.57; changing exterior sampling to 91² keeps the signs. No rigorous continuous-BZ error bound or material knob is asserted.

The exact static response contains occupied→empty virtual transitions. Stage7c's mixed channel coefficients split by LL/RL/LR/RR as (-.03345,-.04790,-.05333,+.04233) at u=.4, t=5 meV, while the positive channel-even coefficient splits as (-.51251,+.20517,+.71253,-.22617). These yield **quantitative response-sector dominance inequalities**; valley chirality and Rashba sign alone do NOT prove B>0 or C<0. All low/remote virtual transitions must be preserved.

**Remaining:** evaluate spin/Zeeman/velocity/same/cross coupling chains in the unperturbed valley resolvent and establish a true model-specific matrix-element sign condition with error bounds. Keep eight-band response until any smaller dressed-operator theory passes the spectral test. Prior stages and frozen paper reproduction are unchanged.

Source module: bsf_nmee/valley_vertex_symmetry.py. New tests: tests/test_stage7d_chiral_vertex.py. Complete independent raw code, scans, and numerical evidence live in the separately delivered offline Stage7d archive.
