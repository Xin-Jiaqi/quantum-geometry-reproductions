# Stage 6: sliding-order selection of the NMEE tensor

## Scope and headline

This work uses the same eight-band synthetic sliding Kane–Mele bilayer as stages 2–5. **All numerical responses refer to the zero-temperature, gapped, equilibrium/static mixed spin–momentum curvature protocol.** They are not a reproduction of the unsettled driven-DC intrinsic NMEE response, a demonstrated sliding-ferroelectric ground state, or a material prediction.

**Main distinction:** The **tensor trace** is constrained to start at *third order* in the deviation from a C3- and horizontal-mirror-symmetric reference stacking, while the **traceless/nematic tensor** can begin at *first order*. This is a general symmetry selection rule, not an assertion that a chosen orbital channel's whole response must scale with either power.

## Symmetry proof: an observable distinction

Fix a reference homobilayer stacking invariant under a threefold rotation about `z` and a horizontal layer-exchanging mirror `m_z`. Let the two-component Cartesian relative displacement from this reference be `d`, and define the in-plane 2x2 tensor

\[
M_{ij}(\mathbf d)=\Gamma_{i z;j}(\mathbf d),\qquad i,j=x,y.
\]

Spin `i` and magnetic field `z` are axial, whereas the electric field index `j` is polar. Under `m_z = diag(1,1,-1)`, `M` is odd (in-plane spin changes sign but `B_z` and in-plane electric field do not). In the bilayer `m_z` also exchanges layers, mapping the relative displacement `d` to `-d`:

\[
M(-\mathbf d)=-M(\mathbf d). \tag{A}
\]

For a proper threefold `C3` with in-plane Cartesian rotation `R`:

\[
M(R\mathbf d)=R M(\mathbf d) R^{T}. \tag{B}
\]

Both conditions are verified on the **independently recomputed Hamiltonians**, not imposed on the numerical output. Let

\[
s=\tfrac12(M_{xx}+M_{yy}),\quad a=\tfrac12(M_{xy}-M_{yx}),\quad N=\tfrac12(M+M^T)-s I_2.
\]

The scalars `s` and `a` are invariant under proper C3 rotations and odd under `d -> -d`. Assuming the insulating response is analytic near `d=0`, no C3-invariant polynomial of first degree exists. The lowest-degree odd C3-invariant polynomials are the two cubic harmonics. With `z = d_x + i d_y`:

\[
s(z,\bar z)=c_1\operatorname{Re}(z^3)+c_2\operatorname{Im}(z^3)+O(|z|^5),\quad
 a(z,\bar z)=d_1\operatorname{Re}(z^3)+d_2\operatorname{Im}(z^3)+O(|z|^5).
\]

For the traceless symmetric sector define `n=(N_xx + i N_xy)`. It transforms under a planar rotation of angle `phi` as `n -> e^(2i phi)n`. Under **C3 only**, `e^(2i 2pi/3)=e^(-i 2pi/3)` and an *allowed linear term* is

\[
n(z,\bar z)=c\bar z+O(|z|^3).
\]

The coefficient `c` may vanish accidentally or by **additional** group constraints. In particular, C3+mz do **not** force the whole component `Gamma_xz;x` to be cubic: its `N_xx` part can be linear.

An exact global consequence of (A),(B) is that the tensor **trace** changes sign when the sliding displacement is rotated by 60 degrees (C3 + inversion in the displacement plane):

\[
\operatorname{Tr}M(R_{60}\mathbf d)=-\operatorname{Tr}M(\mathbf d).
\]

This is a *response/geometry* relation; it is not evidence that these configurations are stable ferroelectric partners.

## Selection of orbital-channel cross derivatives

At **fixed orbital basis and stacking** write `H = H0 + s V_same + c V_cross`. Let `C_{alpha gamma;beta} = d_s d_c Gamma |_(s=c=0)`. A unitary point-group operation that preserves the reference Hamiltonian and the **same stacking** may act on each perturbation with 1D characters `chi_s(g)` and `chi_c(g)`, without mixing the channels. Then the mixed derivative is symmetry-allowed only in the tensor representation satisfying

\[
C = \chi_s(g)\chi_c(g)\,T_{R_g} C,\qquad
(T_R C)_{\alpha\gamma;\beta}=R_{\alpha a}R_{\gamma b}R_{\beta c}C_{ab;c}.
\]

The projector is `P C = (1/|G|) sum_g chi_s(g)chi_c(g) T_(R_g) C`. This assumes the group and characters are genuine; the code does **not** infer orbital-channel characters from a point-group label. If channels mix, a full orbital-channel representation is required instead; in that case this 1D-character projector must not be used.

At zero sliding displacement, direct Hamiltonian tests confirm that both **same-sublattice** and **cross-sublattice** hopping perturbations in this toy model are individually `m_z`-even; `Gamma_xz;x` is `m_z`-odd. Consequently each quadratic coefficient `A`, `B`, and `C` for that tensor component must vanish **at this high-symmetry stacking**. When sliding breaks `m_z`, the response and its cross derivative can turn on continuously. This is a rigorous symmetry reason for the *onset* of the response, but not for the sign or its magnitude away from the symmetry point.

## Frozen independently evaluated numerical data

For the synthetic model with `t_perp=5 meV` and `d(u)=u(0.23,0.41)` in fractional coordinates, the integrator independently evaluates both in-plane spin and electric indices, with `B_z`. The plotted quantities are the **full static** tensor trace and the Frobenius norm of the **traceless symmetric odd-under-cross-hopping-sign** part.

| Sliding u | Trace M | (Trace M)/u^3 | Mixed-channel nematic norm | Nematic norm/u |
|---:|---:|---:|---:|---:|
| 0.025 | 2.5483e-9 | 1.63094e-4 | 2.8321e-7 | 1.13286e-5 |
| 0.050 | 2.0300e-8 | 1.62404e-4 | 5.6390e-7 | 1.12780e-5 |
| 0.100 | 1.6063e-7 | 1.60634e-4 | 1.1078e-6 | 1.10778e-5 |
| 0.200 | 1.2322e-6 | 1.54025e-4 | 2.0616e-6 | 1.03081e-5 |
| 0.400 | 8.3717e-6 | 1.30808e-4 | 3.0686e-6 | 7.67142e-6 |

Numbers are rounded for this report; exact data and all mesh settings are in `results/sliding_stage6_symmetry.json`. At u=0.05 vs u=0.10, the observed trace ratio is ~7.91 (leading cubic prediction: 8), whereas the mixed-channel nematic norm ratio is ~1.96 (leading linear prediction: 2). Departures at larger shifts are expected; the selection rule constrains **leading order near the symmetric reference** and the trace's exact symmetry parity, not a global power law.

Rotate the Cartesian sliding displacement at `u=0.1`: the trace ratio at +60 degrees is about `-0.99978`; +120 degrees yields `+0.99982`; +180 degrees yields `-0.999999999`. The +120-degree independently computed 2x2 matrix differs from `R M R^T` by at most `2.65e-11` in the calculation's response units. This is a direct Hamiltonian/response covariance check.

Compared with a denser valley quadrature, absolute trace differences at u=0.05/0.1/0.2 are `2.94e-11`, `5.56e-11`, and `8.62e-11`, respectively. Reducing the magnetic-field finite-difference step to 0.125 T changes the maximal tensor component by at most `~5.54e-11` in those test cases. A fully rigorous continuum BZ error bound is **not** claimed.

## Physical interpretation and caveats

1. **Generic, symmetry-enforced:** `M(-d)=-M(d)`; `Tr M` is cubic or higher near a C3-symmetric mirror stacking, and changes sign under a 60-degree rotation of `d`; the nematic tensor *may* appear linearly.
2. **Model-specific, observed:** `Tr M` is nonzero and approximately cubic; the mixed-channel nematic magnitude is approximately linear. The leading interference contribution is sensitive to the declared `A/B` hopping basis. A small even-channel nematic response does not imply symmetry-forbidden linear terms for all models.
3. **Not proven:** Stable ferroelectric pairs, path energetics, universal sign of the mixed derivative, material-scale magnitude, and agreement between static equilibrium and the original paper's driven-DC NMEE.

**Next disciplined check:** derive a compact low-energy `K/K'` orbital-interference model explaining the **coefficients**, not just the allowed powers, and then test robustness against symmetry-preserving changes to monolayer SOC and hopping ranges. This remains at tight-binding level; no DFT request is necessary yet.
