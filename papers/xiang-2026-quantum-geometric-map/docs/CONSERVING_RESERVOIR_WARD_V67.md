# v67 | Two genuinely semi-infinite reservoirs, differentiated scattering Ward, projected E/B contacts

**Scope:** A new **microscopically specified, noninteracting, coherent and elastically charge-conserving** two-reservoir scattering oracle. **It is a *different transport observable* from v66:** the differential **two-terminal Landauer charge conductance** at infinitesimal *reservoir voltage* and its derivatives with respect to an *independently controlled static device electric gate* and magnetic Peierls field. It is **not** the original Xiang PRB orbital Hall coefficient due to a spatially homogeneous time-dependent electric field.

## Exact Hamiltonian and reservoir definition

Reuse v65/v66's explicit finite four-orbital Peierls Hamiltonian `H(B)`, `M(B)=d_B H`, noncommuting supplied Hermitian projected polarizations `P_a(B), P_b(B)`, `d_B P_a`, `d_B P_b`. A *separate* electrostatic gate parameter E enters the device as `H_D(B,E)=H(B)-E P_b(B)`. It is **not** the same as the electrochemical potential difference `V_L-V_R` that drives the two-terminal current.

The two semi-infinite one-dimensional single-channel leads have explicit microscopic Hamiltonians `h_00=eps0`, `h_01=t` and retarded surface Green function

\[
 g_s^r(z)=\frac{z-\sqrt{z^2-4t^2}}{2t^2},\quad z=\mu-\epsilon_0+i\eta,
\]
with `Im g_s<0` in the propagating band. This solves the lead's quadratic surface Dyson equation. `eta=1e-10` is a numerical retarded regulator, **not a fitted scalar-gamma collision/bath**; the eta=0 clean propagating limit is separately checked.

Left and right microscopic coupling vectors in the symmetric reference are `c_L=(0.48+0.07i,0.39-0.06i,0,0)`, `c_R=(0,0,0.52-0.09i,0)`. Their phase interference means contact magnetic derivative terms are nonzero after B-dependent magnetic gauge transformations. In a B-dependent basis `U(B)`, physical leads must transform as `c_a(B)=U(B)c_a`, `c'_a(B)=U'(B)c_a`, not remain artificially fixed. The 1D leads are field independent in this **chosen physical model**; there is magnetic flux through finite sample loops, and no claimed realistic homogeneous B field in arbitrary extended metallic electrodes.

\[
\begin{aligned}
 \Sigma^r_a&=g_s^r c_a c_a^\dagger, &\Gamma_a&=i(\Sigma_a^r-\Sigma_a^{r\dagger}),\\
 \partial_B\Sigma_a^r &=g_s^r (c_a' c_a^\dagger+c_a c_a'^\dagger),
 & \partial_B\Gamma_a&=i(\Sigma_a'-\Sigma_a'^\dagger),\\
 G^r&=[\mu-H_D-\Sigma_L-\Sigma_R]^{-1},
 & G_B^r&=G^r(M-E\,\partial_B P_b+\Sigma_L'+\Sigma_R')G^r.
\end{aligned}
\]

The measured **reservoir-bias linear-response conductance** at `T_lead=0` is `G_2term=(q²/h) * T`, where `T=Tr[Gamma_L G^r Gamma_R G^a]`. Lead-boundary current conservation is independently checked by the full two-channel scattering matrix `S=I-iW†G^rW`, `W=(sqrt(-2 Im g_s)c_L,sqrt(-2 Im g_s)c_R)`. It must obey `S†S=I`, `|S_RL|²=T` and `I_L+I_R=0` for the same energy. Unlike a phenomenological gamma, these channels are explicitly microscopic semi-infinite chains.

## Mixed magnetic/electric gating derivative, and why contacts matter

`∂_B T` differentiates all of `Gamma_L`, `Gamma_R`, `G^r` and `G^a`. The gate derivative uses `G_E=-G P_b G`. The mixed derivative `∂_B∂_E T` also uses

\[
G_{EB}=-G_B P_b G-G(\partial_B P_b)G-G P_b G_B.
\]

Eight differentiated Caroli terms account for the product's four factors. This is a **magneto-gating coefficient of two-terminal conductance**, not a direct mixed orbital Hall transport coefficient. A separately implemented 16-corner fourth-order E/B finite difference checks the entire result, avoiding the analytic derivative implementation.

**Additional differentiated Ward identities:** with `G_B`, contact linewidth derivatives, and `W_B`, verify the magnetic derivative of the *spectral identity* `i(G-G†)=G(Gamma_L+Gamma_R)G†` and of the scattering unitarity `S†S=I`. Independently evaluate `d_B T=2Re(S_RL* d_B S_RL)` and compare with differentiated Caroli trace. An uncontacted fourth orbital maintains local charge continuity, and in the noncommuting orbital frame the **physical site projector is also B-dependent** (`Q(B)=U Q0 U†`); this projector's derivative is necessary for the local magnetic continuity Ward test.

## Independently frozen dimensionless benchmark

Reference: `mu=.13`, `B=.23`, signed `q/hbar=-.84`, static device gate `E=.04`, microscopic lead onsite 0, hopping 1.1, numerical regulator `1e-10`. Explicit finite open graph and projected dipoles from v66.

| Physical observable | Symmetric | Landau | Noncommuting orbital frame |
| --- | ---: | ---: | ---: |
| Landauer transmission T | 0.900626086558 | same within 2e-15 | same within 2e-15 |
| ∂B T | -0.180487543399 | same within 2e-15 | same within 2e-15 |
| ∂E_gate T | -0.557089765326 | same within 1e-14 | same within 1e-14 |
| ∂B∂E_gate T | +0.884057478696 | same within 1e-14 | same within 1e-14 |
| Physical internal Tr(i[H,P_a] rho_bias) | -0.365610142747 | same | same |
| ∂B internal-current diagnostic | +0.090535511480 | same | same |

The scattering unitarity residual is `~1e-15`, differentiated scattering Ward and differentiated optical identity `~1e-14` or smaller, and uncontacted local charge continuity and its magnetic derivative are `~1e-16`. The independent 16-corner mixed finite-difference matches the analytic result within `4e-10` in all three representations, with `h=0.002`. Caroli and independent scattering amplitude's ∂B T agree within `1e-12`.

**Failing negative controls:** hold lead couplings fixed under a noncommuting orbital frame and T changes from `0.90062608656` to `0.89317054195`. Drop magnetic derivatives of the orbital-frame lead selfenergies and mixed ∂B∂E_gate T incorrectly becomes `0.61437900958`, rather than `0.88405747870`. Drop `d_B P_b` and it incorrectly becomes `0.78429578275`. In Landau gauge the corresponding wrong values are `0.87532437425` and `0.90311061288`. Holding the electrical gate origin fixed while shifting `P_b→P_b+cI` produces a spurious result; **shift both reservoir chemical potential and lead onsite energy by -E*c** and exact invariance is restored. This tests a consistent electrostatic reference, not an apparent orbital Hall coefficient.

## Scientific acceptance boundary and next steps

**VERIFIED:** specified two-channel elastic scattering, semi-infinite chain retarded self energies, exact optical/scattering Ward, local continuity and their magnetic derivatives, contact/position-frame covariance, independent finite gate/B derivative.

**CONDITIONAL:** finite four-orbital Peierls open-device model with supplied projected off-diagonal orbital dipoles; coherent noninteracting reservoirs; device-only static electric gate. A physical material needs field-complete Wannier position operators, lead electrostatics, and full covariant current and scattering vertices.

**UNRESOLVED:** **same-observable** homogeneous electric-field-driven orbital mixed E/B current versus conserving reservoir Keldysh (including displaced polarization and boundary-current vertices), any conserving *interacting* collision kernel and bulk/dc limits. The original Xiang PRB disputed Supplemental Material and any omitted three-band channels remain unverified. The Xiang source-paper science status stays **PARTIAL**, unchanged.

**Reproduction:** from QGMap root, `PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -p 'test*.py' -v` and `PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/reproduce_v67.py`; resulting numbers in `results/frozen/v67_conserving_reservoir_ward.json`. v67 has **45 new physical science tests + one independent frozen regeneration test**, plus prior v60-v66 independent 219 = **265 standalone**. Full GitHub repository test count is a separate CI gate.
