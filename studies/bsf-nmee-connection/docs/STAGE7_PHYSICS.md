# Stage 7: localized valley theory of orbital interference

## 1. Question and protocol

Do K and K' account for the observed sliding-induced static NMEE, and can a momentum-local expansion reproduce not only its sign but the **competition of orbital tunneling channels**? This work continues the stage-2 to stage-6 **toy eight-band sliding Kane–Mele** Hamiltonian, with in-plane displacement `tau(u)=u*(0.23,0.41)` in fractional primitive lattice coordinates. The unperturbed layers carry opposite Rashba SOC and are related by the horizontal layer-exchanging mirror. This is a **T=0, gapped, equilibrium static finite-Zeeman** coefficient, **not** an established driven-DC intrinsic NMEE, material calculation, or physically stable ferroelectric switch.

At finite Zeeman parameter Bz, define the occupied-to-empty spin-electric mixed curvature (in the established program convention)

\[
\chi_{ij}(B_z;\mathbf k)=-2\operatorname{Im}\sum_{v\in\mathrm{occ},c\in\mathrm{emp}}\frac{\langle v|\sigma_i|c\rangle\langle c|v_j|v\rangle}{(\epsilon_v-\epsilon_c)^2},\quad
\Gamma_{iz;j}=\int_{\mathrm{BZ}}\frac{d^2k}{(2\pi)^2}\,\partial_{B_z}\chi_{ij}.
\]

The finite-B derivative is evaluated with field ±0.25 T (and a ±0.125 T spot check); each finite-field Hamiltonian is diagonalized independently with the same frozen insulating occupation. The response is not a Fermi-surface Boltzmann quantity. The approximate k·p Hamiltonian uses **exactly the same** spin matrices and Zeeman convention as the full model.

## 2. Microscopic valley expansion

For each valley `nu=K,K'`, let `q=k-K_nu`. We compute

\[
H_\nu^{(2)}(q,u)=H_\nu(0,u)+q_i V_i^{\nu}(u)+\tfrac12 q_iq_j W_{ij}^{\nu}(u),\quad
v_j^{(2)}(q,u)=V_j^\nu(u)+q_iW_{ij}^\nu(u).
\]

`V_i` is the *actual* velocity derivative of the microscopic 8-band Hamiltonian. `W_ij` is a symmetric finite difference of the same velocity, and its finite-difference-step dependence was checked. In particular, the original microscopic interlayer hopping is generated from real-space displacements `d` with phase `e^{ik·d}`: `∂k_i T=Σ_d i d_i t_d e^{ik·d}`, `∂k_i∂k_j T= -Σ_d d_i d_j t_d e^{ik·d}`. The orbital embedding and the velocity therefore transform together. A momentum-linear k·p truncation omits these curvature terms in both H and v.

The valley model remains **8 bands**; it has not eliminated remote bands via Schrieffer–Wolff or extracted universal small-parameter-only closed-form response coefficients. Taylor data is local to each valley and not a globally periodic Hamiltonian.

## 3. Direct low-energy numerical verification

Four independent settings `u∈{0.1,0.4}`, `t_perp∈{5,20} meV` were tested with two nonoverlapping disks and an independent BZ calculation. At valley radius 0.5 nm^-1, the **signed** sum over both disks equals 99.8%–100.6% of the full BZ `Γ_xz;x`. At radius 0.75 nm^-1, the respective signed ratio differs from unity by at most **0.13%** among the four settings. A ratio slightly larger than 1 reflects weak opposite-sign exterior contributions, not >100% physical probability. This is **not** the integral of the absolute value of the k-resolved response.

The quadratic k·p prediction for the *disk-local* Γ has ≤0.0081% relative error at radii ≤0.75 nm^-1 over this four-setting sample. The linear k·p model deviates roughly 2%–10% on the tested disks; retaining quadratic dispersion and velocity terms is important for quantitative agreement. See `sliding_stage7_valley_kp.json` and `sliding_stage7_capture.json`.

## 4. Orbital-channel physics survives in the valleys

In the fixed A/B orbital basis, split the tunneling as `t_s V_AA/BB + t_c V_AB/BA` and invert the sign of the cross channel, while keeping all other Hamiltonian parameters fixed. The even and cross-sign-odd response coefficients are defined as

\[
C_{\rm mix}=\frac{\Gamma(+t,+t)-\Gamma(+t,-t)}{2t^2},\quad
C_{\rm even}=\frac{\Gamma(+t,+t)+\Gamma(+t,-t)}{2t^2}.
\]

They approach the weak-t² mixed and even coefficients under analyticity; at finite `t=5 meV`, higher-order terms remain. For a two-channel perturbation the odd component is a **controlled Hamiltonian comparison**, not a uniquely gauge-invariant separation of a response observable.

| Sliding u | Full-BZ mixed | Valley k·p mixed | Full-BZ even | Valley k·p even |
|---:|---:|---:|---:|---:|
| 0.10 | -0.031057 | -0.031093 | +0.003542 | +0.003542 |
| 0.20 | -0.058545 | -0.058611 | +0.026666 | +0.026670 |
| 0.40 | -0.092261 | -0.092347 | +0.179001 | +0.179026 |

All values are `Γ/t²` in the original program units divided by eV²; both valleys are summed. At these u, the negative mixed term competes with the faster-growing positive even term, explaining a sliding-induced sign reversal. **No general requirement fixes the sign of either coefficient in another material.**

The original full-BZ t=20 meV sign-reversal root was ≈0.1974; the quadratic valley model predicts `u_c≈0.19751` (disk radius 0.75 nm^-1, finite-B step 0.25 T). The root differs by O(10^-4) in u, and residuals from changing integration density/field difference should not be converted into fictitious more precise digits. The local model also recovers the stage-6 near-symmetric-stack `Tr Γ ∼ u³` and the nematic sector `∼u`; the former is symmetry-enforced as a minimum allowed odd C3 invariant, the latter is a symmetry-allowed linear term.

## 5. Parameter neighborhood and limits

Using full-BZ integration at `t_perp=5 meV`, changing Rashba SOC by ±20%, intrinsic SOC by ±20%, or tunneling decay length by ±12.5%, **all seven tested parameter settings** retain `Γ(u=0.1)<0<Γ(u=0.4)` with positive μ-to-band gap. But the response amplitudes vary substantially. In the most weakly gapped cases the minimum sampled μ distance drops to about 2.6–4.1 meV, so responses become significantly larger/more sensitive. A limited seven-point scan is not a proof of model universality or robustness to all perturbations.

## 6. What we actually learned and what is missing

Established in this **synthetic model**:

1. The valley neighborhoods, rather than a large uncompensated background, account for the major **signed** integrated static response.
2. A second-order valley Taylor expansion independently recovers the local response and its orbital-channel competition without fitting the response itself.
3. Sliding reversal follows competition between two opposing coefficients, not a gap closing or a K-versus-K' cancellation.
4. The response change persists for a modest symmetry-preserving parameter neighborhood, although amplitudes are nonuniversal.
5. C3 plus layer-exchange mirror dictates the leading sliding powers of the trace and allows a linear traceless response, but does not fix signs of microscopic coefficients.

Still unresolved: fully analytic 2/4-band valley downfolding and an explicit low-energy closed formula for the competition coefficients, a basis-independent characterization of the orbital-channel contribution, actual BSF energetics and polarization, complete Wannier spin/position matrix input, and the author NMEE equilibrium static versus driven-DC definition. **No PRL-level novelty claim or real-material prediction has been established.**

## Next narrow physics task

Use symmetry-adapted valley eigenstates to identify a minimal subspace and perform a **controlled Schrieffer–Wolff projection**, explicitly bounding discarded-band corrections, before claiming an analytic low-energy theory. Avoid attempting generic DFT until the relevant valley orbitals and symmetry representations are identified.
