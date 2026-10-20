# BSF × NMEE | Stage 2: Microscopic-shift test of interlayer response

**Scope:** independently computed *synthetic model* predictions, not DFT/Wannier or evidence of actual ferroelectric stability. All numbers below refer to one **zero-temperature, gapped, static interband** NMEE response protocol, the central finite-Zeeman-field derivative of an E-linear spin susceptibility. This is deliberately distinct from the PRL's unresolved driven-DC final expression.

## 1. Microscopic Hamiltonian, one recipe for all shifts

A Bloch Hamiltonian (A/B sublattice × spin × top/bottom) is constructed as

\[
H(k,\tau)=\begin{pmatrix}h_{\rm KM}(k,\lambda_R)&T(k,\tau)\\T^\dagger(k,\tau)&s_z h_{\rm KM}(k,\lambda_R)s_z\end{pmatrix},
\]

where `s_z` is the spin Pauli matrix acting on monolayer 4x4 states; therefore the two isolated layers have opposite Rashba coupling and opposite local `Gamma_xz;x`. For orbital labels `a,b` in A/B,

\[
T_{ab}(k,\tau)=I_{\rm spin}\sum_{R\in\Lambda}t_\perp
\exp[-|R+\rho_b-\rho_a+\tau|^2/\ell^2]
\exp[i k\cdot(R+\rho_b-\rho_a+\tau)].
\]

The sum is truncated at integer Bravais translations `-5..+5`, with `ell=0.16 nm`; verified invariance under shifting tau by one Bravais vector to below `3e-17 eV`. The *physical displacement* also appears in the analytic velocity `dH/dk`. A model with `m_z` exchanges top/bottom and flips the in-plane spin while leaving `B_z` and `E_x` even. Its partner is **recomputed** with `tau -> -tau`, not obtained by assigning or fitting an NMEE tensor.

Default model: Kane–Mele hopping `t=0.85 eV`, intrinsic SOC `10 meV`, layer-opposite Rashba SOC `20 meV`, A/B staggered potential `10 meV`, bond length `0.142 nm`, Zeeman factor `g=10`; `tau=(0.23,0.41)` in primitive fractional coordinates; `mu=0.02541381265 eV`, layer potential difference **zero**, `t_perp=28 meV`.

## 2. Static susceptibility definition and numerical safeguards

With the zero-field occupied projector fixed, for each small static `B_z` calculate the *interband* linear spin-electric susceptibility `chi_{x;x}(B_z)` with the published project's `chi_k` routine, and form `[chi(+B)-chi(-B)]/(2B)`. The observable `S_x` is split as `P_top S_x`, `P_bottom S_x`; total equals their sum. This matches the existing **static diagnostic protocol**. `B=0.25 T` and `0.125 T` are compared; the code rejects B-induced Fermi-gap closure. Uniform BZ grids are demonstrably inadequate near the valleys, even when they appear fine. A two-disk K/K' polar mesh with finer resolution near q=0 plus BZ-exterior integration is used. The code checks Hermiticity, velocity finite difference, reciprocal-space gauge covariance, and interface translation periodicity.

The model's *sampled* middle-band gap on the reference quadrature is `~0.03419 eV`, with `mu` inside it; this is not a rigorous global minimum proof over continuous k.

## 3. Frozen baseline results

All values use `(nm V T)^-1` under the NMEE project's 2D Pauli-spin and Zeeman conventions; not material observables.

| Numerical setting | Total static Gamma_xz;x |
|:--|--:|
| Exterior 81²; valley 90×120 per cone, R=1.5 nm^-1 | +2.03111485 × 10^-4 |
| Exterior 141²; valley 140×180 per cone, R=1.5 | +2.03110469 × 10^-4 |
| Exterior 151²; valley 140×180 per cone, R=2.0 | +2.03110519 × 10^-4 |

Absolute variation `~1.02e-9`, relative `~5e-6`, under independent quadrature refinement. This verifies numerical stability **within the chosen model and integration scheme**, not model realism or scientific novelty.

At the coarse-but-verified resolution, with `B=0.25 T`:

| Control model | Total Gamma_xz;x |
|:--|--:|
| `t_perp=0` | 6.4 × 10^-18 |
| `t_perp=2.5 meV` | ~5.70 × 10^-6 (separate exploratory run) |
| `t_perp=5 meV` | +2.19739 × 10^-5 |
| `t_perp=10 meV` | +7.62759 × 10^-5 |
| `t_perp=20 meV` | +1.82891 × 10^-4 |
| `t_perp=28 meV` | +2.03111 × 10^-4 |
| `t_perp=40 meV` | +1.44461 × 10^-4 |
| `t_perp=28 meV`, A↔B cross-sublattice hoppings **removed** | +2.83083 × 10^-8 |
| `t_perp=28 meV`, Rashba SOC = 0 | ~10^-17 |
| mirror partner `tau -> -tau` | -2.03111 × 10^-4 |

The effect is roughly quadratic at very weak tunneling; it is **not globally monotone** and should not be summarized as `Gamma ∝ t_perp²` over the full range. Top and bottom responses at `t_perp=28 meV` are approximately `-0.00835857` and `+0.00856168`, respectively: the net is the small **imperfect cancellation** of two sizeable local responses. The same-sublattice-only vs all-hopping comparison changes the microscopic hopping matrix and its norm; it establishes *orbital-channel sensitivity* in this model but not a universal exact AB-only selection rule.

## 4. Critical symmetry checks

- `m_z` Hamiltonian and `dH/dk` covariance under `U = layer_swap ⊗ I_sublattice ⊗ s_z`: errors `<3e-17` in native units.
- `tau -> tau + a1` invariance of H and dH/dk with five-shell cutoff: `<3e-17`.
- `k -> k+G` Bloch covariance under the diagonal position-phase unitary: `<2e-15`.
- Total-vs-layer sum in computed response: floating-point accuracy.
- Mirror partner total response sign reverses; each target-layer response is the negative of its **opposite** source layer.
- `B_step=0.125 T` vs `0.25 T` changes total response by `~2.28e-9` (`~1.1e-5` relative).
- At high-symmetry `tau=0`, the default fixed `mu` is *not* within the local direct bandgap at K. Such a point cannot be evaluated with the same insulating-chemical-potential contract without reexamining occupations; it is excluded rather than force-fitted.

## 5. Scientific interpretation and next question

**Validated within the model:** symmetry-related sliding partners have opposite net `Gamma_xz;x`; the isolated monolayers' large opposite local contributions cancel when decoupled; switching on a real-space displacement-dependent interlayer hopping lifts their exact cancellation; the effect depends dramatically on A/B orbital mixing and vanishes without Rashba SOC.

**Hypothesis, not yet established:** the response is associated with interlayer reconstruction of Zeeman/spin quantum geometry, rather than being explainable as a generic field-induced redistribution of layer populations. Further work requires resolving valley/energy contributions, a gauge-invariant microscopic channel decomposition, and tests on distinct model families. This is not yet a PRL-scale general mechanism or a demonstrated sliding-ferroelectric material effect. No claims about energetics, switchability, barriers, Berry connection from real Wannier functions, or device magnitude are made.

**Important provenance:** Earlier dialogue mentioned a `5.798×10^-5` static signal for an informal second-stage model, but its original runnable source and model specification were not present among the finalized stage-1 files. The present **independently specified and archived** model instead predicts `2.03110×10^-4`. These values are not independent convergence estimates of the same parameterized Hamiltonian, and they must not be compared as such.
