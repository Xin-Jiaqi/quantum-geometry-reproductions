# IGMC theory and response protocol · v0.1

Source: PRL main Eqs. (2)-(17), SM sections I-IV. These are **paper conventions** except where clearly marked as independent check.

## External perturbation and measured quantity

- Single **oscillatory, spatially uniform magnetic field**, `B(t)=B cos(ωt)`, enters only through `H_1=-g μ_B B^β(t) σ_β` (the paper uses Pauli matrices σ rather than σ/2 in this coupling).
- Linear *charge current* `J_α(t)`; current response of order `B`; clean intrinsic limit without phenomenological `τ` in the quoted C and D expressions.
- `e=ℏ=1` in derivation, restored factors in parts of main model example. Specify whether ω denotes angular frequency (rad/s) or cyclic frequency (Hz); article Fig. 2 caption says `10^12 Hz` while the response is written `cos ωt`.
- `f_n = f(ε_n, μ,T)`, `ε_mn = ε_m-ε_n`, `r_nm^α = i⟨u_n|∂_{kα}u_m⟩`, `σ_mn^β=⟨u_m|σ_β|u_n⟩`, `z_nm^αβ = r_nm^α σ_mn^β`, `Q_nm=Re z_nm`, `Z_nm=-2 Im z_nm`.
- Interband connection uses the paper's `r_nm` convention, fixed phases are not required for the gauge-invariant product `r_nm σ_mn`. For a differentiable 2-band Hamiltonian, `r_nm^α=i⟨n|∂_α H|m⟩/(ε_m-ε_n)` for `n != m`; this is an **independent computational identity**, invalid at degeneracy.

## Paper-defined low-frequency coefficients

`J_C^α(t) = -g μ_B σ_C^{αβ} B^β cos(ωt)` with

`σ_C^{αβ}=∑_{n,m≠n}∫[dk] f_n Z_nm^{αβ}= -∑_n∫[dk](∂f_n/∂ε_n) v_n^α σ_nn^β` (paper Eqs. (11), (13)).

`J_D^α(t) = -g μ_B σ_D^{αβ} B^β sin(ωt)` with

`σ_D^{αβ}= ω∑_{n,m≠n}∫[dk] f_n · 2 Q_nm^{αβ}/(ε_m-ε_n)` (paper Eq. (12)).

Displacement is the derivative of an induced AC polarization, not a finite DC transport current. For a gapped T=0 insulator, C vanishes whereas D can be nonzero at finite frequency and tends to zero linearly as ω→0. The paper also includes an `O(ω²)` Zeeman curvature correction to polarization, dropped in final D.

## Independently computed AC Kubo protocol

For nondegenerate eigenstates, define **complex** `K_ab(Ω) = ∑_{nm}∫[dk] (f_n-f_m) v_nm^α σ_mn^β / (Ω-(ε_m-ε_n)+iη)`, where Ω denotes **photon energy** `ℏω` in the model's energy units (our code uses ℏ=1). Compare the nonresonant low-frequency expansion:

`K_ab(Ω)= σ_C^{αβ} + i σ_D^{αβ} + O(Ω²)` when η→0 after keeping Ω below gap. Physical current is `Re[-g μ_B B K_ab(Ω)e^{-iωt}]`. This complex reconstruction is not interchangeable with a static dc coefficient. Check real O(Ω²) and imaginary O(Ω³) corrections, or define alternate finite-η regulator explicitly.

The SM Green's-function calculation and Liouville theory provide a second derivation for `K`; a numerical implementation using *the same formula twice* is not independent. Numerical C from a Fermi-surface derivative and D from a finite-field induced polarization will be preferred as independent controls.

## Model boundaries

- Rashba Eq. (15): **continuum**, not a periodic BZ model. Interband gap closes at `k=0`, geometric kernels are singular. Do not use uniform cartesian coarse grids as a convergence proof.
- Dirac Eq. (17): periodic square lattice, two nondegenerate bands for the published `β>0`, well-defined finite BZ integral `k_x,k_y ∈ [-π/a,π/a)`. At `β=0`, near gap closure, switch to a special integrator and flag failed nondegenerate limit.
- For arbitrary Wannier models, physical spin matrices must be available in same gauge. This stage implements only *both paper models*.

## Symmetry audit (not an accepted correction yet)

Main Table I lists `Z` as **T-even** and `Q` as **T-odd**. Accordingly `σ_C` and `σ_D/ω` have different behavior in a T-invariant solid: a Rashba model is T invariant with nonzero C but zero D. Yet main Eq. (14) prints a common `η_T=-1` on antiunitary symmetry for `σ^{C/D}`, while SM uses different `eV²` vs `aeV²` Jahn symbols. **OPEN**: verify whether Eq. (14) omitted the C/D-dependent antiunitary factor or whether a distinct finite-frequency convention resolves it. Preserve main Eq. (14) unchanged in the source map; do not silently treat it as the same physical transformation for both channels.
