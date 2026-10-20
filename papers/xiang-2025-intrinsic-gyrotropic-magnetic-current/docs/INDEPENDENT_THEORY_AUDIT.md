# Independent theory audit · round 1 (2026-10-19)

These are independent checks against the uploaded PRL Eqs. (5)-(17), not an assertion that the entire SM derivation has been audited.

## A. Two-band gauge-stable computational expression: PASS (nondegenerate scope)

Use `r_nm^α=i⟨n|∂_α H|m⟩/(ε_m−ε_n)` and `z_nm^αβ=r_nm^α σ_mn^β`. This directly avoids finite differences of arbitrarily phased eigenvectors. `Q=Re z`, `Z=−2Im z` under the PRL convention. Tested by unitary basis rotation of the **Hamiltonian, velocity and Pauli spin together**. The implementation raises at degeneracies instead of using a finite nondegenerate denominator.

## B. Rashba model Eq. (16): PASS for continuum model normalization

`H(k)=α k² I + λ_R (k_y σ_x−k_x σ_y)` with `α=ℏ²/(2m)=3.80998212/(1.1)=3.46362011 eV Å²`. Lower band minimum `−λ²/(4α)≈−0.0259843 eV`. At `μ≥0`, `σ^C_xy(e/ℏ removed)=λ/(4πα)=0.0137851385 Å^-1`, independent of μ. At `−λ²/(4α)≤μ<0`, it is `sqrt(λ²+4αμ)/(4πα)`. Direct polar momentum integration of the independent `z=rσ` eigenvector matrix elements reproduces this to ~0.1% at representative chemical potentials. `σ_C` vanishes below the band bottom.

Independent checks are **not** only repeating Eq. (11): (i) closed-form dispersion and Fermi radii, (ii) numerical interband geometric sum, (iii) explicit spin-texture/magnetic-point-group parity.

## C. Dirac-lattice finite-ω versus low-ω: PASS for gapped model

`d=(v sin(ak_y)/a,−v sin(ak_x)/a, β+2M[2−cos(ak_x)−cos(ak_y)]/a²)`, `H=d·σ`. Γ direct gap = `2β=0.5 eV`.

Compute separately `C=Σ f_n Z_nm`, `D/Ω=Σ f_n 2 Q_nm/(ε_m−ε_n)`, and complex full Kubo `K(Ω)=Σ (f_n−f_m) v_nm σ_mn/[Ω−(ε_m−ε_n)]` without broadening, with `Ω=ℏω` (eV), `Ω<gap`.

At `N=101²`, `μ=0 eV`, `Ω=0.001 eV`:

- `σ^C_xy` vanishes to grid precision as `N→∞` (reported separately below).
- `σ^D_xx/Ω≈0.015234078` in the code's **natural-unit** 2D response conventions (not `A m^-1 T^-1`).
- `Im K_xx(Ω)=1.5234188588e−5` versus `σ^D_xx=1.5234078036e−5`, difference `1.11e−10` due to permitted `O(Ω³)` contribution.
- Doubling Ω doubles paper's low-frequency D; full Kubo deviation follows higher-order frequency behavior.

## D. Fermi-surface alternative to occupied-band curvature: PASS, selected model points

At T=200 K, compare paper Eq. (11) `∑ f_n Z_nm` to distinct Eq. (13) `−∑_n f'_n v_n^x σ_nn^y`, using the *analytically differentiated Hamiltonian*. At `μ=0.85 eV`, `N=201²`, the absolute discrepancy is `7.98×10⁻¹¹` natural 2D units. This is a distinct mathematical integrand with a distinct Fermi-surface sampling requirement.

## E. Potential internal inconsistency in paper Eq. (14): OPEN, high priority

Printed main Eq. (14) assigns antiunitary factor `η_T=-1` to both `σ^C` and `σ^D`. Taking `R=I` and time reversal `T` would enforce `σ_C = -σ_C`, hence **zero** in the nonmagnetic Rashba model. But (i) the main text and Fig. 1 explicitly report nonzero `σ^C_xy` for T-invariant Rashba (MPG `4mm1'`), (ii) Table I states Zeeman curvature Z is **T-even** while metric Q is **T-odd**, and (iii) SM §V employs distinct `eV²` and `aeV²` for the two tensors. Numerically, for Rashba, explicitly checked `σ_y H(k)^* σ_y=H(−k)`, nonzero angular average Z_xy and zero angular average Q_xx.

**Working hypothesis (independent, NOT attributed to authors)**: at this level, the antiunitary sign should be channel dependent (C even, D odd); a more complete AC tensor transformation may involve frequency reversal and complex conjugation. Do not silently repair Eq. (14); the model test is evidence of the mismatch that must be reconciled before magnetic-point-group exhaustive automation.

## Not yet performed

A fresh full wavepacket Lagrangian derivation of SM Eqs. (12)-(34), fresh lesser-Green's-function derivation of SM Eqs. (35)-(47), exhaustive independent magnetic-group tensor classification for SM Tables I-III, absolute laboratory unit conversion and direct figure pixel/digitized comparison. These remain explicit acceptance tasks.
