# Round 7: exact Fermi-contour partition of the lattice Dirac IGMC

2026-10-20 · paper model Eq. (17) · NO publisher data, broadening, or fitted amplitude

## Definition and energy surface

Keep exactly the previous paper inputs `β=0.25 eV`, `v_F=0.825 eV Å`, `M=0.1 eV Å²`, `a=1 Å`, `g=2`, zero temperature, driving frequency **f=1 THz cyclic Hz**, so `ℏω=0.0041356676943898555 eV`. The Hamiltonian `H=d·σ` has `ε±=±E(k)`. Hence the occupation difference in the interband Kubo/geometric response is

\[
f_- -f_+ = \Theta(E(\mathbf k)-|\mu|).
\]

Unlike a uniform k mesh, the integration evaluates the occupied region *without any sampling of the Fermi step*. With `x=a k_x`, `y=a k_y`, `u=cos x`, and `B=2M/a²`, `A(y)=β+B(2−cos y)`, `w=v_F/a`,

\[
E^2(x,y)=(B^2-w^2)u^2-2ABu+A^2+w^2(1+\sin^2 y).
\]

For every fixed y the roots of `E²=µ²` are therefore obtained by a quadratic equation in u, **including two different crossings at the same y** (e.g. `µ=1 eV, y=0.5`). Map roots to `x=arccos u`, partition `[0,π]` into occupied/unoccupied intervals, and integrate the analytic two-band Pauli kernel only on occupied intervals using 64-point Gauss. In the y direction, first split at every root of `E(0,y)=|µ|`, `E(π,y)=|µ|`, and `E(x_stationary,y)=|µ|`; Gauss integration on each segment no longer crosses a Fermi-topology kink.

By kx↔−kx, ky↔−ky even parity for **C_xy** and **D_xx** of Eq. (17), integrate one BZ quadrant and restore `4/(2π a)^2`. This parity has explicit tests. This calculation does **not** require eigenvectors or a phenomenological relaxation term.

## Second mathematical route for C (distinct Fermi-line formula)

The paper's Eq. (13) is, at T=0,

\[
C_{xy}=-\sum_n\int\frac{d^2k}{(2π)^2}f'(\epsilon_n)v^x_n s^y_n
=\int\frac{d^2k}{(2π)^2}\delta(E-|\mu|)(\partial_{k_x}E)n_y.
\]

For each fixed y, replace the x integral with an **actual delta-function contour integral**, summing over *all* Fermi crossings with `sign(∂_{k_x}E)n_y`. This does not evaluate the occupied-sea Zeeman curvature `Z`; it is an independent implementation of the identity `∑_m Z_nm=∂_k σ_nn` and integration by parts. Both approaches share the analytic Fermi-boundary location, so this is an independent response-integrand / 1D integration cross-check, not an independent dispersion model.

## Metal convergence and units

Data in `results/raw/round7_metal_convergence.csv`; the primary magnitude is the published plotting convention `−g µ_B σ`, SI `A/(m T)`; the raw C/D kernels have Å-based units. All numerical estimates below are **not** original figure-digitization values.

| µ (eV) | C_xy [A/(m T)] | D_xx [A/(m T)] | status |
|---|---:|---:|---|
| 0 | 0 within numerical precision | +0.017754 | insulator; Round6 consistency |
| ±0.80 | ≈−20.47 | ≈−0.05675 | Fermi boundary integrated |
| ±0.95 | ≈−30.15 | ≈−0.06303 | Fermi boundary integrated |
| ±1.00 | **−33.48623534** | **−0.06417844743** | selected metallic benchmark |
| ±1.05 | ≈−34.72 | ≈−0.06372 | metallic peak region |
| ±1.10 | ≈−30.81 | ≈−0.05449 | Fermi boundary integrated |
| ±1.20 | ≈−17.23 | ≈−0.02942 | Fermi boundary integrated |

For the tested `µ = 0, 0.8, 0.95, 1.0, 1.05, 1.1, 1.2 eV`, the change between outer orders `160` and `224` is bounded by

- `max |Δ C_xy| = 1.063 × 10^-9 Å^-1`
- `max |Δ (D_xx/ℏω)| = 3.436 × 10^-10 Å^-1 eV^-1`
- `max |C(Eq11) - C(Eq13)| = 5.35 × 10^-10 Å^-1` against independently integrated contour delta term (448-point Gauss outer per topology segment).

These are **resolution comparison residuals**, not strict mathematical upper bounds on integration error. `inner order 32 → 64` also changes the fixed-metal integrals by less than `1e−10` in tested cases. Finite-frequency denominator changes `C` and `D` by about `3×10⁻⁶` relative to their leading low-frequency terms at `µ=1 eV`, with frequency kept non-resonant.

Eq.(17) has exact particle–hole symmetry. The numeric results at `µ` and `−µ` coincide exactly (same θ occupancy); no independent prediction of experimentally asymmetric metallic bands is implied.

The scan `|µ|=0..1.45 eV` at listed *discrete points* gives a large conduction amplitude at `|µ|≈1.05 eV` and displacement at `|µ|≈1.01 eV`. These are scanned sample maxima, **not claims of mathematically optimized peak positions**. Unlike old uniform sampling, the Fermi surface is represented exactly for each integration y slice.

## What this does not resolve

In the plotted Fig. 2(d), the paper visually appears near twice the displacement amplitude of the independently confirmed main-text Eq. (12), even after rigorous metal convergence. `µ=±1` corresponds to `−0.06418 A/(m T)` on the **validated** main-text branch; the raw authors' graph data and implementation settings remain unavailable. The ×2 SM Eq. (34) *literal intermediate* must remain a **separate diagnostic**, not a corrected physical output. In-gap amplitude is `+0.01775389`, not the approximate visually read `+0.04`.

Complex-frequency/antiunitary author's Eq. (14) remains a source-interpretation question; Table II group `−42′m′` lacks a printed `D_xy=D_yx` equality. `SM §VI` experimental Faraday-field applicability has not been quantitatively accepted, and no real material/Wannier numerical result is attempted.

## Executable commands

```bash
cd papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
PYTHONPATH=src python scripts/audit_round7.py --no-write
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/verify_release.py
```
