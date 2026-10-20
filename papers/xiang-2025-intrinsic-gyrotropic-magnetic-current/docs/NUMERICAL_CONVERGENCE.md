# Numerical convergence · round 1

All below are measured numbers from the local Python implementation. The source figures are not available as machine-readable numerical arrays. Code convention: model energies eV, momentum Å^-1, e=ℏ=1 for response kernels. These raw coefficients must **not** yet be compared to the published A m^-1 T^-1 axes.

## Rashba continuum Eq. (15)/(16), direct polar quadrature

Angle samples: 96, radial midpoint samples 1200, `kmax=0.5 Å^-1` (covers all Fermi surfaces of the quoted 3 μ points). `μ=−0.015 eV`: analytic `0.00896277841186 Å^-1`; at 64 angular samples the direct result was `0.00895246554892 Å^-1` (absolute `1.03e-5`, relative `0.115%`). `μ=0` and `+0.009 eV`: analytic `0.01378513850357 Å^-1` and direct at 64 angles `0.01379342840130 Å^-1` (absolute `8.29e-6`, relative `0.0601%`). Change of nphi from 64→96 cannot remove radial step error, and the radial mesh is a convergence parameter, not a figure-fit parameter.

Rashba `k=0` is degenerate: heatmaps deliberately mask `|k|<0.023 Å^-1` only for visual readability; the response integral uses an unmasked polar midpoint rule and no arbitrary broadening. A finer radial convergence series remains to be archived before unconditional Fig. 1(c,d) singular-hotspot acceptance.

## Dirac lattice fixed parameters

`β=0.25 eV`, `a=1 Å`, `vF=0.825 eV Å`, `M=0.1 eV Å²`. Periodic midpoint mesh, `k↔−k` matching. Results are `σC_xy` and `σD_xx/Ω` in code natural-unit conventions:

| `N×N` | `μ=0`, C | `μ=0`, D/Ω | `μ=0.85 eV`, C | `μ=0.85 eV`, D/Ω |
|---:|---:|---:|---:|---:|
| 51² | +3.8842e−8 | +0.01523422 | see CSV | see CSV |
| 101² | ≈0 at machine precision | +0.015234078 | −0.08055233 | −0.05036566 |
| 151² | ≈0 | +0.015234078 | −0.08055553 | −0.05038049 |
| 201² | ≈0 | +0.015234078 | −0.0806328138 | −0.0503994566 |
| 301² | ≈0 | +0.015234078 | −0.0807804136 | −0.0504206589 |

For `μ=0.85 eV`, `201²` vs `301²`: conduction absolute diff `1.4760×10^-4` / relative `0.1827%`; displacement-per-Ω absolute diff `2.1202×10^-5` / relative `0.0421%`. Near `μ=0.3 eV`, D changes sign/cancels across k and the relative error becomes large even though its absolute value is small; report absolute error and do not declare such µ converged from relative tolerance alone.

Raw data: `results/raw/dirac_convergence.csv`. Fig2 model image drawn with N=181²; separate better-converged reference points N=301² support the main μ=0,0.85 regimes only. A peak-local adaptive scheme, digitization and SI audit remain open.

## Frequency cross-check

With µ=0, Dirac `N=101²`, Ω=0.001 eV: low-frequency `σD=1.5234078036×10⁻⁵`, exact `Im K=1.5234188588×10⁻⁵`, absolute `1.1055×10⁻¹⁰`. A full frequency/η order-of-limits experiment is still required.

## Round-2 SI convergence (no fitted response amplitude)

For `g=2`, `f=1e12 Hz`, `hbar omega=0.004135667694 eV`, `mu=0`: `N=101,201,301,401` gives in-gap `Dxx=0.01775389321786, 0.01775389321777, 0.01775389321777, 0.01775389321777 A/(m T)`. Full finite-frequency Kubo gives `0.01775609692115`; difference `0.0124%`, tiny compared to Fig2(d) visual value around `0.04`. `mu=1 eV`, N=201 -> `-0.06416676`, N=401 -> `-0.06418607 A/(mT)`. Thus mismatch cannot be attributed to grid or the leading finite-ω correction at these points. Full csv: `results/raw/round2_SI_convergence.csv`. This does not settle non-smooth Fermi-edge regions elsewhere.
