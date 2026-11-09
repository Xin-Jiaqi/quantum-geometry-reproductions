# Ezawa 2025 · Bulk photovoltaic effects in altermagnets

**Source**: Motohiko Ezawa, *Physical Review B* **111**, L201405 (2025), [DOI 10.1103/PhysRevB.111.L201405](https://doi.org/10.1103/PhysRevB.111.L201405). Independent, non-author model reproduction.

## What is reproduced (and what is not)

- **Eq. (9)**: 2D Rashba + d-wave altermagnet Hamiltonian with Néel vector along y and finite B sigma_z gap.
- **Eq. (7)/(11)**: symmetric **magnetic injection current (MIC)** for linearly polarized light, especially sigma^{x;xx}. Two-band optical matrix elements are computed from the quantum metric g_ab=Re(r_{-+}^a r_{+-}^b) and multiplied by the conduction–valence group-velocity difference. The k-space integral includes the T=0 Pauli occupation difference.
- **Eq. (17)/(18)**: sharp-delta resonant-contour results reproduce the analytic expression and plateau with **no fitted normalization**. Tests check J=0, J -> -J, lambda=0 inversion-symmetric limit, gap threshold, projector-based optical elements, k-space quadrature and an independent Cartesian Lorentzian integration.
- **Figures**: band asymmetry, resonant contour, MIC spectrum, magnetic reversal, and k-resolved kernel.
- **Not yet included**: shift current sigma_shift^{x;xy}, electron–hole interactions, true material BZ/Wannier data, finite-temperature occupation, calibrated relaxation-time/thickness SI spectra, or magnetic point-group automation. This package is not an experimental solar-cell efficiency prediction.

## Model and normalization

In model units hbar=e=tau=1, energies in epsilon0, wavevectors in k0, W k0=1:

~~~text
H(k) = a(kx^2+ky^2)I + (-lambda ky, lambda kx+J(kx^2-ky^2), B).sigma
C^{c;ab}(E) = sigma^{c;ab}/[tau e^3/(hbar^2 W k0)]
            = -(1/pi) integral d^2k (f_- - f_+) (partial_{kc}|h|)
                         g_ab delta(2|h|-E)
g_ab = (1/4)[partial_a h.dot.partial_b h / |h|^2
             - (h.dot.partial_a h)(h.dot.partial_b h)/|h|^4]
~~~

Here a=hbar^2/(2M) is dimensionless in epsilon0/k0^2; lambda in epsilon0/k0, J in epsilon0/k0^2; all are numerical quantities scaled by those units. Restoring the paper's SI units requires restoring k0, W, relaxation time tau, optical Fourier-amplitude convention and actual J/lambda **with its length dimension**. The scalar parabolic term affects the Pauli occupation boundary but drops from Delta-v and the geometric factors.

**Critical distinction from nlo-response-workflow**: its maintained injection_current.py computes the antisymmetric r^b r^c-r^c r^b channel (circular/helicity-related NIC). It necessarily vanishes for optical xx. The present code implements the independent **symmetric quantum-metric channel**. Do **not** swap the kernels or add this as an unreviewed change to the production Wannier pipeline. A general multiband adapter would need a distinct magnetic_injection protocol, magnetic symmetry/tensor filtering, optical matrix-element convention checks, actual orbital embedding, and n(k)/-n(k) magnetic-domain data.

## Run

~~~bash
python -m pip install numpy matplotlib
export PYTHONPATH=src
python -m unittest discover -s tests -v
python scripts/reproduce_core.py --outdir results
~~~

reproduce_core.py emits JSON to stdout even with no --outdir, for the parent repository's runner. results/ includes CSV, five PNG and SVG plots, and validation.json. No source PDF or private DFT data is included.

## Key independent discrepancy to retain

Eq. (19) in the published PDF states, in the same normalization,

~~~text
C^{x;yy}_printed = -(J/8lambda) [1 - 2(B/E)^2 + 4(B/E)^4].
~~~

Directly evaluating Eq. (7) with Eq. (9), using both g_yy and the resonant-contour integration, instead gives (verified numerically over multiple B, E, J in the single-contour region)

~~~text
C^{x;yy}_numeric = -(J/8lambda) [1 - 8(B/E)^2 + 16(B/E)^4].
~~~

Example B=0.1, E=0.4, J=0.005, lambda=1: printed -0.000556640625, contour -0.0003515625; relative difference **36.8421%**. A separate Cartesian finite-gamma integral approaches the contour result as gamma->0 and nk->infinity. This is an **OPEN_SOURCE_DISCREPANCY**: no author correction is claimed and the printed expression is preserved in code for audit. It does not affect the main x;xx reproduction of Eq. (17).

## Acceptance and future transfer

Scientific scope is **MODEL_READY (Ezawa continuum, clean limit)**, not **WANNIER_READY**. To transfer to a Wannier engine, use a separate sigma_mic^{c;ab} term proportional to Delta v^c * Re[r^a_nm r^b_mn], preserve Fermi factors and finite scattering tau, and compare the magnetic-domain odd part. Verify magnetic point group and P breaking before running. Only then normalize the 2D sheet output by an explicitly declared effective thickness if a bulk coefficient is required.
