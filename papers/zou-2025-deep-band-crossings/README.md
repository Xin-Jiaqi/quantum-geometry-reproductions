# Zou et al. 2025 — Nonlinear optical effects enhanced by deep band crossings

**Source:** Nianlong Zou *et al.*, *Physical Review Letters* **135**, 226901 (2025), [DOI 10.1103/sqwy-xhdf](https://doi.org/10.1103/sqwy-xhdf). Independent/non-author reimplementation from the published text and Supplemental Material (SM).

**Status: `MODEL_READY` — source-defined three-orbital TB, local SHG leading-singularity cancellation only.** This does *not* mean full SHG, full shift spectrum, GeTe, SOC/HSE or SI-calibrated BPVE are reproduced.

## Verified results

- Reimplements the Hermitian three-orbital tight-binding model printed in SM Sec. IV.5, with exact first/second momentum derivatives, mirror and spinless-time-reversal tests.
- Reproduces the six velocity ratios and band gaps in **SM Table S1** to the paper's displayed precision.
- Confirms the **SM (S10)** covariant-derivative sum rule by independent parallel-transport-gauge finite differences.
- Demonstrates that the leading near-degenerate **SM (S12)** Hermitian connection term scales as $1/q$ and approaches the full S10 result.
- Reconstructs **SM (S6) interband** and **SM (S19) mixed** SHG singular contributions, checking the $2\omega$ and $\omega$ groups independently. The two channels separately diverge as $1/q$; their sum remains bounded and the relative residual vanishes as $q\to0$.
- Verifies optical-index exchange, conduction-band relabeling and random U(1) eigenvector gauge invariance.

**Detailed step-by-step derivation (Chinese):** [DERIVATION_SHG.md](DERIVATION_SHG.md), with the checks of SM (S20)–(S26) and restrictions.

## Quick run

```bash
python -m pip install numpy matplotlib
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/reproduce_core.py --outdir results
```

The module is also callable via the parent repository's `python scripts/core.py run deep-dbc` and `python scripts/core.py test deep-dbc` once registered in `physics/capabilities.json`. The model recipe emits JSON to standard output even when no `--outdir` is supplied. Optional plots require Matplotlib.

### API

```python
from deep_dbc import nodal_parameters, band_data, covariant_derivative, shg_pieces
n = nodal_parameters(1.03)
d = band_data(1e-4, n['ky'] + 0.37e-4, 1.03)
p = shg_pieces(d.energies, d.berry, 0.7*n['epsilon0'], axes=(0,0,1))
print(p.interband, p.mixed_singular, p.residual)
```

`axes=(0,0,1)` represents $\chi^{xxy}$. Energy is measured in the SM TB energy unit and momentum in dimensionless reciprocal lattice units; the common $e^3/(2\hbar^2)$ and Brillouin-zone measure are deliberately **not** included in the local response budget.

## Scope and numerical caveats

1. This is an explicit audit of the **leading DBC-induced $1/q$ terms**. In particular, `mixed_singular` intentionally omits the bounded metric/Berry-curvature/velocity pieces of SM (S7). Thus `interband + mixed_singular` is **not** the full finite SHG susceptibility.
2. The approach to the nodal point is oblique and non-exceptional: $k_x=q$, $k_y=k_{y0}+0.37q$. Optical frequency is off-pole, $\hbar\omega=0.7\epsilon_0$, and $f_v=1,f_c=f_{c'}=0$. Exact degeneracies require a non-Abelian description.
3. No numerical linewidth or broadening is inserted into Berry connections; near-node tests use $q\ge10^{-5}$ to avoid floating-point ill-conditioning.
4. Earlier exploratory finite-BZ shift-current spectra require separate grid-convergence and physical-unit calibration. They are not silently promoted here to a quantitative Fig. 1(d) reproduction.
5. No original publisher PDF, SM PDF, private material data or credentials are bundled.

## Citation and status

Zou, N., Li, H., Ye, M., *et al.* (2025), “Nonlinear Optical Effects Enhanced by Deep Band Crossings,” *Phys. Rev. Lett.* **135**, 226901, DOI: 10.1103/sqwy-xhdf. For the implemented equations see SM (S6), (S10)–(S14), (S19)–(S26) and Table S1. Source-derived statements and independent numerical checks are separated in [DERIVATION_SHG.md](DERIVATION_SHG.md).

