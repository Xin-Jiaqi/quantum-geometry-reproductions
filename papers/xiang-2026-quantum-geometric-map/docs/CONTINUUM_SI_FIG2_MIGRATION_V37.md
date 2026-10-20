# v37: Continuum tilted Dirac, SI Liouville and published Fig. 2(b) audit

**Paper:** Xiang et al., *Quantum Geometric Map of Magnetotransport*, PRB 113 L201406 (2026). **Scientific acceptance remains PARTIAL.**

This is a provenance-preserving migration from the frozen local v30 reproduction archive, not a new theory derivation. Four original source modules and four 2026 validation suites contribute 62 independent model tests.

## Physics protocols

1. **v24 continuum**: `continuum_ti_v24` implements main Eq. (10), with alternative SM Eq. (62) helicity available by explicit setting. (H=\beta k_x I+k_x\sigma_y-k_y\sigma_x), and `bbar` multiplies \(\sigma_x\) with \(\bar B_x=-g\mu_B B_x\). Independent direct eigensystem and analytic spectral projectors reproduce two-field Liouville currents. SM Eq. (20) is compared with a separate radial/angular analytic kernel.
2. **v25 literal SI**: `si_fig2_v25` distinguishes the printed `literal_hbar` formula from an explicit `single_2pi_diagnostic` branch, not an accepted correction. The Hall-voltage calculation additionally assumes the quoted resistance can be interpreted as sheet resistance and a 100 μm width.
3. **v26 SI Liouville**: `physical_si_liouville_v26` constructs the two-band current from matrix elements using energies in joules, momentum inverse meters, and current A/m. A dimensional cross-check against Eq. (11) supports the literal \(\hbar\) normalization, leaving the figure discrepancy.
4. **v27 figure shape**: `fig2_source_audit_v27` fits 20 K and 40 K *approximately raster-extracted* published curves with a Fermi tanh model. The CSV is derived coordinates, **not original author raw data**; no publisher image or PDF is mirrored.

## Commands

From this module after `python -m pip install -e .`:

```bash
python scripts/reproduce_advanced.py --case continuum-dirac-v24
python scripts/reproduce_advanced.py --case fig2-si-v25
python scripts/reproduce_advanced.py --case si-liouville-v26
python scripts/reproduce_advanced.py --case fig2-curve-v27
python scripts/verify_release.py
```

## Reference numbers

- v24: with beta=1/8, mu=0.2, model T=0.035, direct SM Eq. (20) `+0.002489935413104685` per \(\bar B\tau\); independent angular result `+0.002489935413104679`.
- v25: Eq. (11) literal zero-temperature ordinate approximately `11.0208361`; explicit `single_2pi_diagnostic` approximately `1.7540205`; Hall voltage literal conditional estimate near 268 μV versus 42.7 μV for the diagnostic.
- v26: physical SI direct current/E yields sheet conductance approximately `-2.68238e-7 S`; independent literal Eq. (11) magnitude `2.68238e-7 S`.
- v27: digitized plateau amplitude near `1.7495`; effective temperatures near 20.16 K and 39.95 K; pixel-calibration uncertainty must not be ignored.

## Method and release boundaries

These findings identify an unresolved amplitude inconsistency, not proof of its source. Angular versus cyclic Zeeman frequency is a candidate explanation only; no author plotting code is available. The fixed-(eta) Liouville coefficient is not automatically a collision-complete measurable steady-state DC PHE. Omitted orbital three-band and complete Kubo–Středa transport remain open. No author contact was performed.

For the original file-level provenance see `results/frozen/v37_manifest.json` (SHA256 and Git blob SHA1). The manuscript and publisher plots are not included in this GitHub repository.
