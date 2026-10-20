# Callable physical interfaces and fail-closed material input contract

Import from `qgmap`, **not** from the other two independent modules. This module deliberately avoids a speculative universal NMEE/IGMC/magnetotransport solver.

| Request | Available API | Units, boundary |
| --- | --- | --- |
| Dirac `H(k)`, bands | `qgmap.dirac.hamiltonian`, `energies` | energy eV, k Å^-1, `vF_eVA` explicitly supplied |
| Zeeman geometry | `local_zbc_analytic`, `local_zbc_eigensolver` | no k=0 nondegenerate formula |
| Printed Eq11/SM63 | `eq11_kernel`, `eq63_literal_integral` | separate, sign-disagreeing protocols |
| Thermal Dirac | `continuum_temperature_kernel` | ideal unbounded clean Dirac; not a full transport calculation |
| Rashba OHE | `rashba.zero_temperature`, `direct_fermi_ring_sum` | e=ħ=1, constant τ, T=0, outside band-edge/crossing |
| Classical OHE | `orbital.ohe_parabolic_2d` | trivial one-band model, not multiband quantum geometry |
| Two-color density matrix | `liouville.two_color_response` | finite matrices, explicit `V_E,V_B,η_E,η_B,ω_E,ω_B` |
| Magnetic groups | `symmetry.table_i_acceptance` | 11 explicit generator settings, not all 122 MPG |

## Required for material transfer

Do NOT infer a quantitative magnetic response from `wannier90_hr.dat` alone. Stage-1 read-only contracts for other workflows can help validate structural dimensions, SOC and checksum completeness, but physical inputs also include:

- lattice and reciprocal units, dimensionality and occupancy;
- same-gauge Hamiltonian H(k), position/Berry connection `A_a(k)`, velocity `D_a H`;
- projected spin and/or orbital moment `S_i(k)`, `m_i(k)`, and covariant derivatives needed by the target response;
- specified field-frequency/collision/thermalization protocol and near-degenerate subspace treatment.

**Without these, real-material magnetotransport is BLOCKED.** No mock HR-only real-material result is claimed. The protocol distinction between `pre-equilibrated-B` and `adiabatic-B` is not an implementation detail; it can change observables.

## Cross-paper reuse

NMEE, IGMC and PRB have different perturbations and observables. Reuse common linear-algebra, integration, tensor symmetry and provenance helpers only after agreeing on observables, units, gauge and limiting procedure. Never copy a response coefficient across papers because the API name sounds similar.

## Newly migrated advanced physical methods (v31–v33)

| Calculable method | Python function | Independent physics gate |
| --- | --- | --- |
| Four-band intrinsic spin MNHE, SM (16–19) | `qgmap.spin_multiband_v5.four_band_audit(N=23)` | sum of four density-matrix currents vs independently evaluated SM Eq.(29) |
| 2D orbital PHE, SM Eq.(56) | `qgmap.orbital_phe_v4.orbital_phe_eq56(N=61)` | literal first/second lines and independent finite-field Boltzmann derivative |
| Strict 2D orbital MNHE SBC | `qgmap.orbital_mnhe_v5.orbital_mnhe_2d(N=61)` | orbital-moment Fermi-surface term vs squared Berry curvature |
| 3D two-band orbital MNHE, SM Eq.(52) | `qgmap.orbital_mnhe_v5.orbital_mnhe_3d(N=17)` | positional-shift Eq.(48)/(49) versus literal Eq.(52), factor-two issue retained |

Each function returns a `protocol` label and normalized model-unit results. The synthetic spin and orbital protocols must **not** be combined into a physical SI coefficient without separately defining the perturbation, magnetic moment units, collision model, dimensionality and band/occupancy convention.

Unified JSON commands (run from this paper module after `pip install -e .`):

```bash
python scripts/reproduce_advanced.py --case spin-four-band
python scripts/reproduce_advanced.py --case orbital-phe
python scripts/reproduce_advanced.py --case orbital-mnhe-2d
python scripts/reproduce_advanced.py --case orbital-mnhe-3d --output results/local-orbital3d.json
```

A real Wannier input pipeline is still **BLOCKED** until consistent magnetic/position operators and response limits are provided. The original v30 SCBA and Ward code paths are not yet published in this repo.
