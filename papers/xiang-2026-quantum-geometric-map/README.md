# Xiang et al. (2026): Quantum Geometric Map of Magnetotransport

**Reference:** L. Xiang, J. Jia, F. Xu and J. Wang, *Physical Review B* **113**, L201406 (2026). [DOI: 10.1103/kd12-4ycj](https://doi.org/10.1103/kd12-4ycj).

**Scientific status: PARTIAL.** Independent reproducible calculations; not an official author implementation.

## GitHub scope versus full offline research

This repository module publishes **directly importable, tested core capabilities** from the v30 research program, not the entire 421-file v30 checkpoint. The complete offline archive (including historical 30-round numerical audits, detailed response branches, and figures) remains separate. The CI here certifies only the particular code checked in, not every historical research result.

Supported:

- Tilted Dirac band spectrum, Hamiltonian conventions (main Eq. 10 vs SM Eq. 62), Zeeman Berry curvature via analytic and gauge-invariant matrix elements, Eq. (11) and literal SM Eq. (63), thermal-step and finite-field Drude diagnostics.
- Rashba orbital Hall response Eq. (9) at `T=0`, with distinct one-ring/two-ring occupancy sectors and independent Fermi-ring calculation.
- Classical orbital Hall effect (spinless parabolic band, finite temperature).
- Finite-dimensional **two-color Liouville resolvent** with explicit `E→B` and `B→E` orderings, independent frequencies, broadenings, and a **separately labelled** equilibrium population correction. This is a method module, not the complete published driven-DC formula.
- SM Table I 11 representative magnetic-group tensor spaces using Neumann operators.
- SCBA self-energy and a noncrossing Bethe-Salpeter ladder using the same scalar-disorder variance.
- Symmetric four-channel spectral Kubo-Greenwood response with separate RA, AR, RR, and AA ladders.
- Finite-(q, omega) resolvent Ward identity with the exact finite-difference longitudinal insertion.
- SCBA self-energy and noncrossing Bethe-Salpeter ladder in an independent periodic two-band regulator.
- Symmetric four-channel spectral Kubo-Greenwood response (`RA+AR-RR-AA`) with a distinct ladder in each channel.
- Finite-((q,\omega)) SCBA resolvent Ward identity using the exact finite-difference longitudinal insertion `[H(k+q/2)-H(k-q/2)]/q`.

## Start

```bash
python -m pip install -e .
python scripts/verify_release.py
```

```python
from qgmap.dirac import Parameters, eq11_kernel, local_zbc_eigensolver
from qgmap.rashba import zero_temperature, direct_fermi_ring_sum
from qgmap.liouville import two_color_response
print(eq11_kernel(1/8))
print(zero_temperature(1.0, 0.4, -0.04))
```

The source and benchmarks use explicit model units; they are not automatically a SI transport or materials solver.

## Scientific limits

- **Fig. 2(b):** temperature dependence verified in the original research; original plotted plateau is about `2π` below the literal SI Eq.(11) result. Do not fix by arbitrary scaling.
- **Steady-state PHE:** adiabatic Liouville pole and equilibrium/static full response differ under some preparation/relaxation protocols.
- **Orbital MNHE/PHE:** SM Eq.(52)/(56) discrepancies remain open. General multiband quantum geometry, omitted three-band terms, SCBA/Středa full transport, and first-principles calculations are **not** established by this published core.
- **Rashba OHE:** model calculation is available, QMQ attribution as a whole-paper statement remains open.
- Degeneracies and response limits must not be treated as covered by nondegenerate expressions.

See [scientific acceptance](docs/SCIENTIFIC_ACCEPTANCE.md), [integration contract](docs/INTEGRATION.md), and [provenance](docs/PROVENANCE.md). Do not redistribute publisher PDF, original figure images, or private logs.

## Migrated advanced-method subset (v31)

The public core now includes independently callable spin-MNHE, orbital PHE/MNHE, SCBA/ladder, four-channel spectral Kubo, and finite-(q, omega) resolvent-Ward method subsets. See the [v31 spin](docs/SPIN_MULTIBAND_MIGRATION_V31.md), [v32 orbital](docs/ORBITAL_MIGRATION_V32.md), and [v34 SCBA/Kubo/Ward](docs/SCBA_KUBO_WARD_MIGRATION_V34.md) protocols. These are model-level methods, not a general degeneracy-robust real-material solver.

### Orbital advanced-method subset (v32)

The migrated `qgmap.orbital_phe_v4` and `qgmap.orbital_mnhe_v5` modules independently check 2D/3D model responses and preserve the literal SM Eq. (52)/(56) and finite-field comparison branches. [Methods and limitations](docs/ORBITAL_MIGRATION_V32.md). They **do not** implement the paper's omitted general multiband three-band orbital terms.


### Disorder/Kubo/Ward subset (v34)

The original v11/v17/v18/v19 SCBA, ladder, four-channel spectral Kubo, and finite-((q,\omega)) resolvent-Ward modules are now migrated with 48 original regression tests. See [the v34 protocol and boundaries](docs/SCBA_KUBO_WARD_MIGRATION_V34.md). This remains a method subset: it is not full Kubo–Středa, physical finite-(eta) continuity, or the paper's driven-DC PHE.


### Exact Kubo and limit-order subset (v35)

The original v20 finite-disorder full-Lehmann Kubo continuity/contact-term benchmark and the separate SCBA external-regulator limit-order scan are now migrated with their 18 original tests. Two new JSON cases are `exact-kubo-v20` and `limit-order-v20`. See [v35 methods, reference numbers and physical limits](docs/EXACT_KUBO_LIMITS_MIGRATION_V35.md). Neither constitutes a complete transverse or mixed-field DC response.
