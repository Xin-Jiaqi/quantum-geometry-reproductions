# v36 migration: second-order EB Liouville, gauge contacts and homogeneous length gauge

**Scientific status: PARTIAL.** This migration imports the **unchanged original** v21, v22 and v23 method source and their 66 original tests from the frozen local v30 research checkpoint. File-level SHA256 provenance is recorded in `results/frozen/v36_source_provenance.json`.

## Protocols are distinct

| Model | Source module | What is independently verified | Boundary |
| --- | --- | --- | --- |
| finite-q static Zeeman derivative | `electromagnetic_eb_v21.py` | exact finite-disorder Kubo central B difference vs Liouville derivative, including derivative of thermal equilibrium density, mixed charge Ward and pure-gauge current contact | finite-q scalar potential at finite frequency, not uniform E or experimental dc PHE |
| finite-q two-color | `two_color_exact_v22.py` | B-then-E and E-then-B contributions, eigenbasis resolvent vs Sylvester solution and direct time-domain evolution | separately switched finite-q scalar potential and uniform onsite Zeeman; eta values are adiabatic regulators, not a scattering time |
| homogeneous E length gauge | `homogeneous_length_gauge_v23.py` | periodic covariant derivative, numerical FFT vs independent algebraic derivatives, SM Eq. (20) leading 1/eta pole and explicit equilibrium population preparation difference | independent lattice model; not the published continuum Dirac Fig. 2, no conserving microscopic collision integral |

## Reproduce the validated method subsets

From `papers/xiang-2026-quantum-geometric-map` after installation:

```bash
python -m pip install -e .
python scripts/reproduce_advanced.py --case eb-finite-q-v21
python scripts/reproduce_advanced.py --case eb-two-color-v22
python scripts/reproduce_advanced.py --case eb-uniform-v23
python scripts/verify_release.py
```

The CLI outputs explicit `protocol`, `scientific_status=PARTIAL`, model/unit conventions, and no high-dimensional density matrices. The original numerical routines retain their complete matrix diagnostics when imported in Python.

## Independent reference values (model units)

- v21: derivative mixed current response: `-0.5013185224459267 - 0.17904312452188173j`; B-centered finite difference differs by approximately `8.0e-8`; exact mixed charge Ward residual around `1e-15` and Zeeman derivative of a pure-gauge current response around `1e-13`.
- v22: B-then-E: `0.559090353576692 + 0.04421258246073823j`; E-then-B: `-0.5671831780252455 - 0.07102534021631202j`; sum: `-0.00809282444855352 - 0.026812757755573707j`; independent Sylvester difference about `2e-15`.
- v23: uniform-E full kernel: `0.03254037287627427 + 0.01059782201791099j`; independent algebraic error about `1.1e-9`; Eq. (20) pole residue: `0.002513459754146514` in model units.

All quoted numbers are reproduced from **local model runs**, not publication SI conductivities.

## Unresolved scientific tasks

1. Match the continuum tilted Dirac Liouville coefficient to a microscopic conserving *physical steady-state* transport protocol. Distinguish re-equilibrated static B from dynamically driven B.
2. Resolve literal Eq. (11) to Fig. 2(b) amplitude discrepancy (about 2π) without inserting a fitted scaling factor.
3. General orbital three-band QMQ/BCQ, full Kubo–Středa/Bastin and complete real-material same-gauge spin/orbital magnetic-moment operators remain unresolved.

A successful software CI run demonstrates repeatable **method-level** benchmarks only, not a full-paper scientific PASS.