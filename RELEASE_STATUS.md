# Quantum geometry releases (2026-10-20)

| Paper | Remote module | Software gate | Scientific status |
| --- | --- | --- | --- |
| Jia 2026 PRL NMEE | `papers/jia-2026-nonlin-magnetoelectric-edelstein` | `scripts/verify_release.py` | PARTIAL |
| Xiang 2025 PRL IGMC | `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current` | `scripts/verify_round9.py` | PARTIAL |
| Xiang, Jin & Wang 2025 PRL spin quantum geometry | `papers/xiang-2025-spin-quantum-geometry` | 11 paper-core tests, Python 3.11/3.12 CI, JSON smoke | PARTIAL; curated Dirac/SM subset; metal driven DC and paper absolute factors OPEN |
| Xiang 2026 PRB magnetotransport | `papers/xiang-2026-quantum-geometric-map` | `scripts/verify_release.py` | PARTIAL; GitHub **curated core + migrated nondegenerate spin/orbital method subset** |

For the PRB, the complete 30-round original **offline** v30 checkpoint is separate from the independently readable/callable core published in this repository. It includes additional model branches, scans, plots, and scientific logs; those are **not** claimed to have been uploaded here. The core implements selected independently checked methods and leaves absent advanced interfaces explicitly unimplemented.

Open physics:
- Spin QGT PRL 135: Fig.1(a) ~2× axis, SM Eq.55/56 and SM Eq.69/70 algebra, finite-eta DC interpretation; original multiband/PT methods remain offline for later migration.
- NMEE: static-versus-driven density-matrix response and CuMnAs PT-degenerate DC.
- IGMC: source Fig.2(d) normalization and response convention.
- PRB: Eq.(11) vs Fig.2(b) ~2π SI amplitude, physical steady-state PHE collision protocol, SM Eq.(52)/(56) internal algebra, multiband orbital response and Rashba QMQ attribution.

Do not infer a full-paper science PASS from the software CI result.

## 2026-10-20 advanced method migration

The PRB module additionally contains reusable four-band spin-MNHE (SM 16–19 vs 29), orbital PHE Eq.56, orbital MNHE 2D SBC and 3D Eq.49/52 diagnostics. These are callable via `python scripts/reproduce_advanced.py --case ...` and covered by independently scoped GitHub CI. The disagreement branches remain OPEN and are never re-labelled as full-paper scientific acceptance. The remaining v30 offline exact finite-system Kubo continuity/contact benchmark, limit-order scans, full Kubo–Středa terms and historical raw figures remain outside this remote subset.


## 2026-10-20 SCBA/Kubo/Ward migration (v34)

Migrated the original v11/v17/v18/v19 kernels and 48 regression tests: matched SCBA self-energy/ladder, RA control, four-channel symmetric spectral Kubo, and finite-(q,omega) resolvent Ward. Latest method scope: `papers/xiang-2026-quantum-geometric-map/docs/SCBA_KUBO_WARD_MIGRATION_V34.md`. Full Kubo–Středa, real-frequency physical continuity, and original paper's driven-DC PHE remain open.


## 2026-10-20 v35 exact Kubo / SCBA regulator migration

The original v20 exact finite-disorder real-frequency longitudinal Kubo continuity/contact tests and independent SCBA regulator-order scan are now in the PRB module, with 18 original tests and JSON entry points. Their PASS status is a mathematical/method validation, not acceptance of the original paper's transverse steady-state PHE. See `docs/EXACT_KUBO_LIMITS_MIGRATION_V35.md`.

## 2026-10-20 v36 second-order Liouville migration

Migrated unchanged original v21–v23 physics modules and 66 historical tests, plus frozen JSON numbers, source-file SHA256 manifest and three callable JSON cases. Finite-q scalar electric and homogeneous length-gauge electric remain distinct protocols. Original continuum steady-state DC PHE and complete Kubo–Středa remain OPEN.
