# Quantum Geometry Reproductions

Independent, **paper-by-paper** reproductions of quantum geometry and magnetotransport. Each paper has an isolated Python package, response protocol, tests, and scientific acceptance status. **CI PASS is not the same as full scientific reproduction.** These are not official author implementations.

## Paper index

| Paper | Available independent capability | Scientific status |
| --- | --- | --- |
| [Jia et al., PRL **136**, 226302 (2026): NMEE](papers/jia-2026-nonlin-magnetoelectric-edelstein/) | Kane–Mele/CuMnAs, spin response, printed vs full-static expressions, non-Abelian diagnostics, Wannier input preflight | **PARTIAL:** equilibrium vs driven DC remains open |
| [Xiang et al., PRL **134**, 116301 (2025): IGMC](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/) | Rashba/Dirac, Zeeman geometry, finite-frequency response, 122 magnetic point groups | **PARTIAL:** Fig. 2(d) amplitude and some conventions remain open |
| [Xiang, Jin & Wang, PRL **135**, 146303 (2025): Spin Transport Revealed by Spin Quantum Geometry](papers/xiang-2025-spin-quantum-geometry/) | Two-valley Dirac NDSC/reduced intrinsic Eq.13, independent SM integrals, finite-eta N-band Liouville and PT projector/U(2) method subset | **PARTIAL:** Fig. 1(a) absolute factor, metallic DC protocol and SM algebra remain open; initial curated core only |
| [Xiang et al., PRB **113**, L201406 (2026): Quantum Geometric Map of Magnetotransport](papers/xiang-2026-quantum-geometric-map/) | Dirac/Rashba, Liouville, spin/orbital responses, SCBA + ladder, four-channel spectral Kubo, finite-(q, omega) Ward, **exact finite-disorder linear Hall + twist Chern (v41)** | **PARTIAL:** full Kubo-Streda, physical dc continuity and real-material multiband acceptance remain open; this is a curated subset of the v30 archive |

The two Xiang papers are **different physics**: 2025 IGMC concerns magnetic-field-driven gyrotropic AC current; 2026 PRB concerns spin/orbital magnetotransport under electric and magnetic perturbations. Do not merge them simply because they share authors and quantum geometry terminology.

## Layout

```text
papers/
  jia-2026-nonlin-magnetoelectric-edelstein/         # NMEE, existing frozen module
  xiang-2025-spin-quantum-geometry/                    # Xiang PRL 135 (2025), curated Dirac paper core
  xiang-2025-intrinsic-gyrotropic-magnetic-current/   # IGMC, existing frozen module
  xiang-2026-quantum-geometric-map/                   # PRB magnetotransport, independent
studies/                                               # Cross-paper research
.github/workflows/                                     # Path-scoped CI for each module
```

## Run

Use Python 3.11; create a separate virtual environment per paper.

```bash
cd papers/jia-2026-nonlin-magnetoelectric-edelstein
python -m pip install -e .
python scripts/verify_release.py
```

```bash
cd papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
python -m pip install -e .
python scripts/verify_round9.py
```

```bash
cd papers/xiang-2026-quantum-geometric-map
python -m pip install -e .
python scripts/verify_release.py
```


```bash
cd papers/xiang-2025-spin-quantum-geometry
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/reproduce.py
```

See each module's README and [PRB physical integration contract](papers/xiang-2026-quantum-geometric-map/docs/INTEGRATION.md). The published PRB core uses its own import namespace `qgmap` and never modifies NMEE or IGMC frozen benchmarks.

## Reuse without overstating physics

- Choose the physical protocol first: intrinsic vs extrinsic, static equilibrium vs driven DC vs AC, magnetic-field preparation, temperature, relaxation, and degeneracy.
- Preserve source equations, units, signs, and convergence evidence. An independently generated model figure does not equal original publication data.
- Never silently replace the printed formula with an independently preferred correction. Label branches and keep disputed factors visible.
- For a real Wannier-material response, `wannier90_hr.dat` alone is **not sufficient**. Same-gauge projected positions/Berry connection, velocities, spin/magnetic moment operators, derivatives where necessary, occupations, and validated collision assumptions are required.
- Publisher PDF/SM, unlicensed figures, private research archives and unpublished author communications are not mirrored here.

## Workflows and licensing

See [.github/workflows](.github/workflows). Workflows are scoped to their own module. The existing [IGMC module's BSD-3-Clause license](papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/LICENSE) does not automatically license the entire repository. No root license is asserted.

See [RELEASE_STATUS.md](RELEASE_STATUS.md) for the cross-paper release status and outstanding scientific gates.


The PRB module now also includes an independently tested SCBA/ladder, four-channel spectral Kubo and finite-q, finite-ω resolvent-Ward method subset; see [PRB integration notes](papers/xiang-2026-quantum-geometric-map/docs/SCBA_KUBO_WARD_MIGRATION_V34.md). These do not constitute full Kubo–Středa or scientific acceptance of the original dc PHE.

The PRB module now also exposes the independently benchmarked v20 full-eigenstate longitudinal Kubo continuity/contact checks and SCBA eta limit-order diagnostics. These are *not* a replacement for a full electromagnetic Kubo–Středa calculation.

The PRB module also exposes v21–v23 second-order EB response methods (finite-q exact mixed Kubo, dynamic two-color orderings and homogeneous length-gauge Liouville with the SM Eq.20 pole). All are independently verified synthetic-model subsets, not full steady-state PHE acceptance.

The PRB reproduction also contains the original v24–v27 continuum Dirac, direct-SI Liouville and Fig. 2(b) source-curve diagnostics. The unresolved factor of approximately \(2\pi\) between Eq. (11) and the plotted ordinate is preserved as an OPEN issue, not silently corrected; see [v37 method notes](papers/xiang-2026-quantum-geometric-map/docs/CONTINUUM_SI_FIG2_MIGRATION_V37.md).

The PRB module now includes an independent **clean antisymmetric Hall Fermi-sea** benchmark against integer Chern numbers; it is separate from v18's symmetric dissipative Kubo and from the original paper's EB response. See [v38 scientific scope](papers/xiang-2026-quantum-geometric-map/docs/V38_CLEAN_STREDA_AUDIT.md).

The PRB module also includes the v39 energy-dependent SCBA derivative and RR density Ward benchmark, with separate finite-torus exact disorder comparisons. This is not a full disordered Hall conductivity or original mixed-field EB response. [v39 method limits](papers/xiang-2026-quantum-geometric-map/docs/SCBA_ENERGY_DERIVATIVE_V39.md).

The PRB module now provides a v40 same-disorder RA/AR/RR/AA **current-vertex** calculation, plus a Hall-zero check proving that the antisymmetric RA correlator alone is **not** a physical Hall conductivity. [Scientific scope](papers/xiang-2026-quantum-geometric-map/docs/CURRENT_VERTEX_V40.md).

The PRB v41 module adds the **full exact finite-torus linear Hall Lehmann response** and an independent disorder-compatible occupied-subspace **twist Chern** reference. The PRB CI tested 291 cases at the PR head. It remains **PARTIAL**: neither a conserving disorder-averaged Středa I/II decomposition nor the original nonlinear EB PHE is complete. [v41 physical protocol](papers/xiang-2026-quantum-geometric-map/docs/EXACT_HALL_V41.md).

The PRB v42 module adds the **exact finite-sample linear Hall Smrčka–Středa I(a)+I(b)+II decomposition** in the **same** v41 disordered QWZ torus and current convention. Independent inverse-matrix current traces, direct energy integration, v38 clean-Bloch sea, and the full Lehmann limit were verified; 10 new tests brought remote PRB GitHub Actions to **301/301 PASS**, run [39530869244](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39530869244). The standard surface/sea contributions cancel strongly and are *not* uniquely physical mechanisms. **Still PARTIAL:** conserving disorder-averaged Středa with matched self-energy/vertices, physical metallic DC, and the published second-order EB response remain open. See [v42 protocol](papers/xiang-2026-quantum-geometric-map/docs/EXACT_STREDA_V42.md).

The PRB **v43** study now checks the **matched noncrossing SCBA Středa I(a) *surface* current vertex only**, versus v42's independent frozen-disorder, twist-averaged finite-torus surface correlator using **the same W and bare eta**. At W=.12, mu=1.1, eta=.12, L4/24 seeds/7x7 twists: bare SCBA -1.27807, one-ladder -1.43032, exact -1.50028±.05502 e²/h, **a mismatch remains**. GitHub Actions [39533046760](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39533046760) passed **312/312** PRB tests (11 new). [v43 scope](papers/xiang-2026-quantum-geometric-map/docs/SCBA_SURFACE_V43.md). **Do not promote to complete conserving disorder-averaged Středa Hall or published EB response; overall PARTIAL.**
