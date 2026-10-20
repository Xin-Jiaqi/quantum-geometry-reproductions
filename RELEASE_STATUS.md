# Quantum geometry releases (2026-10-20)

| Paper | Remote module | Software gate | Scientific status |
| --- | --- | --- | --- |
| Jia 2026 PRL NMEE | \`papers/jia-2026-nonlin-magnetoelectric-edelstein\` | \`scripts/verify_release.py\` | PARTIAL |
| Xiang 2025 PRL IGMC | \`papers/xiang-2025-intrinsic-gyrotropic-magnetic-current\` | \`scripts/verify_round9.py\` | PARTIAL |
| Xiang 2026 PRB magnetotransport | \`papers/xiang-2026-quantum-geometric-map\` | \`scripts/verify_release.py\` | PARTIAL; GitHub **core-code subset** |

For the PRB, the complete 30-round original **offline** v30 checkpoint is separate from the independently readable/callable core published in this repository. It includes additional model branches, scans, plots, and scientific logs; those are **not** claimed to have been uploaded here. The core implements selected independently checked methods and leaves absent advanced interfaces explicitly unimplemented.

Open physics:
- NMEE: static-versus-driven density-matrix response and CuMnAs PT-degenerate DC.
- IGMC: source Fig.2(d) normalization and response convention.
- PRB: Eq.(11) vs Fig.2(b) ~2π SI amplitude, physical steady-state PHE collision protocol, SM Eq.(52)/(56) internal algebra, multiband orbital response and Rashba QMQ attribution.

Do not infer a full-paper science PASS from the software CI result.
