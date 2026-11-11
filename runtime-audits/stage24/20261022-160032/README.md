# Stage 24 — Local Runtime Audit (PR #66 review checklist)

Repository: Xin-Jiaqi/quantum-geometry-reproductions
Branch: audit/stage24-20261022-160032 (PR #67 target: main; no merge; main untouched)

Contents:
- STAGE24_LOCAL_RUNTIME_AUDIT.md        — full scientific report (Chinese/English tables)
- STAGE24_SCIENTIFIC_REPORT.json        — machine-readable results
- STAGE24_HALDANE_SYMMETRY_GATE.json    — C3/C3h/m_z Hamiltonian gate results
- STAGE24_NLO_RUN_ONCE.json             — one-shot NLO shift-current run record
- STAGE24_SHIFT_CURRENT_COMPARISON.json — independent re-implementation vs NLO
- STAGE24_HALDANE_PACKAGE.json          — synthetic package mapping
- stage24_reimpl_raw_band_sum.npy       — re-implementation raw weights (8x401)
- stage24_photon_eV.npy                 — photon grid (401, 0..2 eV)
- spectra/                              — NLO run output (spectrum.csv, run-manifest.json, summary.txt)
- scripts/                              — reproducible Stage 24 scripts
- logs/                                 — command logs (redacted)

Verdict: PASS_WITH_UNRESOLVED_CONVENTIONS. Science repos untouched; dirty trees preserved.
