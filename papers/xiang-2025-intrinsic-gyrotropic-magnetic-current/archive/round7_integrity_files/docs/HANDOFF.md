# IGMC handoff · round 7 · 2026-10-20

Status **PARTIAL**. New metallic C/D Fermi-contour and independent Eq13 Fermi-line benchmarks accepted in `docs/ROUND7_METALLIC_FERMI_CONTOUR.md`; final scope matrix in `docs/FINAL_ACCEPTANCE_MATRIX.md`. Original Fig2(d) amplitude and Eq14 source convention remain OPEN.

Absolute path: `/mnt/data/igmc_reproduction_round7/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`. Repo-relative path: `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`.

```bash
cd /mnt/data/igmc_reproduction_round7/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
PYTHONPATH=src python scripts/audit_round7.py --no-write
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/verify_release.py
```

Inputs: `parameters_round7.json`; source `src/igmc_lab/fermi_contour.py`; outputs `results/raw/round7_metal_*.csv`, `results/raw/round7_metal_audit.json`, `results/figures/round7_*`; checkpoint `results/frozen/round7_metal_reference_v1.json`. Older frozen standards unchanged; `results/frozen/ROUND6_RELEASE_MANIFEST.json` preserves previous manifest bytes.

Key numbers at ±1 eV: C_xy `-0.11883230822771439 Å^-1` = `-33.48623534315167 A/(m T)`; D/omega `-0.05506958188347791 Å^-1 eV^-1`; D `-0.06417844742657552 A/(m T)`. 160→224 outer order maximal tested C change `1.062324506e-9`, D change `3.435038087e-10`; independent Eq13 FS maximum mismatch `5.349639898e-10`. Wave frequency cyclic Hz f=1e12, ℏω=0.004135667694 eV. No g/omega adjustment, no empiric factor 2.

**No author email. No GitHub remote push.** Continue only after deciding whether to seek original Fig2(d) source/physics clarification or stage an explicitly PARTIAL code release. Toolkit 0.18.0 pinned source, native library still skipped; tests show skipped reason.

---

# IGMC reproduction handoff · round 6 · 2026-10-20

**Scientific status:** `PARTIAL`, `MAIN_EQ12_INDEPENDENT_MULTIROUTE_PASS`, `FULL_COMPLEX_TR_QUADRATURE_PASS`, `FIG2D_AUTHOR_PLOT_OPEN`. SM Table I/II source audit from round 5 kept unchanged (243/244 exact; one `-42'm'` independent-equality omission).

**Execution directory:** `/mnt/data/igmc_reproduction_round6/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/` (repo-relative `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`).

```bash
cd /mnt/data/igmc_reproduction_round6/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
PYTHONPATH=src python scripts/audit_round6.py --no-write
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/verify_release.py
```

**New physically independent work:** `src/igmc_lab/d_vector_independent.py`; `scripts/audit_round6.py`; `tests/test_round6_dvector.py`; explanation `docs/ROUND6_DVECTOR_FREQUENCY_AUDIT.md`; data `results/raw/round6_*`; frozen `results/frozen/round6_dvector_reference_v1.json`.

**Reproducibility constants:** `f=1e12` cyclic Hz, `ℏω=0.0041356676943898555 eV`, `g=2`, Dirac model caption constants. `Gauss128² D/ω=0.015234078015630925`, `midpoint301²=0.015234078035703156` in `Å⁻¹ eV⁻¹`. Restored D= `0.01775389319438161 A/(m T)`; source figure visually close to ~`0.04`, not original numerical data. Complex TR and old spectral Kubo residuals `1.73e−17` and `5.56e−17 Å⁻¹`.

**No author email / no remote push.** The author-oriented physics questions are only drafted in `docs/AUTHOR_DISCUSSION_BRIEF.md`; do not send without user review. Do not overwrite any earlier frozen values. If continuation is needed, first prioritize resolving original Fig2(d) source data/explicit normalization; then metallic occupation-boundary convergence; attempt external live-toolkit API smoke test only if installed.

---

# IGMC reproduction handoff · round 5 · 2026-10-20

**Paper status:** `PARTIAL`. **New source-level result:** `244 source vs derived C/D tensor subspaces: 243 PASS, 1 AUTHOR_TABLE_CONSTRAINT_OPEN`. Table I all 122 group cases match; Table II 121 group cases match, the lone source group `-42'm'` lacks a published equality `D_xy=D_yx` for its two otherwise-allowed entries. See `docs/ROUND5_SM_TABLE_SOURCE_AUDIT.md`; this does **not** change SM Table III allowed classification.

**Executable source directory:** `/mnt/data/igmc_reproduction_round5/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/` (repo relative: `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`).

**Current strict commands:**

```bash
cd /mnt/data/igmc_reproduction_round5/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
PYTHONPATH=src python scripts/audit_round5.py
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/verify_release.py
```

**Input and output locations:** `src/igmc_lab/paper_table_source.py` manually transcribes pp.8–10, `data/toolkit_mpg_generators_v018.txt` pins upstream v0.18.0 from `551a0daf...`, `results/raw/round5_source_table_244.csv`, `results/raw/round5_source_table_audit.json`, `results/frozen/round5_source_table_reference_v1.json`. SM photo/PDF is external to module and not redistributed.

**Keep unresolved:** Eq.(14) paper-printed antiunitary sign vs channel-resolved response, SM Eq.(34) mid-derivation factor, original author Fig.2(d) physical SI amplitude. Earlier frozen baselines must not be mutated. Live original toolkit API unavailable without pin-install. Do not contact authors or push code without a later clear publication decision.

---

## Previous Round 4 and earlier handoff (preserved)

**Round-4 v2 addendum (2026-10-19):** 244/244 full magnetic-group character averages `d=(1/|M|)Σ η^θ det(R)(tr R)^2` match independent 9D Neumann/SVD tensor dimensions. CSV now contains `C_char_dim`, `D_char_dim`, plus full canonical tensor shapes JSON. The new frozen checkpoint is `results/frozen/round4_mpg_reference_v2.json`. 40 unit tests: 39 pass and one intentional external native-toolkit skip. Run `PYTHONPATH=src python scripts/audit_round4.py` and `PYTHONPATH=src python scripts/verify_release.py` from the **paper directory**. Full paper remains PARTIAL; no remote push or author email.

# IGMC reproducibility handoff · latest round 4 · 2026-10-19

**Global status:** `PARTIAL`; **new gate:** `122_MPG_CLOSED; 244_C_D_TENSOR_SUBSPACES_CHECKED; SM_SOURCE_TABLES_PARTIAL; FIG2D_OPEN`.

**Repository-relative path:** `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`

**Local source path:** `/mnt/data/igmc_reproduction_round4/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`

**External upstream:** `Xin-Jiaqi/group-theory-operations-toolkit`, v0.18.0, commit `551a0daf8d4ba5a5f24b2811694bfe438616de48`. Read via authenticated GitHub connector, not locally cloned (network unavailable).

## Run next time

```bash
cd /mnt/data/igmc_reproduction_round4/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
PYTHONPATH=src python scripts/audit_round4.py
PYTHONPATH=src python scripts/verify_release.py
# In another development environment only after toolkit v0.18.0 installed at pinned commit:
PYTHONPATH=src python scripts/audit_round4.py --require-live-toolkit
```

**Results and evidence:** `docs/ROUND4_TOOLKIT_SYMMETRY_AUDIT.md`, `results/raw/round4_mpg_122.csv`, `results/raw/round4_mpg_audit.json`, `results/raw/round4_mpg_tensor_shapes.json`, `results/frozen/round4_mpg_reference_v2.json`. There are **58 C-allowed MPG and 58 D-allowed MPG** among 122. For C: longitudinal-only 24, transverse-only 18, mixed 16; for D: 23, 19, 16. `4mm1′` C 1/D 0, `4m′m′` C 1/D 2; monoclinic y/z conventional-setting test passes. `paper_printed_eq14` differs from channel-resolved C dimension for 58 groups.

**Preservation:** older rounds 1–3 references untouched and verified; old scientific question Fig.2(d) remains open. Generated CSV/JSON contain **model calculations**, not publisher data. Original SM PDF and figure crops kept outside public paper subtree. No GitHub commit/CI run or author email attempted.

**Next priorities:** (1) run full 244 native-toolkit API checks in a network-enabled/cached toolkit environment; (2) finish every symbolic SM Table I/II tensor equality in the source convention; (3) obtain publisher source Fig.2(d) raw author data or prepare question for user approval; (4) keep distinct `paper_printed_eq14` and `channel_resolved_symmetry` until complex-frequency conventions settle.

---

## Previous Round-3 handoff (preserved verbatim)

# IGMC reproducibility handoff · latest round 3 · 2026-10-19

**Status:** `MAIN_EQ12_DISPLACEMENT_INDEPENDENT_PASS; SM_EQ34_FACTOR2_SOURCE_ISSUE; FIG2_D_ABSOLUTE_OPEN; FULL_PAPER_PARTIAL`.

**Repository-relative path:** `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`.

**Local source path:** `/mnt/data/igmc_reproduction_round3/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`.

```bash
cd /path/to/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
python -m pip install -e .
python scripts/audit_round3.py
python scripts/verify_release.py
```

**New science verified:**

- Main Eq.(12) `σ_D/ω=∫∑fn 2 Re(r_nmσ_mn)/(ε_m−ε_n)` under `H1=-gμ_B B σ` agrees with finite-field Wilson-loop polarization derivative. N301²: Wilson `0.015233103766...`; main `0.015234078036...`, relative error `6.40e-5`. Shift parameter `b=1e-4 eV`; smaller/larger b checked; source model exact.
- Independently solved **orbital-basis** Liouville 4×4, frequency 1 THz `ℏω=0.004135667694 eV`, N101²: imaginary current kernel `6.3010904648973e-5 Å^-1`, spectral Kubo `6.3010904649138e-5 Å^-1`, absolute difference `1.65e-16 Å^-1`.
- **Localization of discrepancy**: in SM p.6 Eq.(34), after a cosine expression using the full `ω±` bracket and `Re z`, the intermediate next line defines `Q=2Re z` without compensating 1/2, while main Eq.(5) has `Q=Re z`. The direct preceding and following lines cannot literally all agree. Factor-2 *diagnostic* gives plotted magnitudes closer to Fig2(d), but no evidence which factor was used for the published graph. Report/source details: `docs/ROUND3_INDEPENDENT_DISPLACEMENT_AUDIT.md`.

**Before further publication:** Fig2(d) remains unmatched on a **physically validated branch**. Do not override response factor, alter `g`, reinterpret source frequency contrary to SM VI, or mark Fig2(d) PASS. Original Eq(14) common antiunitary sign open as author-intent issue, but clean channel-specific symmetry already independently tested.

**Independent checks:** `src/igmc_lab/independent_checks.py`, `tests/test_round3.py`, `scripts/audit_round3.py`, `parameters_round3.json`, `results/raw/round3_*`, `results/figures/round3_independent_factor_audit.png`. Frozen round1/round2 references preserved unchanged, new round3 frozen entries added separately.

**Next possible work:** Search for original author Fig2 numerical outputs/code; if none exists, develop a concise non-accusatory author question with p6 Eq34 first-line algebra + independent two-route evidence for user review (do not send). Then expand full SM MPG tables and finite-ω metallic response after this gate.

**GitHub state:** no push, PR, or CI run for this paper. NMEE upstream untouched. No author email sent.

**Third-party comparison:** See `docs/GITHUB_PRIOR_ART_AUDIT.md` for WannierNLQG formula definitions (agree with **main** Q/Z) and checked upstream NMEE architecture, with exact blob SHAs. Third-party GPL code **not copied**, and no verified author original Fig2(d) code found in targeted public search.
