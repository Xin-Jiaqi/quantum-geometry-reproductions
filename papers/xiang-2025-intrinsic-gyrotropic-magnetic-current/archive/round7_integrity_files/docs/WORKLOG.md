# Round 5 checkpoint · 2026-10-20

- Read SM pp. 8–10 printed Tables I/II directly (not source inferred from existing SVD output): 11 rows C + 12 rows D, each source group and tensor pattern transcribed in `src/igmc_lab/paper_table_source.py`. 58+58 source-supported MPG names, 244 combinations including each group's forbidden channel.
- Ran independent 9×9 tensor-projector comparison with toolkit v0.18.0 pinned MPG generators (commit `551a0daf8d4ba5a5f24b2811694bfe438616de48`): **243/244** exactly agree. SM Table II p.9 `-42'm'` printed 2 independent entries `D_xy,D_yx`; analytic `S4` + `C2xT` require `D_xy=D_yx`, i.e. **1**, with no change to Table III longitudinal/transverse allowed classification. See `docs/ROUND5_SM_TABLE_SOURCE_AUDIT.md`.
- Preserved explicit SM monoclinic unique-y -> toolkit unique-z coordinate transformation. Recorded anomaly as OPEN_SOURCE_TABLE_CONSTRAINT; did not silently correct the original source transcription.
- Added `scripts/audit_round5.py`, `tests/test_round5_source_tables.py`, `results/raw/round5_source_table_244.csv`, `results/raw/round5_source_table_audit.json`, `results/frozen/round5_source_table_reference_v1.json`. New tests include generator-level analytic tensor check and source matrix subspace equality.
- Existing `results/frozen/*` round1–4 SHA-256 are separately rechecked and unchanged. Source PDF and SM page image crops remain outside publication subtree; no remote push or author contact.
- Next: original author Fig2(d) data/protocol check, exact finite-frequency antiunitary interpretation of printed Eq.(14), and native toolkit API live integration when independently available.

## Previous worklog (preserved)

# Work log · 2026-10-19 · Round 1

## Confirmed work

1. Extracted and inspected the supplied 7-page main PRL and 11-page SM, with figure-image inspection of PRL page 4. Inventoried 8 main figure panels, Eqs. 1–17, SM §§I–VI and Tables I-III; SM has no additional numerical figures. Full row-by-row verification of the large MPG tables remains open.
2. Inspected repository `Xin-Jiaqi/quantum-geometry-reproductions` remote `main` at SHA `78a34431bd4a165484ac3f714786b98120132f81`, including prior paper module docs, scripts, tests, frozen reference and CI path filter. **Read only**; neither main nor any other old file changed.
3. Independently implemented source-defined two-band Zeeman QGT using eigenstate matrix elements with analytic velocity matrices, plus full finite-frequency unbroadened Kubo. Fail-closed degenerate bands.
4. Created Fig1 Rashba (a-d) and Fig2 lattice Dirac (a-d) model diagnostics. For Fig2(d), displayed D/Ω so the unconfirmed Hz-versus-rad/s conversion cannot be hidden. Plotted physical -gμB prefactors are **not yet calibrated**.
5. Verified analytic Rashba Eq16, gap/full-frequency Kubo, Fermi-surface identity Eq13, symmetry tests, mesh convergence; isolated substantive printed Eq14 T parity inconsistency.
6. `python scripts/verify_release.py` gives `PASS_MODEL_OFFLINE`, **11** passed unit tests, 6 frozen model values. This verifies only current subset, **not** full PRL method or published SI amplitudes.

## Exact rerun commands

From the paper folder (`papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`):

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/reproduce_all.py
python scripts/verify_release.py
```

Observed runtime: Python 3.13.5, NumPy 2.3.5, Matplotlib 3.10.8. No network, GPU, author code, real-material input or figure-digitized values used.

## Files / evidence

- `parameters.json`, `src/igmc_lab/core.py`, `tests/test_models.py`.
- `results/raw/fig1_rashba.json`, `results/raw/fig2_dirac.json`, `results/raw/dirac_convergence.csv`.
- `results/figures/fig1_rashba_model.png`, `results/figures/fig2_dirac_model.png`.
- `results/frozen/model_reference.json` checked by `scripts/verify_release.py`.
- `docs/{SOURCE_MAP,THEORY_PROTOCOLS,FIGURE_ACCEPTANCE,INDEPENDENT_THEORY_AUDIT,OPEN_PHYSICS_QUESTIONS,NUMERICAL_CONVERGENCE}.md`.

## Next checkpoint (do not silently mark done)

Address Eq14 time-reversal/antiunitary transformation with exact frequency protocol; independently derive C/D from wavepacket or spectral polarization, redo quantitative SI conversion with `ω` ambiguity recorded, improve μ-local peak integration, then digitize figures with error estimates, and complete SM Tables I-III symmetry generation. Afterwards integrate CI on a fresh GitHub branch, preserving existing NMEE frozen baselines.

## 2026-10-19 · Round 2 (incremental, round-1 baseline untouched)

- Reviewed original main page 3 (Eq.14), main page 4 (Fig.1/Fig.2, model parameters), SM §V Tables I/II, and SM §VI photon-energy statement.
- Added explicit TR-related Dirac partner through `mass_sign=-1` (flip **beta and M**, not just beta). Exact identity `H_TR(-k)=sigma_y H(k)^* sigma_y`; k-local Z even, Q odd (residual <1e-15). Clean Kubo real even / imaginary odd under TR. No finite-band fake damping.
- Added `src/igmc_lab/tensor_symmetry.py` and two named protocols `paper_printed_eq14`, `channel_resolved_symmetry`. Check 4mm1′, 4m′m′, P against source SM: C dimensions 1/1/0; D dimensions 0/2/0. Exhaustive 122 group table remains open.
- Source SM §VI uses ~4 meV for THz: `f=1e12 cycles/s`, `hbar omega=0.004135667694 eV`; counterfactual rad/s gives 0.000658211957 eV.
- Restored SI with `e*g*mu_B/hbar` and Å^-1→m^-1; for `g=2`, factor `281.794032638 A/(m T)` per Å^-1. Rashba Fig1(b) plateau `3.884569769 A/(mT)` vs image ~4; Fig2(b) C at mu=1eV ~-33.49 vs image ~-35.
- New unresolved: paper Fig2(d) visual plateau ~+0.04 A/(mT) vs source-formula result +0.017753893 A/(mT), grid-independent from N=101² to 401² and within 0.013% of exact-ω Kubo. Same roughly 2x mismatch near mu=1eV. No fitted correction applied.
- All **23 unit tests passed** (11 old + 12 new). Executed `scripts/audit_round2.py` and stored SI CSV/JSON and figure. Updated physical evidence and source provenance; old six frozen model baselines not changed. No GitHub branch, commit, CI run or author email performed.

Exact reproducibility commands (from this paper directory, after editable installation):

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/audit_round2.py
python scripts/verify_release.py
```

All generated files from this round are under `results/raw/round2_*`, `results/figures/round2_TR_and_SI.png`, `docs/ROUND2_TIME_REVERSAL_AND_UNITS.md`, `REPORT_ROUND2.md`. Original first-round files remain available separately in `/mnt/data/igmc_reproduction_round1/`; round-one original manifest is archived unchanged at `results/frozen/ROUND1_RELEASE_MANIFEST.json` within this round-two snapshot.

Additional independent direct SI test: restored main Eq.(16) `e*g*mu_B*lambda_R*m/(2π*hbar^3)` in kilogram/joule/metre units, independent of Brillouin-zone `1/Å` integral converter. **23 tests pass.** A separate local `venv --system-site-packages` environment also installs the editable package without network and passes the full verifier; this is a clean interpreter/package installation with host NumPy/Matplotlib, **not** a zero-system-package dependency-resolution proof. Round-one archived manifest SHA256 matches the original archived source exactly.

## Round 3 · 2026-10-19 · independent displacement prefactor

- Retained rounds 1 and 2 original snapshots under `/mnt/data/igmc_reproduction_round1`, `/mnt/data/igmc_reproduction_round2`, plus the two original ZIPs. New work is separate `/mnt/data/igmc_reproduction_round3/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/`.
- Read SM p.6 Eq.(34) directly from the supplied PDF image; found a local `Q=2Re(z)` versus main Eq.(5) `Q=Re(z)` and an intermediate cosine coefficient doubled relative to the preceding first line. Do **not** assert author plot source.
- Added `src/igmc_lab/independent_checks.py`: finite-field occupied-band Wilson Berry phase and direct orbital Liouvillian 4×4 solves.
- At 301² / `b=1e-4 eV`, Wilson `0.01523310376637059` versus main `0.015234078035703156` [Å⁻¹/eV], rel. error `−6.39533e-5`. At 101² / 1 THz, direct Liouville Im `6.301090464897298e-5` versus spectral Kubo `6.301090464913824e-5` [Å⁻¹], difference `1.65e-16`.
- Added five tests in `tests/test_round3.py` and `scripts/audit_round3.py`; frozen source parameters are in `parameters_round3.json` and raw results in `results/raw/round3_*`. Model status `PASS`; Fig2(d) published ordinate `PARTIAL`.
- Verifier command `PYTHONPATH=src python scripts/verify_release.py`. Do not modify older frozen baselines, and do not publish full source PDFs or publisher figure crops.


## Round 4 · 2026-10-19 · full 122 magnetic point group symmetry

- Inspected the user's `Xin-Jiaqi/group-theory-operations-toolkit` GitHub repo on `main`, v0.18.0 commit `551a0daf8d4ba5a5f24b2811694bfe438616de48`. Retrieved actual magnetic point-group registry (blob `406309c...`) and group operations (blob `1c43c66...`); frozen only required generators with source SHA provenance under `data/toolkit_mpg_generators_v018.txt`. Source repo **not modified**.
- Implemented independent 9-dimensional Neumann tensor solver with complete `122 MPG × 2 channels = 244` tensor-space audits, each MPG generator closure checked against registry group order. `C` allowed: 58 (24 longitudinal, 18 transverse, 16 both); `D` allowed: 58 (23, 19, 16), and 64 forbidden for each.
- Checked SM p8 Table I monoclinic `2`/`m` unique-y setting vs toolkit z setting via explicit y↔z permutation, and 4mm1′ / 4m′m′ / 3m1′ / 6mm1′ tensor forms. `paper_printed_eq14` yields different C dimensions for 58 groups under separately declared real-amplitude protocol.
- New files: `src/igmc_lab/toolkit_bridge.py`, `scripts/audit_round4.py`, `tests/test_round4_mpg_toolkit.py`, `data/toolkit_mpg_generators_v018.txt`, `docs/ROUND4_TOOLKIT_SYMMETRY_AUDIT.md`, `results/raw/round4_mpg_122.csv` + `round4_mpg_audit.json` and `round4_mpg_tensor_shapes.json`, `results/frozen/round4_mpg_reference_v2.json`.
- Local `git clone` could not run because DNS/network unavailable; external live `magnetic_tensor_basis` API check explicitly skipped; use `python scripts/audit_round4.py --require-live-toolkit` in pinned installed toolkit environment.
- Rerun: `PYTHONPATH=src python scripts/audit_round4.py` then `PYTHONPATH=src python scripts/verify_release.py`; old frozen model/physical results preserved. Still PARTIAL due to original Fig.2(d) and original printed Eq14 author convention. No remote push, no author email.

- Additional independent complete-group character-average route verifies 244/244 SVD dimensions. Round4 checkpoint split into read-only `round4_mpg_reference_v1.json` and new `round4_mpg_reference_v2.json` to preserve incremental frozen evidence. The canonical 3×3 tensor basis of each MPG/sector is saved in `results/raw/round4_mpg_tensor_shapes.json` and hashed in v2.

## Round 4 v2 · 2026-10-19 · independent character crosscheck

- GitHub native source checked again: toolkit main commit `551a0daf8d4ba5a5f24b2811694bfe438616de48`, registry Git blob `406309c...`, operations Git blob `1c43c66...` unchanged; toolkit is **not installed** in this offline container.
- Added `all_group_operations()` and `character_dimension()`: exhaustive `(R, θ)` closure and 9D representation character `η^θ det(R)(tr R)^2`, independent of Neumann/SVD. All 244 C/D groups pass.
- Added an independent 122-group character audit of `paper_printed_eq14`'s C parity variant. Symmetry protocol difference persists for 58 groups; not a claim that complex-frequency author formula is universally invalid.
- Updated `results/raw/round4_mpg_122.csv` (`C_char_dim`,`D_char_dim`) and original canonical tensor shapes output; new frozen v2 checksum `results/frozen/round4_mpg_reference_v2.json` preserves rounds 1–3 baseline outputs.
- Commands: `PYTHONPATH=src python scripts/audit_round4.py`; `PYTHONPATH=src python -m unittest discover -s tests -v`; `PYTHONPATH=src python scripts/verify_release.py`.
- Status: complete 122 MPG mathematical tensor-dimension map PASS, symbolic SM chart comparison PARTIAL, native toolkit installed API comparison BLOCKED here, Fig 2(d) ABSOLUTE OPEN.

## Round 6 · 2026-10-20

- 新建 `src/igmc_lab/d_vector_independent.py`，独立 Pauli-trace 推导实现，不重用前轮本征矢/几何/谱 Kubo 内核。
- 新建 `scripts/audit_round6.py`（`--no-write` 验证）、`tests/test_round6_dvector.py`；验收数据在 `results/raw/round6_*`，`docs/ROUND6_DVECTOR_FREQUENCY_AUDIT.md` 给出完整公式及独立数值证据。
- 已测 12 个随机构型、12 个完整有限频率复响应。`Gauss128² D/ω=0.015234078015630925`, `midpoint301²=0.015234078035703156`。`μ=0` 主公式 SI `0.01775389319438161 A/(m·T)`；未经授权没有从图形读数拟合系数。
- 推导证实 `K_TR(ω)=K_H(ω)*`，C/D 实/虚分量具有不同 TR 符号；Eq14 原文统一符号仍需作者说明。
- 第 6 轮新的冻结数据在 `results/frozen/round6_dvector_reference_v1.json`；旧冻结文件原样保留，源文件 SHA256 从初始 Round-5 manifest 对照。
- 未解决：Fig2(d) 原作者绘图数值溯源、磁金属费米面严格自适应收敛、Eq14 复杂响应约定、Table II `-42'm'` 参数关系。没有发表图、发送邮件或 GitHub 推送。


## Round 7 · 2026-10-20 · metallic Fermi-line acceptance

- Copied Round6 root unchanged to `/mnt/data/igmc_reproduction_round7/`; kept original and frozen model references intact.
- Added analytic Fermi-contour quadratic solver, y topology split, 2D Gauss occupied-sea response and independent 1D Eq13 delta contour integration in `src/igmc_lab/fermi_contour.py`. Input fixed in `parameters_round7.json`.
- Ran `PYTHONPATH=src python scripts/audit_round7.py` from paper root: C absolute outer-order residual 1.0623e−9 Å−1, D/ω residual 3.4350e−10 Å−1/eV, Eq13 vs Eq11 residual 5.3496e−10 Å−1; µ=±1 eV C=−33.48623534, D=−0.06417844743 A/(m T). Observed two Fermi crossings at fixed ky near 1 eV, so naive one-root integration would be wrong.
- Added `tests/test_round7_metal_fermi.py` with explicit topology, symmetry, independent Eq13 and convergence cases; no code imported from previous spectral Kubo in the metal integrator.
- Retained BLOCKED author original Fig2d graph numeric source; no arbitrary normalization, email, or GitHub upload.
