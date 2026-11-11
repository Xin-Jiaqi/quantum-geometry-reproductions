# AgF₂ B1 partial relaxation: independent coordinate audit (2026-10-22)

**Status: READ-ONLY VERIFIED AGAINST GITHUB POSCAR COORDINATES; NOT FORCE CONVERGENCE; NO NEW VASP EXECUTION.**

## 1. Evidence boundary

The coordinates were independently parsed from two publicly committed files:

- Initial main: agf2_preflight/stage2_geometry/candidates/B1_g3.0/POSCAR
- PR #38: agf2_preflight/stage2_real_dft/b1_wsl_relax/results/POSCAR.relaxed (ionic-step-15 snapshot, explicitly unconverged).

The accompanying pure-stdlib [coordinate audit script](audit_agf2_b1_positions.py) implements x/y minimum-image displacement and distance without wrapping through z vacuum. Fixed atom order from the C2DB-derived B1 construction: top Ag ions 1–2; bottom Ag ions 3–4; top F ions 5–8; bottom F ions 9–12; inner F top 5/8; inner F bottom 9/12. These are 1-based indices.

Independent SHA256 checks of the exact POSCAR bytes: initial = 2d323c3fe2353ef984330b3b99fc24a75a85337883bbaf36d08b20ed45e6d60a; step 15 = d49c82807fc2309575e5d3757904ab78caab348aef10e895ac2327a904ac2d68. Both match the recorded file hashes.

**This verifies internal geometrical consistency of the uploaded POSCARs only.** It does not independently verify private OUTCAR ionic-step provenance, forces or total energy. Those remain reported by the WSL runner.

## 2. Independent coordinate-derived quantities

| Quantity | Initial | Step 15 | Difference |
| --- | ---: | ---: | ---: |
| Top Ag z (Å) | 15.953906 | 15.490255 | −0.463651 |
| Bottom Ag z (Å) | 10.984635 | 11.445990 | +0.461355 |
| Ag-plane separation (Å) | 4.969271 | 4.044265 | −0.925006 |
| Top inner F z (Å) | 14.969271 | 14.582070 | −0.387201 |
| Bottom inner F z (Å) | 11.969271 | 12.370195 | +0.400924 |
| **Vertical inner F-plane gap (Å)** | **3.000000** | **2.211875** | **−0.788125** |
| Nearest interlayer F–F 3D distance (Å) | 3.680562 | 3.052256 | −0.628305 |
| Top Ag mean Δy (Å) | — | −0.02635 | — |
| Bottom Ag mean Δy (Å) | — | −0.03150 | — |
| **Ag-layer relative Δy (Å)** | — | **+0.00515** | **+0.000876 of b** |

### Important correction A: vertical F-plane gap

The historical B1_RELAX_RUN_RECEIPT.json records final_geometry.vertical_F_plane_gap_A = **2.687**. The two inner F z planes in the *actual committed step-15 POSCAR* instead give:

\[
g_{F,z}=z_\mathrm{top,innerF}-z_\mathrm{bottom,innerF}
=14.582070-12.370195=\boxed{2.211875\ \mathrm{\AA}}.
\]

The **3D nearest interlayer F–F distance** is a distinct quantity. The top ion 5 / bottom ion 9 (and mirrored ion 8 / 12) have minimum-image separation approximately (1.876575, 0.949930, 2.211870) Å, whose 3D norm is **3.052256 Å**, in agreement with that separate receipt field. The old value 2.687 must be marked as **inconsistent with this explicit plane definition**, pending the original source's explanation or correction. Do not quietly replace a historical numerical receipt.

### Important correction B: lateral slide

The old report uses −0.3400 → −0.3391 in fractional y but associates this with ~0.03 Å of sliding. For b = 5.879665746 Å, the independently computed Ag-layer *relative* displacement is 0.00514848 Å = 0.00087564 b. The individual layers move by −0.02635 Å and −0.03150 Å in y, which nearly cancel. Thus ~0.03 Å approximates each layer's *absolute* translation, **not** the layer-relative sliding. This is an Ag-sublattice registry proxy; F internal distortions mean there is not one unique shift parameter for a distorted bilayer.

### Important correction C: convergence vs geometric compression

The step-15 maximal force is still source-reported as **0.497 eV/Å**, exceeding the EDIFFG=−0.02 eV/Å convergence threshold by ~24.9 times. A geometric energy descent plus a substantial z compression **does not establish a local minimum or equilibrium interface distance**. No Hessian eigenmode has been computed, so “CG soft-mode oscillation” is an interpretation requiring the raw line-search record rather than an established physical instability.

## 3. Physical interpretation from these coordinates

- The Ag planes approach by ~0.925 Å. The inner F planes approach by ~0.788 Å, so the two layers also undergo internal distortions rather than translating as rigid slabs.
- Along this particular ionic path the two Ag layers have only ~0.005 Å of relative y drift. This is **not** an independent verification of positive sliding curvature or a finite energy barrier: a high-symmetry starting path may exhibit small net y force even at a saddle point.
- The step-14 reported nearest F–F 2.644 Å and subsequent rebound may reflect line-search trial/corrector behavior, a steep structural energy landscape, or both. These possibilities cannot be separated with the public step summary alone.
- The magnetic order and preliminary gap observations in the existing PR are *runner-reported electronic results*, not recomputed by this geometry-only audit.

## 4. Next DeepSeek action: strictly read-only

1. Read the already-existing **private WSL OUTCAR/vasprun.xml/OSZICAR**, align CG trial/corrector events (especially trialstep, gam, g(F), ort) to the ionic forces, energies and coordinate sets at steps 6/13/14/15.
2. Explain where 2.687 Å came from; if a *different physical definition* was used, identify its indices and definition. Otherwise add a dated correction in a new audit without suppressing provenance.
3. Determine which force-bearing atoms and Cartesian components drove the rise in max|F| and trial-step overshoot. Without a dynamical matrix, do not label a true unstable eigenmode.
4. Verify the step-15 snapshot against the **last completed** VASP force block and check STOPCAR/CONTCAR behavior before selecting *any* new restart geometry.
5. Push sanitized read-only output back to **this PR #38 branch**. Do **not** run VASP, submit PBS, change parameters, or begin B2/NEB/SOC/Wannier/IGMC.

See the existing [postrun scientific review and detailed task](B1_POSTRUN_INDEPENDENT_REVIEW_AND_READONLY_TASK_20261022.md). After evidence is available, separately decide whether any single continuation is justified and where its initial geometry should come from.

## 5. Reproduce without VASP

Run the following script with exactly two existing paths:

python audit_agf2_b1_positions.py [original B1_g3.0 POSCAR] [step-15 POSCAR.relaxed]

The script prints JSON and does not create or modify files. It uses only Python stdlib and the known Ag4F8 orthorhombic atom ordering. It does not work unchanged for skewed cells or arbitrary bilayers.

**Method reference:** [VASP Structure optimization](https://vasp.at/wiki/Structure_optimization), [POTIM](https://vasp.at/wiki/POTIM). The documented conjugate-gradient trial and corrector steps must not be treated as separate potential-energy minima.
