# Stage 8h | Independent raw-output science review of the real Japan PBS WSe2 glide pilot

**Date:** 2026-10-21 · **Mode:** READ-ONLY (zero new PBS jobs) · **Site:** Tohoku IMR `cmspbs1`
**Upstream Git SHA:** `52f126345ee3d7a89d3b6c4ea199a359f35bd094` (inputs frozen for this audit)
**Reviewer code:** `hpc-stage8g/review_phase1.py` SHA256 `9b62b53a0540d1fde241c4b98ac465ed0cfaac56d1fcaeaff7096fa717796532` (with `onsite_pbs.py` `3a6bb12a…` and `material-screening/ideal_stack_symmetry.py` `4f314936…`, byte-exact on HPC)
**Raw evidence:** actual remote `OUTCAR`/`CONTCAR` read directly on the HPC, not pasted summaries.

## 1. Verified raw-output facts (independent of the harness report)

PBS accounting (`qstat -x`): all three jobs `Exit_status = 0`, walltime 17.5 / 12.7 / 6.5 min.
Final SCF per job converged below EDIFF=1e−7 eV (last DAV dE in OSZICAR: s0 −0.88e−7, s0.25 +0.84e−7, s0.5 +0.65e−7 eV). `reached required accuracy` present once in each OUTCAR; final force tables parse cleanly; cell fixed (error 0.0 Å), xy frozen (drift ~1.6e−14 Å), only z relaxed.

| job | PBS ID | final TOTEN (eV) | ΔE (meV/cell) | max \|Fz\| (eV/Å) | ionic steps | net upper-layer Fx,Fy (eV/Å) |
|---|---:|---:|---:|---:|---:|---:|
| path_A_s0.00 | 490621 | −45.664806 | 0 | 0.006172 | 16 | (0, 0) |
| path_A_s0.25 | 490622 | −45.651172 | +13.633 | 0.006209 | 13 | (−0.0316, 0.0183) |
| path_A_s0.50 | 490623 | −45.640689 | +24.117 | 0.009974 | 6 | (0, 0) exactly |

The harness-reported 13.634 / 24.117 meV/cell and the second finite difference −3.150 meV/cell are reproduced from raw OUTCAR TOTEN (sanitized evidence: `STAGE8H_AUDIT_RAW_20261021.json`).

## 2. Group-theory gate on the ACTUAL relaxed six-atom CONTCARs (tol 0.001 Å, robust to 1e−5 Å)

Whole-molecule enumeration (2D lattice modulo translation, origin-independent z) of the six final atoms:

| job | surviving ops | point group / space motif | total tensor dim | layer-staggered dim |
|---|---|---|---|---|
| s=0.00 | 6 (E, 2×C3, 3×m_v; all zflip=+1) | **C3v** | 1 | 1 |
| s=0.25 | 2 (E, m_v) | **C_s** | 2 | 2 |
| s=0.50 | 4 (E, **glide {m_z\|t}**, 2× vertical C2′/screw with zflip) | **order-4 nonsymmorphic, glide present** | **0** | **2** |

**The highest-priority test PASSED: the DFT-relaxed s=0.50 CONTCAR retains the layer-exchange horizontal glide mirror.** The z-only ISYM=0 relaxations preserved the full atomic-motif symmetry to better than 1e−5 Å at every point (s0: C3v, s0.25: C_s, s0.5: glide order 4, verified by tolerance sweep 1e−3 → 1e−5 Å).

## 3. Sliding-stationarity diagnostic (net upper-layer in-plane force)

Raw final force tables, summed over the upper layer (W2 + Se5 + Se6):

- s=0.00: net (Fx,Fy) = (0, 0) — AB-like stacking stationary.
- s=0.25: net (−0.0316, +0.0183) eV/Å — off-stationary; gradient points back toward the AB minimum (consistent with E(s=0.25) > E(0)).
- s=0.50: net (0, 0) **exactly** — the glide midpoint is a stationary point for rigid in-plane sliding, as required by the retained glide. Individual frozen-atom Fx/Fy are nonzero and cancel layer-wise; the force pattern is exactly mirror-symmetric between layers (Fx,Fy equal, Fz opposite).

## 4. Symmetry selection rules at the validated geometries (ideal vs actual)

- s=0.00 (C3v): uniform static (equilibrium) NMEE tensor restricted to one scalar component.
- s=0.25 (C_s): 2D restricted (nematic-type) total static response.
- s=0.50 (glide order 4): **uniform total-spin static NMEE FORBIDDEN (total tensor dim 0); layer-staggered static response ALLOWED with dim 2.** This is a symmetry selection rule on the geometry actually relaxed by DFT; it does not establish a real WSe2 NMEE magnitude.

## 5. Numerical issues found and minimal next compute

1. **s=0.50 max|Fz| = 0.00997 eV/Å ≈ EDIFFG = −0.01 (marginal).** The W atoms carry ±0.00997 eV/Å. Tighten EDIFFG to −0.005 (or −0.003) and re-relax ONLY s=0.50, or continue the relaxation from the current CONTCAR; verify ΔE stability at meV level.
2. **k-mesh 15×15×1 is initial, not production-converged.** Repeat ONLY s=0 and s=0.5 at 21×21×1 (identical PBE-D3/POTCAR/ENCUT 520/dipole protocol) to test the +24.1 meV/cell rise. 
3. **Discrete second difference −3.15 meV/cell is not a barrier curvature** and cannot locate the continuous-path saddle; 3 constrained points ≠ CI-NEB.
4. **Mirror-partner energies not independently run** (s=0.75, s=1.0 on arm A). With the glide/C3v symmetry now verified on the actual relaxed geometry, E(s)=E(1−s) is symmetry-guaranteed along arm A, but a numerical check of one mirror pair is cheap and recommended.

Proposed minimal next wave (requires explicit approval, NOT authorized by this read-only phase): 2–4 jobs = {s=0, s=0.5 @ 21×21×1} ∪ {s=0.5 re-relax EDIFFG −0.005} ∪ optionally {s=0.75} for the mirror check.

## 6. Explicit non-claims

No CI-NEB switching barrier, no ferroelectric bistability, no SOC full-BZ insulating gap, no real WSe2 static NMEE value, no equivalence to the driven-DC Edelstein result. `NEB_barrier = null`, `material_static_NMEE = null` in the evidence JSON. Remaining 12 generated inputs were NOT submitted.

## 7. Security / provenance

Original OUTCAR/CONTCAR stay on the HPC; only sanitized hashes, energies and derived symmetry quantities are versioned here. No POTCAR, credentials or private site settings in this repository.
