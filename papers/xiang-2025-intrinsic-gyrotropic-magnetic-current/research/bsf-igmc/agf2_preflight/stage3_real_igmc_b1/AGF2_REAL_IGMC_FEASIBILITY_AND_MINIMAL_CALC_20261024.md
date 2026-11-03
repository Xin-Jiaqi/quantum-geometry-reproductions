# AgF₂ × BSF × IGMC — next stage: real-material response feasibility and the ONE minimal B1 calculation

**2026-10-24 | Stage 3 (proposal round) | READ-ONLY + PLAN ONLY | PR builds on the archived PR #38 methodology**

Scope of this round: inventory of existing AgF₂ evidence, a light reproducible
Python analysis (glide symmetry, cost model, Wannier manifold estimate), a
**single most-informative calculation proposal** (not authorized yet), and
explicit scientific decision criteria. **No VASP/PBS/Wannier/SCF and no new
AgF₂ material calculation was run** (`no_compute_executed: true`).

---

## 0. The open scientific question

**Does real AgF₂ B1 possess a nonzero, physically meaningful, and
layer-sliding-switchable displacement IGMC \(D_{xx}\)?**

Everything validated so far is *methodology and screening*:

- PR #38 (archived): strict Wannier90 H/S/A import algebra (V = U_dis U,
  S^W = V† S^B V, r.dat final gauge, centers = diag A(0)), the multiband
  Zeeman response gate (direct low-ω Dxx vs Wilson slope, 0.41% agreement at
  41²), and the Wilson-winding fix — proven on **synthetic** four-band /
  six-band fixtures only.
- SOC screen (RUN_RECEIPT.json, PASS_SCREEN): B1 stays a compensated Type-I
  AFM insulator under SOC (|S_Ag| = 0.568 μB along x, net 0.038 μB;
  occupation gap 1.3616 eV), the glide {m_x|(1/2,1/2,0)} is exact to
  1.28e-5 Å (re-verified here), and Ag develops orbital moments |L|≈0.142 μB
  with a common +z component — a diagnostic hint of an orbital channel, not
  an orbital-IGMC value.
- The channel-resolved protocol (η_C = +1, η_D = −1) allows **D_xx under
  m'_x under both conventions** — allowedness ≠ a nonzero value.

**The material value \(D_{xx}(\omega\to 0)\) of real B1 has never been
computed.** That single number — with its stability and sign — is the
gate to the whole "sliding-switchable IGMC in AgF₂" story.

## 1. Step 1 — what the existing SOC calculation can and cannot provide

### Reusable (already on disk, hashes recorded)
| Artifact | Value | Use |
|---|---|---|
| B1 screening geometry (step-11 ionic) | POSCAR.soc_input, glide exact to 1.28e-5 Å (re-verified) | the SCF is a single point on this structure |
| INCAR/KPOINTS/POTCAR hashes | byte-identical to the WSL run | exact reproducibility of the new run |
| Converged SOC electronic record | OUTCAR/vasprun.xml hashes, TOTEN −35.55001581 eV, 64 SCF iterations, EDIFF=1e-6 | density/energy provenance; **not** k-resolved operators |
| Band/gap diagnostics | occupation gap 1.3616 eV, HOMO −5.7746 / LUMO −4.4130, clean 0/1 occupations | frozen-window design for Wannier90 |
| Moments | Ag ±0.568 μB (x), Ag y-canting ±0.019, net mz −0.038, Lz +0.070 | magnetic-state provenance; **not** the Zeeman operator |
| Cost calibration | 8 ranks × 4.66 h = 37.3 rank-h; 262 s/iteration; 907 MB/rank | cost model for the next job |

### NOT reusable for a response (the real gap)
1. **WAVECAR was not saved** (`LWAVE = .FALSE.`): there are no Bloch
   wavefunctions, hence no k-resolved matrix elements of anything.
2. **Local moments are not a Zeeman operator**: MAGMOM/MOMENTS.csv are
   site-projected densities; the response needs
   \(\,S_{mn,j}(\mathbf k)=\langle\psi_{m\mathbf k}|\hat\sigma_j|\psi_{n\mathbf k}\rangle\,\)
   in the same gauge as H(k) — a matrix, not a scalar per site. Using the
   moments instead would silently manufacture the physics we must measure.
3. **No position/Berry-connection matrices** (r.dat final gauge), no
   amn/mmn/spn/eig, no Wannier interpolation, no multiband Dxx.
4. B2 is untouched: no second structural minimum, no switching path, no
   partner SCF.

Conclusion of Step 1: **the operator chain is entirely missing; the only
missing ingredient is one wavefunction-preserving SOC SCF plus a Wannier90
projection.** Everything downstream (parsers, same-gauge algebra, response
gate, Wilson cross-check, decision rules) already exists and passed tests in
PR #38.

## 2. Step 2 — the single most informative minimal calculation

**One job: "B1 SOC Wannier-producing single-point"** (B1 only; no B2, no
NEB, no MAE, no k-convergence series, no different magnetic states).
It is a *strict superset* of the validated SOC screen: identical inputs,
identical magnetic seed and functional, with the wavefunctions and charge
density **kept** and a Wannier90 projection on the same mesh.

### 2.1 Inputs (all already byte-hashed; only tag changes are listed)
- POSCAR: byte-identical `POSCAR.soc_input` (step-11 geometry).
- KPOINTS: 12×8×1 Gamma (identical; sufficient for W90 interpolation on
  the slab and for the response mesh).
- INCAR: identical tag set to `INCAR.soc` **except**:
  - `LWAVE  = .TRUE.`  (write WAVECAR — the wavefunction record)
  - `LCHARG = .TRUE.`  (write CHGCAR for the density record)
  - `EDIFF  = 1E-7`    (one more digit than the screen: Wannier amn/mmn/spn
    are sensitive to the wavefunction tail; 64 iterations → ~72–90
    expected, still far below NELM=160)
  - everything else unchanged: GGA=PE, ENCUT=500, LDAU (Ueff=4 eV on
    Ag-4d), LSORBIT/LNONCOLLINEAR, SAXIS=0 0 1, MAGMOM x-axis Type-I seed,
    ISYM=−1, ISMEAR=0 SIGMA=0.05, NSW=0, IBRION=−1, ISIF=2, LREAL=.FALSE.,
    NCORE=4.
- Wannier90 interface (preflight-gated, see 2.4): preferred two-step
  `WAVECAR → vasp2wannier90 → wannier90.x`; fallback `IWANNIER90=.TRUE.`
  if the private vasp_ncl build has the wannier90 library.
- seedname.win (spinor mode): `spinors = .true.`, `saxis = 0 0 1`,
  `num_wann = 48` (Ag-4d 20 + F-2p 24 + Ag-5s 4; spinor count), frozen
  window ≈ [−7.0, −3.0] eV covering HOMO −5.77 and LUMO −4.41 and the
  first ~1.4 eV of conduction, `dis_win_max ≈ +1 eV`,
  `write_hr/write_rmn/write_u/write_v/write_spn = .true.`,
  `use_ws_distance = .false.`, `write_ndegen_applied = .true.`.
  (Final window values are set from the actual band structure inside the
  run; the 48-WF default gives 24 occupied / 24 empty Wannier states for
  clean interband pairs.)

### 2.2 Expected cost and resource cap (hard)
- SCF: the same cell/mesh/functional as the 4.66 h screen; EDIFF 1e-7 adds
  a small tail → **estimate 5–6 h wall at 8 ranks ≈ 40–48 rank-hours**;
  W90 projection (48 WFs × ~88 spinor bands × 96 k-points) ≈ 0.2–0.5 h.
- **Cap: 6 h wall / 48 rank-hours, 8 MPI ranks, ≤ 1 node.** WAVECAR disk
  ≈ 2–4 GB, CHGCAR ≈ 1–2 GB — checked against free space before launch.
- **Hard-stop independent of monitoring**: the job is launched under a
  local wall-clock watchdog (`timeout --kill-after=60 21600` on the VASP
  process plus a background monitor that writes a STOPCAR-style abort and
  kills the MPI job at 5 h 50 m if convergence is not reached), so a lost
  SSH/chat connection cannot extend the run. This directly addresses the
  4.66 h vs 4 h overshoot recorded for the screen: the cap is enforced by
  the local process, not by the executor's attention. **No retry, no
  parameter change, no second launch** (single continuous authorized run;
  identical to the screen's discipline).

### 2.3 Success / failure criteria (explicit, physics-first)
Success requires **all** of:
1. **SCF converges** to EDIFF=1e-7 within the 6 h cap (expected 72–90
   iterations); WAVECAR/CHGCAR saved; hashes recorded; TOTEN ≈
   −35.55 eV (within ~1 meV of the screen) — the same physical state.
2. **Wannier interpolation fidelity**: the Wannier bands reproduce the
   direct SOC eigenvalues within **10 meV** in the frozen window around the
   gap, and the Wannier gap ≈ 1.36 eV within 50 meV (gap and topology of
   the insulator must survive the projection).
3. **Same-gauge operator chain** passes the archived PR #38 acceptance:
   V†V−I ≤ 1e-8 per k; H(k) = V† diag(eig) V equals the Wannier-Fourier
   H(R); spn Hermitian on a k sample; S^W(k) = V† S^B V; r.dat read as the
   final gauge; **centers = diag A(0) counted once**; dimensionless Pauli
   spin norms ≤ 1.
4. **Response evaluation on the real manifold**: the multiband gate's
   `direct_dxx_omega` generalized from the fixed 4-band indices to the
   actual 24 occupied × 24 empty Wannier subspace; low-ω slope extracted;
   **Wilson cross-check** (embedded per-ky Wilson rows with the A_x
   transporter and spin insertion, `wilson_slope`) agrees with the direct
   value within the documented mesh-limited tolerance (synthetic
   reference: 0.41% at 41²; expectation ≤ ~5% at 12×8×1, sub-sampled).
5. **Nonzero, stable, physical magnitude**: |Dxx/ω| ≥ 10× the numerical
   floor of the same computation; insensitive to (a) Wannier spread
   tolerance, (b) mesh subsampling 12×8×1 → 8×6×1 (cheap post-hoc on the
   Wannier Fourier), (c) the direct-vs-Wilson protocol choice. The
   quantitative unit/scale bar (vs the paper's \(\sigma^D\) convention) is
   fixed in the post-processing with the Eq. (14) protocol note already
   archived in PR #38 (η_D = −1, channel-resolved).

Failure (any one): SCF not converged within cap → **fail, no retry**;
Wannier fit > 10 meV or gap destroyed → fail (window/spread redesign is a
*new* job needing new authorization); same-gauge checks fail → fail (the
operator chain is not credible); response below the floor or direct-vs-
Wilson disagreement → the value is not credible. In every failure the
output is a diagnosis, not a tuned fit.

### 2.4 Preflight gates (no compute, done before requesting authorization)
- Confirm `wannier90.x` (or the vasp_ncl W90 library) exists on the WSL
  machine, which interface (`vasp2wannier90` vs `IWANNIER90`), and that the
  **spn** (spin-matrix) output path is available in the noncollinear build
  — this is the make-or-break risk for "same-gauge physical spin".
- Confirm disk headroom and that no prior dirs are touched.
- Confirm the PR #38 parsers/gate import unchanged (they do not need
  changes for the projection; only the response-index generalization and
  the unit conversion are added — small, tested on the synthetic fixtures
  first).

## 3. Step 3 — scientific decision criteria (pre-committed)

| Outcome of the ONE run | Decision |
|---|---|
| **GO**: Dxx/ω significantly nonzero (≥10× floor), stable (mesh/spread/protocol-insensitive), Wannier-gap intact | **Next stage = B1↔B2 response reversal and switching**: B2 = glide image geometry, one partner SCF at identical cost/cap; test the connector fingerprint Dxx(B2) ≈ −Dxx(B1) vs the AHE-like same-sign failure mode; then (separately authorized) a minimal NEB/sliding barrier assessment. |
| **PIVOT**: Dxx ≈ 0 or cancels/unstable at credible precision | Physical analysis first (spin-channel strength at the gap, orbital-channel contamination, window content), then decide AgF₂ demotion; screen one alternative magnetic sliding ferroelectric with the same one-job template. |
| **BLOCK**: operator chain unreliable (spn unavailable, same-gauge impossible, W90 interface broken) | Report the **single minimal blocker** (e.g., "vasp_ncl cannot export spin matrices; need build with wannier90 lib") and stop. Do not build larger tools. |

## 4. What this round delivers (this PR)

- `AGF2_REAL_IGMC_FEASIBILITY_AND_MINIMAL_CALC_20261024.md` (this file).
- `AGF2_MINIMAL_WANNIER_IGMC_JOB_SPEC_20261024.md` — the one-job spec.
- `scripts/agf2_real_igmc_feasibility.py` +
  `AGF2_REAL_IGMC_FEASIBILITY_EVIDENCE.json` — reproducible read-only
  analysis (glide residual, cost model, manifold/window estimates).

No new material compute was performed or authorized by this round.
