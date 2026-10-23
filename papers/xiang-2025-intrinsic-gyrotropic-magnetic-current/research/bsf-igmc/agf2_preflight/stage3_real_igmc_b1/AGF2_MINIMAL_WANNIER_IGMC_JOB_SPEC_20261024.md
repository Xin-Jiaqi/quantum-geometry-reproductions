# AGF2_MINIMAL_WANNIER_IGMC_JOB_SPEC — the single most-informative calculation

**2026-10-24 | NOT AUTHORIZED — waits for explicit user approval | B1 only**

## One-line proposal

Re-run the validated B1 SOC single-point **exactly as the screen, but with
LWAVE/LCHARG kept and a Wannier90 spinor projection on the same 12×8×1
mesh**, so that real-material H(R), S(R), A(R) (same gauge) become
available and the first quantitative displacement IGMC \(D_{xx}/\omega\) of
AgF₂ B1 can be computed with the archived PR #38 gate and cross-checked
with its Wilson protocol.

## Why this job (and only this)

Every downstream piece already exists and passed tests (import parsers,
same-gauge algebra, response gate, Wilson winding, validate). The **only**
missing object is the real-material operator chain, which requires exactly
one wavefunction-preserving SOC SCF + one W90 projection. B2, NEB, MAE and
k-convergence series are deliberately deferred: they cost more and answer
nothing until B1's Dxx is known to be nonzero and credible.

## Inputs (byte-reproducible)

| file | source | change vs the archived screen |
|---|---|---|
| POSCAR | `POSCAR.soc_input` (sha256 b4f1aa46…; step-11 geometry) | none |
| KPOINTS | 12×8×1 Gamma | none |
| INCAR | `INCAR.soc` (sha256 1ab493df…) | `LWAVE=.TRUE.`, `LCHARG=.TRUE.`, `EDIFF=1E-7` only |
| POTCAR | hash-only 8c0da262… (content never on GitHub) | none |
| seedname.win | new (synthetic-only provenance; see §Wannier) | new |

### Wannier (spinor, projection-mode)
- `spinors = .true.`, `saxis = 0 0 1` (matches VASP SAXIS),
  `num_wann = 48` (Ag-4d 20 + F-2p 24 + Ag-5s 4, spinor count),
  `dis_win_min/max` from the actual band structure (guess [−7.0, −3.0] eV,
  covering HOMO −5.7746 and LUMO −4.4130), `write_hr = .true.`,
  `write_rmn = .true.`, `write_u = .true.`, `write_v = .true.`,
  `write_spn = .true.`, `write_ndegen_applied = .true.`,
  `use_ws_distance = .false.`.
- W90 interface: preflight-gated — (a) preferred two-step
  `WAVECAR → vasp2wannier90 → wannier90.x`, or (b) `IWANNIER90=.TRUE.`
  if the private vasp_ncl build links the wannier90 library. The spn
  (spin-matrix) export path must be confirmed before launch: **without it
  the "same-gauge physical spin" requirement cannot be met and the job
  should not start.**

## Resource envelope (hard caps)

- 8 MPI ranks, 1 node, OMP=1; wall cap **6 h** (≈ 40–48 rank-hours;
  screen calibration: 262 s/iteration, 907 MB/rank, 64 iterations at 1e-6).
- Disk: WAVECAR ≈ 2–4 GB + CHGCAR ≈ 1–2 GB; preflight free-space check.
- **Hard-stop independent of monitoring**: local watchdog
  `timeout --kill-after=60 21600` + background monitor that aborts the MPI
  job at 5 h 50 m wall if EDIFF is not reached; a lost SSH/chat connection
  cannot extend the run (direct fix for the recorded 4.66 h/4 h overshoot).
- Single continuous run: **no retry, no parameter change, no second
  launch**. NELM=160 is the SCF hard stop inside VASP.

## Success / failure criteria

Success (ALL):
1. SCF converged EDIFF=1e-7 within cap; WAVECAR/CHGCAR saved + hashed;
   TOTEN within ~1 meV of −35.55001581 eV (same physical state); moments
   reproduce ±0.568 μB Type-I.
2. Wannier bands fit the direct eigenvalues ≤ 10 meV in the frozen window;
   Wannier gap ≈ 1.36 eV ± 50 meV (insulator survives the projection).
3. Same-gauge chain passes the archived checks: V†V−I ≤ 1e-8, H = V†diag(ε)V
   = Fourier H(R), spn Hermitian on a k sample, S^W = V†S^B V, r.dat final
   gauge, centers = diag A(0) once, Pauli norm ≤ 1.
4. Dxx/ω from the generalized multiband gate ≥ 10× numerical floor, stable
   under mesh subsampling (12×8×1 → 8×6×1) and Wannier spread tolerance;
   direct vs Wilson slope agree within the documented mesh-limited tolerance
   (synthetic reference 0.41%; expectation ≤ ~5%).
5. Value recorded with the Eq. (14) channel-resolved protocol (η_D = −1,
   archived in PR #38) and explicit unit conversion vs the paper's σ^D.

Failure (ANY → diagnose, do not tune, no retry without fresh authorization):
- no EDIFF convergence within 6 h (SCF tail pathology);
- Wannier fit > 10 meV or gap collapse (window/spread design issue);
- any same-gauge check fails (operator chain not credible);
- Dxx at/under the numerical floor, or direct vs Wilson disagreement beyond
  tolerance (the response is not credible — this is a legitimate scientific
  NO-GO datum, not a bug).

## What the authorized run is FORBIDDEN to do

B2 SCF, NEB, MAE/orientation sweeps, k-mesh convergence series, different
magnetic orders, orbital-L coupling beyond the recorded LORBMOM diagnostic,
or any expansion of the 48-WF window without a new authorization.
