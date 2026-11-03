# PR #106 | Independent physics-first review before computational approval
Date: 2026-10-23. **Verdict: NO_GO_AS_WRITTEN / REVISE_READ_ONLY.** No new material computation is authorized or requested here.

The scientific priority is right: measure actual AgF2 B1 insulating spin-Zeeman displacement IGMC Dxx before B2, NEB or further symmetry-only work. However the proposed wavefunction/Wannier job risks constructing the wrong low-energy Hilbert space, and should not be started yet.

## 1. Material occupied space and proposed Wannier count conflict

The plan's OWN JSON records NELECT=100 and **50 occupied SOC spinor bands**, but the proposed 48-Wannier Hamiltonian is interpreted as "24 occupied + 24 empty". That occupation assumption is physically inconsistent. A 48-band Hamiltonian cannot represent all 50 occupied spinor states of the proposed electron-count model. One could intentionally downfold to a frontier-band model, but that is a *partial response model*; omitting deep occupied bands in a Berry/geometric polarization response requires independent justification and a controlled omitted-band test.

The trial orbitals 4 Ag × five d =20; 8 F × three p =24; and 4 Ag × one s =4 add to **48 spatial orbitals**. For ordinary Wannier90 spinor projections containing both up and down by default, the list gives **96 spinor Wannier functions**, not 48. An explicit spin-selected construction would need its own defensible basis. The original Wannier90 spinor documentation and Fe SOC tutorial confirm up/down projection doubling:
https://wannier90.readthedocs.io/en/latest/tutorials/tutorial_17/
https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/projections/

If a full 96-WF manifold is chosen, actual NBANDS must be at least 96 with suitable empty bands; the previous receipt has not pinned actual NBANDS. VASP noncollinear defaults are different and often doubled; inspect the existing private OUTCAR read-only, not an estimate:
https://vasp.at/wiki/NBANDS

**Hard requirement before run**: show actual NELECT, NBANDS and per-k spinor eigenvalues from the previously archived SOC output; explicitly count the bands in the outer and frozen energy windows at EVERY k. Propose either (A) a full occupied-space basis, or (B) an openly labeled low-energy model with a quantitative test that discarded occupied bands are negligible. Do not infer any frozen-band count from only HOMO/LUMO.

## 2. Expected 6-hour cost not justified

The prior SOC result required 4.658 h to reach 1E-6 (64 iterations × ~262 seconds). The new proposal assumes 72–90 iterations at EDIFF=1E-7: extrapolation is 5.24–6.55 h **for SCF alone**, before full spinor Wannier projections/localization and possibly a larger NBANDS/WF count. A 6h hard wall clock is a valid future safety limit but **not a credible promise of complete one-job response delivery**. Do not authorize based on the estimate as written. Rank0 memory ≈907 MB is not total process memory. No discretionary overrun is acceptable; a lost monitor connection must not disable local hard termination.

## 3. Source of real spin operator must be proven before the clock starts

Official VASP tags are LWANNIER90, LWANNIER90_RUN and **LWRITE_SPN** for compatible VASP2WANNIER90-enabled noncollinear builds. The proposed fallback "IWANNIER90" is not a documented VASP option. Merely finding wannier90.x does not prove vasp_ncl supports same-source formatted spn, mmn, amn and spinor projections. Verify the installed binary/build and workflow READ-ONLY in advance:
https://vasp.at/wiki/LWANNIER90
https://vasp.at/wiki/LWRITE_SPN
https://vasp.at/wiki/LWANNIER90_RUN

## 4. Definition of the physical observable and B2 domain

- For the 2D supercell (vacuum c ≈26.94 Å), report vacuum-independent **sheet response** and SI units alongside Dxx/omega, plus a plausible magnetic AC drive and Faraday-electric-background comparison. Ten-times the internal numerical floor means detectable in code, not experimentally useful.
- After reliable Wannier interpolation, test Dxx on **increasingly dense 2D k meshes** and show signed k/band cancellations. 12×8 → 8×6 is not a useful convergence demonstration.
- B1's internal glide {m_x|(1/2,1/2,0)} maps B1 onto **itself**. It does not make an independent opposite FE B2. Source-defined B1/B2 switching uses the layer slide and Type-I **m_z connector**, in the same laboratory axes. Even a symmetry-image B2 does not prove structural bistability or a slide barrier.
- The first decisive physical number is total interpretable real-material spin-Zeeman Dxx (or an explicitly labeled truncated-subspace estimate) with its sign, band/k-resolved contributions, units and uncertainty. An allowed tensor component is not proof of a substantial response.

## Required single read-only revision, before asking for computational approval

1. Archive-source NELECT/NBANDS/eigenvalue census and correct physical spinor projection count (not spatial-orbital count), with an explicit full-occupied versus controlled downfolding decision.
2. Verify availability of vasp_ncl + actual VASP-Wannier90 spinor matrix export on THIS build without running VASP, and establish the authoritative source for H/S/A.
3. Re-budget under that revised NBANDS/WF/window workload with an independent hard watchdog and an honest PARTIAL outcome if the cap expires; do not quietly alter inputs or exceed authority.
4. Define Dxx/omega normalization as a sheet response, mesh/cancellation checks and the distinct B2 magnetic connector; defer B2 calculation and claims of switching.

**Final judgement: REVISE_READ_ONLY; DO NOT AUTHORIZE THE 48-WF/6h JOB AS WRITTEN.** This is not a negative material result. The next scientific objective remains one interpretable B1 Dxx calculation, but only after the subspace and response definition are physically defensible.

No VASP/PBS/Wannier computation was run by this review. The previous archived PR #38 is unchanged.
