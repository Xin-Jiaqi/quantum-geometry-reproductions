# AgF₂ Stage-2 · Readiness checklist for the first qsub batch (advisory only)

Status: **NOT_RUN / no jobs submitted**. This checklist only *suggests*
preconditions, costs and failure criteria; nothing here authorizes a
submission.

## 1. Preconditions (all must be true before any qsub)
- [ ] User explicitly approves: compute-budget (core-hours/walltime), queue
      (`CP_001` candidate), PBS shape (32 cores / 120–128 gb / ≤24 h),
      VASP cpu build (6.6.1_cpu candidate) and PAW set
      (private PBE.64 Ag 02Apr2005 + F 08Apr2002) with `potcar_lock_sha256`
      `8c0da262849433cad092946af66c1b584d457e46a70d4c1eb768ea68305f0fff`.
- [ ] POTCAR assembled privately on the cluster (element order Ag → F);
      its sha256 matches the recorded combined hash; **never uploaded**.
- [ ] POSCAR/INCAR/KPOINTS hashes match the protocol draft (§2 of
      INPUT_PROTOCOL_DRAFT.md); species order `Ag F` matches POTCAR order.
- [ ] Site oneAPI/MPI loader path and `mpirun -np 32 -ppn 32 -hostfile
      $PBS_NODEFILE` invocation verified with the user (site-proven pattern).
- [ ] IVDW D3 variant (zero vs BJ) confirmed against the chosen VASP build.
- [ ] ISMEAR/SIGMA/ISYM/LDIPOL/NCORE finalized and entered in the
      convergence record (they are operator choices, not paper parameters).

## 2. Suggested minimal first test (before the three relaxations)
1. **Single point** (NSW=0, same settings) on B1_g3.0 (identical to the
   Stage-1 frozen B1 geometry): verify SCF stability, MAGMOM seed does not
   collapse, POTCAR order correct, `EDIFF=1E-6` reached.
2. If SCF OK: run the three ion relaxations (ISIF=2) B1_g2.6 / g3.0 / g3.4
   with identical settings; walltime 24 h each.
3. Then 14×9×1 / 16×10×1 convergence checks on any converged minimum.

## 3. Estimated cost (rough planning figure, not a quote)
- 12 atoms, 500 eV, 12×8×1, spin-polarized non-SOC; ~2–5 h walltime per
  relaxation on 32 cores → order 64–160 core-hours per case; three cases +
  single-point + convergence checks ≈ **a few hundred core-hours** in total.
  Queue occupancy at check time was high (CP_001 ~470 jobs), so walltime
  estimates are pessimistic; re-check queue before submission.

## 4. Failure criteria (stop and review)
- SCF fails to converge (NELM exceeded, `EDIFF` not reached).
- Magnetic order collapses (Ag local moments → ~0 or ferromagnetic
  alignment) — stop; re-examine MAGMOM seed and ISYM.
- Unphysical short interlayer contact (min F–F distance below ~2.2 Å) or
  ionic coordinates drifting out of cell.
- Relaxation slides into another stacking or shows a high-symmetry
  minimum — record; do **not** interpret mirror-equivalent energies as a
  ferroelectric double well (needs a verified energy double well + polarization).
- Forces not converged within NSW for any case.
- `OUTCAR/vasprun.xml` hashes stored privately; raw large files never
  uploaded to the public repo.

## 5. After the pilot (future, requires separate approval)
14×9×1 / 16×10×1 convergence, B2 partner construction via unitary m_z with
axial sign-flipped in-plane moments, sliding-coordinate scan, then G3 items
(SOC, Type-I vs Type-II, Berry polarization, NEB, Wannier, IGMC) only with
separate scientific and resource clearance.
