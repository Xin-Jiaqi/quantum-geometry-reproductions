# AgF2 BSF × IGMC: Stage 2 geometry and minimum DFT pilot
**Status: GEOMETRY-PREPARED / DFT-NOT-RUN.** Separate follow-up research; frozen Xiang PRL reproduction source, numbers and tests are untouched.

## Physics first
The issue is not whether two symmetry-related pictures can be drawn. We need to establish that C2DB-derived AgF2 bilayers possess **locally stable sliding structures** under an independently converged magnetic DFT protocol. Two exactly symmetry-related structures can still relax to a single high-symmetry minimum; symmetry-related electronic responses alone are NOT proof of a ferroelectric double well.

For Type-I, the magnetic B1/B2 connection is unitary {m_001 | (0.5,0.16)}. Physical spin is AXIAL: S' = det(R) R S, so m_z reverses the mapped in-plane S_x. The non-SOC collinear input corresponds to opposite interlayer Neel vectors, not a claim about an SOC easy axis.

## What is included
- Six **unrelaxed, constructed** POSCARs: B1 and B2, each with nominal vertical F-F interface-plane gaps of **2.6, 3.0, 3.4 Angstroms**. The original 3.0 Angstrom is a Stage-1 assumption; none of these gaps is a fitted or literature equilibrium spacing.
- The same C2DB-derived lateral cell, 20 Angstrom vacuum and source provenance. Check the independent B1/B2 mirror partners and 12-atom Ag4F8 stoichiometry with `python stage2_geometry/verify_candidates.py` from the `agf2_preflight` directory. No VASP inputs or calculation outputs.
- `CALCULATION_GATE.md` and `EXECUTOR_BRIEF.md` give the pilot and handoff contract. All results remain **NOT_RUN** until physical evidence is supplied.

## Minimum-cost physical campaign
1. **Gate 0, before any job:** licensed VASP, PAW Ag/F identities and checksum kept private, machine/account and approved compute budget. Freeze PBE + Ag 4d Dudarev U_eff=4 eV + D3 variant, cutoff baseline 500 eV, vacuum 20 Angstrom, electronic tolerance 1e-6 eV and force threshold 0.01 eV/Angstrom from Zhu et al.'s SM. kmesh, SIGMA, MAGMOM, dipole, optimizer and D3 variant must be explicitly chosen, not misrepresented as author inputs.
2. **Pilot:** three B1 initial gaps at 0% strain, Type-I magnetic seed, **non-SOC collinear magnetic relaxation** with in-plane cell fixed and relative sliding **not artificially constrained**. Example Ag moments in POSCAR order: [-m,+m,+m,-m] with m a recorded starting guess, F zero. Non-SOC spin direction is not a magnetocrystalline anisotropy result. Use inexpensive mesh 12x8x1 for screening and check convergence using 14x9x1 / 16x10x1 before scientific comparisons.
3. **Only if plausible:** transform final B1 atomic positions and AXIAL magnetic moments by the connector to form a B2 partner. Independently compare energies and transformed forces at equal settings; then test release of sliding coordinate. Numerical equal energies verify symmetry but not two distinct minima.
4. **Do NOT jump to** SOC anisotropy, competing Type-II magnetic energy (the reported 0.26 meV/f.u. is very small), Berry polarizations, NEB, Wannier or IGMC without separate evidence and approval.

## Run the deterministic preparation check
```bash
cd papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight
python stage2_geometry/verify_candidates.py
```
This is a geometry/provenance test, **not a DFT test**. Neither script nor CI has access to VASP/POTCAR/Slurm.

## Scientific status
`source=CONSTRUCTED_C2DB`, `physics_approved_by_user=false`, `data_origin=NOT_RUN`, `potcar_lock_sha256=null`. No job has been submitted. Upload DFT evidence in a *new* `stage2_real_dft` research directory, preserving these frozen initial conditions.
