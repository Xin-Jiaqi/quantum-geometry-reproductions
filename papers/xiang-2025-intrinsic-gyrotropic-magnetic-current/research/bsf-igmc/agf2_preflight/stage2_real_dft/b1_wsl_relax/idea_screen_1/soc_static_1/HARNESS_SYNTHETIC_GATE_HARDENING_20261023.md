# DeepSeek Harness: synthetic multiband IGMC test integrity repair (zero DFT)

**2026-10-23. Read-only original research evidence, lightweight local Python patch/testing, no VASP/PBS/Wannier/SCF.**

## Goal

Independently correct and re-run the scientific guards for the existing PR #38 synthetic four-band response test (commit adcd75e). **Do not re-implement the two-band Dirac or previously solved Wilson test.** Keep its existing main positive result unless new tests actually contradict it. This task is **quality hardening**, not fresh material evidence.

Source of truth: https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38. Continue the existing branch `research/agf2-b1-wsl-relax-20261021`. Read `POST_ADCD75E_INDEPENDENT_CODE_AUDIT_20261023.md` and new `scripts/audit_source_gate_invariants.py` before making changes.

### 1. Mandatory scientific fix: NC4 was not spin-conserving

The old `build_no_soc()` used `H=I⊗h2(k)` and `S=I⊗σ`, i.e. both h2 and sigma act on the same factor. For h2(0)=0.6 sigma_z, the physical x spin commutator is **||[H,Sx]||₂=1.2 eV**, independently reproduced (not zero). Thus the old `D≈0` outcome is **not evidence for the stated spin-conserved zero theorem**.

Fix the control consistently as **H=h2(k)⊗I_spin, S_j=I_orb⊗σ_j**. Verify the commutators for all S_j at many k points: max norm <1e-10 eV. Recompute both the direct Dxx/omega and the finite-field Wilson derivative, and insist both vanish within a justified numerical tolerance. Add an explicit unittest/assertion that rejects a wrong tensor-factor arrangement.

The supplemental `audit_source_gate_invariants.py` computes the prior defect and illustrates corrected tensor ordering; compare results against source. **Do not claim this helper was executed by ChatGPT; it was added for your local execution.**

### 2. Mandatory correction: NC2 value-to-label inversion

The code and JSON for **v2_spinZ_pos** say `response_cov=-0.025357180978` and `response_bare=+0.034293747531`. The Markdown says the opposite. Repair the explanatory table, retaining a dated erratum to the original commit. The absolute mismatch 0.0596509285 is correct, but **the sign assignment matters physically**. Distinguish v2_positive and v2_spinZ_pos.

### 3. Mandatory missing consistency gates

- **Indirect/global band gap**: current gap_scan finds min_k (ε₃(k)−ε₂(k)) only. Also calculate min_k ε₃(k)−max_k ε₂(k) for b=0, ±1e-4, ±1e-3, over a suitably converged synthetic k grid. Reject `INSULATOR_PASS` if the *indirect* gap is not positive, even if the direct gap is positive. If invalid, create a separately named synthetic fixture and compare controls transparently; do not silently alter old results.
- **Wilson phase branch**: present fail_closed merely checks |det W| of an **unembedded b=0** Wilson at N=21. That does NOT prove the embedded Wilson phase used for Dxx/omega has no branch jump. Evaluate the **embedded** Wilson eigen/determinant phases at ±b and neighboring transverse ky, check principal-branch jumps and finite-b continuity (or use a documented periodic unwrap/branch method). A Chern/winding case must be explicitly handled or flagged, not silently smoothed away.
- Convert the two-band anchor, b-step convergence, direct/Wilson discrepancy, joint gauge invariance and negative controls into assertions with nonzero exit codes on failure. Avoid purely self-reported success booleans; run a deliberate injected-fault test to show CI rejects wrong physics.
- The ad hoc orbital-Liouville `+0.034304` vs direct `+0.034294` check is **not in the committed script**. Either commit a genuinely independent, reproducible finite-frequency density-matrix solver for this exact fixture, with exact units/quadrature, or explicitly mark it `DEVELOPMENT_ONLY_NOT_ACCEPTANCE_EVIDENCE`. Do not add a fake number or call a spin-Zeeman density-matrix computation an orbital-field IGMC.
- Confirm that the reused Cui pinned interface is content-identical after the provenance header (verified independently except the final trailing newline), and keep upstream SHA pinned. Do not update/merge the independent Cui paper branch.

### 4. Resource and scope limits

Use only Python/NumPy, existing installed dependencies and **≤15 minutes of local ordinary CPU** for this round. No install/recompile/unlimited mesh, no new VASP, qsub/PBS, Wannier run, AGF2 DFT, B2, NEB, or material response integration. Stop on unresolved physics or quota and report PARTIAL. No new third-party network services or leaked private pseudopotentials.

### 5. Deliverables (all directly in same PR #38)

1. Correct the existing synthetic test code and associated gate JSON and Markdown **with explicit dated erratum/provenance**; preserve history and do not rewrite old source DFT records.
2. A new `MULTIBAND_GATE_HARDENING_20261023.md`, detailing NC4 error/fix, NC2 label, global-gap values, phase branch check, Liouville status, test commands and numerical evidence.
3. A machine-readable QA JSON or test output with failures exercised and exact code hash/commit.
4. At least one fast regression test that fails if H and physical S act on the same tensor factor, and one that fails for negative indirect gap or a deliberate Wilson-phase jump.
5. One short PR Conversation comment: PASS/CONDITIONAL/PARTIAL, whether the corrected multiband low-omega Dxx result survives, and the **one** next step toward authenticated material operators.

### Exit criterion

**Only call full synthetic-operator gate PASS** after the corrected spin-conserving zero control, indirect gap, embedded Wilson phase continuity and all failure assertions actually pass. Passing means you may develop a genuine source-authenticated H/S/A file importer; **never** call it material IGMC or permission to run expensive DFT/Wannier. Continue to treat pinned SM Eq.14 one-row disagreement as a separate source-availability caveat.

The user's goal is to find a high-information, genuinely new BSF × IGMC physical story without spending time or resources on repeat calculations. All intermediate information stays in GitHub; no manual agent-to-agent copy-paste.
