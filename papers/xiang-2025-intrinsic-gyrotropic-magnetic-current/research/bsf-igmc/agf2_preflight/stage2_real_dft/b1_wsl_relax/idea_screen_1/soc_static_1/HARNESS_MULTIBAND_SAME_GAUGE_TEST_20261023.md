# DeepSeek Harness: bounded four-band same-gauge IGMC falsification

**2026-10-23. Task: synthetic Python/NumPy validation only. No VASP, PBS, Wannier, material SCF or material IGMC.**

## Goal

Use an insulating four-band SOC spinor toy Hamiltonian with nontrivial physical-spin S and position/Berry-connection A to test the low-frequency spin-Zeeman displacement response with **two independent techniques**:

- direct occupied–empty interband Zeeman-geometric Dxx/omega;
- rank-2 occupied **non-Abelian finite-field Wilson polarization derivative** for H-b Sx with the same position embedding.

This test must also fail deliberately when spin operators are in a different gauge or the position connection is omitted.

## Avoid doing previous work again

- Existing original Xiang reproduction already has a two-band Dirac Eq.(12)/finite-field Wilson and Liouville check in papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/scripts/audit_round3.py and frozen results/raw/round3_wilson_convergence.csv.
- Existing BSF eight-band hopping toy has a four-occupied-state non-Abelian Wilson check, documented in research/bsf-igmc/docs/MODEL_AND_VALIDATION.md and docs/SCIENTIFIC_ASSESSMENT.md.
- Read the new independent review POST_9097F42_INDEPENDENT_REVIEW_AND_MULTIBAND_GATE_20261023.md in this same PR before implementing.

## Reuse existing cross-project synthetic operator infrastructure

The **same repository at pinned commit cdac6418166466b5243e73d8eb291d6dea8d4772** already contains:

- papers/cui-2026-curvature-flux-zeeman-qgt/ROUND26_SAME_GAUGE_WANNIER_OPERATORS.md
- papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/wannier_operator_contract.py
- papers/cui-2026-curvature-flux-zeeman-qgt/scripts/audit_wannier_operator_contract.py
- associated synthetic-bundle tests

That work has an **in-memory** common-basis H(R), physical Pauli S(R), full A(R) contract, analytic/covariant derivatives, and positive/negative gauge controls. It has **NO** real Wannier file import or AgF2 material operators. This component may not be present on PR #38's branch. Pin its SHA and use as a reference; do not merge/rebase unrelated branches or overwrite existing files.

## Scientific and execution requirements

1. Check the chosen four-band fixture has a global insulating gap at two occupied bands over all sampled k and at both small b signs. If source fixture is gapless, make a separately named, explicit synthetic modification; do not alter or claim real material data.
2. Use covariant DkH = partial_k H - i[A,H] and physical spin, with consistent real-space Bloch/Wannier origin, BZ measure, and Zeeman energy b sign. In the finite-field Wilson-loop computation **include A/embedding in parallel transport**, not just neighboring eigenvector overlaps.
3. Compare direct multiband Dxx/omega to Wilson slope, with separate mesh and b-step convergence. Use small meshes 21², 31², 41² where feasible. No fitted scalar factor, no silent e/hbar or 2pi changes.
4. Jointly rotate H/S/A by the same constant unitary: responses must stay invariant. Intentionally rotate H/A but not S; omit A; substitute orbital pseudospin for physical spin. These **negative controls must produce a detectable incorrect result**. Include spin-conserving, SOC-free zero-displacement control.
5. Fail closed on band crossings, ambiguous finite-field Wilson branches, inconsistent Wannier centers, or a non-gapped half-filling state. Distinguish numerical covariance from evidence of a material response.
6. Execute only light Python/NumPy work, with total ordinary CPU time target under 15 minutes. If slower or blocked, write PARTIAL and stop. **Do not start any VASP, qsub, PBS, Wannier, B2, SOC or other new material job.**
7. Keep source paper's unresolved '-42 prime m prime' supplemental row separate; do not label it resolved without the actual primary SM.

## GitHub output contract

Push to the EXISTING PR #38 branch research/agf2-b1-wsl-relax-20261021:

- idea_screen_1/soc_static_1/MULTIBAND_SAME_GAUGE_IGMC_GATE.md
- idea_screen_1/soc_static_1/MULTIBAND_SAME_GAUGE_IGMC_GATE.json
- idea_screen_1/soc_static_1/scripts/test_multiband_igmc_gate.py, plus minimal synthetic fixtures if necessary.

Report source SHA, equation and sign conventions, gap, b and N grids, direct and Wilson values and relative discrepancy, gauge invariance residual, negative-control errors, any unresolved units, and PASS/PARTIAL/FAIL. Test executable must print reproducible output; claims of PASS require actual successful execution.

Update **the same PR #38 Conversation** with an evidence-based result and the **one** recommended next step. Do not ask the researcher to copy outputs between DeepSeek and ChatGPT.

## Go criterion

GO for **developing an authenticated material operator interface**, not GO for costly DFT/Wannier, requires a finite-field/direct multiband match plus joint-gauge invariance and failed negative controls. PARTIAL if gauge/units/source operator remains uncertain; FAIL if methods fundamentally disagree. No outcome is by itself proof of nonzero real AgF2 IGMC.

Prior 4-hour VASP authorization has expired and had a documented overrun. No new run is authorized.
