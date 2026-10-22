> **SUPERSEDED AS A PLAN (2026-10-22):** User approved the read-only science review. The executable instruction is now [BD-MNPS3-L2-R2P0-READONLY.md](BD-MNPS3-L2-R2P0-READONLY.md), indexed from [CURRENT.md](CURRENT.md). **Zero compute submission authority.** This older draft is preserved as history.

# DRAFT R2-P0 | MnPS3 L2 site-resolved Néel operator and SCF failure audit

**Prepared:** 2026-10-22. **Status: READ-ONLY SCIENTIFIC PLAN / NO DFT SUBMISSION PERMITTED.**  
**Target user workflow:** local DeepSeek harness reads this GitHub file and uses already authorized local private evidence. It must not treat a published GitHub file as permission to spend compute allocation or restore a historical R1 job. \`tasks/CURRENT.md\` must continue showing **R1 CLOSED / no active compute** until a distinct new task is deliberately activated.

## Why R2-P0 exists

Our [Gate6 exact 32-band toy decomposition](../docs/BD_GATE6_EXCHANGE_SOC_MAGNETIC_TEXTURE_AUDIT_20261022.md) shows most layer-transfer cancellation goes into **electronic exchange torque**, but that toy Hamiltonian contains *layer-A-type* exchange \(m\tau_z\sigma_z\), whereas the actual source MnPS3 has a within-layer two-Mn Néel pattern. Thus the toy layer-odd spin response is **not** an accepted MnPS3 magnetic-order response, and a correct four-Mn projector is the next scientific necessity. R1 produced a **source-based L2 moment map and accepted geometric symmetry** yet no converged DFT, so this phase is read-only and must not imply physical response results.

## R2-P0 action for an existing local harness, no user log copy/paste

1. **Read immutable R1 evidence locally** from \`BD-MNPS3-L2-P1-R1\`, including private preflight manifests and source structures. The remote private GitHub archive only exposes a sanitized \`result.json\`/\`results_summary.md\`; do not assume a file mentioned there is downloadable through the public repo. Source article is Sun et al., *Nature Communications* **17**, 3101 (2026), Eq.13 \`[E || M_z M_y]\` for **L2 path 1**, not a toy magnetic connection.
2. Extract **every physical Mn site index, layer T/B, honeycomb A/B sublattice, position, collinear sign \(\eta_{\ell a}\) and local moment reference** for \(\mathbf P_1^+(\mathbf L_2)\) and \(\mathbf P_1^-(\mathbf L_2)\). Verify true site permutation under the source Eq.13 with any required spin basis change, preserving ordering and case-sensitive site metadata; do not infer the order from merely \`MAGMOM=-1,1,1,-1\` unless atom numbering is independently checked.
3. Design projector/vertex **definitions**, not computed matrices:
   \[
   S_{\rm layer-odd}^{\alpha}=\sum_{\ell a} l_\ell P_{\ell a}s^\alpha P_{\ell a},\qquad
   S_{\rm N\acute eel}^{\alpha}=\sum_{\ell a}\eta_{\ell a}P_{\ell a}s^\alpha P_{\ell a}.
   \]
   Identify real \(L_2\) normal modes and define the electronic generalized force \(F_{\theta_\nu}=-\partial H_{\rm ex}/\partial\theta_\nu\). Distinguish \(i[H_\perp,S_{\rm layer-odd}]\) from the torque exerted **on the local Mn moments**, and do not equate either with a bulk spin Hall current without a transport geometry.
4. Independently derive the **SOC magnetic point symmetry** and which Cartesian \(E_\beta\to F_{\theta_\nu}\) or \(E_\beta\to S_{\ell}^\alpha\) **linear** tensor components are allowed or forbidden in the two source-correct states. The paper Eq.13 operates in a **nonrelativistic spin-space group**; do not apply its pseudo-spin-label reversal unmodified to a SOC Kubo tensor. Identify which projector/orbital embedding terms a would-be Wannier implementation needs.
5. Write **one private, source-hashed result** \`bd_mnps3_r2p0_operator_contract.json\` plus a short human-readable report with \(L_2\) site tables and magnetic symmetry verdict; avoid private paths, raw coordinates and unpublished matrices in public output. Classify \`OPERATOR_SPEC_READY\` or \`BLOCKED_SITE_MAP\`.
6. Separately review R1 *existing* VASP logs **read-only**, compare the two frozen INCARs semantically, confirm which fields differ intentionally, and propose **at most two genuinely controlled SCF stabilization pilots** with one intervention at a time (mixing **or** algorithm/damping, not all simultaneously), preserving the published functional/U/D3/cutoff/k mesh and the two magnetic structures. Do not run or schedule them; include estimated walltime/memory and a no-go if the prior 2h cap is scientifically insufficient. Mark \`SCF_RESTART_REQUIRES_NEW_USER_APPROVAL\`.
7. Produce a **sanitized, minimal pointer** suitable for updating the public B+D research tree with the operator definition and whether the site-map gate passed, after source/disclosure review. Detailed input coordinates, exact hash evidence, raw logs, POTCAR and full operator matrices stay private.

## Stop conditions

- The R1 source site map or original MnPS3 structure mapping is unavailable or ambiguous: \`BLOCKED_SITE_MAP\` — **no guessed site index**.
- Geometry passes but SOC magnetic point group cannot be verified: \`BLOCKED_MAGNETIC_SYMMETRY\`.
- Wannier \`hr.dat\` alone is available but no spin/position/site operators: \`BLOCKED_RESPONSE_VERTICES\`; no fake local spin response.
- R1 VASP lacks numerical convergence, and approval for new compute is absent: **NEVER SUBMIT**. This file and R1 closure give no implicit approval.

## Output summary

\`\`\`text
Task: BD-MNPS3-L2-R2P0-SITE-OPERATORS
Mission: READ_ONLY
R1: CLOSED; NUMERIC_CONVERGED=false
L2 +/− site map: READY / BLOCKED + source hashes
Magnetic connection Eq13 under real SOC: VERIFIED / PROVISIONAL / BLOCKED
Actual physical Néel vs layer projection: defined with correct Mn labels / BLOCKED
Allowable frequency-resolved response component: <source-derived> / UNDETERMINED
SCF diagnostic: one-variable-at-a-time pilot plan / BLOCKED
New jobs authorized or submitted: NONE
Next decision: choose if separate bounded R2 COMPUTE task is scientifically justified
\`\`\`

**Novelty note:** Nothing about the exact electronic Ward identity, the exchange-vs-SOC partition, the source's pseudo-spin connection, or a van der Waals magnetic field-induced torque is automatically new. Compare [2025 vdW layer SOT](https://doi.org/10.1038/s41467-025-63966-2), [2026 THz MnPS3 symmetry modulation](https://doi.org/10.1038/s41467-026-78519-4), [2026 octupole-driven AFM torque](https://doi.org/10.1038/s41586-026-11172-5), [2025 NSOT](https://doi.org/10.1038/s41467-025-63171-1), and the MnPS3 FE-spin literature before proposing a paper claim.
