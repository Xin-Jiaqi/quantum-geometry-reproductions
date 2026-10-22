# BD-MNPS3-L2-R2P0 | ACTIVE READ-ONLY SCIENCE AUDIT

**Approved in conversation:** 2026-10-22. **Mission:** 0-VASP, 0-Slurm, 0-WSL/IMR/para jobs, 0 new compute allocation, **private evidence review only**.  
**Actual runtime status:** \`APPROVED_FOR_LOCAL_HARNESS_INTAKE / NOT_CONFIRMED_STARTED\`. Creating/updating this GitHub task does not remotely trigger the DeepSeek/Codex harness.  
**Historical R1:** \`BD-MNPS3-L2-P1-R1\` remains **CLOSED** after two scheduler TIMEOUTs with no converged first SCF. Do **not** re-submit the R1 pair.

## Research goal

Confirm whether our proposed **site-/mode-resolved electrical exchange driving force** can even be defined for the published **actual four-Mn L2 AFM** in sliding polar bilayer MnPS3 before more expensive SCF, Wannier or finite-frequency Kubo calculations. Independently inspect **existing** R1 failure evidence and propose a bounded future SCF stabilization experiment, not a new job. A reliable \`BLOCKED\` is a legitimate result.

### Public sources to fetch autonomously

- [Primary MnPS3 source](https://doi.org/10.1038/s41467-026-69635-2), Sun et al., *Nat. Commun.* 17, 3101 (2026): states \`P_1^+(L_2)\` and \`P_1^-(L_2)\`; **Eq.13 \`[E || M_z M_y]\`** is the **nonrelativistic spin-space connection**. Do not replace by L1 Eq.11 or a toy \`m_z T\`.
- [Gate 6 science and exact exchange/SOC split](../docs/BD_GATE6_EXCHANGE_SOC_MAGNETIC_TEXTURE_AUDIT_20261022.md).
- [2026-10-22 freshly verified competitor/equation-level novelty gates](../docs/BD_PEER_UPDATE_20261022_SOT_THz_OCTUPOLE.md). Avoid unearned novelty claims.
- [Four-Mn site spin operator algebra test](../scripts/mnps3_l2_site_operator_gate.py). **Run \`--selftest\` locally first**. Its synthetic site labels are not real MnPS3 results.
- [Earlier R2P0 checklist](BD-MNPS3-R2P0-SITE-OPERATORS-DRAFT.md). This ACTIVE file supersedes its planning status. [R1 public closeout](../docs/BD_MNPS3_R1_CLOSEOUT_20261022.md).

### Private evidence (fetch from user's existing local workspace, not by copying files into chat)

- Local R1 campaign, including \`BD-MNPS3-L2-P1-R1\` / \`preflight.json\`, \`preflight.md\`, two exact POSCARs and source construction transforms, full SHA256, nonpublic site indexing/magnetic seeds.
- Private GitHub summary: \`Xin-Jiaqi/wannauto-research-archive/research/20261021/bd-mnps3-l2-p1-r1/result.json\` and \`results_summary.md\`. **Do not confuse narrative references to \`evidence/preflight.json\` with a file actually present in the private GitHub branch**; it can exist locally even if remote archival copy is absent.
- Two original R1 VASP input decks and actual OSZICAR/OUTCAR/logs, **read-only**. Never transmit licensed POTCAR, private cluster paths or login information into public GitHub or general model context.

## Ordered tasks — no user cut-and-paste

**1. Reconstruct exact source-backed four-Mn L2 mapping.** For **each** polar daughter, determine for every Mn index: chemical site reference, T/B layer, A/B honeycomb sublattice, collinear magnetic sign and axis, identity under the published spatial Eq.13 pairing, and source+POSCAR SHA256. Do not infer four-Mn site order from \`MAGMOM=-1 1 1 -1\` alone. Record all four-moment physical maps; if source files cannot be located, mark \`BLOCKED_SITE_MAP\`.

**2. Run the published four-site algebra preflight.** Create a **private** JSON matching the script fields:
\`\`\`json
{
  "source_doi": "10.1038/s41467-026-69635-2",
  "poscar_plus_sha256": "<private full 64-hex digest>",
  "poscar_minus_sha256": "<private full 64-hex digest>",
  "sites_plus": [
    {"id":"<real Mn site>", "layer":"T", "sublattice":"A", "eta":1}
  ],
  "sites_minus": [
    {"id":"<real Mn site>", "layer":"B", "sublattice":"A", "eta":1}
  ],
  "source_index_to_target": [0,1,2,3]
}
\`\`\`
**Schema example ONLY**, not physical MnPS3 mapping: replace *all arrays/IDs/signs/permutation* with four genuine entries per daughter, input full SHA digests and independently verified site matching. No made-up labels. Run the script \`--selftest\` before real data. A passing script means **ideal spin/site projector algebra** only, not a SOC magnetic-group PASS.

**3. Explicit nonrelativistic-vs-SOC symmetry gate.** The paper uses independent spin-space \`E\` in Eq.13, while the *spatial* product \`M_z M_y=C_{2x}\` (in conventional axes) acts on axial spin as \`(s_x,s_y,s_z)→(s_x,-s_y,-s_z)\` in a relativistic spinor representation. Compare \`spin_E_nonrelativistic\` and \`spin_C2x_axial_SOC_candidate\` from the source-backed site mapping; reconstruct **actual SOC magnetic point symmetry** on each daughter with the true magnetization axis, allowing combined time reversal if required. **Do not certify a physical Kubo response tensor merely from Eq.13.**

**4. Define the physically coupled mode, not only the layer spin.** Distinguish \`S_layer-odd^alpha\` (weights T:+1/B:−1) and \`S_Mn-Neel^alpha\` (four-site ordered \`eta_{la}\`). Define independent local-moment collective coordinates \`theta_nu\` and **electronic generalized force** \`F_nu=-dH_ex/dtheta_nu\`. Enumerate what actual same-gauge spin/site/position/velocity vertices the future calculation would need. If unsupported, \`BLOCKED_RESPONSE_VERTICES\`. Return allowed candidate component(s) of \`E_beta→F_nu(omega)\` **only after** actual magnetic/SOC symmetry.

**5. Read-only R1 SCF stabilization diagnosis.** Explain **individually** the oscillatory positive-state first SCF and the diverging negative-state first SCF, frozen common physical settings, differences between the two INCARs, and the 2 h walltime/16-rank efficiency. Propose at most two narrow pilot arms, each changing only one algorithm/mixing aspect, preserving PBE+U4/D3/ENCUT500/9×9×1 and magnetic state. Do not prescribe 64 ranks blindly, misdiagnose SCF as a geometry/FE issue, or submit either pilot. Record proposed budget and whether new approval would be required (**YES**).

**6. Compare real response novelty, not formula labels.** Specifically match the proposed field-induced exchange force against 2025 vdW layer SOT, 2025 insulator Néel SOT, 2026 MnPS3 THz phonon symmetry dynamics, 2026 spin photocurrents and 2026 octupole spin-transfer torque. Identify the *closest same-observable competitor* and a discriminating measurable prediction. If none, state \`NOVELTY_UNRESOLVED\`.

**7. Archive once, without additional user copying.** Output one PRIVATE \`BD-MNPS3-L2-R2P0-SITE-OPERATORS/result.json\` and short \`results_summary.md\` in the existing private work/archive scheme, with git refs, hashes, true site map, algebra reports, symmetry and SCF analysis, status and blockers. Follow existing \`wannauto-research-archive/MAINTENANCE.md\` branch/PR/secret-scan conventions. Public output may contain **only** sanitized acceptance level, next scientific gate and journal references. Do not overwrite closed R1 truth or auto-update \`CURRENT.md\` to \`COMPLETE\` without actual return evidence.

## Required response to user (one short final message)

\`\`\`text
BD-MNPS3-L2-R2P0 | READ_ONLY
R1: CLOSED / 0 SCF converged
Site map: PASSED_SOURCE_VERIFIED / BLOCKED_SITE_MAP
Four-site operator algebra: PASS / FAIL / NOT_RUN
Actual SOC magnetic group: VERIFIED / PROVISIONAL / BLOCKED
Candidate mode driving force: DEFINED / BLOCKED
R1 SCF failure diagnosis: COMPLETE / LIMITED
Real new compute jobs: ZERO
Private result evidence: <relative record name and commit/hash>
Next decision: one precise physics or approval blocker; NO automatic reruns
\`\`\`

### Scientific and permission hard stops

Zero Slurm/para/WSL/Japan submissions; no \`--execute\`, \`sbatch\`, \`mpirun\`, VASP invocation, modification of historical R1 INCAR/POSCAR or new remotely charged compute. **R2-P0 approval covers source reading, programmatic algebra checks and reporting only.** If current private files are not accessible by the harness, return \`BLOCKED\` with a concrete path-free provenance gap; **do not invent site maps**. The latest public URL is merely an instruction pointer, not proof that the local harness has started.
