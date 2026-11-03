# DeepSeek Harness: zero-VASP IGMC response-theory gate after AgF₂ SOC screening

**Version: 2026-10-23; execution type: READ_ONLY_RESEARCH_AND_LOW_COST_PYTHON; ABSOLUTELY NO NEW VASP/PBS/WANNIER RUNS.**

## 0. Objective

Turn the verified AgF₂ SOC B1 screening into a firm next-stage **go/no-go** on physical IGMC without spending more DFT resources. Use the user's connected WSL and GitHub source archives only for read-only file review; run existing Python symmetry/model unit tests where practical, but do not initiate DFT, use HPC, write any VASP inputs, or make extra resource reservations.

**One and only coordination surface:** https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38  
**Branch:** research/agf2-b1-wsl-relax-20261021. Read latest head; do not assume fee4e62 is still the head. Push additions to this SAME PR, preserving historical receipts and review notes.

Start by reading these new files in the existing PR:
1. \`soc_static_1/B1_SOC_MAGNETIC_AND_RESPONSE_FEASIBILITY_AUDIT_20261023.md\`;
2. \`soc_static_1/B1_RESPONSE_OPERATOR_FEASIBILITY_MATRIX_20261023.md\`;
3. \`soc_static_1/audit_soc_mxprime.py\`;
4. older \`idea_screen_1/B1_IGMC_TENSOR_ERRATUM_20261022.md\` and \`B1_BSF_IGMC_NOVELTY_STORY_GATE_20261022.md\`;
5. existing frozen Xiang model \`papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/docs/SOURCE_MAP.md\`, \`src/igmc_lab/tensor_symmetry.py\`, \`core.py\`, \`independent_checks.py\`.

Important: user's PR #38 source SOC B1 is a valid **electronic screen**, not IGMC or FE discovery. The previous runner's statement that common \(M_z\) and Ag \(L_z\) canting necessarily breaks \(m'_x\) has been independently disproved at printed site-moment precision; use the antiunitary **time-odd axial** mapping instead. The SOC run also violated the 4h authorization, which is documented and must not be normalized or repeated.

## 1. Phase A: reproduce the symmetry audit and correct the physical claim

1. Run using the committed data only: \`python audit_soc_mxprime.py MOMENTS.csv\`. Optionally check the structural glide again with the committed B1 POSCAR.
2. Site pairs are Ag (1,2),(3,4) and F (5,8),(6,7),(9,12),(10,11). For physical \(\,m'_x=\{m_x|\tfrac12,\tfrac12\}\mathcal T\,\), both spin and orbital magnetic vectors obey \((M_x,M_y,M_z)\to(-M_x,+M_y,+M_z)\) between paired sites. Compare all 12 projected vectors at actual printed precision. **Site mapping and magnetic density symmetrization are not the same test.**
3. Record paired residual norms/components and projection uncertainties, including F moments. Explain why uniform \(M_z\) is symmetry-allowed and why \(\,L_z\approx0.07\,\mu_B\,\) is not orbital IGMC.
4. If original private OUTCAR is accessible read-only, confirm moment-coordinate convention (\`SAXIS=0 0 1\`), units and \`LORBMOM\` site/orbital projection; **no file upload** and no new SCF. Do not attempt to artificially force exact symmetry.

## 2. Phase B: source-authoritative Eq. (14) audit (MAIN TASK)

1. Locate the **already-owned** Xiang et al. 2025 PRL **7-page article PDF and 11-page Supplemental Material PDF**, whose hashes are pinned in \`docs/SOURCE_MAP.md\`. These are private archive copies; confirm hashes, read locally, **do not upload publisher PDFs or full SM text**.
2. Transcribe only a **short permitted equation definition** and your own derivation for main Eq. (9)-(14), plus the symmetry tables I–III from SM as independent claims, with page/equation anchors and without reproducing long passages. Identify whether tensors refer to the coefficient of \(\mathbf B\), \(\dot{\mathbf B}\), in-phase/out-of-phase AC field, complex \(\omega\)-resolved causal susceptibility, or real field-quadrature tensors.
3. Compare the actual **printed** Eq. (14) antiunitary sign and SM P/T tables against the code's two explicit protocols \`paper_printed_eq14\` and \`channel_resolved_symmetry\`. Do NOT silently fix the paper; both paths need source-grounded evidence. Evaluate how antiunitary \(i\mapsto -i\), \(\omega\mapsto-\omega\), complex conjugation, and \(\mathrm{Re}/\mathrm{Im}\) or \(\omega\)-odd displacement slope affect the tensor parity.
4. Re-derive the symmetry constraints of \(m'_x\) for C and D from the *actual paper-defined observables*, then the partner signs under unitary \(m_z\) and antiunitary \(m'_z\), with **the same real Cartesian axes**. Test also a pure time-reversal-invariant Rashba model and an explicitly time-reversal-breaking Dirac insulator as independent sanity examples from existing repository models.
5. **Do not declare Eq. (14) wrong** unless source-plus-independent-physical derivation shows a convention-independent contradiction. If ambiguous, return \`BLOCKED_SOURCE_CONVENTION\` and show the two predictions side by side. Pay particular attention to conduction response vanishing in an insulator but displacement slope remaining allowed.

## 3. Phase C: actual material operator readiness (NO computing)

Read the repository and any already connected/available related materials code (if relevant, W90Flow) **read-only**. Inventory whether there is a verified route to:

- SOC spinor Wannier \(\,H_{nm}(\mathbf k)\,\) including correct Ag/F bands and fixed geometry;
- full **same-gauge** interband physical-spin/Zeeman matrices and the correct \(g\mu_B\) prefactors;
- a covariant velocity or position/Berry matrix with intracell orbital embedding, not only HR;
- a multiband non-Abelian occupied Wilson polarization finite-field cross-check;
- separately justified **orbital magnetic-field coupling** (not PAW \`LORBMOM\`);
- clear finite-frequency \(\sigma^D\), \(\sigma^C\) units and Maxwell-Faraday backgrounds.

For each item assign \`EXISTS_AND_TESTED / EXISTS_UNVERIFIED / MISSING / NOT_ACCESSIBLE\`, pointing to exact GitHub paths and tests. Do not conflate previously tested *two-band toy model* with a real-material spinor export interface.

If any item is missing, write the **shortest feasible implementation route** and one falsification/unit test. Do not code a big new framework, automatically install packages, or run Wannier.

## 4. Phase D: research decision, useful even if blocked

Return exactly one prioritized science story, one real blocker and one low-compute validation recommendation:

- **Primary:** Type-I unitary connector could switch a nonzero insulating displacement \(\,D_{xx}\,\) while preserving weak \(M_z\)/orbital \(L_z\) and AHE parity. This must be real-material tested, not just read off symmetry.
- **Alternative:** strained Type-II antiunitary connector predicts even \(\,D_{xx}\,\), odd \(M_z\)/AHE. Note that Type-II has a **different strain/physical branch** in Zhu et al.; not an automatically stable zero-strain phase.
- **Stop conditions:** physically selected tensor conventions invalidate the story; wrong magnetic group; spin-only D response vanishes for actual SOC material; material operator extraction impractical, no measured amplitude; no well-defined domain pair.

Remember \(\,C\,\) is a clean Fermi-surface term and vanishes for in-gap insulating chemical potential; do not suggest simultaneous bulk in-gap conduction-IGMC and displacement-IGMC current fingerprints.

## 5. Required GitHub deliverables (append-only, not fabricated)

Add under \`idea_screen_1/soc_static_1/\`:

1. \`SOC_MXPRIME_INDEPENDENT_RECHECK.json\` (from committed CSV and optional nonsecret private-log cross-check), with 12-site pairs and max vector mismatch.
2. \`IGMC_EQ14_SOURCE_CONVENTION_AUDIT.md\` (primary page/eq anchors, both conventions and causal/real quadrature transformation, clear conclusion \`PASS / PARTIAL / BLOCKED\`).
3. \`IGMC_RESPONSE_OPERATOR_FEASIBILITY_MATRIX.md\` (existing modules, missing same-gauge pieces and concrete tests).
4. \`BSF_IGMC_GO_NO_GO_AFTER_SOC.md\` (one-page physical judgement and exactly one recommended next task, *not executed*).
5. Optional tiny unit tests for the mirror-moment operation or toy-level parity **only** if they reuse current code and require no HPC/big dependencies.

Keep PDF/SM copies, POTCAR, WAVECAR, CHGCAR, large raw OUTCAR/vasprun, private machine paths, credentials, personal identifiers and other non-redistributable sources **out of GitHub**. The original SOC receipt and scientific notes must not be rewritten to hide the budget overrun or old interpretation; new dated corrections only.

Add a short PR Conversation comment with the four results and evidence level, so ChatGPT can independently fetch them without the researcher relaying DeepSeek's output. Never tell the researcher to copy intermediate messages between agents.

## 6. Strict stop boundary

**No VASP, qsub, PBS, B2 structure optimization, MAE, Wannier execution, material IGMC integration, or other paid compute.** No reusing approved SOC hours to justify a new task. Python read-only parsing and unit tests are allowed.

If primary PDFs or related source code are unavailable, mark \`BLOCKED_SOURCE_UNAVAILABLE\` or \`PARTIAL_SOURCE_REVIEW\` precisely. **Deliver the remaining independent symmetry and operator matrix review anyway**, rather than proposing an automatic VASP job.

The goal is high information gain before any new compute authorization.

### One-line launch

“Open PR #38 and fully execute \`B1_SOC_RESPONSE_THEORY_READONLY_HARNESS_TASK_20261023.md\` using GitHub and the already-held local paper/SM archive. This is a read-only research/audit task; no VASP or PBS or Wannier jobs. Push all results to the same PR.”
