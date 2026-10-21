# DeepSeek Harness execution: MnPS3 L2 polar pair, R1

**Version:** 2026-10-21-R1 | **Status:** READY FOR LOCAL AGENT INTAKE, **NOT YET RECEIVED OR SUBMITTED**.  
**Owner:** local DeepSeek/Codex harness running on user's own machine.  
**Study:** independent B+D research, NOT original Xiang 2025 PRL reproduction.
**Primary real-material target (one pair only):** bilayer MnPS₃ \(\mathbf P_1^+(\mathbf L_2)\) and \(\mathbf P_1^-(\mathbf L_2)\), published path 1.
**Priority:** REAL STRUCTURES + SMALL PAIRED DFT BASELINE, not more toy models, not early Wannier/response calculations.

## Input sources and exact science

1. [Real-material research/physics R0](../docs/BD_MNPS3_REAL_MATERIAL_PREFLIGHT_20261021.md), which includes already-reviewed scientific ambiguities and peers.
2. [R0 source/magnetism preflight](BD-MNPS3-L2-P1-STRUCTURE-R0.md).
3. [GPT–local agent–WSL/HPC responsibility and acceptance rules](../docs/BD_RESEARCH_TO_COMPUTE_HANDOFF_20261021.md).
4. Primary original scientific paper: [Sun et al., *Nature Communications* 17, 3101 (2026)](https://doi.org/10.1038/s41467-026-69635-2); read **original supplementary information**, particularly the source reference monolayer and full layer/moment map, published \(\mathbf L_2\) and path 1.
5. Primary source's **L₂ path-1 connection is [E || M_z M_y], author Eq. (13)** in its spin-space notation. Do **not** replace it by the toy's \(m_z T\) or the source's L₁ Eq. (11). Distinguish nonrelativistic spin-group from actual **SOC magnetic space/point group**. 
6. Closest competing papers: same team's *Science Advances* 2026 six-state spin-texture/NLH study, other group *npj 2D Materials and Applications* 2026 MnPS₃ stacking energy, *ACS Nano* 2026 FE-spin photocurrents, Tang et al. 2025 insulating NSOT. Do not claim generic FE-switched spin readout as new.

### Source protocol to be independently confirmed, not blindly copied

The Nat. Commun. paper reports VASP/PAW/PBE+D3, Mn-d Dudarev \(U_{\mathrm eff}=4\) eV, plane-wave cutoff 500 eV, k mesh 9×9×1, vacuum >15 Å, total energy convergence 1e-6 eV and force convergence 0.01 eV/Å. **These are reference values, NOT a complete validated submission input.** Validate exact POTCAR order/functional, D3 type, MAGMOM site mapping, ISPIN/noncollinear+SOC policy and layer orientation in original SI, and correct VASP template/launcher environment for each host. Avoid ad hoc altered physical settings that prevent paired comparisons.

## Mission / boundaries

You are the **executor and independent numerical QA agent**. **Proceed autonomously within this one bounded job mission:** fetch authoritative structures, perform source/structure/magnetic checks, prepare dry-run, and (if ALL scientific, technical and permission gates are met) **submit and monitor only TWO small paired baseline calculations** on authorized WSL or Japanese HPC. Do not ask user to paste routine logs, paths, command output or screenshots; inspect the local workspace, existing credentials/host profiles and repositories yourself, and maintain state on disk. If a hard blocker requires human decision, ask **one consolidated question** reporting the exact blocker and a recommended single action. Never report a job as submitted before the scheduler/process returns a real ID.

Existing code roles:
- \`extension-to-BSF\`: geometry generation only, and its default clone-both-layers logic DOES NOT perform the required paper's upper-layer mirror/rotation automatically.
- \`materials-structure-core\`, \`batch-symmetry-checker\`, \`group-theory-operations-toolkit\`: structures/symmetry; spglib alone cannot identify an AFM SOC magnetic group.
- Private \`wannier-workflow\`: stage-1 source/strict native verification; **do not launch the full stage1-coldstart chain in R1**, nor change Wannier fitting rankings.
- Private \`wannauto-research-archive\`: historical execution evidence and **private** current-state/evidence handoff; previous RUNNING records are NOT live.
- \`nlo-response-workflow\` is only certified for shift current/SHG; **do not claim it solves layer spin torque**.
- \`local-codex-bridge\` separates \`local\` and \`hpc\` instances; never modify bridge secrets/config automatically, and do not infer Japan access just because the hpc profile exists.

## Sequential autonomous execution — stop on exact failures

### R1.0 Retrieve and identify original structure **without guessing**
- Pull the Nature Communications article and SI/source files into PRIVATE local scratch space, with URL/DOI, copyright notices and SHA256; check whether fully relaxed \(\mathbf P_1^\pm(\mathbf L_2)\) structures are **actually downloadable** (CIF, POSCAR, repository/Zenodo etc.).
- If found, do not retype coordinates. Otherwise obtain a trustworthy monolayer structure with source ID and construct only explicitly **unrelaxed** candidate bilayers using the paper's monolayer orientation mapping \(\{M_z|\tau_o\}\), vertical separation and registry geometry. Apply magnetic \(\mathbf L_2\) mapping only as explicitly supported by SI.
- Strictly verify stoichiometry, lattice row/column convention, species order, layer indexing, all Mn local moments and exact pair operator; record fractional \(\delta\), vertical separation, vacuum, magnetic axis and source coordinate uncertainty. If the SI does not provide enough data to define both physical states and moment mappings unambiguously, report **BLOCKED_STRUCTURE** and stop before any DFT.
- Do not substitute a MnPS₃ AA stacking or the wrong \(\mathbf L_1\) family simply because a POSCAR is available.

### R1.1 Cheap geometry/symmetry gate — **no compute cost**
- Parse structures with the existing structure packages, verify minimum interlayer distances, cells, checksums and no conflicts. Test pairing under the **published** operation to tolerances, including magnetic site mapping and spin orientation, rather than only a visual or point-group label.
- Confirm both structures are source-backed or clearly marked \`GENERATED_UNRELAXED\`. For generated unrelaxed states, ensure at least a consistent starting geometry, not imaginary author-relaxed results.
- Identify the actual collinear magnetic order and check which linear \(\chi_{S_{\ell}^{\alpha},E_\beta}\) components could be permitted for a later SOC study; mark **PROVISIONAL** until the self-consistent SOC magnetism is established.
- Write a frozen \`preflight.json\` + \`preflight.md\`: source URLs/hashes, two POSCAR SHA256, moment site map, transformation, pass/fail reasons and actual code commits.

### R1.2 Platform/evidence gate — **read-only dry-run first**
- Inspect WSL local CPU/RAM, available VASP binaries and license, pseudopotentials, MPI, known existing Python/wannier tools and pre-existing signed VASP launch template; do not create or print credentials. Prefer **WSL** if an actually available, approved short baseline can complete there.
- If WSL lacks essential resources, inspect existing Japanese HPC SSH/scheduler profile, project/partition/account, available allocation, software/module and working submit template **without guessing partition names or requesting new paid capacity**. HPC choice requires an already configured, authorized allocation and scheduler policy permitting these runs.
- Pin executable/build, POTCAR component hashes (private), actual Git SHA(s), scientific configuration, exact units and magnetic configuration; compare both inputs byte-for-byte except intended structural/magnetic differences.
- Use current project's dry-run/check/audit utilities. Before real execution check no duplicate matching job or run folder exists. Never invoke **\`--execute\` on a pipeline that chains Band/Wannier/optical jobs**; R1 is specifically just short baseline electronic-structure runs. If neither host has authorized, genuinely runnable VASP, report **BLOCKED_PLATFORM** and stop.

### R1.3 Bounded material baseline — **submit only if all gates pass**
- For a paper-relaxed author pair: submit TWO matched short **static SCF** tests (one per \(\pm P_z\) daughter) in already-approved run templates.
- For **source-supported but GENERATED_UNRELAXED** candidates: use the pre-existing trusted workflow to run at most TWO matched short ionic-relaxation tests (one per daughter), as the starting geometry is not relaxed; use a bounded, independently tested relaxation input. **Do not run both relax + additional static jobs automatically.**
- Resource/budget guard: **max 2 submitted jobs total** in this campaign; **max 2h wall time per job** (or stricter existing queue limits), using an existing allocation with **no new external spending**. No job arrays, NEB, broad stacking grid, unconstrained retries, DFT+U parameter sweeps, high-cost SOC/noncollinear production loops, multi-GPU training or Wannier stage in R1.
- If these bounds do not fit the real host/job science, the scheduler cannot meet them, or preflight cannot support a meaningful calculation, return \`READY_AWAITING_RESOURCE_APPROVAL\` and a single compact proposed alternative rather than selecting an arbitrary configuration. **Do not submit if approval, license or profile is ambiguous.**
- Submission can use WSL process launcher or \`sbatch\`/other **verified real Japanese scheduler** (never invent a cluster command). Capture actual PID/job IDs, work directory private path, input SHA256, executable version, scheduler/account/partition label (sanitize as needed), time budget, and completion state. No false \`RUNNING\`.

### R1.4 Recovery and scientific QA
- Follow running jobs through established supervisor or scheduler query. Do **not** invent an automation cron/daemon; stop this invocation when the local DeepSeek/Codex environment is about to end, but persist enough state for manual/agent resume. Use one campaign directory and one immutable input fingerprint; dedupe rather than resubmit.
- For completed SCF: check native electronic convergence, \`E_0\`, band/gap if actually available, magnetic moments and atom order, per-state equal protocol, and structural energy differences only when comparable. For completed relax: check ionic + electronic convergence, final forces, lattice, layer separation, magnetic moments, collapse of one polar state into another, spontaneous spin/layer changes; both families must be compared under same U, dispersion and pseudopotential settings.
- A process/scheduler exit code 0 is \`RUN_FINISHED\`, not automatically \`NUMERIC_CONVERGED\`, still less \`FE_VERIFIED\`. Do not claim polarization or switch barrier without proper calculations; no SOC layer-spin response or complete torque observable is expected from R1.
- Failed jobs: preserve actual error provenance and only propose corrective steps with evidence; **zero blind retries** in R1. Do not silently modify INCAR, k mesh, Hubbard U or magnetic seed to chase a result.

### R1.5 Record and return without repeated copy-paste
- Keep complete author SI, POSCAR, VASP/Slurm job scripts, licensed POTCAR, WAVECAR/CHGCAR, remote paths, unpublished calculations and full large data **only in private WSL/Japan local storage**.
- Record signed sanitized run summary and provenance in **private** \`wannauto-research-archive\` \`research/YYYYMMDD/\`, following \`MAINTENANCE.md\` (\`CURRENT_STATUS.md\`, \`MAINTENANCE_WORKLOG.md\`, archive scanner, source/commit versions). This is **only after GitHub branch/PR rules and disclosure review**; otherwise leave a local patch/handoff. Do not push raw private evidence into a public repo.
- Produce a concise sanitized \`results_summary.md\` and \`result.json\` accessible to the user: exact task ID \`BD-MNPS3-L2-P1-R1\`, stages/stop point, source DOI/URLs, frozen input+code hashes, true host/job IDs if available, RUN_FINISHED vs NUMERIC_CONVERGED, gap/energy/moments if measured, and next action; failures are results.
- If configured to push **public** metadata through an approved GitHub workflow, only append status/evidence class and science uncertainties to this B+D research subtree, leaving confidential source/output paths out. Do not overwrite this task prescription, original paper reproductions, benchmark truths, public material structures, or successful older acceptance.

## Important scientific warning

Our prior toy model's \(m_zT\) switch relation and \(\omega\)-dependent 32-band spin/torque cancellation are **NOT** proofs about the specific MnPS₃ \(\mathbf L_2\) switching pair. The authors' nonrelativistic spin-space Eq.13 is not automatically the exact SOC spin/torque response connection. Choose the actual \(\chi_{S,E}\) component only after structure/magnetic symmetry are accepted. Subgap AC layer torque and spin-density work is **future R2**, requiring separately validated same-gauge spin, layer/position operators and a complete finite-\(\omega\) Kubo oracle. **Do not compute fake "spin Hall/NMEE" observables from hr.dat alone.**

## Return format (do not require user to collate logs)

\`\`\`text
BD-MNPS3-L2-P1-R1 | STATUS: PREPARED / BLOCKED_* / SUBMITTED / RUNNING / FINISHED / FAILED
Source structures: VERIFIED_SOURCE / SOURCE_RECONSTRUCTED_UNRELAXED / BLOCKED
Magnetic L2 mapping: PASS / FAIL / PROVISIONAL
Host: WSL / JAPAN_HPC / NONE (what verification actually happened)
Input hashes and code SHAs: ...
Job 1, job 2: actual identifiers or NOT_SUBMITTED
Native convergence: YES / NO / NOT_FINISHED / NOT_APPLICABLE
Comparable pair energies and magnetic order: only if truly available
Private evidence pointer: hash/manifest ID (no secret path)
Action required from user: NONE or ONE precise blocker
Next scientific gate: REAL_R2_ONLY_AFTER_REVIEW
\`\`\`

**Instruction on activation:** This versioned file should be invoked through the stable [CURRENT.md](CURRENT.md) task entry. Simply cloning this public repository or reading this document **does not start or submit anything**. An agent with actual user-local file/SSH/scheduler access must carry out the steps and report its own outcome.
