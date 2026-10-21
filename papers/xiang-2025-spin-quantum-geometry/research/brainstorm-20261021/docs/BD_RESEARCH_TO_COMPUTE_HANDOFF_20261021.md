# B+D | GPT → GitHub → DeepSeek harness → WSL/HPC scientific handoff

**Research workflow contract, 2026-10-21. Proposed operating procedure, NOT proof of a connected/executed automation.**

This handoff belongs to **independent B+D magnetic stacking research** inspired by Xiang/Jin/Wang, *PRL* 135, 146303 (2025). Original paper reproduction lives outside this \`research/\` subtree. Nothing in this protocol authorizes editing published-paper acceptance numbers, changing thresholds mid-run, remote HPC submission without owner approval, disclosing private structures, or confusing a model PASS with real-material acceptance.

## 1. Roles and repositories: no new central software package

| Responsibility | Existing owner/repo | What it demonstrably provides | Boundary |
| --- | --- | --- | --- |
| Physics ideas, toy models, formulas, targeted peer audit, scientific go/no-go | [\`quantum-geometry-reproductions\`](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions), **this research subtree** | Dated independent reports/scripts and original-paper isolation; separate [physics capability registry](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/tree/main/physics) distinguishes VERIFIED primitives from UNRESOLVED physical coefficients | No original-paper source code, fitting to expectations, or claim of new effect from tensor symmetry |
| Structural/magnetic symmetry and stacking connection operators | [\`group-theory-operations-toolkit\`](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit) | 80 layer groups, 122 magnetic point groups, 528 magnetic layer groups, LEE/IMHE/optical response tensor screens and connection framework | Structural space/point group alone does not determine AFM ordering; verify magnetic representation and response definition |
| Geometry and physically embedded structure identity | [\`materials-structure-core\`](https://github.com/Xin-Jiaqi/materials-structure-core) | StructureRecord, row-lattice, SHA256/origin, 2D minimum image | 0.x integration candidate; not bonding, energy or FE proof |
| Public benchmark/parent selection | [\`materials-structure-benchmark\`](https://github.com/Xin-Jiaqi/materials-structure-benchmark) | Public provenance, structure IDs, source hashes, 2742 monolayers and 59 bulk-parent references | Existing materials mostly \`not-evaluated\`, NOT validated magnetic sliding FE |
| Produce homobilayer candidate structures | [\`extension-to-BSF\`](https://github.com/Xin-Jiaqi/extension-to-BSF) | Deterministic POSCAR+manifest, shift scans, geometry collision tests, immutable digests | Pure structural engineering alpha; no VASP or FE stability calculation |
| Heterobilayer match, if specifically needed | [\`heterojunction\`](https://github.com/Xin-Jiaqi/heterojunction) | Enumerated supercells, controlled strain/gap and manifests | Does not include relaxations or energetic ranking |
| Independent nonmagnetic crystallographic check | [\`batch-symmetry-checker\`](https://github.com/Xin-Jiaqi/batch-symmetry-checker) | CIF/POSCAR parsing, multi-tolerance spglib reports, hashes | Slab 3D group and 2D magnetic layer group are NOT interchangeable |
| Real SCF/Band / candidate Wannier fitting (private) | [\`wannier-workflow\`](https://github.com/Xin-Jiaqi/wannier-workflow) | Input/POTCAR locks, **dry-run**, bounded attempts, resume/supervise, independent band-error-ranked top-K, frozen \`package-manifest.json\` | Not automatically qualified to compute layer spin QG or full interlayer torque from \`wannier90_hr.dat\` alone |
| NLO shift current/SHG (private) | [\`nlo-response-workflow\`](https://github.com/Xin-Jiaqi/nlo-response-workflow) | Stage-2 manifest/hash handoff, SOC/nocc/symmetry preflight, 2D sheet response | Currently **shift current / SHG** only; do not call it a verified spin-response solver |
| Run evidence, platform divergence, signoff (private) | [\`wannauto-research-archive\`](https://github.com/Xin-Jiaqi/wannauto-research-archive) | Dated evidence, \`CURRENT_STATUS.md\`, \`MAINTENANCE.md\`, runtime snapshots, strict convergence corrections | Historical RUNNING is not live; private data only, no credentials or licensed raw files |
| Local agent control | [\`local-codex-bridge\`](https://github.com/Xin-Jiaqi/local-codex-bridge) | Persistent Codex/DeepSeek thread; distinct local, hpc, maintenance instances; least privilege and task cwd guards | **Only macOS arm64/Codex 0.147.0 tested per repo**; existence of hpc profile does not establish WSL/Japan scheduler connectivity, authorization, or run status |

**Division of labor:** GPT designs/reviews physics, verifies evidence presented, publishes *publicly safe* dated idea/model/peer reports here; **local DeepSeek harness** prepares dry-runs and controls job execution/queue interrogation/acceptance on **WSL or Japanese HPC**; those platforms are the truth for raw calculation outputs. The user remains the scientific/cost/credential approval authority. GitHub records reproducible claims and hashes; it is **not** a job scheduler or the place for confidential raw DFT/Wannier matrices.

**Conflict avoidance:** one active owner for each file/run. Treat GPT edits to paper-inspired \`research/\` notes as independent. DeepSeek owns local/HPC run state and private evidence. Do not let two agents push different versions of the same manifest to default branch; for substantial compute-code changes use a reviewed PR/isolated branch and pin the commit deployed on each platform. Preserve old failures/corrections and refer to them with dated status rather than erasing.

## 2. Mandatory sequence and go/no-go decisions

| Gate | GPT-scientific question | DeepSeek / actual platform checks | Required outcome |
| --- | --- | --- | --- |
| **0: idea & peers** | Distinct mechanism vs Xiang spin QG, insulating NSOT, sliding FE, layer Edelstein, quantum geometric multipoles | None; no compute requested | Specific falsifiable observable, tensor rank/order, rival literature equation, no-go criteria. \`IDEA\` only |
| **1: toy physics** | Correct magnetic connection, \(\pm P\) symmetry, orbital embedding, SOC source, full layer spin/torque response, Ward identity | Small local test may run in the independent model | \`MODEL_PASS\` or \`NO_GO\` with repo SHA, script and independent checks, NOT real-mat success |
| **2: physical candidate preflight** | Has **published evidence** for magnetic order and a genuine \(\pm P_z\) sliding pair? Does its actual magnetic subgroup permit the chosen observable? | Verify source POSCAR/orientation, moments, structure ID/hashes, two candidates/connection map, layer indices, collisions, SOC/functional/U policy | \`CANDIDATE_READY\` only if evidence supports it. Existing benchmark "not-evaluated" does not count |
| **3: cheap structural DFT** | Can published geometry/polar order be independently supported? | WSL small dry-run then actual relax/energy pair only after explicit approval; check both \(\pm\) states, nonpolar competitors, convergence, spin ordering, vacuum, mapping, polarization, path/barrier | \`STRUCTURE_DFT_VERIFIED\` only with all specified controls; otherwise \`BLOCKED/FAILED\` |
| **4: full electronic/Wannier** | Are near-gap spin-orbit bands and subspaces reproducible? | Frozen SCF/Band → top-K; strict native convergence and same-gauge spin/layer/position/velocity operator exports; no unauthorized parameter search | \`MODEL_ACCEPTED\` if independent band/vertex thresholds passed. **Band fit alone is insufficient for spin QG** |
| **5: physical response** | Full complex \(\chi_{S_\ell,E}(\omega)\), \(\chi_{\mathcal T_\ell,E}(\omega)\), proper layer spin continuum, material-specific selection rule vs prior work | Compute only on accepted same-gauge, correct-occupation model; frequency/damping/k convergence, Ward/continuity, \(\pm P\), \(m\to-m\), \(t_\perp\to0\), SOC-off controls; compare WSL and HPC where feasible | \`PHYSICS_ACCEPTED\` only with frozen manifest, physical units, independent oracle and honest novelty review |

**Special B+D test suggested by Gate 5:** the toy can show interlayer electronic spin-transfer torque \(\sim10^{-4}\) with much smaller net \(\dot S_-^x\); this is *not* enough to label a real AFM switching torque. Require an exchange-coupled **magnetic order** dynamic readout and distinguish it from conventional insulating Néel torque.

## 3. Minimal task handoff (copy into private campaign; placeholders, not submitted jobs)

Use the same manually reviewed task identity on both sides. This is a **human-readable example**, not an implemented scheduler/config schema:

\`\`\`yaml
task_id: BD-YYYYMMDD-candidate01-gate3
status: PROPOSED                 # never RUNNING without fresh scheduler evidence
research_note: research/brainstorm-20261021/docs/BD_GATE5_TORQUE_WARD_20261021.md
github_repository: Xin-Jiaqi/quantum-geometry-reproductions
github_commit: "<freeze actual git SHA>"
target_material: "<private candidate ID or published source DOI>"
input_structure_sha256: "<hash of actual POSCAR>"
paired_structures: ["<state_A hash>", "<state_B hash>"]
magnetic_order: "<layer-resolved moments + connection operation>"
scientific_hypothesis: "nonzero interlayer torque spectral channel; independent of known NSOT?"
reference_equation: "<closest literature formula/equation or outstanding audit>"
null_tests: [soc_off, interlayer_off, neel_reverse, polar_partner]
platform: "<WSL or approved Japanese HPC partition>"
calculation_protocol: "<frozen functional, U, SOC, kmesh, vacuum, relaxation>"
code_commit: "<specific source tree actually deployed>"
scheduler_job_id: null            # fill only after confirmed submission
approval: PENDING                 # scientific + time/budget clearance
acceptance_manifest_sha256: null
raw_evidence_location: "<private WSL/HPC store, never public>"
result_status: NOT_SUBMITTED
\`\`\`

An internal scheduler recipe (e.g. Slurm \`sbatch\`, VASP MPI invocation, module names, walltime and allocations) is **platform specific** and belongs in the **private execution package** after checking actual local installation/account. The template **does not create or submit jobs**. No credentials, private endpoints, licensed POTCAR, WAVECAR/CHGCAR or unpublished Wannier/spin matrices are to be copied to the public B+D research directory.

## 4. Acceptance ≠ successful exit status

The private acceptance report for every run must distinguish:

1. **SETUP_OK:** input structure and hashes, environment/dependencies, executable/queue submission, job ID, and correct scientific protocol.
2. **RUN_FINISHED:** scheduler finished normally; no inference yet about physics quality.
3. **NUMERIC_CONVERGED:** verified *native* SCF/ionic/Wannier convergence, not merely output-file presence.
4. **REPRESENTATION_ACCEPTED:** DFT↔Wannier bands, gaps, occupied count, symmetry, spin gauge, and required position/spin/layer operator matrices all independently validated.
5. **PROPERTY_ACCEPTED:** energy/polarization or spin/torque spectrum physically passes paired-state, k/eta and null-control tests in physical units.
6. **NOVELTY_REVIEWED:** outcome compared equation-by-equation against actual published peer mechanisms. This is a science judgment, **not an automatic output classification**.

If any gate fails, mark \`BLOCKED\` or \`FAILED\`, keep the data/attempt hashes, and **never** silently change the convergence thresholds, output signs, frozen wannier rank or choice of reference observable to make the run pass. Re-runs are new or explicitly resumed, linked campaigns.

## 5. Return packet for GPT review

The local harness may return **only** a compact, sanitized \`handoff.md\` or \`results_summary.json\`, containing task ID, source and code commit, input hashes, scheduler completion evidence, explicit statuses by gate, two magnetic/polar state mapping, normalized energies/polarization, fit/response convergence, independent failures and any authorized public figures. Raw complete structures/operators from unpublished work should remain in the private archive, with a private hash/location reference. GPT can then update **this research folder** with model-vs-material status clearly separated and review the candidate's scientific interpretation; it must not claim to have accessed HPC when given only a historic summary.

## 6. Immediate plan — no remote job requested yet

1. **Freeze Gate 5 toy acceptance** in [Gate 5 report](BD_GATE5_TORQUE_WARD_20261021.md) and standalone script. It is a response-mechanism *screening signal*, not a discovery.
2. **Preselect one published real magnetic sliding-FE pair** with known structural file, layer magnetic order, polar switching path and credible source rights (e.g. inspect published MnPS3 bilayer evidence; do not assume it passes). Use existing structure benchmark/core and BSF generator, but mark the candidate \`not-yet-evaluated\` until proven.
3. **Before any DFT**, write a two-state, one-observable **private** preflight task: same moment orientation, explicit \(\pm P\), symmetry allowed/forbidden response, no unusual external field assumptions, and a cheap reference conventional response.
4. **Only after approval** ask the DeepSeek \`hpc\` / local execution instance to run WSL/Japan small tests, not a broad blind high-throughput search.

**Workflow assessment:** the repositories support the requested division of responsibility *as an architecture*. GitHub connector access is confirmed, but **this conversation has not connected to the local DeepSeek bridge, WSL, or any HPC queue**, so no real run, scheduler status or data acceptance is asserted. See the private AutoWann archive for its separate historical platform provenance; don't reuse a historical \`RUNNING\` label.
