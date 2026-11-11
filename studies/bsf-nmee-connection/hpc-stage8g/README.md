# Stage 8g: autonomous Japan PBS WSe2 pilot, real-run report and independent raw-output science gate

**Updated 2026-10-21:** The user's local DeepSeek harness **reports three actual completed, converged Japan IMR PBS VASP pilot jobs**, not jobs executed by ChatGPT. Numeric input/output summaries: [real run report JSON](RESULTS_HARNESS_REPORTED_20261021.json) and [Stage8h scientific review](STAGE8H_SCIENCE_REVIEW_20261021.md). ChatGPT has **not independently opened the remote raw OUTCAR/CONTCAR**: status is `HARNESS_REPORTED_CONVERGED_REVIEWED__RAW_RECHECK_PENDING`. New [read-only zero-job audit](review_phase1.py) and [Stage8h autonomous task](STAGE8H_READONLY_TASK.json) target this gap. No new PBS jobs are authorized by the review.

The research repository itself now contains every NON-SECRET input needed for a local DeepSeek agent to perform the first computational structural gate: original BiDB 2685 WSe2 source structure; deterministic 15-configuration input builder; one dry-run/submit-on-PBS utility; cryptographic input locks, bounded three-job pilot; a scientific output checker; machine-readable TASK_CONTRACT.json and a full DEEPSEEK_HARNESS_HANDOFF.md. No ZIP extraction or repeated user copy/paste is required.

### Agent start instruction

Open and execute the research-branch file
DEEPSEEK_HARNESS_HANDOFF.md in this folder, with TASK_CONTRACT.json as its
binding operational contract. Use already authorized local WSL/SSH and Japan
PBS resources yourself. Prepare all 15 non-secret WSe2 AB→BA inputs from the
SHA-verified BiDB row2685 source in this Git branch. Preflight the REAL
scheduler, license, local W/Se PAW POTCAR and numerical setup, then submit at
most the three selected A-arm jobs if all gates pass. Monitor, collect and
science-audit outputs, keep resumable job IDs, and report authentic status
to the user. Do not ask the user to manually transfer files or run commands.
Do not modify the frozen NMEE paper reproduction or upload secrets/POTCAR.

### Distinct technical assertions

1. WSe2 input is a BiDB PBE-mono-layer + interlayer z-scan geometry, **not** an author-authenticated fully relaxed bilayer.
2. First-wave jobs are only path_A_s0.00, path_A_s0.25 and path_A_s0.50; 15 generated for deterministic future use, **not approved for immediate automatic batch execution**.
3. Initial non-SOC PBE+D3(BJ) z-only relaxations probe constrained total energy along an AB→derived-BA displacement, **not** a converged switching barrier or genuine CI-NEB.
4. A lower energy at the ideal glide midpoint would still not prove hidden-layer static NMEE; its symmetry merely ALLOWS such a layer-staggered response while forbidding the uniform total-spin static tensor.
5. The user has previously executed Japan PBS jobs, but current Japan/WSL permissions, queue status, VASP binary, licensing and pseudopotential details have NOT been re-probed by ChatGPT.
6. The agent may only report finite energies derived from real OUTCARs, with forces, CONTCAR and immutable input versions; it must distinguish raw completion from reviewed scientific convergence.

### Provenance and separation

- Source paper: Jia et al. PRL 136, 226302 (2026), original reproduction frozen.
- Original research: BSF×NMEE Stages 7a–8f in studies/bsf-nmee-connection; Stage8f physical glide mechanism is separate from current execution plumbing.
- Computational handoff: this Stage8g module and tests only.
- Public code may contain safe, derived geometries and test evidence; locally licensed VASP POTCAR, passwords, private SSH setup, private job output and credentials must **never** be committed.
- **No VASP/QE, relaxation, NEB or real-material NMEE job has been run by this ChatGPT session.** Future progress must be evidenced by the independent local harness and actual scheduler outputs.
