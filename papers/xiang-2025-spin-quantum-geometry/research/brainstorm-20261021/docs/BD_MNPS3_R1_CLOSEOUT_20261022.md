# MnPS3 bilayer B+D real-material R1 | closed without numerical convergence

**Recorded:** 2026-10-22. **Status:** R1_CLOSED_BY_USER; **RUN_FINISHED / NUMERIC_CONVERGED=false / FE_VERIFIED=false**.  
**Scope:** sanitized status from the user's locally executed DeepSeek harness, cross-checked against the **private** \`wannauto-research-archive\` results summary, machine-readable result and current-status/worklog on its authenticated GitHub branch. **The research assistant did not operate the HPC scheduler or independently inspect raw VASP outputs in this turn.**

## Mission and original scientific input

- Real candidate: bilayer MnPS3, paired \`P_1^+(L_2)\` and \`P_1^-(L_2)\` path-1 states from [Sun et al., Nature Communications 17, 3101 (2026)](https://doi.org/10.1038/s41467-026-69635-2). This is **independent B+D research**, not Xiang et al. PRL source reproduction.
- Both atomic models originated from **source-anchored/generated unrelaxed candidates**, **not downloaded and verified original author-relaxed POSCARs**. Literature geometry and magnetic-site mapping were checked and documented privately.
- Preflight results in private evidence: the L2 magnetic ordering/site mapping passes for both daughter structures, published Eq.13 mirror/connection site check gives an RMS displacement of **0.0000 Å**, and the structural/magnetic screening reports the expected **Cm** spatial family and plausible chemical bonds. These **geometric/symmetry checks did not depend on the unsuccessful SCF**.
- Passing a geometry/magnetic symmetry check is **not** equivalent to verifying self-consistent magnetism, energetic FE stability, a Berry-phase polarization or an actual spin/torque response.

## Executed computational experiment and scientific outcome

1. The local harness evaluated WSL and another Japanese compute route, but selected an **explicitly user-approved existing para Slurm allocation** after technical access/resource checks. A separate existing-host preflight verified the chosen DFT+U/dispersion settings. **No unsupported claim of Japan HPC execution**.
2. The authorized campaign submitted **exactly two matched paired magnetic ionic-relaxation jobs** to Slurm with a two-hour wall limit each. Both jobs reached **TIMEOUT / scheduler-terminated** at the hard walltime, with the scheduler reporting \`exit 0:0\`. Here \`RUN_FINISHED=true\` means **the job lifecycle ended**, **not a successful VASP relaxation**. Job identifiers, exact input fingerprints, platform paths and detailed logs remain in the **private** evidence archive.
3. Neither run obtained a converged first electronic SCF solution: one displayed charge-mixing oscillation across many diagonalization iterations, the other a strong SCF divergence. **Zero ionic steps were completed**. There are **no qualified paired total energies, local magnetic moments, band gaps, self-consistent relaxed geometries, FE barriers, polarizations or spin spectra**.
4. Frozen input provenance is present in the private R1 record. R1 did not silently retry, change U/functional thresholds or extend allocations.

**Acceptance matrix**

| Gate | Outcome | Interpretation |
| --- | --- | --- |
| Published-state / exact-site / L2 mapping preflight | **PASS** | Symmetry/geometric setup valid **for the specified source-derived input models only** |
| Actual authorized scheduler execution (two bounded attempts) | **RUN_FINISHED=true** | Both reached walltime limit; \`exit 0:0\` is **not** electronic or ionic convergence |
| Native electronic self-consistency | **FAIL** | Oscillation / divergence during first SCF |
| Ionic relaxation | **NOT ACHIEVED** | Zero accepted ionic steps |
| Verified FE bistability / inter-state energy / layer spin or torque properties | **NOT VERIFIED** | No trustworthy numerical result to use in physics conclusions |
| User decision | **R1 CLOSED** | Explicitly accept the negative numerical finding as the R1 stopping point; **no additional jobs are authorized** |

## Scientific interpretation / useful negative result

The original B+D hypothesis remains **UNTESTED ON REAL MATERIAL**, not falsified or confirmed: the two self-consistent daughter Hamiltonians and their spin-resolved response operators were never obtained. R1 is a valuable **computing-protocol and source-preflight diagnostic**, not a physical NO-GO theorem.

The evidence suggests numerical stabilization of the existing frozen DFT settings should be separately benchmarked if/when R2 is approved. The private archive records options such as an alternative electronic algorithm or modified mixing/damping and a resource scaling check, with the published PBE+U/D3/structural physics unchanged. **These are proposals only, not actions taken and not permission to retry.** No 64-rank run, extra 2-hour cycle, Wannier workflow, spin current or layer-torque calculation was executed under this closeout.

## Provenance and records policy

**Private source of truth:** \`Xin-Jiaqi/wannauto-research-archive\`, dated R1 evidence record \`research/20261021/bd-mnps3-l2-p1-r1/\`, plus current-status/worklog update. This public note deliberately omits raw/hyper-specific private computational inputs, scheduler IDs, remote paths, licensed POTCAR and detailed OUTCAR, using only verified, non-sensitive scientific status from the user-authorized record.

**Public research control plane:** [single-link task pointer](../tasks/CURRENT.md) should now read **R1 CLOSED / NO ACTIVE COMPUTATION**. Its former [R1 execution prescription](../tasks/BD-MNPS3-L2-P1-EXECUTE-R1.md) is historical only. Do **not** automatically re-submit this historical task from a still-cached DeepSeek instruction.

**No original-paper acceptance change:** all B+D toy analyses remain labeled *independent model tests*, and published PRL reproduction state is not modified. A new real-material run would require a **new R2 task ID**, scientific and resource approval, preserved R1 fingerprints, an explicit convergence-fix A/B protocol and separate native numerical acceptance.
