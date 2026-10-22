**Stage8l update 2026-10-22:** New [scientific integrity review](STAGE8L_EVIDENCE_RECONCILIATION_AND_EFFICIENT_NEXT_GATE.md) detected 7 nondegenerate Stage8k key_result fields incorrectly labeled band35 but storing actual band36 energies/spins. Raw reader has been guarded and the code now proposes a [zero-PBS-job pseudo occupied-empty spin vertex gate](stage8l_pseudo_interband_gate.py) with [autonomous agent handoff](STAGE8L_DEEPSEEK_ZERO_JOB_HANDOFF.md). Prior historical evidence remains unchanged; physical PAW spin, momentum vertex and static WSe2 NMEE are still uncomputed. VASP official ALGO=None+LORBIT=11 WAVECAR postprocessing is preferred as a *separately permissioned future job* over blindly repeating two full SOC SCFs.

# Stage 8i: two-job WSe2 spinor-SOC feasibility test

**Updated 2026-10-22: Two real Japan spinor-SOC VASP SCFs have been performed by the local DeepSeek harness, with usable sampled EIGENVAL data but MPI exit=255 and truncated PROCAR. [Stage8j independent corrective science review](STAGE8J_INDEPENDENT_VALLEY_SPIN_ERRATUM_20261022.md) finds misidentified Kprime and spin-splitting partners in the original report. Verified sampled positive gaps remain, but true non-Ising spin-mixing and static NMEE are NOT established. The original task files below remain useful as historical execution provenance.** This stage is independent original BSF × NMEE research, NOT a Jia et al. PRL reproduction or a quantitative NMEE calculation.

**Direct harness entrypoint:** [DEEPSEEK_HARNESS_HANDOFF.md](DEEPSEEK_HARNESS_HANDOFF.md). Machine contract: [TASK_CONTRACT.json](TASK_CONTRACT.json).

- Exact raw Stage8h `CONTCAR` SHAs authenticated before source files may be reused; *not* inferred from ideal structure or re-created from BiDB.
- `prepare_soc.py` produces only the real AB and glide midpoint 3R-WSe2 nonmagnetic spinor-SOC SCFs (21×21×1 sampled BZ, no ionic relaxation).
- `onsite_soc.py` executes validated LOCAL-on-Japan PBS dry run, at-most-two-job submission, read-only status and output audit. Site must have an actual authorized SOC-capable VASP binary, often `vasp_ncl`. Previous Stage8g `vasp_std` does not by itself prove LSORBIT support.
- `eigenval_gate.py` reports sampled 2D indirect and direct gap with assumptions and warning flags, **not** a rigorous continuum-BZ gap or nonzero spin transition matrix elements. Actual remote `PROCAR` spin projections are a diagnostic, not the cross-gap Kubo response.
- No POTCAR, WAVECAR, CHGCAR, private outputs or credentials are included, and none should be pushed to GitHub.

**Most important physics:** Stage9 exact spin-conservation null test means large Ising W SOC does not guarantee an in-plane transverse equilibrium/static NMEE. The two calculations ask whether a quantifiable microscopic spin-mixing channel is plausible enough to justify a later independent operator-complete Wannier test. Neither the static total response nor the hidden layer-staggered response has yet been calculated for real WSe2.

No new job has been submitted from this repository update. The local harness execution report must retain provenance and distinguish real completed SCFs from mere input preparation.

## Stage8k isolated-band spinor prescreen (zero jobs)

- [Stage8k physics review](STAGE8K_NONDEGENERATE_SPIN_SCIENCE_REVIEW.md) and [frozen candidate JSON](STAGE8K_NONDEGENERATE_SPIN_CANDIDATES_20261022.json) identify three nondegenerate AB spinor pseudo-spin targets from already-archived data; the midpoint contrast is **not yet physical** without PAW operator/validated WAVECAR parser.
- [Zero-job local DeepSeek self-contained assignment](STAGE8K_DEEPSEEK_READONLY_HANDOFF.md), [candidate selector](stage8k_existing_spin_prescreen.py) and [offline regression](../tests/test_stage8k_nondegenerate_spin_prescreen.py).
- No additional VASP, Wannier, spinor Kubo or actual static NMEE has been run by this Stage8k addition.

## External 2026-10-22 result and independent erratum

- [Original sanitized SCF / sampled-gap audit](STAGE8I_AUDIT_RAW_20261022.json); [spinor texture evidence](STAGE8I_SPIN_TEXTURE_20261022.json) and [original valley diagnostic](STAGE8I_ISING_SPLIT_20261022.json) retained as historical outputs.
- [Original DeepSeek scientific report with prominent independent erratum](STAGE8I_REVIEW_CONCLUSIONS_20261022.md).
- [Stage8j independent reciprocal-space/spin-pair corrections](STAGE8J_INDEPENDENT_VALLEY_SPIN_ERRATUM_20261022.md), backed by [standard-library evidence rechecker](stage8j_evidence_check.py) and [4 regression tests](../tests/test_stage8j_valley_spin_erratum.py).
- This correction is **analysis of already committed evidence, ZERO NEW VASP jobs**. Read existing real EIGENVAL at K and its true time-reversal partner, validate WAVECAR interpretation and nondegenerate off-axis spinor mixing before authorizing heavy operator-complete Wannier/NMEE work.
