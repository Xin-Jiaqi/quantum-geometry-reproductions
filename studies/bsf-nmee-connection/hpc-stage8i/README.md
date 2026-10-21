# Stage 8i: two-job WSe2 spinor-SOC feasibility test

**Ready for local DeepSeek harness execution, NOT YET submitted.** This stage is independent original BSF × NMEE research, NOT a Jia et al. PRL reproduction or a quantitative NMEE calculation.

**Direct harness entrypoint:** [DEEPSEEK_HARNESS_HANDOFF.md](DEEPSEEK_HARNESS_HANDOFF.md). Machine contract: [TASK_CONTRACT.json](TASK_CONTRACT.json).

- Exact raw Stage8h `CONTCAR` SHAs authenticated before source files may be reused; *not* inferred from ideal structure or re-created from BiDB.
- `prepare_soc.py` produces only the real AB and glide midpoint 3R-WSe2 nonmagnetic spinor-SOC SCFs (21×21×1 sampled BZ, no ionic relaxation).
- `onsite_soc.py` executes validated LOCAL-on-Japan PBS dry run, at-most-two-job submission, read-only status and output audit. Site must have an actual authorized SOC-capable VASP binary, often `vasp_ncl`. Previous Stage8g `vasp_std` does not by itself prove LSORBIT support.
- `eigenval_gate.py` reports sampled 2D indirect and direct gap with assumptions and warning flags, **not** a rigorous continuum-BZ gap or nonzero spin transition matrix elements. Actual remote `PROCAR` spin projections are a diagnostic, not the cross-gap Kubo response.
- No POTCAR, WAVECAR, CHGCAR, private outputs or credentials are included, and none should be pushed to GitHub.

**Most important physics:** Stage9 exact spin-conservation null test means large Ising W SOC does not guarantee an in-plane transverse equilibrium/static NMEE. The two calculations ask whether a quantifiable microscopic spin-mixing channel is plausible enough to justify a later independent operator-complete Wannier test. Neither the static total response nor the hidden layer-staggered response has yet been calculated for real WSe2.

No new job has been submitted from this repository update. The local harness execution report must retain provenance and distinguish real completed SCFs from mere input preparation.
