# Repository collaboration rules

- User preference (2026-11-09): do not enable GitHub Actions email notifications.
  Report relevant check outcomes in the conversation or project worklog instead.
- Research archives, figures and worklogs should not trigger unrelated model
  regressions. Preserve checks for model code, tests and numerical test fixtures;
  do not disable validation or security checks merely to suppress email.
- A skipped workflow is not a passed workflow. Keep local and remote validation
  claims separate, and batch routine result uploads when practical.

## Layer NMEE continuation

- Start at `studies/bsf-nmee-connection/docs/manuscript/worklog.md`: read its continuation entry and current physical status before running anything. Dated sections are history, not automatic next-step instructions.
- The canonical Obsidian worklog and this repository mirror must stay synchronized; preserve the existing manuscript and avoid fragmented new user-facing notes. Raw results and scripts belong in the study result archives.
- User request 2026-11-10: finish the bounded KM + MoSSe/InSe/MoS2 criterion comparison, then stop numerical expansion for a discussion of physics, manuscript, and next directions. A universal criterion is a goal to test, not a result to assume.

- Updated 2026-11-10: the user accepted the representative BiDB screening direction, but authorized only a first-batch proposal and GitHub/worklog archival. Read the latest worklog entry and `results/screening/bidb-pilot-plan-20261110`. Do not submit DFT or a batch until the user approves proceeding after discussion.

- Execution authorization 2026-11-10: user approved full-direction kernel completion and 7269 validation, then the bounded 7268 pilot only after acceptance; KM/existing-DFT analysis proceeds concurrently. Estimate pure runtime after benchmarking. Respect the two-job project cap and Japan compute-node-only science. Other candidates remain planning, not blanket dispatch.

- Structural exception approved 2026-11-10: after 7268 A0 failed the frozen 0.10 eV/Å force gate, user explicitly approved one zero-field z-only relaxation with fixed cell/xy, at most 12 ionic steps and 110 minutes on 4 CPUs. No automatic retry; require force <0.03 eV/Å and bounded displacement, then revalidate zero-field response and new matched isolated references before gate release. Keep raw-geometry result and failure.

- Reporting preference 2026-11-10: start with compact submitted-job/status list, screening progress with an explicit denominator, and expected acceptance clock time in Asia/Shanghai. Then discuss material-physics contribution, counterexamples, manuscript-worthy additions, and remaining PRL evidence gaps. Do not equate jobs finished with candidate acceptance.

- Structural/physics language 2026-11-10: use toolkit BSF/Enhanced BSF terminology; specify coordinate basis/origin, monolayer and bilayer point groups, lower-to-upper Seitz stacking operator and registry, and internal layer-preserving/exchanging operations. Database IDs are secondary indexes. Give induced-spin directions (±x, ±y, or explicitly rotated KM x′) for stated E/B signs, and separate gate-induced amplitude changes from spin direction. Use 支持/不支持/尚无法判断 for hypotheses, not 命中. Preserve distinction between a construction operator and a symmetry of the complete bilayer; old protocols and symmetry-inferred components are labeled.
