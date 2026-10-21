# AgF₂ Stage-2 · G0 preflight reports (non-sensitive)

Status: **G0 done, NOT_RUN, BLOCKED_PROTOCOL, no jobs submitted.**

This directory contains the public, non-sensitive G0 reports produced by the
local DeepSeek Harness run on 2026-10-21 (read-only site probe of the Japan
PBS cluster via the researcher's existing private SSH setup + prepared input
protocols):

- `SITE_G0_REPORT.json` — site facts: PBS 2024.1.2, queues and resource
  maxima, observed VASP cpu builds, MPI modules, private PAW identities and
  SHA-256 references (no POTCAR content, no hostnames, no accounts).
- `INPUT_PROTOCOL_DRAFT.md` — three B1 input sets (POSCAR/INCAR/KPOINTS/PBS
  template) summary + SHA-256, literature-frozen vs operator-choice vs missing
  parameters, Type-I MAGMOM ordering with per-atom mapping.
- `READY_FOR_QSUB_CHECKLIST.md` — advisory preconditions, suggested minimal
  first test, cost estimate, failure criteria. Suggests only; does not run.
- `G0_RECEIPT.json` — `status=BLOCKED_PROTOCOL`, `data_origin=NOT_RUN`,
  `submitted_job_ids=[]`, `physics_approved_by_user=false`,
  `actual_jobs_submitted=0`.

Private artifacts (SSH config, credentials, POTCAR files/content, personal
cluster paths, the full private input workdir) are NOT committed anywhere.
