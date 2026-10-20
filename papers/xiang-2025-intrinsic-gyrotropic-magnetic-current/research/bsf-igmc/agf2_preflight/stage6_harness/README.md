# AgF₂ BSF × IGMC: Stage 6, Harness execution handoff

**Research extension only** under `papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage6_harness/`. The Xiang PRL reproduction and frozen benchmarks are untouched.

**Current status: COMPUTATION NOT STARTED / SOURCE STRUCTURES AND AUTHOR DFT RECIPE NOT VERIFIED.** This package is a *real DFT task contract and scientific audit interface*, not a simulated calculation or a claim that AgF₂ has a computed IGMC.

## First action for Harness

Read [`HARNESS_TASK.md`](HARNESS_TASK.md). Search the official Zhu et al. Supplemental Material and associated author-provided resources for the actual monolayer/bilayer coordinates and computational protocol. The printed relative translations `(0.5,-0.16)` and `(0.5,+0.16)` are **not** POSCAR files. If source-quality structures or computational details remain unavailable, submit a machine-readable `BLOCKED_SOURCE` report and request review; do not create fictitious structures or select an undocumented Ag Hubbard U.

## Target physics and priority

The *insulating* low-frequency spin-Zeeman current is the time derivative of a field-induced electronic polarization, `j_i^D(t) = d[alpha_ij^{PB} B_j(t)]/dt`. Our symmetry prediction concerns `D_xx` **and** `D_yy`: for the reported Type-I `m_z` connection the two domain coefficients are odd, while for a **separate 4% strained** Type-II `m'_z` campaign they are even. These are **inter-domain relations conditional on the magnetic mapping**, not evidence that either tensor component is numerically nonzero.

Prioritize zero-strain Type-I `B1/B2` as the initial campaign. Do not compare zero-strain Type-I against 4%-strained Type-II as if they were two branches of the same unstrained sample. The published `0.26 meV/f.u.` Type-I/II splitting requires unusually tight independent convergence and sensitivity checks before ranking magnetic states.

## Auditable handoff (no raw private DFT files in Git)

Local private run directory (outside this public repository):

```text
<private_campaign_root>/
  source/              authenticated structure and method sources
  inputs/              B1 / B2 POSCAR, frozen templates, PAW lock (POTCAR stays private)
  runs/                DFT output, relaxation, SOC, bands, Wannier
  receipts/            one verified receipt per stage plus final `stage6_receipt.json`
```

A **schema example only**, without synthetic material results, is supplied at [`stage6_receipt.template.json`](stage6_receipt.template.json). Run:

```bash
python stage6_harness/audit_receipt.py --receipt /private/receipts/stage6_receipt.json
python -m unittest discover -s stage6_harness/tests -v
```

Commands assume the working directory is `agf2_preflight/`. `audit_receipt.py` reads a local receipt and, when a directory `--root` is given, verifies local file SHA256; it never opens a VASP job or uploads raw data. The output is a staged status and response-parity diagnostic, not an automatic certification of first-principles physics.

Publish only the **reviewed, sanitized** receipt and convergence tables after consulting the user. Private VASP outputs, PAW files, MAGMOM-sensitive unpublished models, and licensed raw inputs must not be pushed automatically. Do not overwrite earlier Stage 2–5 frozen values.

## Baseline sources and blockers

- Zhu et al., *PRL* **135**, 056801 (2025), DOI `10.1103/dmzg-ck2t` (AgF₂ BSF and magnetic order).
- Xiang et al., *PRL* **134**, 116301 (2025), DOI `10.1103/PhysRevLett.134.116301` (spin-Zeeman IGMC, still PARTIAL on published Fig. 2(d) ordinate).
- Existing Stage 4 `agf2_preflight.py` is the first *nominal* input lock check; passing it does not validate structural provenance or magnetic symmetry.
- Existing `wannier-workflow` Stage-1 exports **do not** include all projected spin and full position operators required by quantitative material IGMC. Stage-6 explicitly gates their extension and validation.
