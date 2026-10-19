# Release gates: paper-model code and scientific status

## Offline one-command acceptance

From this paper directory:

```bash
python -m pip install -e .
python scripts/verify_release.py
```

The command fails if an inventory file changes without updating `results/RELEASE_MANIFEST.json`, if a frozen model result drifts, or if the unit tests fail. Five frozen outputs are recomputed, including the Kane–Mele 451×451 reference. It requires **no private research files**. Its numerical tolerance is a reproducibility tolerance, not the claimed uncertainty of published experimental or theoretical data.

This verifier **must not** turn the physics status to PASS. Source-level unit tests, source-file integrity, and a frozen calculation are independent of the unresolved intrinsic finite-T driven-DC response protocol.

## Scope of model tests

| Subsystem | Accepted | Not accepted / explicit limit |
|---|---|---|
| CuMnAs Eq. (5), Fig. 1(b) | PT doublets, energies and Dirac-node structure | exact digital reconstruction of the original figure |
| CuMnAs intrinsic Fig. 1(c) | staggered symmetry, static/finite-frequency checks, separate BGK response | paper-defined degenerate finite-T driven-DC intrinsic coefficient; fixed-grid Dirac limit |
| CuMnAs extrinsic Fig. 1(d) | published Eq. (2) Fermi-surface projector kernel, tau scaling, Neel reversal | unique real-material absolute lattice normalization, microscopic disorder response |
| Kane–Mele Fig. 2 | tight-binding kernel, finite insulating NMEE, 451² frozen reference, K/K' hotspot diagnostics | all published curves agree quantitatively; static Eq. (9) equals published DC formula |
| W90Flow Stage-1 handoff | signed package preflight and fail-closed operator validation | end-to-end numerical NMEE on a real material with HR alone |

## Preservation

Historical rounds and private notebooks remain offline in the research archive; this public folder is the **curated release**, not the full archive. Reproduction scripts never delete or edit frozen outputs unless the user explicitly redirects their output.

## Scientific correspondence

One physical-definition question is prepared in `docs/AUTHOR_DISCUSSION_BRIEF.md`. Sending an email, changing the theory conventions, or marking Fig. 1(c) as paper-validated requires a subsequent independent decision.
