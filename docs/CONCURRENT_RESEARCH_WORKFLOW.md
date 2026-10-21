# Concurrent Development Contract — Avoid Cross-Paper Conflicts

The repository is being modified simultaneously by other research threads. Parallel work is **safe with conditions**, not automatically race-free.

## Path ownership

| Work | Write scope |
| --- | --- |
| Jia NMEE | papers/jia-2026-nonlin-magnetoelectric-edelstein/** |
| Xiang IGMC | papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/** |
| Xiang Spin QGT | papers/xiang-2025-spin-quantum-geometry/** |
| Xiang PRB map | papers/xiang-2026-quantum-geometric-map/** |
| Independent new science | studies/**, one separate study per directory |
| Handoff / unified interface v1 | RESEARCH_HANDOFF.md, docs/PHYSICS_INTERFACE_SPEC.md, docs/CONCURRENT_RESEARCH_WORKFLOW.md, integration/physics_interface/**, .github/workflows/physics-interface-handoff.yml |
| Shared root indexes | README.md, RELEASE_STATUS.md; **not touched by this task**, to avoid other writers |

This task never modifies existing papers/**, study modules, frozen numeric baselines or their CI.

## Git workflow

1. New purpose-specific branch from an observed main SHA. Never force-push another thread's branch.
2. One PR with a clearly bounded diff. Confirm changed paths before merge.
3. Scope CI to new interface documents and tests; leave each paper's existing CI path rules alone.
4. Before merging, check updated main SHA, mergeability, PR diff and successful Python 3.11/3.12 scoped CI. GitHub normal merge should integrate concurrent *disjoint* paths; conflicts must be resolved rather than overwritten.
5. Each new scientific study records immutable full repo commit SHA and, if used, Group Theory Toolkit full SHA in its study_manifest.json. **Do not pin a mutable branch name alone.**
6. When another paper's entrypoint or scientific method changes, revise the registry via a *separate PR* after new paper-level tests pass. Registry is a snapshot, not an automatically synchronized scientific truth.
7. A new adapter never silently rewrites normalization, tau, tensor axes, printed/derived branches or incomplete paper source results.

## CI and scientific claims

The metadata validator checks the declared paper directories, evidence documentation and target script files exist; it **does not run paper calculations**. It intentionally blocks a missing/renamed CLI. The interface design's evaluate raises NotImplementedError until a separately reviewed adapter exists.

Distinguish explicitly:
- SOFTWARE_CI_PASS: the test suite ran successfully.
- MODEL_OR_METHOD_VALIDATED: named mathematical/numerical subset accepted.
- PAPER_PARTIAL: not all source figure, DC physical or real-material gates completed.
- NEW_STUDY: new independent hypothesis and numerical validation; not attributable to the original paper by default.

## Security, licensing and reproducibility

Publisher PDF/SM, restricted author assets and private notes must not be copied to public GitHub without distribution rights. Do not assume the repository has a universal root license: inspect each folder's licenses, borrowed tables and provenance. For Wannier/DFT workflows, fail-closed if same-gauge projected spin/position/magnetic operators or material parameters are absent.

## Before merging a handoff PR

- [ ] Only the new handoff/interface/workflow paths are changed.
- [ ] Four paper directories remain untouched and independently runnable.
- [ ] Registry references four real checked-in paper entrypoints and evidence files.
- [ ] Physical response protocols are distinct and marked PARTIAL/method-only.
- [ ] No cross-paper numerical solver is falsely advertised.
- [ ] New interface tests and both CI environments pass.
- [ ] Latest PR merge is performed without discarding concurrent changes.
