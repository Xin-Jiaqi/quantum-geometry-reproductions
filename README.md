# Quantum Geometry Reproductions

论文级量子几何理论与数值复现资料库。**每篇论文一个独立目录**，各自保留公式约定、计算程序、可执行示例、回归测试与未解决的科学问题。后续工作可以引用已经验证的论文模块，而不修改其冻结基准。

## Papers

| Paper | Main physics | Method status |
|---|---|---|
| [Jia et al., PRL 136, 226302 (2026)](papers/jia-2026-nonlin-magnetoelectric-edelstein/) | Intrinsic and extrinsic nonlinear magnetoelectric Edelstein effects (NMEE) | Tight-binding method reproduction and workflow integration tests; scientific DC-limit question documented |

未来的 quantum geometry 工作新增到 `papers/<first-author>-<year>-<short-title>/`，不要直接在第一篇论文目录中添加不相关的新模型。

## Reproducibility principles

1. Every numeric result states its response **protocol**, parameters, BZ mesh, physical units, and caveats.
2. No unverified use of `wannier90_hr.dat` as a complete quantum-geometric input; projected spin/position operators must be provided explicitly for a quantitative material claim.
3. Distinguish an independently verified *static equilibrium* coefficient from the paper's *driven dc* coefficient when the two differ.
4. Reproduction results are not the original authors' data. No publisher PDFs, original figure crops, private material data or computational secrets are committed.
5. Results in this repository are research implementations and **not** evidence that every figure has been pixel-perfectly reproduced.

## Quick start

```sh
cd papers/jia-2026-nonlin-magnetoelectric-edelstein
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/reproduce.py --case fig2-printed-vs-static --mesh 451
python scripts/reproduce.py --case fig1d-printed --mesh 151 --temperature 100 --mu 0.03 --theta 0.15
```

See the individual paper directory for physics, normalization, and Stage-1 Wannier package handoff details.

## Status and license

This is an independently developed research reproduction, **not** an official implementation by the authors. This repository is public, but no explicit reuse license has been selected. Public access does not by itself authorize reuse. Each paper folder states its validation scope; software tests do not imply full scientific agreement.

## Offline reproducibility gate

From the Jia et al. paper directory run `python scripts/verify_release.py`. This checks SHA256 integrity, five frozen model outputs, and the full regression suite. See [`RELEASE_GATES.md`](papers/jia-2026-nonlin-magnetoelectric-edelstein/docs/RELEASE_GATES.md).
