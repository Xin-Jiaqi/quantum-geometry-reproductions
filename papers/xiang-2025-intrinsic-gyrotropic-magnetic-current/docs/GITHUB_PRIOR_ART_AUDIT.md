# Public GitHub prior-art check (2026-10-19)

This is outside-source corroboration, not data taken from the user-supplied paper. No third-party source copied into this package.

## Existing personal project architecture

`Xin-Jiaqi/quantum-geometry-reproductions`, default branch `main`, inspected at commit `78a34431bd4a165484ac3f714786b98120132f81`. Its published `jia-2026-nonlin-magnetoelectric-edelstein` paper module supplies repository organization, `docs/PHYSICS_PROTOCOLS.md`, `scripts/verify_release.py`, `tests/`, `integration/stage1_handoff.py`, strict SHA256/frozen regression philosophy, and forked **named-response** conventions. New IGMC work is strictly separate and has not been pushed to this repository. Existing NMEE scientific blockers must not be silently claimed solved or changed.

## Separate public implementation of Zeeman interband geometry

A GitHub content search for Xiang et al. (2025) and the Zeeman-geometry terms located **a third-party numerical package with explicit formula provenance**:

- Repository: `https://github.com/ZhuochengLu/WannierNLQG`, Julia package. Source project is GPL-2.0-only; do not incorporate code into this Python module without checking the license and architecture.
- `theory/QuantumGeometry/ZeemanInterbandQuantumMetric.md`, observed blob SHA `e7d4e1b10182c0d8c2a6bf225e500fe66a4156ea`: defines `g_Z,nm=1/2 Re(A_nm+B_nm)` with `A=r_nm σ_mn` and `B=r_mn σ_nm`. This equals **`Re(z)`** since `B=z*`, consistent with main Xiang Eq. (5).
- `theory/QuantumGeometry/ZeemanInterbandBerryCurvature.md`, observed blob SHA `a3259f0621d9c3b997f4a3daa84a80156e230cd7`: defines `Ω_Z=Re[i(A−B)]` consistent with **`−2Im(z)`** in main Eq. (6).

The package lists Zeeman interband geometric quantities as **k-slice** features. These Markdown source definitions are independently documented conventions; they **do not compute the paper's Fig. 2(d) displacement current**, do not settle the SM Eq. (34) normalization, and are not evidence of the original authors' plotting inputs. No copied code, synthetic data or model normalization was taken from this repository.

## Author numerical source availability

A targeted GitHub content search (full title, `116301 Zeeman`, Zeeman QGT Fig2), plus public-index search, did not identify a verifiable **author-published Fig.2(d) raw dataset or original numerical script** in the inspected results. This does **not** prove that no such data exist. Main Eq. (12) and the SM Eq. (34) local discrepancy remain separate and explicit until source-level numerical evidence is obtained.
