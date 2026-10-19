# Current release status, 2026-10-19

## Purpose

Create one expandable quantum-geometry **paper reproduction** repository, one folder per paper. Current first paper: Jia et al., *Nonlinear Magnetoelectric Edelstein Effect*, PRL 136, 226302 (2026).

## Tested and included

- 30 automated tests (20 inherited and 10 new): model band kernels, PT symmetries, Eq(2) tau & length scaling, static-versus-printed protocol split, import regression, signed Stage-1 handoff preflight.
- Recalculated 451×451 Kane–Mele: `printed=-0.010739723280418214`, `full=-0.009548371489696384`, `diagonal=+0.0011913517907195352` in 2D convention `(nm T V)^-1`.
- CuMnAs generated examples `fig1b`, `fig1c-static-diagnostic`, `fig1d-printed`, `fig1-bgk` with fixed model parameters and protocol names; no claim of absolute paper curve identity.
- Exact inverse-pair k mesh used in the NEW CuMnAs CLI for symmetry-protected response checks; legacy modules retained for historic reproducibility.
- `package-manifest.json` SHA256/quality/NOCC/SOC/dimension preflight wired to the existing Stage-1 → Stage-2 export contract.

## Blockers still present

1. CuMnAs Fig.1(c) *paper-defined finite-temperature driven-dc intrinsic* coefficient in exact PT doublets is not established from independent final formula; static and BGK paths are labelled diagnostics, not promoted.
2. PRL printed intrinsic formula and SM full-static density-matrix term have a documented physics protocol discrepancy. Both implementations are provided, neither falsely overwritten.
3. No real projected Wannier spin or position matrices were provided, so no real material NMEE can be calculated by this reproduction alone.
4. No clean new GitHub repository exists yet under the connected account; the GitHub integration in this session supports file/branch commits but does not provide repository creation. Publication requires a repository to exist first.

The source release excludes the 33 internal research ZIPs, original publisher PDF/figures and any private input data. Those retain their separate offline archive.
