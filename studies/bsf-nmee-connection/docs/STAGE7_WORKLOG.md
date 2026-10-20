# Stage 7 reproducibility / handoff log

- Starting point: frozen stage-6 BSF + NMEE snapshot; no edits to paper reproduction under `papers/jia-2026...` and no changes to stages 1–6 files.
- Protocol: noninteracting, T=0, gap-protected equilibrium/static spin–momentum mixed curvature; finite Bz difference. No material or author-driven DC claims.
- Added `valley_kp.py`: from **identical** microscopic H, derive H(K/K'), v_i, and symmetric Hessian W_ij; compute static interband response for full and Taylor models.
- Added stage-7 scripts to compare local k·p order, whole BZ signed capture, channel-mixed derivatives, response zero, C3 sliding scaling, and small material-parameter perturbations.
- Main numerical result: two valleys radius 0.75 nm^-1 recover signed total Γ_xz;x to within 0.13% among four tested u/t combinations; quadratic valley model errors ≤0.0081% for valley-local response.
- Mixed coefficient at t=5 meV, u=.1: BZ -0.0310567; k·p -0.0310929. At u=.4: BZ -0.0922609; k·p -0.0923471. The positive even contribution increases from +0.00354 to +0.1790.
- Predicted valley-model reversal root u≈0.19751 vs previous full-BZ u≈0.1974. C3 trace cubic and traceless/nematic linear responses both recovered within local k·p.
- Modest parameter screen: baseline, Rashba ±20%, intrinsic SOC ±20%, hopping range ±12.5%; all 7 settings preserve endpoint reversal with gapped occupations. No universal sign of coefficients asserted.
- Cross-checks: velocity really differentiates the Taylor Hamiltonian; Hessian finite-difference step stable; full-model kernel in new implementation agrees with stage-6 response module on the same points; ±B and μ gap checked.
- Caveat: k·p is still 8 bands and coefficients come from derivatives of full microscopic H, not closed analytic 2-band Schrieffer–Wolff projection. Signed capture is not absolute integrand capture; local Taylor truncation is not globally periodic.
- Source/provenance reference: Stage6 GitHub main `3947a6fbaeeda3056d5ea38385f743386690a0e6`; prior stage6 offline archive remains unchanged.
