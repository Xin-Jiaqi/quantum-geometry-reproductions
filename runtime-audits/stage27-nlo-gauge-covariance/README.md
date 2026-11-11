# Stage 27 — NLO shift-current: source-gauge vs cell-gauge covariant geometry audit

**Independent, third-party-reproducible runtime audit.** Executable oracle +
raw JSON + fixture + redacted log for the G1–G8 gauge-covariance checks
described in the NLO review (PR [Xin-Jiaqi/nlo-response-workflow#3](https://github.com/Xin-Jiaqi/nlo-response-workflow/pull/3)).
All files here were produced by the exact scripts and fixture in this
directory; `SHA256SUMS` pins their digests.

> Scientific status: **fixture-scope evidence, not a universal theorem.** The
> measured constant `sigma/(-raw) ≈ 3.1406` is a finite-difference estimate
> **approaching π** under the exact conventions below for ONE spinless Haldane
> fixture; it is NOT an analytic proof and NOT a universal FDR equality. The
> earlier `-π/2` claim remains withdrawn.

## Frozen conventions (identical across every comparison)

| Quantity | Value / convention |
| --- | --- |
| H_src | analytic 2-band spinless Haldane, M=0.4, φ=π/2, t1=1.0, t2=0.2, NN=(3,2), NNN=(3,2) (d-vector form, eV) |
| Gauge | U(k)=diag(1, e^{−ik·δ1}), δ1=(1,0) Å; H_cell(k)=U†H_src U |
| τ (centers) | τ_a = i(∂_a U)U† in Cartesian-k units; verified against `fixture/wannier90-wcenter.dat` (G0b, dev 1.7e-9) |
| Covariant derivative | D_a X = ∂_a X + i[X, τ_a] (commutator i[X,τ_a], NOT i[τ_a,X]) |
| k mesh | 5×5×1 primary; 9×9×1, 13×13×1 convergence (reduced, Γ-centered) |
| η_sp | 0.05 eV (Lorentzian delta η_sp/π/((ω−ε)²+η_sp²)) |
| η_band | 0.005 eV primary; scan {0.05, 0.02, 0.01, 0.005, 0.002} |
| Components | (111)=xxx, (122)=xyy, (212)=yxy |
| Field symmetrization | C_sym = ½(C + C^{bc↔}); σ from Im C_sym with ordered pair n<m, f_nm=f_n−f_m |
| Projector FD | central differences, h_d1=1e-4 (first/covariant), h_d2=3e-3 (second) |
| Raw assembly | exact mirror of NLO `_compute_block` (r=v·f, iabc+iacb, occdiff, 1/Nk) |
| Strong-point mask | |raw| > 0.02 for residual reporting |

## Reproduce

```bash
# numpy required (any recent version); no other dependencies
python3 stage27_gauge_covariance_audit.py <output-root> <fixture-dir>
# e.g.
python3 stage27_gauge_covariance_audit.py ./out ./fixture
# writes <output-root>/STAGE27_GAUGE_COVARIANCE_AUDIT.json
```

The three `*.py` modules must stay in the same directory (the audit imports
`stage26_oracle` and `stage26_second_derivative_validation`, which implement
independent HR/centers readers, Fourier kernels and the projector machinery;
they contain no NLO source-code imports).

## Acceptance record (this fixture, 2026-10-23)

| Check | Meaning | Value |
| --- | --- | --- |
| G0 | H_cell analytic == HR-kernel H(k) | 7.9e-16 |
| G0b | τ from U == wcenter differences | 1.7e-9 |
| G1 | D_a H_cell == U†(∂_a H_src)U | 9.6e-9 |
| G2 | D_a P_cell == U†(∂_a P_src)U | 1.2e-8 |
| G3 | σ(P_src, ∂) == σ(P_cell, D) — gauge covariance | 7.0e-5 |
| G4 | σ(P_cell, ∂) ≠ σ(P_cell, D) — plain ∂P_cell is wrong object | 1.38 |
| G5 | σ(P_cell, D) == −π·R_dressed (NLO production raw) | 3.2e-4, const 3.1406 |
| G6 | σ(P_cell, ∂) == −π·R_plain (gauge-dependent pair, kept) | 6.3e-4 |
| G7 | R_dressed == R_src | 1.4e-15 |
| G8 | σ(P_src, ∂) == −π·R_src | 3.5e-4, const 3.1407 |
| G5b | η_band stability (max 0.0059 @0.05, else ≤6.2e-4); k-grid 5→9→13: 3.2e-4→2.0e-4→1.7e-4 | bounded / convergent |

Pairing table in the JSON documents ALL σ × raw combinations; only the
matched pairs (src/cov ↔ R_src/R_dressed; plain ↔ R_plain) are consistent
(ratio ≈ 1.000), mixed pairs give ≈0.53 or ≈2.16.

## Evidence provenance

- NLO fix under audit: commit `299f20e` (Stage 26 cartesian slot fix) on
  branch `fix/shift-current-cartesian-derivative-stage26`; Stage 27 branch
  head `1eb4cf0` (PR #3).
- Fixture: synthetic spinless Haldane package (stage-24 model, 2×2 cell-gauge
  HR with Wannier centers), same data used by the NLO CLI regression runs.
- Raw output: `evidence/STAGE27_GAUGE_COVARIANCE_AUDIT.json` (identical values
  reproduced in this repository's layout; deterministic).
- Redacted log: `evidence/stage27_gauge_covariance_audit.log.txt` (machine
  paths removed).

## Scope / unresolved

Applies to the tested spinless Haldane fixture and the agreed covariant
connection only. Not covered: literal FDR SM Eq. (75) normalization/sign
(original text unavailable, `SKIPPED_NO_FDR_ARCHIVE`; the projector form is an
independent reconstruction of arXiv:2409.16358 Eqs. (3)+(6)), absolute SI
prefactors, 2D sheet normalization (Lz×0.1), bath/scattering protocol, and
DFT/W90Flow material certification. The NLO historical legacy spectra and
Stage 26/27 numerical results are untouched.
