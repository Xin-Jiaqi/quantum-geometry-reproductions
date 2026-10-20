# Stage 4 reproducibility log

- Existing published PRL reproduction, stage 1, stage 2 and stage 3 files remain unchanged.
- Added `bsf_nmee/response_mechanism.py` to implement the algebraically equivalent occupied/empty interband expression for the *static* mixed spin/momentum curvature derivative.
- Added three explicitly controlled hopping variants: full, same-only, cross-only, plus a physical relative AB/BA hopping-sign counterfactual. Independent Hamiltonians are diagonalized each time; no fitted response scaling.
- Tested model-specific nonadditivity \(\Gamma[T_s+T_c]-\Gamma[T_s]-\Gamma[T_c]+\Gamma_0\) and exact layer-gauge evenness in global tunneling sign.
- Recorded zero-root movement from \(u\approx0.1973\) at 20 meV to \(u\approx0.2676\) at 10 meV, both with nonzero sampled insulating gaps. Unreleased exploratory SOC/decay variations are not used to establish a universal claim.
- Compared new mixed-curvature expression to the original finite-B susceptibility and the denser valley quadrature.
- Documented limitations of matrix-element/denominator separation under possible degeneracy. Total response remains the core observable.
- `results/sliding_stage4_mechanism.json` is the frozen quick-gate calculation with two parameter sets; run with `--quick`. The full optional extended scan may be expensive and is not a frozen result.
- Planned literature/authorship status: no author contact yet for BSF extension; no experimental/DFT claims. Retain NMEE static vs driven-DC distinction.
