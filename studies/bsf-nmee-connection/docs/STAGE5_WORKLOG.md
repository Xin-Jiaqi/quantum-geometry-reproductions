# Stage-5 worklog

- Audited all stage-4 physical claims without changing frozen PRL/stages-1–4 source.
- Defined two independent hopping energy parameters `s,c`, and tested all four sign choices. Layer-sign covariance and mixed derivative agree to numerical precision.
- Extracted `C_sc=lim_{t→0}(Gamma(t,t)-Gamma(t,-t))/(2t²)` from three small t values, as a controlled mixed Taylor derivative.
- **Scientific correction:** negative `C_sc` at u=.10, .20, .40; the finite-coupling stage-4 inclusion-exclusion `Gamma_mix` contains more than odd-channel interference. Total sign reversal arises from competition with an increasing positive even-channel response.
- At fixed u=.20, found a tunnel-strength-induced sign crossing near 19.63385 meV, without sampled Fermi-gap closure; tested a finer valley quadrature and smaller Zeeman step.
- Restored untouched stage-4 `response_mechanism.py`: early development had added optional channel parameters there, which invalidated the stage-4 frozen manifest test. Created independent `channel_perturbation.py` instead and verified it matches stage-4 kernel at the default amplitudes to ~1e-18 in response units.
- Coarse circular-patch BZ grid area normalization is approximate; documented measured quadrature area deficit instead of silently declaring it exact.
- No actual DFT/Wannier calculation, ferroic stability or author-verified driven DC response is asserted.
