# v70 | Multi-frequency Lehmann–Keldysh ordered-limit audit

Date: 2026-10-22. Original Xiang PRB source audit: **PARTIAL**. The v70 observable is the **finite-device projected polarization current**, NOT the paper's bulk orbital magnetotransport tensor.

## Protocol

Explicit full-Hermitian finite lead Hamiltonians (v69) and independent semi-infinite Keldysh lesser-Green benchmark (v68). Same four-orbital device Hamiltonian, off-diagonal projected positions, magnetic frame, device electric coupling, electric current contact and device-only observable. Three frequencies: ω=0.20, 0.29, 0.40. At each frequency, evaluate lengths 250 and 350 with spectral η=0.02, 500 with η=0.01, and 1000 with η=0.005. The paired L/η path is **not** a mathematical extrapolation, and η is a transition-spectrum regulator, **not** a physical reservoir relaxation rate.

## Numerical evidence (model units)

| ω | Semi-infinite Keldysh σ | Finite-lead L=1000, η=.005 σ | Absolute error | Relative error |
| ---: | --- | --- | ---: | ---: |
| 0.20 | -0.002705541 + 0.001756620i | -0.002675252 + 0.001636795i | 1.236×10^-4 | 3.8314% |
| 0.29 | -0.005292650 - 0.001194257i | -0.005179062 - 0.001373655i | 2.123×10^-4 | 3.9135% |
| 0.40 | -0.006430101 - 0.005990459i | -0.006278300 - 0.006127154i | 2.043×10^-4 | 2.3245% |

For ω=0.20, the error **increases** from 7.77388×10^-4 at L=250 to 7.93198×10^-4 at L=350, before decreasing to 4.73456×10^-4 at L=500 and 1.23594×10^-4 at L=1000. This falsifies an assumption of pointwise monotonicity. Likewise the L=500, η=.01 boundary-cap switch at ω=.29 from 0 to ±1.1 shifts the complex σ by ≈2.18×10^-5. The cap effect is smaller than the L=500 continuum error (3.02×10^-4), but is **not zero**. Sign-related cap differences are below 2×10^-8 in this particular model. No boundary independence theorem is claimed.

### Interpretation

- **VERIFIED:** finite Hermitian Lehmann spectral results, independently constructed continuum Keldysh comparator, three-frequency convergence *trend*, strong finite-size resonance negative control and distal boundary-cap sensitivity.
- **CONDITIONAL:** joint L/η path, arbitrary supplied projected orbital dipoles, ideal leads, device-local E field and finite-device projected polarization current.
- **UNRESOLVED:** controlled η→0/L→∞ continuum extrapolation with systematic uncertainty; density/physical finite q→0 homogeneous bulk orbital response; source SM Eq20/31 sign mismatch, SM Eq49/52 factor-two source discrepancy, m≠l three-band omitted sectors, full magnetic Wannier embedding, published quantitative curve closure. Original manuscript status **PARTIAL**.

No raw author source figures were used to fit normalization. No original-paper COMPLETE claim is made.

## Progress assessment (planning estimates, not numerical physics claims)

| Workstream | Estimate | Completed | Actual blockers |
| --- | --- | --- | --- |
| GitHub engineering/reproducibility | ≈90–95% | Modular source, frozen artifacts, PR CI, merged-commit verification, versioned handoffs | v70 last remote CI; ongoing maintenance/release gates |
| Reusable finite-model physics capability | ≈75–80% | Peierls H, supplied off-diagonal P, seven magnetic/electric contacts, Ward, semi-infinite elastic Keldysh and independent finite-lead Lehmann | General material-ready Wannier field embedding, bulk observable API, robust limit-error bars |
| **Complete original Xiang PRB reproduction** | **≈30–40%** | Named source formula/model subproblems independently checked | Original PHE/MNHE source sign/factor and omitted three-band terms, full tensor/figure quantitative closure, source-identical relaxation prescription and actual bulk orbital current |

Percentages are rough project management guidance, **not** the percentage of equations proved and not inferred from version/test counts. The bottleneck is source-specific physical meaning, not adding more finite-device regressions.

## Next target v71

**Priority 1:** A source-equation-to-observable map distinguishing original PHE (Zeeman spin channel), orbital MNHE and the finite projected (i[H,P]) device current. Explicitly derive the limiting and contact terms, identify which cannot be identified without extra assumptions. Quantify L/η uncertainty at multiple paths and keep ω/q/L limits separate.

**Priority 2:** A concrete end-to-end material-input contract: full (H_{nm}(mathbf k)), velocity and spin, position/Wannier connection and magnetic positional shift with gauge-covariant derivatives; reject incomplete inputs instead of filling dipoles with orbital centers.

**Priority 3:** Preserve parallel literal source branches for disputed signs/coefficients, create independent three-band controls, then compare published panels only with SI and source normalization documented.

v70 artifacts: `scripts/reproduce_v70_ordered_limits.py`, `tests/test_ordered_limits_v70.py`, `results/frozen/v70_ordered_limits.json`.
