# v70 → v71 Handoff | 2026-10-22

## Scope and evidence

Original article: L. Xiang et al., *Quantum Geometric Map of Magnetotransport*, PRB 113 L201406 (2026). **Scientific status: PARTIAL**. This work is independent and neither author-endorsed nor a confirmation of the published bulk orbital MNHE.

v70 computes the same finite-device polarization-current observable in two physically distinct but comparable *implementations*: v69 explicit Hermitian device+finite leads Lehmann sum vs v68 continuum semi-infinite elastic reservoir Keldysh. Unlike v66 scalar-γ, neither method substitutes scalar phenomenological dephasing for the lead self-energy. The finite Lehmann response has spectral η>0 and requires coordinated L→∞/η→0. Do not treat η as physical reservoir scattering.

## Results at three frequencies

| ω | v68 semi-infinite σ | v70 L1000/η.005 σ | Relative difference |
|--:|---|---|---:|
| .20 | −.002705540717+.001756620087i | −.002675251809+.001636795325i | 3.8314% |
| .29 | −.005292649616−.001194257399i | −.005179061976−.001373654831i | 3.9135% |
| .40 | −.006430100592−.005990458960i | −.006278299919−.006127153776i | 2.3245% |

Four stages L=250(.02),350(.02),500(.01),1000(.005). At ω=.20 error *increases* from 7.77388e-4 to 7.93198e-4 for L250→350: nonmonotonic finite-size spectral resonance. Later stages decrease error; none has proved a full continuum extrapolation. At ω=.29 and L500, end_cap 0 vs ±1.1 changes σ by ≈2.18e-5; it is small relative to continuum finite-size residual but not zero. No particular extrapolated value or error bar is asserted.

## Scientific acceptance

**VERIFIED:** full-Hermitian finite-lead v69 Lehmann implementation, independent v68 elastic Keldysh comparator, three-frequency numerical checks, nonmonotonic and far-boundary negative controls, prior exact finite Lehmann Ward, static isothermal derivative and source-independent magnetic frame covariance.

**CONDITIONAL:** prescribed finite-length/spectral-resolution path; ideal noninteracting leads, externally supplied projected Pa/Pb; electric drive only on device; reported device polarization-current observable rather than bulk uniform-field orbital current.

**UNRESOLVED:** original supplemental sign differences (SM Eq20 vs Eq31 and main/SM helicity), orbital Eq49→52 factor-2 and omitted m≠l three-band terms, original Fig2(b) SI plateau discrepancy, complete magnetic Wannier position operator and response, thermodynamic finite-q / ω / disorder order of limits, and source-paper *physical* tensor/figures. **DO NOT promote original paper beyond PARTIAL** merely because of 926-ish passing integrated tests.

## Management progress (rough, non-physical estimates)

- GitHub + reproducible engineering: **90–95%** of target workflow. Multiple historical PRs and remote CI runs are available; v70 PR and postmerge CI should be recorded here once actually completed.
- Reusable finite-model physics capability: **75–80%**. Magnetically covariant finite Peierls device, supplied nonlocal orbital dipoles, all position/electric/magnetic contact derivatives, Ward, semi-infinite leads, two-energy Keldysh, independent Hermitian Lehmann all implemented with tests. Missing complete material-driven Wannier and bulk-observable API.
- Original Xiang article *complete source-level reproduction*: **30–40%** planning estimate. This is the bottleneck and NOT inferred from software volume; source equations, physical channels and quantitative published plots require separate closure.

## Next scientific work v71

1. **Source-observable audit before more code**: distinguish published PHE (spin Zeeman), orbital OHE and orbital MNHE from finite projected device J=i[H-EPb,Pa]; derive candidate finite-(q,ω) gauge-invariant bulk tensor and explain when device/transport can or cannot equal it.
2. **Source literal alternative branches**: SM Eq20/31 sign, Eq49/52 factor-2, independent 3-band offdiagonal interference; do not use fitted factor correction.
3. **Wannier input contract**: full magnetic field dependent orbital embedding (positions and intercell connections), spin and velocity, filling, correct time-reversal/point group and SI normalization. Fail closed for missing operator data.
4. **Numerical limit error bar**: multiple η(L) trajectories at three ω values, end caps and finite-size-spectral drift; only then consider controlled extrapolation.

v70 code `scripts/reproduce_v70_ordered_limits.py`, tests `tests/test_ordered_limits_v70.py` and `tests/test_ordered_limits_v70_frozen.py`, frozen file `results/frozen/v70_ordered_limits.json`; science report `docs/MULTIFREQUENCY_ORDERED_LIMITS_V70.md`. PR/CI SHA should be appended only after actual check.
