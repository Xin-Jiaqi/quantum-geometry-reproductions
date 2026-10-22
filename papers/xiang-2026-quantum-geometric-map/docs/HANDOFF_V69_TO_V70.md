# v69 -> v70 Research Handoff | 2026-10-22

Repository: Xin-Jiaqi/quantum-geometry-reproductions. Latest baseline at branch start: main 1662b926a4552aa479e9a53ef63a73a573d049d2. Xiang PRB source-paper status remains **PARTIAL**.

## What v69 adds

An independent finite-Hamiltonian Lehmann implementation. Build the same four-orbital Peierls device plus two explicit finite tight-binding chains and diagonalize the **entire Hermitian device+lead Hamiltonian**. Compute the AC susceptibility by an eigenstate-pair sum; do not call the semi-infinite surface self-energy or v68 Keldysh integral to construct this response. Use the same device-only projected polarization current J_a(E)=i[H_device-E P_b,P_a] as v68, including J_a^E=-i[P_b,P_a]. Lead coupling vectors match the v68 self-energy convention Sigma^r=g_surface v v†; device-to-lead matrix elements are v, not t*v.

## Key numeric acceptance

At mu=.13,T=.37,B=.23,q/hbar=-.84,omega=.29,tlead=1.1:
- v68 continuum sigma = -0.005292649615898 - 0.001194257399095 i.
- Explicit finite-lead Lehmann, L=1000, eta=.005: -0.005179061975941 - 0.001373654830605 i, absolute difference 2.1233e-4 (~4% of continuum response magnitude).
- Joint ordered-limit path (L,eta)=(250,.02),(350,.02),(500,.01),(1000,.005) reduces discrepancy 9.269e-4 -> 8.622e-4 -> 3.016e-4 -> 2.123e-4. **This is convergence toward the continuum, not yet zero residual.**
- v68 dB sigma = +0.007070053900131 - 0.004146290848383 i. Independent five-point derivative of the finite-lead Lehmann at L=350, eta=.02 is +0.007025618379405 - 0.004280408212469 i (difference ~1.4e-4).
- Finite-system two-energy Ward residual 7.8e-17; total AC charge exactly zero and device charge is balanced by opposite explicit lead charge.
- Static isothermal Lehmann including the diagonal -f'(epsilon_n) P_nn term gives 0.002996376693954; independent full E re-diagonalization FD gives 0.002996376693709 (absolute error 4.2e-13).
- 45 new physical tests + one quick frozen numeric regeneration test. v60-v68 baseline 306 plus v69 46 = **352 local tests**. Existing tests were executed in bounded module batches because the complete suite can exceed the current local command time window; GitHub full CI remains the release gate.

## Project progress assessment (management estimate, not a physics metric)

1. **GitHub/reproducible engineering pipeline: ~90% of this project's current target.** The v60-v68 sequence has repeated PR + remote CI releases; v69 is the next scoped PR. Remaining work is this PR's final remote acceptance, release log and handoff archive. This measures pipeline maturity, not source-paper correctness.
2. **Finite-model physical capability integrated into the repo: ~70-80%.** The repo now has finite graph, magnetic embedding, supplied projected position contacts, scalar-gamma AC, semi-infinite elastic-lead scattering, two-energy Keldysh response, finite-Hamiltonian Lehmann, Ward/gauge/finite-difference checks. It remains a controlled model stack with externally supplied orbital dipoles, not a material-ready complete Wannier workflow.
3. **Complete original Xiang PRB reproduction: ~30-40%.** This estimate is deliberately much lower. v69 provides an independent method and a controlled ~4% finite-L/eta comparison, but does not derive the original bulk transport tensor from the source's equations. Source SM sign/factor/three-band issues, physical orbital-current observable, real magnetic Wannier embedding, and bulk/limit-order normalization remain unresolved.

Percentages above are subjective project planning estimates. Do not equate code/test volume with fraction of paper physics reproduced.

## Scientific boundaries

- **VERIFIED:** finite closed-Hamiltonian Lehmann identity, finite-system Ward and charge partition, static isothermal derivative versus full E re-diagonalization, frame covariance and independent finite-B derivative.
- **CONDITIONAL:** finite L, spectral eta, device-local E field, supplied Pa/Pb, ideal noninteracting leads and device polarization current.
- **NOT YET CONVERGED:** AC continuum limit has a residual 2.1e-4. Increase L and decrease eta jointly; avoid eta->0 at fixed L. Add fit/extrapolation and stability/error bars.
- **UNRESOLVED:** source paper's intrinsic bulk orbital E/B response, full magnetic Wannier embedding, the exact source SM sign/factor/three-band questions.

## v70 priority

1. Continue the controlled joint L/eta limit and quantify extrapolation uncertainty, including different end caps and frequency values.
2. Derive the exact source observable (bulk orbital-current / homogeneous-field tensor) from the paper equations and show how it maps, or fails to map, to the finite device polarization-current observable.
3. Only then implement the source-specific Kubo/Keldysh expression, with separate unit tests for disputed signs, factors and three-band channels. Do not fit corrections to force agreement.

## Source paths

- src/qgmap/finite_lead_lehmann_v69.py
- tests/test_finite_lead_lehmann_v69.py (45 gates)
- tests/test_finite_lead_lehmann_v69_frozen.py (quick frozen numeric gate)
- scripts/reproduce_v69.py
- results/frozen/v69_finite_lead_lehmann.json
- docs/FINITE_LEAD_LEHMANN_ORACLE_V69.md

## GitHub release verification (2026-10-22)

- [PR #62](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/62) merged by squash; final PR HEAD `40c752f36b2cbfcd86105592a9e4f3a42a8dcdeb`; merge commit `be8bb8344b0c1149b7f2a9944c88202ba19f649e`.
- Final PR-head [Actions 39722605718](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39722605718): **890/890 PASS**, `Ran 890 tests in 89.488s`, completed/success; overall source-paper status **PARTIAL**.
- Independent post-merge main [Actions 39722881262](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39722881262): **890/890 PASS**, `Ran 890 tests in 87.704s`, completed/success on exact merge commit `be8bb8344b0c1149b7f2a9944c88202ba19f649e`.
- Clean-extraction archive: 83 SHA256 entries, all matched; v69 module/tests **46/46 PASS** after fresh extraction. v60-v68 baseline **306/306 PASS** in bounded module batches; cumulative local result **352/352**.
- Archive is an independently tested standalone snapshot. Exact byte identity with every remote Git blob is not asserted unless each blob hash is explicitly matched.

Original paper status remains **PARTIAL** until source-equation sign/factor/three-band questions, full magnetic Wannier embedding and bulk orbital-current correspondence are resolved.
