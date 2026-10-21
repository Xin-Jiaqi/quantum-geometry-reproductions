# v43: Matched SCBA Hall surface vertex versus exact finite-disorder ensemble

**Scope: PARTIAL.** This is a *linear charge Hall* Středa **I(a) RA Fermi-surface channel**, not the complete Hall I+II conductivity and not Xiang et al., PRB 113 L201406 (2026)'s second-order mixed EB PHE. A green software CI is not acceptance of the complete paper.

## Exactly matched microscopic inputs

Both methods use the same two-orbital QWZ Hamiltonian, chemical potential μ, bare positive Green regulator η, charge current derived from H(k+A), and scalar Gaussian onsite disorder with **variance W**. Each exact real-space sample has Vᵣ distributed as Normal(0,√W), identical in the two onsite orbitals, frozen while boundary twist changes. The SCBA equation is

\[
G^R(E,k)=[E+i\eta-H(k)-\Sigma^R(E)]^{-1},\qquad
\Sigma^R(E)=W\langle G^R(E,k)\rangle_k,\quad G^A=(G^R)^\dagger.
\]

Using *the same* W (not a fitted relaxation time), solve the RA Bethe–Salpeter vertex:

\[
\Gamma_y^{RA}(k)=v_y(k)+C_y^{RA},
\qquad C_y^{RA}=W\langle G^R(v_y+C_y^{RA})G^A\rangle_k.
\]

The evaluated surface observable in e²/h is

\[
\sigma_H^{I(a)}=\frac12\operatorname{Re}\langle\operatorname{Tr}[
v_xG^R\Gamma_y^{RA}G^A-
v_yG^R\Gamma_x^{RA}G^A]\rangle_k.
\]

Only **one** current vertex is dressed: dressing *both* legs would generally double-count ladder diagrams. Independently using the AR ladder on the left current leg reproduces the same result. The SCBA RA density resolvent Ward condition is checked independently:

\[
\Gamma_0^{RA}=I+W\langle G^R\Gamma_0^{RA}G^A\rangle
=I-\frac{\Sigma^R-\Sigma^A}{2i\eta}.
\]

This is **not** by itself a proof of physical finite-scattering-rate continuity or the full Hall response.

## Independent exact reference

For each frozen realization and twist, compute the *dense inverse-matrix real-space* Středa I(a) trace from the v42 exact finite-torus module, with the same μ and η. Twist average inside each sample; average **independent disorder seeds** separately. The standard error is calculated across **seeds**, not from pretending that twist points are independent random samples.

**Frozen primary protocol:** m=−1, μ=1.1, η=.12, SCBA mesh 35×35, exact L=4, 7×7 twists, 24 seeds (1947–1970).

| Variance W | SCBA bare I(a) | SCBA one ladder I(a) | Exact finite-sample mean I(a) | Seed SEM |
|---:|---:|---:|---:|---:|
| 0 | −2.011033 | −2.011033 | −2.008081 | 0 |
| .03 | −1.736138 | −1.810174 | −1.854201 | .043972 |
| .12 | −1.278070 | −1.430324 | −1.500275 | .055025 |

In these specific finite samples a single ladder improves agreement over bare SCBA, but **does not close it**. At W=.12 the residual is roughly +0.070 e²/h or 1.27 sample SEM; this SEM is not an estimate of systematic finite-size or missing-diagram error.

**Further NONacceptance scaling diagnostic**, W=.12, μ=1.1, η=.12, SCBA 49×49: exact L=3 (96 seeds, 9×9 twists) gives −1.56588±.03112, L=4 (96 seeds, 9×9 twists) −1.52948±.02251, L=5 (64 seeds, 7×7 twists) −1.54349±.02190, versus ladder −1.43032. The visible mismatch has not been shown to disappear with L. Neglected crossed diagrams, finite-size effects, and regulator/DC protocol remain possible explanations, **not established diagnoses**.

## Checks and use

11 independent new test cases verify the earlier v40 frozen RA coefficients, left/right BSE equality, a separately iterated Neumann ladder versus the four-Pauli solver, W=0 normalization against twist-averaged exact clean reference, direct v42 current trace, weak/moderate disorder sampling, mesh offset, RA Ward/BSE residuals, proper seed-based uncertainty, and invalid-input rejection.

To reproduce from the PRB paper directory:

    python -m pip install -e .
    python -m unittest tests.test_scba_surface_v43 -v
    python scripts/reproduce_scba_surface_v43.py --W 0 --seeds 24
    python scripts/reproduce_scba_surface_v43.py --W .03 --seeds 24
    python scripts/reproduce_scba_surface_v43.py --W .12 --seeds 24

## What v43 does NOT establish

It does not implement a **conserving disorder-averaged sea II**, the necessary electromagnetic derivatives and contact terms, or cross-diagram corrections. Combining v39's *bare-current energy-dependent sea* with this I(a) result is **not** yet justified as a complete Středa response. It does not address a controlled thermodynamic and ω,q,η order of limits, the paper's original second-order EB response, or first-principles material operators.

**Next v44:** derive and independently test the matching sea II and an alternative magnetic-field Středa density route, then compare **the full averaged** charge Hall result against v42 exact I+II with systematic seed/size/twist/convergence analysis. Original Fig. 2(b) approx. 2π amplitude issue stays OPEN.
