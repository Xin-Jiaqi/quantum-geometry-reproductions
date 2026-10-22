# v69 | Independent finite-lead Lehmann oracle and ordered-limit audit

Date: 2026-10-22. Parent source: v68 exact two-energy semi-infinite-lead Keldysh response. Overall Xiang PRB source reproduction remains **PARTIAL**.

## Why v69 is a genuinely different implementation

v68 obtains the finite-device response from energy-dependent semi-infinite lead self-energies and a two-energy lesser-Green-function integral. v69 instead constructs one explicit, finite, Hermitian Hamiltonian consisting of the same four-orbital Peierls device and two finite nearest-neighbor tight-binding chains, diagonalizes the **entire device+lead Hamiltonian**, and computes the response from an eigenstate-pair Lehmann sum. No surface self-energy or Keldysh integral is used in this response calculation.

The device-endpoint coupling vectors are the same physical vectors v_alpha for which v68 has Sigma_alpha^r=g_surface^r v_alpha v_alpha^dagger. The internal lead hopping is t=1.1; the device endpoint matrix element is v_alpha, **not** t v_alpha. The field acts only on the device through H(E,B)=H(B)-E P_b(B), while the reservoir chain Hamiltonians and chemical potential remain fixed. The observable intentionally matches v68's DEVICE projected polarization current J_a(E)=i[H_device(B)-E P_b(B),P_a(B)], embedded as zero on the explicit lead sites. It is neither terminal Landauer current nor yet the source paper's bulk orbital Hall tensor.

For finite frequency with e^(-i omega t) convention and finite spectral resolution eta>0, the independent density-response matrix in the eigenbasis of the full closed Hamiltonian is

\[
(\rho_E)_{mn}=\frac{(f_m-f_n)(P_b)_{mn}}{\omega+\epsilon_n-\epsilon_m+i\eta}.
\]

Then sigma_ab^(L,eta)=Tr(J_a^0 rho_E)+Tr(J_a^E rho_0), including J_a^E=-i[P_b,P_a]. This expression is obtained directly from the Lehmann sum and is not formed by reusing the v68 Green-function integral.

## Numerical evidence

Parameters: (mu,T,B,q/hbar,omega,t)=(0.13,0.37,0.23,-0.84,0.29,1.1), finite four-orbital device, two leads of length L, dimensionless hbar=1.

The v68 converged continuum reference is

\[
\sigma_\mathrm{Keldysh}=-0.005292649615898-0.001194257399095i,
\quad
\partial_B\sigma_\mathrm{Keldysh}=0.007070053900131-0.004146290848383i.
\]

| Explicit lead length L | Transition resolution eta | Finite-Hamiltonian Lehmann sigma | Absolute difference from v68 |
|---:|---:|---|---:|
| 250 | 0.020 | -0.004480632176 - 0.001641176558 i | 9.269e-4 |
| 350 | 0.020 | -0.004515423359 - 0.001567582191 i | 8.622e-4 |
| 500 | 0.010 | -0.005064862978 - 0.001391880742 i | 3.016e-4 |
| 1000 | 0.005 | -0.005179061976 - 0.001373654831 i | 2.123e-4 |

Jointly increasing L and decreasing eta improves agreement. The production point is within 2.5e-4 absolute, about 4% of the continuum response magnitude, but is **not exact convergence**. It is explicitly a bounded finite-L/eta agreement, not evidence that the bulk response or all limit orders are settled. Eta is a regulator for the discrete closed-system transition spectrum; it is not a reservoir relaxation rate. Finite-size resonances mean that L and eta must be scaled together, and eta must not be sent to zero at fixed L.

The five-point B derivative of the independent finite Hamiltonian at L=350, eta=0.02 is 0.007025618379405 - 0.004280408212469 i, within 1.5e-4 of the v68 mixed derivative. Three magnetic/orbital representations (symmetric, Landau and noncommuting orbital) give the same v69 result to 2e-11 at L=60, eta=0.03.

## New exact Ward and static result

The finite closed Hamiltonian obeys its own exact two-energy Ward identity

\[
(\omega+i\eta)\rho_E-[H,\rho_E]+[P_b,\rho_0]=0,
\]

with residual 7.8e-17. The total induced charge is zero; the device's induced charge -0.2134428574 - 0.3643610823 i is balanced by the explicit finite-lead charge +0.2134428574 + 0.3643610823 i. This is an independent finite-Hamiltonian identity, distinct from v68's open-reservoir collision/injection Ward.

At omega=0, the static **isothermal** response requires the diagonal -f'(epsilon_n)(P_b)_nn term in addition to off-diagonal Lehmann terms. Including it yields sigma_static=0.002996376693954; independently re-diagonalizing the complete finite Hamiltonian at E=+/-1e-4 gives 0.002996376693709, absolute error 4.2e-13. Omitting the diagonal occupation response silently interchanges the isolated dynamical and equilibrated static limits and fails this gate.

## Acceptance and boundaries

- **VERIFIED:** independent full-Hamiltonian Lehmann formula; exact closed-system Ward; explicit device/lead charge partition; isothermal static Frechet derivative versus full E re-diagonalization; frame covariance; five-point B derivative.
- **CONDITIONAL:** finite L and nonzero eta; supplied projected Px,Py; field applied to device only; noninteracting ideal chain leads; observable is device polarization current.
- **NOT YET CONVERGED:** production AC Lehmann vs v68 continuum remains a 2.1e-4 absolute discrepancy. A denser joint (L,eta) sequence and a principled extrapolation remain warranted.
- **UNRESOLVED:** real-material complete magnetic Wannier embedding, thermodynamic bulk homogeneous-field orbital current, relation of device-current observable to the source's transport tensor, and original source SM sign/factor/three-band questions. The original Xiang paper stays **PARTIAL**.

## Reproducibility

- Source: src/qgmap/finite_lead_lehmann_v69.py
- 45 physical gates: tests/test_finite_lead_lehmann_v69.py
- independent quick frozen-data regeneration: tests/test_finite_lead_lehmann_v69_frozen.py
- benchmark: scripts/reproduce_v69.py; default writes production L=1000, eta=0.005; --quick omits this long diagonalization.
- frozen numeric evidence: results/frozen/v69_finite_lead_lehmann.json
- run OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m unittest discover -s tests -p 'test_finite_lead_lehmann_v69*.py' -v.

The v69 science gates are not a claim that the source article is completely reproduced. The next step is tighter ordered-limit control and then an explicit derivation/implementation of the source's **bulk** orbital-current observable, not just further finite-device tests.
