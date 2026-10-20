# Stage 7e: exact weak-hopping quantum-geometry coefficients and two distinct velocity mechanisms

**Only for the frozen synthetic gapped eight-band sliding Kane–Mele model, zero-temperature equilibrium STATIC finite-Zeeman response.** This is not the disputed driven-DC NMEE, DFT, a stable ferroelectric, or a universal material sign claim. Prior Stage 1–7d code and results are preserved.

## A physical result beyond a numerical sign decomposition

In the declared A/B orbital embedding, write \(H=H_0+sV_s+cV_c\) and \(v_x=\partial_{k_x}H=v_0+sW_s+cW_c\). Total-layer sign parity enforces even total hopping order. A gapped, analytic expansion gives
\[
\Gamma_{xz;x}(u;s,c)=A(u)s^2+B(u)c^2+C(u)sc+O(t_\perp^4).
\]

**NEW:** The coefficients \(A,B,C\) are now evaluated **directly at zero interlayer hopping** from exact occupied-projector derivatives, not obtained by extrapolating finite hopping signals. Moreover, \(C=C^{[v_0]}+C^{[W_s]}+C^{[W_c]}\) identifies *two different physical routes*:
- **Positive \(B\)**: the two cross-sublattice hoppings \(V_cV_c\) dress Bloch eigenstates; the original intralayer \(v_0\) reads out almost all of the resulting mixed Berry curvature.
- **Negative \(C\)**: the explicit interlayer **cross-sublattice velocity** \(W_c\) interferes with the same-sublattice hopping \(V_s\); band dressing from \(v_0\) reinforces that negative term, while \(W_s\) partially cancels it.

At \(u=0.4\), \(K/K'\) radius \(0.75\,\mathrm{nm}^{-1}\), 26 radial Gauss × 72 angular per valley:

| Contribution | B (AB/BA-only) | C (same/cross interference) |
|:--|--:|--:|
| Intralayer velocity \(v_0\) after state dressing | +0.168010297235 | -0.036101925385 |
| Same-sublattice explicit velocity \(W_s\) | — | +0.009206881356 |
| Cross-sublattice explicit velocity \(W_c\) | +0.000222532524 (both equivalent copies) | -0.064657013479 |
| **Total** | **+0.168232829759** | **-0.091552057507** |

Thus 99.8677% of \(B\) arises from the \(v_0\) path, whereas the \(W_c\) contribution amounts to 70.6232% of the magnitude of \(C\). For \(A\), the \(v_0\) piece is numerically zero at second order, leaving a very small positive same-hopping velocity contribution. These partitions are **defined in the fixed original orbital embedding** and are not individually arbitrary \(k\)-dependent Wannier-gauge invariant; the complete coefficients and response are.

## Exact derivation and how to reproduce it

For a source \(\lambda S_x\) and a local \(q_x\) perturbation,
\[
\chi_{x;x}=i\operatorname{Tr}\{P[\partial_\lambda P,\partial_{q_x}P]\},
\quad \Gamma_{xz;x}=\int (d^2k/(2\pi)^2)\partial_{B_z}\chi_{x;x}.
\]
Use five *formal independent* source variables \(\{s,c,B_z,\lambda,q_x\}\), with \(sq_x W_s\) and \(cq_xW_c\) explicitly included in the local Hamiltonian. The 32 projector Taylor coefficients \(P_M\) are generated recursively by \(P^2=P\) and \([H,P]=0\). For every nonempty multi-index \(M\), let
\[
T_M=\sum_{\varnothing\neq L\subsetneq M}P_LP_{M\setminus L},\quad
R_M=-\sum_{\varnothing\neq L\subseteq M}[H_L,P_{M\setminus L}].
\]
In the eigenbasis of \(H_0\), the occupied–occupied block of \(P_M\) equals \(-T_M\), the empty–empty block equals \(+T_M\), and the cross-gap matrix elements obey
\[
(P_M)_{ve}=(R_M)_{ve}/(\epsilon_v-\epsilon_e).
\]
**Only denominators across the occupied/empty gap** occur, avoiding naked nondegenerate denominators inside either subspace. Extract \(s c B_z\) for \(C\), or use two independent bookkeeping copies \(s_1,s_2\) / \(c_1,c_2\) and divide by 2 for \(A/B\). This is an **exact algorithm for weak-hopping coefficients**, not an elementary 2-band closed formula or sign theorem.

| Sliding u | A | B | C |
|---:|---:|---:|---:|
| 0.0 | ≈0 | ≈0 | ≈0 |
| 0.1 | +8.5293213e-6 | +0.00314181250 | -0.0299881248 |
| 0.2 | +1.6663865e-5 | +0.0239747881 | -0.0568899842 |
| 0.4 | +3.0249521e-5 | +0.168232830 | -0.0915520575 |

The weak-coupling upper positive root of \(A+C\rho+B\rho^2=0\) at \(u=.4\), \(\rho=c/s\), is **0.5438679** (the second small root is 0.0003306). Stage 7d's independently varied channel ratio approaches the same root as \(s\to0\).

## Independent validations and limitations

1. Direct full-eight-band finite-hopping curvature at \(s=0.625\) meV agrees with zero-hopping projector derivatives to ~0.0107% (A), 0.0198% (B), 0.0171% (C) at \(u=.4\); residuals decrease with \(s^2\).
2. Random **k-independent** eight-dimensional unitary changes of every Hamiltonian and response vertex leave the local coefficient invariant to about \(8.4\times10^{-12}\) in program units.
3. The \(u=0\) symmetry point has vanishing \(A,B,C\), and both valleys contribute with the **same sign** at nonzero slide.
4. Denser 40×96 valley quadrature at \(u=.4\) reproduces the 26×72 numerical coefficients. This is convergence testing, **not** a mathematical full-BZ error bound.

The velocity path decomposition depends on the declared orbital embedding; a \(k\)-dependent basis requires Berry-connection-covariant velocities. The sign mechanism is toy-model specific and does not resolve the static vs driven-DC PRL discrepancy. All remote virtual states remain in the eight-band projector; a response-preserving four-band SW reduction has **not** been constructed.

**Next physics question:** resolve the projector derivative into source-resolvent loops \(\{S_x,S_z,V_s,V_c,W_c,v_0\}\) and determine the smallest explicit valley matrix-element products responsible for the signs, without silently discarding remote occupied or empty states.

**Files:** \`bsf_nmee/projector_jet_response.py\`, \`scripts/run_stage7e_analytic_coefficients.py\`, \`tests/test_stage7e_projector_jet.py\`. Complete standalone independent code and data are in the separately delivered Stage 7e archive.
