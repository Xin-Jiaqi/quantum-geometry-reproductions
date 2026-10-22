# Round 15: microscopic antiunitary finite-frequency Zeeman response

**Status:** verified model-method subset for three microscopic antiunitaries,
NOT a complete 31-MPG dynamical classification; full Cui et al. paper PARTIAL.

Source: Cui et al., arXiv:2604.09725v3, Eq.(16)-(19), Table I and Appendix B.
This is an independent synthetic Bloch-model validation, no author material data.

For A=U K with k'=Rk and U sigma_b* U^dag=lambda_b sigma_b,
the true-spin Zeeman tensor T_nm^{ab}=r_nm^a sigma_mn^b obeys:
  T_nm^{ab}(Rk)=-R_ac lambda_b [T_nm^{cb}(k)]*.

With the existing causal_response.py conventions (H_B=-kappa B.sigma,
H_E=+E.r, z=omega_drive+i eta), the exact interband frequency identities are:
  chi_JB(Rk,z)=R chi_JB(k,-z) diag(lambda)
  alpha_ME(Rk,z)=-diag(lambda) alpha_ME(k,-z) R^T.

This result uses complete complex z. For positive eta the -z argument is
**advanced**; replacing it with -omega+i eta is incorrect. Negative-control
differences for the three models are O(0.009-0.063), not roundoff errors.

| Group | k mapping | spin lambda (x,y,z) | CLEAN sector rule |
| --- | --- | --- | --- |
| 1' (TR) | -I | (-,-,-) | All paired odd-frequency J/B terms cancel |
| 2' (C2zT) | I | (+,+,-) | Bx,By: even; Bz: odd |
| m' (mxT) | diag(+,-) | (-,+,+) | xx/yy in-plane odd; xy/yx even |

The rules are verified by (a) analytic spinor-transformed H and exact dH/dk,
(b) spectral pairwise QGT transform, (c) full complex-frequency band Kubo,
(d) independent orbital-basis Liouville, (e) four-band spin x orbital
hybridization and (f) full 2D BZ cancellations/nonzero allowed entries.
The local complex-frequency maximum residual is O(1e-15).

Selected clean BZ (24x24) Im chi_yz^(odd) for C2zT is -0.00105258575;
Im chi_yy^(odd) for m_xT is 0.00062927689 (natural units). Unweighted TR
odd-frequency BZ response cancels to O(1e-18) already at 16x16. These are
illustrative model numbers and are NOT precision-converged material predictions.

Eq.(16) source Fourier frequency = -omega_drive; its O(omega_source)
OmegaA/gN and O(omega_source^2) OmegaN/gA terms obey these exact
microscopic parity constraints. The 31-group static Jahn projector of Round14
is still CONDITIONAL for other groups; the paper's six Table I discrepancies
are retained as open convention/observable questions, NOT adjudicated as
author mistakes.

### Run

    cd papers/cui-2026-curvature-flux-zeeman-qgt
    python -m pip install -e .
    python -m unittest discover -s tests -v
    python scripts/audit_antiunitary.py

20 new tests and 80 preceding ones = 100 tests. AutoWann operators/real Wannier
spins, collision vertices, non-Abelian Kramers subspaces, SI, all 31 direct
microscopic models and macroscopic observable KME remain **UNVERIFIED**.
