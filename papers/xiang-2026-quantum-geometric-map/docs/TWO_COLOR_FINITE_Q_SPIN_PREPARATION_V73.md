# v73 | Independent two-color finite-q spin response and magnetic equilibrium preparation

Date 2026-10-22. Source Xiang et al., PRB 113 L201406 (2026). Original source scientific status **PARTIAL**. Primary-source conventions referenced from previous SHA256-pinned source audits v49/v53, not from a fresh reading of the publisher PDFs. This v73 work is an independently verified *method precursor*, **not original SM Eq15's selected electric 1/eta_a pole or source-identical DC PHE**.

## Precisely what is calculated

Starting from the v72 periodic 2-band Bloch model with spin-Zeeman coupling \(H(k,B_0)=H_0(k)+B_0\sigma_x\), introduce a scalar-potential \(\phi_q\) at finite commensurate q along y, \(E_y=-iq\phi_q\). The exact intermomentum density vertex is \(W=U(k+q)^\dagger U(k)\) in a momentum-dependent Bloch frame. The longitudinal Ward vertex is \(J_y=[H(k+q)-H(k)]/q\). Transverse \(J_x\) is a separately **declared symmetrized endpoint velocity vertex** and is NOT identified with the full transport/magnetization-subtracted author observable.

Two independent retarded frequencies are kept: \(w_E=\omega_E+i\eta_E\), \(w_B=\omega_B+i\eta_B\), output \(w_\Sigma=w_E+w_B\). For \(\mathcal L_w^{+,-}R=wR-H_+R+RH_-\) and \(F_\pm=f_T(H_\pm-\mu)\), compute

\[
R_E=\mathcal L_{w_E}^{-1}(WF_- - F_+W),\qquad R_{B,\pm}=\mathcal L_{w_B,\pm}^{-1}([Z_\pm,F_\pm]).
\]
\[
R_{B\rightarrow E}=\mathcal L_{w_\Sigma}^{-1}(W R_{B,-}-R_{B,+}W),\qquad
R_{E\rightarrow B}=\mathcal L_{w_\Sigma}^{-1}(Z_+R_E-R_EZ_-).
\]

**Both** orders contribute. This is a collisionless finite-q scalar-potential theory, not a substitute for the source's selected tau pole. Individual and total full response satisfy exact periodic charge Ward with \(w_\Sigma\); the source integrated trace cancels. Eigenstate-pair spectral response and independent dense 4×4 orbital Sylvester inversions match to machine precision, with full momentum-dependent U(k) frame covariance.

## Accepted numeric results

Parameters N=24, q=2π/24, \(\omega_E=.27,\eta_E=.045,\omega_B=.12,\eta_B=.03,B_0=.11,\mu=.32,T=.12,\beta=.24,\alpha=.13,m=.32,w=.35\), dimensionless charge=1, hbar=1. All values are \(\chi_{J_x,\rho B}/(-iq)\):

| Contribution | Complex response |
|---|---|
| B then E | +0.001208260956 + 0.009325386045i |
| E then B | −0.043217183392 − 0.019011725386i |
| **Full** | **−0.042008922436 − 0.009686339341i** |

Mixed finite-q charge Ward residual around \(10^{-18}\). Distinct eta_e, eta_b and w_sum strongly affect the result; deleting one order is demonstrably wrong for the specified collisionless protocol.

## Key missing population channel and physical preparation

**The collisionless time-dependent Zeeman linear response at omega_B=0 has no diagonal Fermi occupation derivative** \(f'(\varepsilon_n)Z_{nn}\). A magnetic field that was statically equilibrated *before* the electric perturbation does include this diagonal response, as part of the Fermi-matrix Fréchet derivative. v73 provides an **opt-in, explicitly separate prepared-static-population prescription**, ONLY at exactly omega_B=0. It adds the diagonal term to R_B; it does not pretend this follows from an author's collision integral.

Compare its limit eta_B→0 with the independent v72 analytic d_B of the fully equilibrium-prepared linear electric finite-q response. At fixed eta_E=.045:

| eta_B | Without population: absolute discrepancy | With equilibrium population: absolute discrepancy |
|---:|---:|---:|
| 1e-3 | 0.03169573 | 0.0002249495 |
| 1e-4 | 0.03155487 | 0.0000225133 |
| 1e-5 | 0.03154062 | **0.0000022515** |

Prepared response approaches the independent equilibrium finite-B derivative with first-order eta_B error; collisionless response does not. The fixed finite-q omega_E=0 scan is **not** sufficient to identify a nonzero author tau coefficient: q→0 and the omega/eta_B/eta_E limits have not been controlled as source prescribes.

## Scientific state and next task

- **VERIFIED:** two independent insertion orders, two complex poles, exact finite-q charge Ward, frame covariance, spectral pair-sum vs separate Sylvester, occupation-preparation-sensitive static-B limit, and explicit negative controls.
- **CONDITIONAL:** scalar potential at finite q, noninteracting 2-band toy Bloch model, chosen transverse current, retarded collisionless regulators, *optional* prepared equilibrium B population (not source bath).
- **UNRESOLVED:** original SM Eq15 selected 1/eta_E pole, source Eq20/31 sign and main/SM Eq10/62 helicity, source DC PHE tensor, physical orbital magnetic Wannier position/moment, Eq49/52 coefficient and omitted 3-band terms, proper transport/magnetization current, SI Fig2. Original paper remains **PARTIAL**.

Module src/qgmap/two_color_finite_q_spin_v73.py; 45 physical tests tests/test_two_color_finite_q_spin_v73.py; 1 full independent frozen regeneration tests/test_two_color_finite_q_spin_v73_frozen.py; generator scripts/reproduce_v73.py; numeric reference results/frozen/v73_two_color_finite_q_spin.json. 46/46 local checks PASS. Full GitHub remote CI and merged-main checks are independently required before release claims.
