> **RETRACTED NUMERICAL RESULTS (v47):** The originally released v46 full equilibrium y-electric source had the wrong sign for the sigma_y derivative. All old v46 frozen mixed-field numbers and any comparison of the xy-only response with Eq.(8) are invalid. The code has been corrected; use the fully symmetrized tensor and independent source/finite-B checks in [v47](ORIGINAL_PHE_SYMMETRIC_SOURCE_CORRECTION_V47.md). The text below is preserved solely as a historical protocol record, not accepted numerical physics.

# v46 | Same tilted Dirac Hamiltonian, distinct mixed-field RTA collision targets

**Scientific status: PARTIAL.** This independent protocol test addresses the *original symmetric spin-induced PHE* of Xiang et al., PRB **113**, L201406 (2026), main Eqs. (8), (10), (11) ([arXiv:2510.02661v2](https://arxiv.org/pdf/2510.02661v2)). It is not the v41–v44 antisymmetric QWZ charge Hall extension. Neither v46 collision target is claimed to reproduce the complete Supplemental Material relaxation prescription without line-by-line checking.

## Fixed-tilt spin Zeeman shift and selected RTA definition

In model units vF=e=hbar=1, fixed beta=t_x/vF and b=g*mu_B*B_x/vF,

\[
H_B(k)=\beta k_x I+k_x\sigma_y-(k_y+b)\sigma_x = H_0(k_x,k_y+b),\quad M=\partial_b H=-\sigma_x.
\]

In the exact infinite-plane model, spin B_x can be absorbed by shifting k_y. The *fully B-equilibrated stationary RTA* is explicitly defined by

\[
\rho_{\mathrm{eq},B}=f(H_B),\quad \mathcal L_B X=\eta X+i[H_B,X],\quad
\rho_E(B)=\mathcal L_B^{-1}(\partial_{k_y}\rho_{\mathrm{eq},B}).
\]

The source derivative uses the **complete spectral projector/Fermi divided difference** rather than a selected band-population or off-diagonal term. Differentiating gives two necessary terms:

\[
\left.\partial_b\rho_E\right|_0=
\mathcal L_0^{-1}[\left.\partial_b\partial_{k_y}f(H_B)\right|_0]
+\mathcal L_0^{-1}[-i[M,\rho_E(0)]].
\]

The first term is the magnetically prepared equilibrium/electric source; the second is B acting on the E-induced density. The exact momentum translation demands that the **sum**, integrated on a complete continuum momentum plane, vanish. The current operator is v_x=beta I+sigma_y and has no explicit field derivative in this ideal model.

At beta=0.125, mu=0.30, T=0.055, eta=0.065, cutoff=2.5, radial 96 and angular 144 quadrature, centered b step 3e-5:

| Prepared target (shared model response units) | Value |
|---|---:|
| B adjustment of the equilibrium E source | +0.0762752980553 |
| B acting on pre-existing E density | -0.0762752975979 |
| **Complete mixed derivative** | **+4.57e-10** |
| Independent finite-B difference of full stationary response | -6.37e-11 |
| Difference of two algorithms | 5.21e-10 |

The analogous yx current and symmetric/antisymmetric xy combinations are separately controlled. The near-zero **full** response is a consequence of the protocol and fixed-tilt translation, **not a proof that the published PHE is absent in real TI materials**.

## Alternate *unprepared* magnetic drive (zero-field RTA target frozen)

To expose the hidden collision assumption, keep the relaxation target f(H0) even after B is applied. Both external-field orderings are then explicitly calculated:

\[
\rho_B=\mathcal L_0^{-1}(-i[M,f(H_0)]),\qquad
\rho_E=\mathcal L_0^{-1}(\partial_{k_y}f(H_0)),
\]
\[
\rho_{EB}=\mathcal L_0^{-1}(\underbrace{\partial_{k_y}\rho_B}_{B\,\mathrm{then}\,E}
+\underbrace{-i[M,\rho_E]}_{E\,\mathrm{then}\,B}).
\]

In general rho_B is **not** the derivative of f(H_B): the magnetic field is not allowed to prepare the field-dependent equilibrium occupations. Therefore this well-defined alternative scalar RTA need not satisfy the prepared-response cancellation.

With the same parameters: B-then-E **+0.1100503039423**, E-then-B **-0.0762752975979**, total **+0.0337750063445**. The E-then-B insertion is exactly shared with the prepared calculation, so **the difference comes from magnetic preparation of the electric source**.

The dimensionless low-eta coefficient eta*total is positive and near the formal original Eq.(8) angular coefficient at beta=1/8, but **not equal**; this discrepancy must not be adjusted by fitting a factor. For T=.02, eta=.0008, beta=.125: unprepared coefficient 0.00248344390, v45 printed Eq.(8) angular reference 0.00250641437 (ratio 0.990835). For beta=.6 the ratio drops to **0.764435**; it is an actual protocol-dependent mismatch, not a rounding artifact. Original v45 printed Eq.(8)/(11) sign issue and v37 Fig.2(b) approximately 2*pi normalization mismatch remain OPEN.

## Provenance and scientific gate

- **Author-derived:** The main text defines a symmetric tau-linear Zeeman PHE (Eq.8), a tilted surface Dirac cone (Eq.10), and a closed step coefficient (Eq.11). The final discussion explicitly calls tau phenomenological and treats microscopic bath coupling as outside paper scope.
- **Independently constructed in v46:** both particular scalar RTA targets, field-order response equations, complete spectral source and commutators, full-B finite difference, and tilt/broadening scans. This does **not** replace or correct the authors' Supplemental Material.
- **Scientific acceptance:** named method-gate PASS only; the author's *exact same* RTA denominator, magnetic equilibrium preparation, order of limits, signed Zeeman coupling and normalization require direct independent agreement. Overall paper PARTIAL.

### Reproduce

```bash
python -m unittest tests.test_dirac_mixed_rta_v46 -v
python scripts/reproduce_dirac_mixed_rta_v46.py --eta .065 --beta .125 --nr 96 --nphi 144 --h 0.00003
```

**Next v47:** verify the *actual* author supplemental collision/resolvent source (do not infer it from general physics or an unrelated RTA), run it unchanged as third protocol against these two limits, then re-evaluate printed Eq.11 and Fig.2(b) without fitting factors.
