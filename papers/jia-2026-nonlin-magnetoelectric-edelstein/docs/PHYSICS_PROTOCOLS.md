# Physical protocol contract (frozen)

This file is a **definition/usage guide**, not a claim that the paper contains an error.

## Published intrinsic (PRL Eq. 1) versus full-static (SM Eq. 9)

For clean, nondegenerate Kane–Mele bands, the code evaluates two- and three-band pieces of the article's reported intrinsic expression. The supplemental density-matrix expression also contains a diagonal-density-matrix term. The stored algebraic identity is

```
SM_Eq9_full_static = PRL_Eq1_printed_two_plus_three + disputed_diagonal_term
```

Under the same sign, spin, lattice and k-integration convention, the 451×451 benchmark `Gamma_xz;x` is:

| Branch | Value | Units |
|---|---:|---|
| Printed Eq. (1) branch | -0.010739723280418214 | `(nm T V)^-1` |
| Full static Eq. (9) branch | -0.009548371489696384 | `(nm T V)^-1` |
| Difference (full − printed) | +0.0011913517907195352 | `(nm T V)^-1` |

The full-static coefficient is independently consistent with finite Zeeman field derivatives in insulating tests. The paper's nonequilibrium dc limit and equilibrium preparation are **not** proven identical; the discrepancy is retained as a scientific question. Do not merge these outputs or promote the diagonal term to a correction to the paper without the author-side response convention.

## CuMnAs intrinsic Fig. 1(c)

The bands are exactly double-degenerate due to PT. The naive *nondegenerate* band-sum is invalid. The code instead includes (i) finite-field static and (ii) finite-frequency PT Liouville and BGK models. The reported `fig1c-static-diagnostic` freezes the zero-field occupations in a static susceptibility derivative; it is NOT independently established to coincide with the author's finite-T driven-dc intrinsic coefficient. Near exact Dirac crossings this clean static limit can be logarithmically singular, so a fixed uniform mesh must not be declared converged without extra tests.

## CuMnAs extrinsic Fig. 1(d)

PRL Eq. (2) contains `f'_n`, a spin-space quantum metric, group velocity, and relaxation time `tau`. The PT-doublet `P_+, P_-` projector trace avoids choosing a basis within a degenerate pair. At fixed spectrum, the result scales as `tau`. At fixed dimensionless CuMnAs model it scales as inverse declared lattice length `a_nm`; both properties are explicitly tested. The numerical sign and exact absolute magnitude require a fixed magnetic perturbation, spin definition, and author length convention.

## Independent phenomenological BGK branch

`fig1-bgk` prescribes a physical sequence: establish thermal equilibrium at finite B, then drive uniform steady E with population relaxation `eta_p` and intergroup dephasing `eta_phi`. Its collision model conserves global charge for the defined case. It is a *comparison protocol*, not a microscopic impurity self-energy or the original article's formal intrinsic/extrinsic split.

## Units and axes

- 2D `Gamma` is magnetic moment density in units of `mu_B/nm²` per electric field in `V/nm` and magnetic field in `T`; hence numerically `1/(nm V T)` under the code's Pauli spin convention.
- Kane–Mele model uses explicit model reciprocal vectors; compare its physical lattice length to the paper before trusting experimental amplitude.
- CuMnAs Eq. (5) momentum is dimensionless. Physical conversion needs **explicit** `a_nm` and `t_eV` (`Model` defaults `a_nm=1`, `t=1 eV` purely as reference scales). Velocity already multiplies `a_nm`, while BZ area divides by `a_nm²`; Eq. (2) scales as `1/a_nm`.
- The finite-temperature BGK path preserves a legacy program normalization and is **not** directly comparable to the printed branch in physical units.
- The unit conversion is conditional on the Zeeman perturbation convention `-g mu_B B sigma` and on measuring spin in Pauli/`mu_B` units, not `hbar/2`.
