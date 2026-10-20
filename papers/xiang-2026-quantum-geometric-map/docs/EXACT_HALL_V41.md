# v41: Exact finite-disorder Hall reference before conserving Středa

**Scientific status: PARTIAL for Xiang et al., PRB 113 L201406.** This is an independent QWZ **linear transverse Hall** gate, NOT the original nonlinear second-order `EB` PHE. It does not establish a conserving SCBA or separate physical Kubo–Středa surface/sea coefficients.

## Observable and response protocol

Use `e=ħ=a=1`, `T=0`, a two-orbital `L×L` torus, chemical potential within a finite sample spectral gap, and the **same frozen** scalar onsite random potential (variance `W` per site, common to both orbitals) at every twist.

```math
H_{QWZ}(k)=\sin k_x\sigma_x+\sin k_y\sigma_y+
(m+\cos k_x+\cos k_y)\sigma_z.
```

Set `A=θ/L`, and construct `H(θ)=F diag_k[H_QWZ(k+θ/L)] F†+V` in a **fixed** zero-twist real-site Fourier basis, with `J_a = ∂H/∂A_a`. Thus velocities and the disorder Hamiltonian share exactly the same microscopic protocol. The exact finite-system antisymmetric transverse response is

```math
\frac{\sigma^H_{xy}(z,\theta)}{e^2/h}
=\frac{2\pi}{L^2}\sum_{n\in occ,m\in empty}
\frac{2\,\mathrm{Im}[(J_x)_{nm}(J_y)_{mn}]}{(E_m-E_n)^2-z^2},
\qquad z=\omega+i\delta.
```

For `z=0`, this is a **full finite-system single-particle Lehmann Hall conductivity**, not a bare or ladder-dressed RA bubble. Independent `hall_retarded_crosscheck` instead computes ordered spectral retarded correlators `χ_xy^R(z)` and `χ_yx^R(z)` and reconstructs `σ_H=−(χ_xy−χ_yx)/(2iz)` with the v38 QWZ sign convention. Here the mixed diamagnetic derivative `∂²H/∂A_x∂A_y` vanishes **identically for this QWZ model** and frozen onsite disorder. It is additionally verified by a four-point numerical derivative. Other models and nonlinear `EB` do NOT inherit zero contact automatically.

## Independent topology check

An occupied-subspace FHS determinant-link Chern invariant is built from all occupied eigenstates on a `Nθ×Nθ` twist grid. It does not use current matrices or Kubo denominators. Large gauge transformations across the `2π` seams are explicitly accounted for by `exp(−2πi r_a/L)`. The independently computed twist-average of the exact Kubo response is checked against this Chern integer. Changing occupied rank is rejected rather than silently treated as topological quantization.

### Frozen local benchmark

`L=4`, `Nθ=7`, `μ=0`, `seed=1947`, `δ=0`:

| m | W | Hall at zero twist (e²/h) | Twist-average (e²/h) | FHS C | Min sampled gap |
|---:|---:|---:|---:|---:|---:|
| −1 | 0 | −1.28238202 | −1.00000002 | −1 | 2.00000 |
| −1 | 0.12 | −1.28761993 | −1.00000002 | −1 | 1.44695 |
| +1 | 0.12 | +1.28761993 | +1.00000002 | +1 | 1.44695 |
| −3 | 0.12 | +0.11908961 | +0.00000001 | 0 | 1.83456 |

At small finite `L`, **fixed-twist** Hall need not be quantized. Twist-averaging recovers the integer in these gapped samples. All numbers are exact finite-sample *model* observables, not a thermodynamic disorder average. The response broadening `δ` is not a phenomenological elastic-scattering time or v39/v40 single-particle `η`.

## Run and science gates

```bash
python -m pip install -e .
python -m unittest tests.test_exact_hall_v41 -v
python scripts/reproduce_exact_hall_v41.py --L 4 --mass -1 --W .12 --seed 1947 --ntheta 7
```

Eight local tests check the four clean Chern phases, clean/disordered Hall/FHS agreement, two finite-difference current derivatives, independent retarded spectral Kubo, mixed current contact, regulator limit and invalid inputs.

## Still OPEN

1. A **single** disorder-averaged, conserving Kubo–Bastin/Středa surface + sea decomposition with matched energy-dependent self-energy, current/density vertices, and appropriate contact/gauge identities. This v41 exact sample handles scattering nonperturbatively, but does NOT separately identify Středa I/II or diagrammatic ladder/crossed contributions.
2. Finite-size/seed-ensemble scaling, metallic Fermi energy, response frequency/bath broadening and thermodynamic/DC order of limits.
3. Original PRB nonlinear `EB` PHE, general multiband orbital terms, and the Fig.2(b) roughly `2π` normalization discrepancy.

**Decision:** v41 establishes a meaningful exact linear Hall **reference target** for subsequent disorder-averaged Středa development. Do not advance original-paper scientific acceptance beyond `PARTIAL`.
