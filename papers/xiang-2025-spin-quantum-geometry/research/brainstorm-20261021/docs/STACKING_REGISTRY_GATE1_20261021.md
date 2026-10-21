# Gate 1: registry-built PT bilayer and negative controls (2026-10-21)

**Scientific status:** SYNTHETIC MODEL GATE PASSED; INTERLAYER-NOVELTY HYPOTHESIS **NOT SUPPORTED** BY THIS TOY MODEL.  
**Scope:** independent B+D physics exploration. **Not** the original Xiang/Jin/Wang PRL reproduction, not a demonstrated sliding ferroelectric, not a first-principles material. Original reproduction's frozen acceptance state untouched.

Executable: [`../scripts/stacking_registry_gate1.py`](../scripts/stacking_registry_gate1.py) (requires NumPy). This file is part of `research/`, not paper CI.

## 1. Idea versus actual test

**Idea:** derive two magnetically related stacking structures from real-space registry, then ask whether a layer-staggered spin transport response can be switched through a genuinely interlayer quantum-geometric mechanism.

**What this first test really checks:** physical-distance **spin-independent** interlayer hoppings, exact PT and C2z relations, rank-two projectors, a diagonal-occupation relaxation-time **Fermi-surface proxy**, and the interlayer spin-transfer/local-torque operator identity. It does **not** compute the full Xiang SM intrinsic or AC response, spin accumulation at a boundary, adiabatic material polarization or any SOC-dependent interlayer hopping.

## 2. Hamiltonian and actual structural registry

Square Bravais layers, one spinful orbital per layer, same-site layer positions offset by `delta=(dx,dy)`. Write orbital-embedded tunneling exactly as

$$
t(\mathbf k,\boldsymbol\delta)=
\sum_{\mathbf R\in\mathbb Z^2} t_0
e^{-|\mathbf R+\boldsymbol\delta|^2/(2\ell^2)}
e^{i\mathbf k\cdot(\mathbf R+\boldsymbol\delta)}.
$$

Implementation truncates `R_x,R_y=-3...3` and differentiates **the physical bond phase** analytically; `ell=.8`, `t_0=.85` unless specified. The embedded basis has BZ sewing under `k_x->k_x+2pi`: `H(k+2pi,delta)=W(delta)H(k,delta)W(delta)^dagger`. It is incorrect to demand literal periodicity when the layer orbitals have displaced centers.

$$
H(\mathbf k,\boldsymbol\delta)=
\varepsilon(\mathbf k,\boldsymbol\delta) I+
\operatorname{Re}t\,\tau_x-
\operatorname{Im}t\,\tau_y+
\alpha\tau_z(\sin k_y\sigma_x-\sin k_x\sigma_y)+m\tau_z\sigma_z,
$$

$$
\varepsilon=-0.42(\cos k_x+0.8\cos k_y)+
\gamma(\boldsymbol\delta)\sin k_x,\quad \alpha=.45,\quad m=.8.
$$

Two distinct model arms:

- **Arm R, actual distance-dependent real interlayer hops only:** `gamma=0`; the only registry-dependent physical input is `t(k,delta)`.
- **Arm C, an *assumed* magnetic/orbital-current intralayer complex hopping:** `gamma=0.65 sin(2pi dx)`. This odd dispersion can be obtained from an oriented imaginary intralayer hopping, but **its registry dependence has not been derived from interlayer hopping or a microscopic orbital loop**. Treat it as a controlled counterexample of a falsely appealing switching result, not as an outcome of the model's structural tunneling.

There is also an optional *deliberately imposed* `gamma ∝ t_0^2` sensitivity mode. Making a term vanish at `t_0=0` by assumption is **not** proof of a new physical interlayer mechanism.

## 3. Magnetic symmetry, projectors, observables and limits

At each `k`, the antiunitary model representation `PT=\tau_x(i\sigma_y)K` obeys `(PT)^2=-1` and preserves H. Therefore each band doublet has an exact, non-Abelian two-dimensional spectral projector, and we trace both states together.

`C2z=i sigma_z`, rotating momentum and structure, gives

$$
H(-\mathbf k,-\boldsymbol\delta)=U_{2z}H(\mathbf k,\boldsymbol\delta)U_{2z}^{\dagger}.
$$

**Crucial:** the mapping `delta->-delta` by `C2z` does **not** reverse `P_z`; these are **not shown to be opposite sliding-ferroelectric states**.

For electron velocity `v_x=∂H/∂k_x` in the orbital-embedded gauge (includes `∂ t/∂k_x`), set the **conventional** spin current `J={v_x, sigma_y}/4` and `J_stag={tau_z,J}/2`. A PT transformation makes `J` odd and `J_stag` even. C2z makes both even as operators while changing the electric-drive direction; consequently the explicitly chosen Fermi-surface response is odd under the configuration switch.

The tested scalar, up to charge and relaxation-time normalization, is

$$
C_O=\int_{BZ}\!\frac{d^2k}{(2\pi)^2}
\sum_{d=\pm}\big[-f'(\varepsilon_d)\big]\,
v_{d,x}\operatorname{Tr}(P_d O),
\qquad v_{d,x}=\frac{\operatorname{Tr}(P_d v_x)}{2}.
$$

This is a **diagonal RTA Drude/occupation response**, not the full equilibrium/interband/spin quantum geometric susceptibility. In the insulator limit, the occupation derivative tends to zero at `T=0`; finite-temperature activated behavior is not an intrinsic insulator response.

For layer-staggered spin density `S_-^x=tau_z sigma_x/2`, check exactly

$$
i[H,S_-^x]=i[H_\perp,S_-^x]+i[H_{\rm local},S_-^x].
$$

A full finite-`q` continuity equation additionally needs spin flux divergences and torque-dipole/boundary terms; this uniform operator identity alone does **not** solve the spin-current-definition problem.

## 4. Tests and observed results (independently run 2026-10-21)

**Verified in Python:**

1. `(PT)^2=-1`, antiunitary H/v/J/Jstag parities, Hermiticity in the model.
2. `C2z` maps `(+dx,+dy,k)` to `(-dx,-dy,-k)` for H/v/J/Jstag.
3. `k_x+2pi` **orbital-embedding sewing** for both H and velocity, rather than a guessed periodic Bloch basis.
4. Uniform layer-spin exchange + local-torque continuity to roundoff.
5. `C_stag(+delta)=-C_stag(-delta)` at the chosen Fermi-surface protocol; global total `C_J=0` to roundoff.
6. Explicit `t_0=0` counterfactual, genuine in-gap chemical potential, and k-mesh convergence.

For `dx=.24,dy=.18; mu=1.3,T=.13`, in Arm C (an **assumed independent** gamma channel), rank-2-projector staggered coefficients:

| interlayer `t_0` | 41×41 | 61×61 |
| ---: | ---: | ---: |
| 0.00 | -0.115783222 | -0.115783219 |
| 0.40 | -0.127612494 | -0.127612481 |
| 0.85 | -0.118605284 | -0.118605077 |
| 1.50 | -0.095696161 | -0.095695990 |

The full switching coefficient remains **nonzero and similar even at `t_0=0`**. Hence our nonzero switch is **not evidence of a necessary interlayer transfer mechanism**. The value's `t_0` dependence is merely a correction to an already nonzero disconnected-layer effect, and a new prediction has not been isolated. Arm R yields zero to floating-point precision for the selected DC response, including nonzero registry `delta`.

At `t_0=.85,delta=(.24,.18)` Arm C, the indirect gap is approximately `0.2329`. Its VBM is `0.2307` and CBM is `0.4618`, so `mu≈0.346` is in the *actual* gap. At `T=.01`, `C_stag≈-6.7e-9`; at `T=.02`, about `-4.7e-6`. **Do not use mu=0 as the insulating reference:** mu=0 lies in the valence band of this particle-hole-asymmetric model.

The spin-independent registry hoppings *genuinely vary* with `delta` (though for `ell=.8` that variation is weak due to the broad real-space range); additional shorter-range `ell=.35` tests change the hopping amplitude substantially, while the Arm R DC response still vanishes in the tested range. Stronger registry dependence alone did not fix the unwanted response cancellation.

The finite `Rmax=3` truncation is small but is **not exactly invariant** under `delta -> delta+integer`; an actual stacking-phase diagram must converge Rmax and impose a periodic registry continuation. The script checks the **k-space embedding seam**, not this separate translation-convergence gate.

## 5. Scientific NO-GO result and why it is useful

- **Proved in this model:** model symmetries/partner maps, rank-2 Fermi-surface proxy and the operator identity.
- **Refuted for this candidate:** a nonzero signed layer-staggered response is sufficient to establish a *new interlayer quantum-geometric mechanism*.
- **Not proved:** a physical FE bistability, layer spin flow conservation in SOC, a measurable edge spin accumulation, a new geometric tensor, or novelty over the 2024–2026 layer-spin-Hall/quantum-metric/NSOT literature.
- **Better next gate:** replace the independent `gamma(delta)` channel by **explicit spinful interlayer hopping derived from orbital symmetry**, and test whether a residual `Delta C=C(t_perp)-C(0)` is nonzero and physically distinguishable in an actually paired stacking configuration. Then compute a *full causal finite-frequency susceptibility* and a spin accumulation/torque observable, not just an RTA Drude number.

## 6. Peer-work boundary

Read with the [peer-work audit](PEER_WORK_AUDIT_20261021.md), especially the 2024 sliding-layer-spin Hall, 2024 extended-space layer multipoles, 2025/2026 magnetic stacking theory, 2025 insulating NSOT, and Shi et al. proper-spin-current discussion. The existence of this model symmetry and of a nonzero `C_stag` is already within familiar ingredients; there is **no originality claim**.

**Decision:** Gate 1 complete as a negative control. Proceed only after finding a microscopic interlayer spinful channel and defining an independently measurable observable. Do not modify paper-specific PRL reproduction code or promote this branch to `PAPER_REPRODUCED`.


## Gate 1b: critical four-band PT constraint (additional audit)

An exhaustive independent test of all 16 Hermitian `tau_a sigma_b` matrices for the **specific** `PT=tau_x(i sigma_y)K` representation gives just six PT-even matrices: `I, tau_x, tau_y, tau_z sigma_x, tau_z sigma_y, tau_z sigma_z`. Thus all direct `tau_x sigma_i` and `tau_y sigma_i` interlayer spinful hopping operators are PT-odd in this four-band representation. **Correction to the preliminary future-work suggestion above:** simply inserting PT-even spinful hopping in this same one-orbital/layer model is impossible. Next consider at least two orbitals per layer (eight-band), with physical orbital-dependent PT action, or a physically justified PT-breaking route. [16-operator numerical proof](../scripts/pt_operator_basis_gate.py). This is representation-specific, not a universal bilayer no-go.
