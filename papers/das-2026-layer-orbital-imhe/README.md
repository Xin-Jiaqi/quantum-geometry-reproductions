# Das & Agarwal 2026 · layer–orbital mixed-response core

**Source:** [arXiv:2604.20249v1](https://arxiv.org/abs/2604.20249v1), *Intrinsic Magnetoelectric Hall Effect from Layer-Orbital Quantum Geometry*. Independent implementation, not author code.

**MODEL_READY for the explicitly declared effective-Hamiltonian protocol below.** This connects the core mechanism and reusable multiband calculation; it does not claim R5G parameter/figure reproduction, the reported 0.05 e²/h magnitude, or a universally valid orbital transport theory. Source disagreements do not define numerical fitting targets.

## What can be called

| Entry | What it computes |
|---|---|
| `Fields` | H₀, layer P, reference orbital M, optional mixed operator C, and their k derivatives in one fixed orbital basis |
| `geometry(*fields.jets())` | All 0, electric, magnetic, mixed coefficients of energies, Berry curvature, interband connection and band velocities |
| `orbital_moment(H,dH,ddH)` | A5 interband reference orbital-moment matrix, A1 diagonal continuation and analytic k derivatives; no spin operator substitution |
| `transport(...)` | Eq. (6) two-term sector, all five A16 terms, Eq. (7) sector, and complete mixed derivative of the stated band-drift functional |
| `exact_transport`, `mixed_difference` | Independent four-corner finite-field diagonalization of the same response functional |
| `projector_curvature` | Independent k-space projector finite difference, without eigenvector phase matching |

All intermediate bands, normalized second-order eigenstates and the two perturbation orderings are retained. Exact/near degeneracies fail closed. The included three-orbital, spinless periodic layered model has nonzero intrinsic response while preserving T. A separate magnetic-order parameter activates the T-odd drift sector. It is **our benchmark**, not a replacement label for the missing R5G setup.

## Protocol and units

At a chosen **background** gate, define independent small probe coordinates u and b:

$$
H(k;u,b)=H_0(k)+uP(k)-bM_0(k)+ubC(k).
$$

M₀ is constructed at that reference and frozen as an orbital-basis operator under the probes. Recomputing M after changing u is a different protocol; its mixed contact must be supplied explicitly in C. Field-dependent P or M, physical minimal coupling and magnetic positional shifts must not be silently inferred from this interface.

Natural units e=ℏ=1; k is dimensionless; the full torus uses explicit weights for d²k/(2π)². The thermal argument is kBT in **model energy**, not kelvin. Results are per u·b; drift coefficients are per τ. These numbers are **not SI conductivities**. Real ℰ/B conversion requires declared physical operator, energy and length units. The two scalar coefficients multiply a separately applied in-plane bias.

The defined full functionals are

$$
\sigma_{yx}^{\rm geom}(u,b)=\int\Omega(u,b)f[\varepsilon(u,b)],\qquad
D_{ab}(u,b)=-\int v_a(u,b)v_b(u,b)f'[\varepsilon(u,b)].
$$

`intrinsic_full` and `band_drift_full` are their mixed derivatives. The latter does **not** include Lorentz deflection, magnetic phase-space drift corrections or scattering vertices. It must not be named the complete physical magnetoconductivity. See [derivation and source audit](PROTOCOL_AND_SOURCE_AUDIT.md).

## Run and reuse

From the repository root:

```bash
python scripts/core.py run imhe --output outputs/imhe.json
python scripts/core.py test imhe
```

Or within this directory (NumPy only):

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/reproduce_core.py --outdir outputs
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -W error -m unittest discover -s tests -v
```

```python
from imhe_lab import mesh, layer_fields, geometry, transport
k, weights = mesh(48)
fields = layer_fields(k, gate=0.23, magnetic_order=0.0)
geometry_jets = geometry(*fields.jets())
result = transport(geometry_jets, weights, mu=0.2, temperature=0.15)
print(result['intrinsic_eq6'], result['intrinsic_full'])
```

To change models, supply `Fields(H,dH,P,dP,M,dM[,C,dC])` with H/P/M/C shape `(Nk,N,N)` and derivatives `(2,Nk,N,N)`. Do not infer layer operators, physical spin or orbital embedding from the matrix dimension. k-dependent basis changes require a separately validated covariant extension; current acceptance covers fixed orthonormal bases and constant orbital-unitary transformations.

## Validation and numerical baseline

`results/validation.json` is a newly calculated benchmark, not publisher data. The recipe recalculates rather than reading this file; `--outdir` also emits `mu_scan.csv`.

At gate=0.23, μ=0.2, kBT=0.15, 32² mesh:

- intrinsic Eq. (6) sector: −1.79624866178×10⁻⁵;
- full A16-structure mixed derivative: −2.40101697118×10⁻⁴;
- the T-symmetric drift response cancels to numerical precision;
- magnetic_order=0.13 gives nonzero full drift xx ≈ 2.67723827×10⁻⁴; reversing that order reverses the drift tensor and preserves the intrinsic coefficient.

The omitted single-field cross terms are significant in this benchmark. They cannot be dropped merely because each factor separately carries one field.

29 regression tests cover normalized eigen-equations, nonzero six-band/contact examples, field-step and k-mesh convergence, projector curvature, orbital-moment derivatives, basis/phase gauges, T and magnetic-domain parities, scalar-gate occupation effects, zero couplings, finite-temperature derivatives and malformed-input rejection. See the audit for the strict limits of these checks.

## Important source checks

1. A4 has the opposite inverse-energy-gap derivative sign under Eq. (1)'s stated P convention. A two-band counterexample is executable.
2. First-order states alone miss second-order ket/bra terms needed for A8-like velocity corrections. The implementation verifies the full eigen-equation, not the truncated product.
3. A16's product-rule structure is retained; Eq. (6) and Eq. (7) are exposed as selected sectors, not relabelled as full mixed derivatives.
4. A7's field-dependent moment convention is not silently identified with the frozen-M₀ model. Full R5G and source-specific magnetic dressing remain unverified.
5. A gapped periodic filled-band control has nonzero local mixed curvature but zero integrated mixed Hall derivative. We do not force a nonquantized insulating plateau. This is a necessary check within the stated smooth periodic effective-Hamiltonian/DC protocol, not a proof about every nonequilibrium boundary protocol.

No publisher PDF/SM, private correspondence, user notes or unpublished material inputs are redistributed. The existing six paper kernels and their frozen numbers are unchanged.
