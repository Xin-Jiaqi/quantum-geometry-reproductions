# V30: finite-η multiband / PT method migration and acceptance

## What has been migrated

A minimal independently tested *method* module, not an independent new material prediction.
Source provenance: offline V16 `src/pt_doublet_liouville_v16.py`, V25 `src/tb_spin_geometry_v25.py`, V26 independent backend checks; V30 implementation is **curated/reimplemented**, not a byte-for-byte copy of the entire V27 offline archive.

- `methods/finite_eta.py`: general Hermitian `H(kx,ky)` in fixed orbital basis, Fermi occupation, resolvent for pη−iω+i[H,⋅], independent eigensystem and orbital-superoperator backends, fourth-order k differentiation and symmetrized (+ω,−ω) second-order density matrix.
- `methods/pt_doublet.py`: exact PT gamma-matrix test model, rank-two subspace projectors and a U(2)×U(2)-invariant block resolvent.
- `tests/test_methods_v30.py`: 12 tests covering N=2…6, six-band three-doublet example, the exact doublet projectors, U(2) covariance, both resolvents, PT symmetry and parity, regulator dependence, frequencies and the **independent offline V16 checkpoint** (k=(0.38,0.71), η=0.1, μ=1.1, T=0.17: staggered spin response 2.0232340477418247, tolerance 4e−8).

## Explicit protocol

```
ρ₁(ω) = − L₁(ω)⁻¹ ∂kx ρ₀
ρ₂(0) = − L₂(0)⁻¹ ∂kx (ρ₁(+ω)+ρ₁(−ω))/2
Lₚ = pη−iω+i[H,⋅]
```

`η>0` is the *adiabatic regulator* in the prescribed calculation. The value of `Tr(O ρ₂)` is a **finite-η full response** with Drude poles, not the paper's universal intrinsic τ⁰ coefficient. The finite cosine field Fourier coefficients, electric charge / dimensional SI units and velocity ℏ factors are not included. Matrix density response is defined at a point, not an already converged BZ transport number. Cannot substitute for a physical collision model or DFT.

## Validation / scope

The V30 implementation was independently compared against offline original V16 for a PT staggered observable and against separate orbital-/band-resolvent calculations, including exact U(2) gauge rotations. The previously reported 2× figure amplitude discrepancy and original source-algebra tensions remain OPEN; neither can be fixed by this adapter. Fig. 1(a) reduced intrinsic Eq.13 lives separately in `paper_eq13.py` and must not be identified with `methods/finite_eta.py`.

**Excluded:** offline V14 full four-source SM Eq.2 term audit, V17–V19 Laurent / multi-doublet non-Abelian decompositions, full band-covariant derivatives, Wannier gauge projected position/spin, real-material DFT, stacking or intercalation. Preserve the older offline checkpoints until these are separately and fully migrated.

## Run

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/reproduce_methods_v30.py
```

Software PASS is not an author-paper physics PASS.
