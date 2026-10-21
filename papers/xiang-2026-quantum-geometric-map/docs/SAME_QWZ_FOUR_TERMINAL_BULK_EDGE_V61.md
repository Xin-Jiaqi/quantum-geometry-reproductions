# v61 | Same-QWZ four-terminal charge Hall edge/bulk comparison

**Scientific outcome:** VERIFIED for the *named* QWZ Hamiltonian, its exact principal-layer continuation, finite-energy Caroli/Landauer-Büttiker calculation, and independently tested bulk Chern number. **CONDITIONAL** as a finite-coherent-device transport observation (finite η, μ, contact geometry, zero T). **UNRESOLVED** as a general thermodynamic disordered/contact-independent dc limit or original Xiang orbital mixed E/B response. Xiang source review remains **PARTIAL**.

## Physically matched Hamiltonian

H(k)=sin(kx) σx+sin(ky) σy+[m+cos(kx)+cos(ky)] σz.

The central Lx×Ly sample has on-site mσz; positive links Tx=(σz−iσx)/2, Ty=(σz−iσy)/2. All four semi-infinite matrix principal-layer leads contain exactly these QWZ blocks. Lead order L,T,R,B uses outward hoppings Tx†,Ty,Tx,Ty†, respectively. This fixes the microscopic tight-binding Hamiltonian, rather than inserting unrelated wideband contacts. The cross remains a *particular* contact geometry, not the periodic bulk.

## Observable and signed convention

Four-terminal Caroli transmissions T_ab=Tr[Γ_a G^r Γ_b G^a] and floating I_T=I_B=0 are solved with V_L=+1/2, V_R=−1/2. Measured R_xy=(V_T−V_B)/I_L, units h/q². Bulk σ_xy is a separately computed spectral/FHS Berry Hall conductivity, units q²/h. These are **not** the same observable. For a clean single-channel chiral plateau, σxx→0 and the Hall resistivity is obtained by conductivity-tensor inversion; this specific lead orientation yields R_xy≈1/σ_xy for C=±1. This is a bounded bulk-edge/transport consistency check, *not* an arbitrary finite-device Hall-voltage equality with a Chern number.

## Reproduction (η=10⁻⁴, μ=0, T=0, exact contacts)

| mass | L | bulk Chern | Rxy/(h/q²) | I_L/(q²/h), unit bias | status |
| --- | ---: | ---: | ---: | ---: | --- |
| −1 | 3 | +1 | +1.000006918 | 0.996911307 | Conditional plateau |
| −1 | 5 | +1 | +0.999999914 | 0.999889021 | Conditional plateau |
| −1 | 7 | +1 | +0.999999997 | 0.999995134 | Conditional plateau |
| +1 | 5 | −1 | −0.999999914 | 0.999889021 | Opposite chirality |
| +3 | 5 | 0 | **undefined** | 1.01277×10⁻¹⁰ | Insulating negative control |

The trivial insulating state gives exponentially small lead transport. Naively dividing its tiny transverse voltage by tiny current gives an arbitrarily large, regulator-dependent **fake** Rxy. The implementation explicitly returns `None` when |I_L| < 10⁻⁷ for a unit bias. This reporting threshold is a documented measurement resolution convention, not a universal physics cutoff.

## Independent acceptance

28 v61 local gates: Fourier reconstruction against independently formulated Bloch H, Hermiticity, lead hopping/coupling orientations, positive Γ, nonnegative and bounded Caroli transmission, independent scattering-matrix S†S and block probabilities, four-contact current conservation and probe solution, ±Chern, mass/sign reversal, size/η/μ controls, finite-contact perturbation, local orbital rephasing of both H and couplings, omission-of-contact-transformation negative control, insulating-current refusal, and invalid-input guards. v60 21 regressions re-run as well.

The retarded lead η is a controlled numerical regulator (no claim a universal η→0 contact limit); model clean zero-temperature four-terminal charge Hall physics is not the original paper's intrinsic orbital E×B, nor a general materials solver.

## Next v62 science gate

Add independent *real-space current-density/local-scattering* probes, multiple genuine voltage-contact widths and independent Kubo current-current conductance at the same filling and thermodynamic order of limits. Only claim model-level equivalence where observables and contact/limits are properly matched. Then return to position-operator-complete orbital magnetic perturbations and conserving bulk mixed E/B.
