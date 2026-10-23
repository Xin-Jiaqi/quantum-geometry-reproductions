# Round26 | Same-Wannier-gauge physical-operator input contract

**Scientific gate: synthetic in-memory operator interface only. FULL CUI PAPER PARTIAL.**

## Gap closed and scope

Round25 established a *static* k-dependent Bloch-basis Ward on an explicitly supplied synthetic spin vertex. Before attempting material Zeeman QGT, the missing upstream gate is a **joint** set of same-gauge real-space Hamiltonian H(R), projected **physical** spin S(R), projected position connection A(R), orbital centers and lattice, not just a coincidentally shaped H and three arbitrary matrices. Round26 provides a fail-closed *in-memory numerical contract* and Fourier bridge into existing multiband geometry. It **does not parse** raw wannier90_hr.dat, spn, amn/mmn, or r.dat, nor reconstruct the physical spin from the source Bloch gauge. Matching metadata IDs alone cannot authenticate an external spin-gauge transformation; this requires an upstream independently auditable provenance trail.

## Conventions and derivation

Two real-space Bravais vectors are rows of a 2×3 matrix in Å. A translation R has two **integer** components and Cartesian vector R_cart=R@lattice. Fractional reciprocal coordinate q is dimensionless. For each of H, S_b and position A_a, use the **same** convention:

\[
X(q)=\sum_R e^{2\pi i R\cdot q} X_R,\quad X_{-R}=X_R^\dagger.
\]
\[
\partial_{k_a}X(q)=\sum_R i R_{\mathrm{cart},a} e^{2\pi i R\cdot q}X_R,\qquad
D_aX=\partial_{k_a}X-i[A_a,X].
\]

H(R) in eV, physical Pauli spin S_b(R) in **dimensionless** units, projected Wannier position connection A_a(R) in Å. The orbital centers **already appear** as diag A_a(R=0), and must not be added again. The projected physical spin may fail to satisfy full Pauli algebra after subspace truncation; a necessary norm bound is used at q=0 rather than imposing exact algebra for all projected models. D_aH is an energy×length vertex, not group velocity (which requires division by ħ). Under a **constant joint** Wannier rotation U, all H/S/A transform as U†XU, and the individual orbital-center diagonal changes accordingly. Round25 covers the more general k-dependent basis algebra; no inhomogeneous k-gauge connection is guessed by this importer.

For a spatially varying SOC spin texture, [H,S_b] can be nonzero. The synthetic four-band example uses true spin-1/2 × two orbitals, an onsite spin-orbit term and a spin-dependent hopping; its nonzero spin torque is a *positive control for physical spin nonconservation*, not a calculation of conserved spin current.

## Independent numerical oracles and negative controls

- Analytic real-space Fourier reconstruction and Cartesian derivatives; periodic reciprocal shifts and Hermiticity.
- **Independent Wilson/parallel-transport finite-difference oracle**: to first order U=exp(-i Δk A_a), compare [U H(k+Δk) U† - U† H(k-Δk) U]/(2Δk), similarly for physical spin, against D_a H and D_a S.
- Constant independent random 4×4 unitary on all operators: full local Zeeman-QGT band-projector tensor must be invariant. Also independently compare with existing multiband spectral and projector primitives.
- Necessary provenance/shape/translation/unit guards; explicit rejection of pseudospin, missing spin evidence, non-Hermitian hoppings, missing -R partner, double-counted/mismatched embedding and spin norm >1.
- **Bad physics that a falsified metadata ID cannot prevent**: rotate H and A but *not* physical spin. The Zeeman QGT differs substantially. A second negative control drops A and computes a raw ∂H in place of the covariant vertex.

At q=(0.183,-0.271), the independent *local* synthetic numerical experiment obtains Wilson D_kH residual 1.51e-10 eV·Å, Wilson D_kS residual 5.83e-16 Å, joint-gauge QGT residual 2.09e-15; WRONG spin gauge QGT difference 1.48e-1 and OMITTED A velocity vertex difference 2.63e-1 eV·Å. The physical spin torque max is 0.44 eV. Values are model diagnostics, **not a measured transport susceptibility**. Freeze: [round26_same_gauge_operator_audit.json](results/frozen/round26_same_gauge_operator_audit.json); reproduce with:
~~~bash
python -m unittest discover -s tests -v
python scripts/audit_wannier_operator_contract.py
~~~

13 new tests in the isolated Cui suite. Actual Python 3.11/3.12 acceptance is determined **only by remote CI**; neither the historical 1105 repository tests nor Round21's 167 count should be relabeled as Cui Round26 coverage.

## Hard exclusions

No authentic Wannier physical-spin matrices or position/cell embedding have been attested for HgI2, AlAgS2, Hg3O3 or any other material. This does not produce non-Abelian metallic transport, finite-q electromagnetic Ward, conserving collisions, Drude response, original Eq.(19) in SI, all original figures or a Table I resolution. A real-material stage must independently verify common source gauge, the *physical* spin provenance, position operator, band/embedding calibration, BZ convergence and normalization. **Original paper remains PARTIAL.**
