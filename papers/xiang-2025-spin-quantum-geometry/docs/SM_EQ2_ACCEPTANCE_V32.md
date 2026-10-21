# V32 — Source-by-source audit: SM Eq.(2) to Eq.(6)

Paper: Xiang, Jin & Wang, PRL 135, 146303 (2025).

Status: PASS for specified *nondegenerate* SM Eq.(2) algebra and finite-eta protocol. PARTIAL for the full paper and physical driven-DC intrinsic response.

## Physics and what was implemented

The code performs exactly four groups of the printed SM Eq.(2):

1. Diagonal second derivative of the Fermi occupations.
2. Covariant momentum derivative of offdiagonal first-order interband coherence.
3. Population-gradient source mixed with interband Berry connection.
4. Intermediate-band coherence products (including true l distinct from m and n).

The independent reference is fixed-orbital density-matrix Liouville recursion, which does not use Eq.(2) termwise formulas. An additional V30 orbital Liouville cross-check provides a third numerical route. Independent original offline V14 benchmark values are preserved in results/REFERENCE_V32.json.

Hamiltonian testbed: 1D smooth N-band nondegenerate H(k)=U(k)diag(e(k))U(k)^dagger, N=2,3,4,5. This is a mathematical testbed, **not a prediction of MnBi2Te4**. Bloch phase gauges are explicitly varied.

## Regulator and physical boundaries

omega1=omega2=0, first-order adiabatic regulator eta, second-order regulator 2*eta, and one ordered electric field channel ExEx. eta is **not automatically a microscopic collision rate**; finite eta is not a generic material intrinsic tau^0.

The N=2,3,4,5 maximum differences between literal Eq.(2) source sum and independent orbital recursion at k=.23, eta=.02 are:

- N2: 1.63653e-9
- N3: 5.21249e-9
- N4: 2.81515e-9
- N5: 7.50075e-9.

Literal Eq.(2)'s offdiagonal eta^0 coefficient **does not equal the printed Eq.(6)** for a varying Fermi occupation. Their difference is:

    delta_C_mn = -i * r_mn * (partial_k f_n - partial_k f_m) / (epsilon_m - epsilon_n)^2

For N4, max extra element 0.0028400173; the frozen synthetic orbital observable's extra response is +0.00066274565 (natural units). The independent triple-band term is nonzero but does not cancel this difference. Freezing the dispersions makes partial_k f zero and the discrepancy vanish.

This observation concerns literal algebra under the stated adiabatic prescription. It **does not** prove the author's finite-frequency/symmetrized/microscopically scattered metallic intrinsic coefficient is incorrect. The publication status remains PARTIAL, and no manuscript numbers were fitted.

## Programmatic verification

From the Xiang paper directory:

    python -m pip install -r requirements.txt
    python -m unittest discover -s tests -v
    python scripts/reproduce_sm_eq2_v32.py

V32 adds 19 source-model tests plus 5 separate V30/frozen cross-method tests. Total paper-module suite: 70 tests. CI uses Python 3.11 and 3.12. The constrained real-material/stacking physics is *not included*.

## Publication and offline provenance

V32 is a curated reimplementation of offline V14 code and source-data checks; detailed original numerical figure extraction, private research notes, and publisher PDFs remain offline. V30 exact PT doublets and V31 generalized isolated-cluster Laurent formulas are **separate algorithms**. Do not substitute this nondegenerate Eq.(2) matrix expression into an exact PT degenerate subspace, where term denominators become ill-defined.

Follow-on objective: freeze accepted *paper model* and reusable methods; only then launch separate studies with the group-theory toolkit.