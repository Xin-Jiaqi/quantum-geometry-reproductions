# Xiang, Jin & Wang (2025): Spin Transport Revealed by Spin Quantum Geometry

Independent **paper-model** subset for PRL **135**, 146303 (2025), DOI `10.1103/14mp-263q`.

## Reproduce without publisher data

```sh
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/reproduce.py
```

This GitHub directory is a **curated V30 paper-model plus method subset**. The complete V1–V27 offline archive (including original published-vector extraction, all 54 V25–V27 candidate tests and advanced V14–V19 derivations) is **not claimed to be uploaded here**. Full on-paper scientific status remains **PARTIAL**; software CI success does not imply article reproduction.

## What is independently executable here

- Main Fig. 1(a): two-valley Dirac NDSC Eq. (9) **and reduced intrinsic Eq. (13)** (Fermi-sea/surface separated, independently checked Bloch matrix oracle), band spin-QGT analytic and matrix checks. Original absolute axis has an unresolved common factor of ~2. No artificial correction is applied.
- Main Fig. 1(b): SM Eq. (67), printed main text, printed SM Eq. (70), literal SM integrand and a separate *static density-matrix* branch. Printed results and independent results are **never silently merged**. Agreement with original publisher curve was checked offline; publisher data are not mirrored.
- SM Fig. 1(a),(c),(d): tilt-Dirac band, spin metric and curvature.
- SM Fig. 1(b): printed Eq. (56), independent 2D integration of Eq. (55), and Eq. (58). The published Eq. (56) does not equal the literal angular integration of Eq. (55); both are reported.

## Not included / not validated

- Reduced Eq. (13) does **not** establish a general intrinsic τ^0 solver: its omitted four-band channels and source DC protocol remain open.
- Paper's uniquely specified **metallic driven-DC intrinsic** `tau^0` response (the static density-matrix branch is not a substitute).
- The V30 **generic finite-η multiband resolvent** and **PT rank-two projector method** are now published in `methods/`. This is a narrowly scoped finite-η method migration, not the full offline band-covariant V14–V19 solver, and not the paper's general intrinsic τ⁰ response.
- MnBi₂Te₄ DFT/Wannier material calculations (out of scope by project decision).
- Stacking, intercalation, and group-theory toolkit integration (deferred to an independent `studies/` module after paper model acceptance).

Scientific discrepancies and limitations are explicitly recorded in `docs/ACCEPTANCE.md`. Nothing here is original authors' source code or a publisher PDF.

## V29 Eq. (13) migration

`paper_eq13.py` retains the reduced printed two-band implementation from the frozen offline V27 project, with independent Bloch-matrix element checks (`kernel_matrix_at_k`). `python scripts/reproduce.py` now prints `fig1a.intrinsic_eq13_reduced_mu_scan` at four chemical potentials; `results/REFERENCE_V29.json` freezes separate Fermi-sea and Fermi-surface branches. No ad hoc factor has been added to reconcile Fig. 1(a) with the published ordinate.

## V30 multiband and PT finite-η method interface

`methods/finite_eta.py` accepts any smooth Hermitian Bloch Hamiltonian `H(kx,ky)` and a same-basis Hermitian observable; it computes the finite-η full second-order density-matrix response with independent band and orbital Liouville inversions. `methods/pt_doublet.py` supplies a synthetic exactly-PT four-band benchmark with projector inversion, U(2) gauge rotations and an independently frozen V16 staggered-spin point. The solver does **not** drop divergent finite-η Drude contributions to imply a verified intrinsic τ⁰ response.

```bash
python scripts/reproduce_methods_v30.py
python -m unittest discover -s tests -v
```

See [`docs/METHODS_ACCEPTANCE_V30.md`](docs/METHODS_ACCEPTANCE_V30.md) for the exact protocol, comparison and exclusions.
