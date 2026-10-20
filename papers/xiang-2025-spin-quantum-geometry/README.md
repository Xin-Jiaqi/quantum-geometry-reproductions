# Xiang, Jin & Wang (2025): Spin Transport Revealed by Spin Quantum Geometry

Independent **paper-model** subset for PRL **135**, 146303 (2025), DOI `10.1103/14mp-263q`.

## Reproduce without publisher data

```sh
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/reproduce.py
```

This GitHub directory is a **curated V28 paper-model core**. The complete V1–V27 offline archive (including independent multiband, non-Abelian Liouville, original published-vector extraction and 54 targeted V25–V27 candidate tests) is deliberately **not claimed to be uploaded here**. Full on-paper scientific status remains **PARTIAL**; software CI success does not imply article reproduction.

## What is independently executable here

- Main Fig. 1(a): two-valley Dirac NDSC Eq. (9), band spin-QGT analytic and independent matrix checks. Original absolute axis has an unresolved common factor of ~2. No artificial correction is applied.
- Main Fig. 1(b): SM Eq. (67), printed main text, printed SM Eq. (70), literal SM integrand and a separate *static density-matrix* branch. Printed results and independent results are **never silently merged**. Agreement with original publisher curve was checked offline; publisher data are not mirrored.
- SM Fig. 1(a),(c),(d): tilt-Dirac band, spin metric and curvature.
- SM Fig. 1(b): printed Eq. (56), independent 2D integration of Eq. (55), and Eq. (58). The published Eq. (56) does not equal the literal angular integration of Eq. (55); both are reported.

## Not included / not validated

- Paper's uniquely specified **metallic driven-DC intrinsic** `tau^0` response (the static density-matrix branch is not a substitute).
- Spinful multiband/`PT` degeneracy solver and generic finite-`eta` Liouville tools: these exist in the offline V27 package but are **not part of this first remote curated core**.
- MnBi₂Te₄ DFT/Wannier material calculations (out of scope by project decision).
- Stacking, intercalation, and group-theory toolkit integration (deferred to an independent `studies/` module after paper model acceptance).

Scientific discrepancies and limitations are explicitly recorded in `docs/ACCEPTANCE.md`. Nothing here is original authors' source code or a publisher PDF.
