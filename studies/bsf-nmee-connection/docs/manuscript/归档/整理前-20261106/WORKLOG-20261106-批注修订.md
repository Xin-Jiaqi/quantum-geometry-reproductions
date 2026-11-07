# WORKLOG — 2026-11-06 annotated manuscript revision

Scope: addressed the comments in the user's Obsidian 正文讨论稿; original preserved. No new DFT, remote job, or integrated response scan.

## Deliverables

- Expanded Chinese results draft; standalone monolayer/bilayer Hamiltonians and static response definition.
- Comment-by-comment response, theoretical definitions, reference figure review, supplementary discussion.
- Three main and three supplementary figures from archived data, with PNG/PDF/SVG and figure_data.json provenance.
- Re-read Zotero NMEE main/SM, disorder, band-edge mixture, LEE, blue-phosphorus interface, and Science interface SM. Private paper renders excluded from this archive.

## Physics carried forward

The decoupled reference resolves the apparently opposite trend: KM lower-layer baseline is negative and the finite interface correction positive; MoSSe projected/embedded lower-layer baseline is positive and the correction negative. This organizes the results but does not yet predict the correction sign. Candidate even-in-hopping expansion and layer self-energy are conditional theoretical starting points, not fitted material laws. Next mechanism question: energy denominators versus field derivatives of spin/electric matrix elements, then identify robust band groups.

Published NMEE two/three-band intrinsic expressions are explicitly distinguished from the full static coefficient actually computed. No claim of PRL novelty or experimental detectability established. MoSSe remains the material reference; changing material deferred.

## Checks

- Existing implementation versus written Bloch and static-response formulas: formula_checks.json; three general k points only.
- KM layer-sum max residual: 1.215e-15; signed-k split residual: 1.067e-15.
- Material A/B layer-sum residual in displayed units: 3.286e-14.
- Occupied-band decomposition differs from direct material response by ~0.05% and ~0.27%; stated in theoretical note, not treated as exact closure.
- All local document and figure links resolved; display-math delimiters paired; no NUL.
- Six figures visually inspected. Static preflight: 12 pass, 2 warnings for draft raster export (no TIFF, 320 not 600 dpi), 0 fail. Vector exports provided.
- Original annotated file SHA256: a14a7cd5397acde40a5b62a87be06a8487c5cc518c361937f39eb3dc5f3b7ff5.

## Reading order

正文讨论稿-修订版.md → 批注回应与机制问题.md → 理论定义与公式核对.md; supplementary figures and reference review are separate.

No remote push in this revision; local research snapshot only.
