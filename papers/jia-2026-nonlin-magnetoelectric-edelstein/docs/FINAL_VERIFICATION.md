# Final reproducibility verification (2026-10-19)

The user's acceptance goal is **reproduction of the physics/methods and a reusable software interface**, not pixel matching to unpublished figure data.

**Current scientific release status: PARTIAL, not complete.** This status is deliberate and must not be replaced by a software-test PASS. The following has been independently rerun using only the released source tree:

| Gate | Evidence | Result |
|---|---|---|
| Software regression suite | `python -m unittest discover -s tests -v` | PASS (30 regression tests; independently rerun) |
| Kane–Mele Eq. (6), clean response | explicit printed vs full-static protocols; 451² reference in `results/` | PASS within declared protocol, not identical to published curves |
| CuMnAs Eq. (5) bands, PT degeneracy | analytic two-doublet dispersion, Dirac positions | PASS |
| CuMnAs Fig. 1(d) Eq. (2) | PT projector / spin-space metric, Fermi-surface derivative, tau, Neel reversal | PASS for formula structure and symmetry; absolute figure magnitude remains conditional |
| CuMnAs Fig. 1(c) | staggered spin response, finite-frequency, Dirac regularization, defined BGK response | **NOT PASS** as the author's finite-T driven-dc intrinsic Eq. (1) response |
| Source references | PDF SHA256 and official DOI / SM links in `references/README.md` | PASS; copyrighted PDFs not mirrored |
| Research note provenance | Verbatim user-supplied PRL Section 3.3 in `docs/RESEARCH_NOTES.md` | PASS |
| Existing Wannier workflow contract | signed Stage-1 package input, fail-closed if necessary physical operators absent | PASS as preflight, **not** full material NMEE |

## The one physics blocker

The PRL/Supplemental Material identifies a split into intrinsic and extrinsic terms. The independently tested **full static** formula contains an additional diagonal density-matrix contribution relative to the final printed two- and three-band intrinsic formula. We have not proven that the authors' nonequilibrium zero-frequency convention is equivalent to the equilibrium/static protocol, especially for **finite-T exactly PT-degenerate CuMnAs**. Therefore, the static Fig. 1(c) output may not be relabeled `paper_defined_intrinsic`.

The paper-defined Fig. 1(c) coefficient and its convergence/normalization form the final scientific acceptance gate. The other outstanding issue is a quantitative, real-operator Wannier integration test; in accordance with the user's scope limitation, no new material study is required to address the source-model PRL reproduction itself.

The independent release gate `python scripts/verify_release.py` additionally checks five frozen model outputs and SHA256 fingerprints. This is a **software/model gate**, not evidence of resolution of the intrinsic DC-limit question.

## What 'complete' would mean

1. An independent, PT-degenerate, finite-temperature evaluation of the paper-defined `intrinsic` response with explicit order of limits and numerical convergence through Dirac regions, *or* a documented derivation showing why the existing operational definition is equivalent.
2. Equation and unit cross-check of intrinsic/extrinsic components with an end-to-end CLI and frozen output metadata, including Fig. 1(c).
3. Once these gates pass, run the existing 30-test suite, integrate the module into the Stage-1 package handoff, and tag a version. Do not mark 100% before (1).
