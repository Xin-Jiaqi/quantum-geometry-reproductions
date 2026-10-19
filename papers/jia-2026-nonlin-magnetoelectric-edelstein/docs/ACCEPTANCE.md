# PRL method-reproduction acceptance (2026-10-19)

**Goal:** physically reproducible Fig. 1/2 tight-binding calculations and a testable pipeline handoff, not author data/pixel-level overlays. The following is the reproducible current status, not an asserted 100% PRL validation.

| Part | Independent model/method | Tests | Acceptance |
|---|---|---|---|
| Fig. 1(b) | Eq. (5) PT-doublet spectra and Dirac nodes | analytic dispersion, band gap | **PASS for model** |
| Fig. 1(c) | staggered PT sublattice response | A+B zero, static Dirac divergence, PT finite-frequency and BGK tests | **PARTIAL**: paper-defined degenerate finite-T driven-dc intrinsic coefficient unresolved |
| Fig. 1(d) | published Eq. (2) projector trace Fermi-surface `tau` response | tau, length, Neel reversal, independent PT tests | **PASS for algebra/mechanism**, numeric absolute paper curve unverified |
| Fig. 2(a) | Kane–Mele Eq. (6) tight-binding | gap and k-resolved bands | **PASS for model physics**, full published path numerical digitization not part of release |
| Fig. 2(b-d) | printed versus SM full-static kernels | 451² baseline; K/K' local integration, SOC criticality and hotspot audits in offline archive | **PASS for internal method**, not fully numerically matched to published figures |
| Real Wannier material | HR and covariance import / stage1 signed package | 2D HR, synthetic nonzero connection, gauge/unitary, contract hash check | **PASS for software input preflight**; material physical operators absent |

## Release gates

- `python -m unittest discover -s tests -v` passes and reports the count.
- All `results/*.json` are generated from fixed commands in README and include protocol names.
- Staggered-response and Neel reversal tests run on symmetric **k / -k** grids; the previous asymmetric offset gave finite-N artifacts.
- `integration/stage1_handoff.py` checks signed files and blocks incomplete quantum geometry data.
- No raw PRL/SM PDF, publisher figure, private data, prior 33 ZIP archives or license assumptions are silently included.
- No claim that a missing material quantity was reconstructed from `hr.dat` alone.

## Remaining scientific issue, not a coding bug

The exact equivalence between SM Eq. (9)'s static density-matrix term and the article's finite-temperature driven dc intrinsic definition is not established. Fig. 1(c)'s PT-degenerate finite-temperature intrinsic response is therefore tagged as a model-specific diagnostic. This is the **single scientific question** reserved for possible correspondence with Jian Wang, not an issue to suppress with fitted broadening or finite grids.
