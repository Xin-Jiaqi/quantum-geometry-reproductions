# Independent follow-up research inspired by Jia et al. (2026)

> **This directory is NOT part of the Jia et al. PRL reproduction or its scientific acceptance evidence.** It is a navigation-only index for original follow-up work carried out separately in `studies/bsf-nmee-connection/`.

## Paper vs follow-up: separate scientific claims

- **Published source:** Jia, Xiang, Qiao, and Wang, *Physical Review Letters* **136**, 226302 (2026), [DOI: 10.1103/r3x2-kb93](https://doi.org/10.1103/r3x2-kb93). Source-paper reproduction and its unresolved *driven-DC versus static* distinction remain in [the frozen paper module](../README.md). Its verification status is **PARTIAL**, not superseded by our follow-up tests.
- **Follow-up:** [BSF × NMEE original research index](../../../studies/bsf-nmee-connection/RESEARCH_INDEX.md). The added synthetic **sliding bilayer Kane–Mele** calculations study an **equilibrium/static**, gapped, zero-temperature Zeeman derivative of mixed spin-electric Berry curvature. This does **not** reproduce the paper's driven-DC response or predict an experimentally realized material.
- **Current status:** Original research **MODEL-VALIDATED / MATERIAL-NOT-TESTED**. Distinguish independent theory/symmetry identities from model-dependent signs and numerical examples.

## Follow-up topics (not claims from the PRL)

| Topic | New original finding | Research source |
|---|---|---|
| Sliding symmetry and valley low-energy expansion | At a C3 / domain-mirror reference, trace begins at cubic displacement order while nematic response may be linear | [Stage 7 / 7f](../../../studies/bsf-nmee-connection/docs/STAGE7F_IRREP_MECHANISM.md) |
| Orbital channel interference | Same/cross hopping drives a negative nematic while cross-only hopping leads to a positive scalar in the tested model; signs are not universal | [Stage 7e](../../../studies/bsf-nmee-connection/docs/STAGE7E_ANALYTIC_PROJECTOR_PHYSICS.md) |
| Remote virtual states | Four-band-only curvature can have incorrect sign; precise eight-band static response contains remote-band compensations | [Stage 7g](../../../studies/bsf-nmee-connection/docs/STAGE7G_SPECTRAL_IRREP_PHYSICS.md) |
| Exact response downfolding | Schur/Feshbach response requires dressed vertices and additional remote/contact contributions | [Stage 7h](../../../studies/bsf-nmee-connection/docs/STAGE7H_RESPONSE_PRESERVING_SCHUR.md) |
| Occupied-bundle mixed quantum geometry | Non-Abelian Wilson loop independently reproduces the gapped static Kubo response on synthetic valley patches | [Stage 7i](../../../studies/bsf-nmee-connection/docs/STAGE7I_GEOMETRY_DFT_BRIDGE.md) |
| DFT/Wannier gauge gate | Under k- or spin-source-dependent orbital reembeddings, **correct physical overlaps and covariant source derivatives** are required; naive overlaps can create a spurious response sign flip | [Stage 7j](../../../studies/bsf-nmee-connection/docs/STAGE7J_WANNIER_GAUGE_PHYSICS.md) |
| Public structure acquisition handoff | Official BiDB `bidb.db` reference and metadata-preserving extractor; pending actual author-relaxed WSe2 AB file, explicitly distinguishes z-scan from full relaxation. | [Stage 8d](../../../studies/bsf-nmee-connection/material-screening/STAGE8D_DOWNLOAD_HANDOFF.md) |
| Public 3R-WSe2 source audit (not NMEE reproduction) | Original relaxed AB/BA coordinate file not verified; database/source provenance and a full-BZ Gamma/Lambda/K caution are documented. Status BLOCKED. | [Stage 8c](../../../studies/bsf-nmee-connection/material-screening/STAGE8C_PUBLIC_SOURCE_AUDIT.md) |
| Ideal TMD stacking full-atom symmetry | Constructed **unrelaxed** 6-atom WSe2 and MoS2 show AA D3h, AB/BA C3v with domain-exchanging mz; **not** an actual CIF/DFT NMEE result | [Stage 8b](../../../studies/bsf-nmee-connection/material-screening/STAGE8B_IDEAL_STRUCTURE_SYMMETRY.md) |
| Literature candidate prescreen | 3R-WSe2 primary and 3R-MoS2 control for later structural/symmetry audit; **no DFT or real-material NMEE calculated** | [Stage 8a prescreen](../../../studies/bsf-nmee-connection/material-screening/STAGE8_CANDIDATE_PRESCREEN.md) |

## Maintenance contract

**Each new research stage** with a validated theory or numerical mechanism gets (i) a standalone source/derivation under `studies/bsf-nmee-connection/`, (ii) a physics scope and independent evidence, (iii) separate regression/CI, and (iv) an update of this research index. No paper figure, original author formula, paper-test baseline, or release hash is altered to support an original-research claim.

Current conclusions are conditional on a synthetic model. Stage 8a is literature-only candidate navigation; no DFT run or real-material NMEE prediction is claimed. The paper reproduction remains frozen and separate.
