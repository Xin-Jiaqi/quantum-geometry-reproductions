# BSF × NMEE: Independent original-research index

**Separate from paper reproduction.** This study builds on the independently reproduced model kernels of [Jia et al., PRL 136 226302 (2026)](../../papers/jia-2026-nonlin-magnetoelectric-edelstein/), but its physical conclusions refer to a NEW synthetic sliding bilayer, not results reported in that article. The published driven-DC NMEE protocol and our gapped 0 K **equilibrium/static** response remain distinct.

Corresponding source-paper page: [Paper-associated independent follow-up index](../../papers/jia-2026-nonlin-magnetoelectric-edelstein/research/README.md). This link is navigation only: no changes to the frozen paper `nmee_lab`, reproduction tests, published-source notes, or historical baseline.

| Stage | Distinct **original** physical or theoretical update | Evidence |
|---|---|---|
| 7a | Exact 8-band K/K' Taylor model resolves sliding sign reversal and same/cross channel interference | [Physics](docs/STAGE7_PHYSICS.md) |
| 7b | Exact middle-four-only static interband response can have the wrong sign; remote virtual bands necessary | [Spectral audit](docs/STAGE7B_SPECTRAL_RESPONSE.md) |
| 7c | Negative same/cross interference competes with positive cross-channel response; neither sign universal | [Channel physics](docs/STAGE7C_CHANNEL_PHYSICS.md) |
| 7d | Gaussian sliding derivative tied to interlayer velocity; chirality at K/K'; Rashba sign reversal | [Valley vertices](docs/STAGE7D_VALLEY_VERTEX_PHYSICS.md) |
| 7e | Weak-hopping coefficients obtained directly from occupied projector derivatives; distinct velocity pathways | [Analytic A/B/C](docs/STAGE7E_ANALYTIC_PROJECTOR_PHYSICS.md) |
| 7f | C3 plus reference mirror constrains scalar trace to cubic or higher, nematic may be linear | [Irrep tensor physics](docs/STAGE7F_IRREP_MECHANISM.md) |
| 7g | Large low/remote virtual-band cancellation can be concealed by apparently converged valley-integrated total | [Radial sector audit](docs/STAGE7G_SPECTRAL_IRREP_PHYSICS.md) |
| 7h | Exact Feshbach downfolding of static response needs energy-dependent spin/velocity vertices and remote contacts | [Schur response](docs/STAGE7H_RESPONSE_PRESERVING_SCHUR.md) |
| 7i | Independent non-Abelian mixed Wilson-loop versus Kubo verification; symmetry and preliminary DFT gates | [Mixed geometry](docs/STAGE7I_GEOMETRY_DFT_BRIDGE.md) |
| **7j** | **NEW: orbital/Wannier gauge must be aligned over *both* momentum and spin-source axes; omitting cross-parameter basis overlap can fabricate even the sign of static NMEE** | [Gauge derivation](docs/STAGE7J_WANNIER_GAUGE_PHYSICS.md) · [Code](bsf_nmee/wannier_gauge_links.py) |

## Scientific status and next decision

- **Validated within synthetic gapped T=0 static protocol:** source-independent occupied-bundle mixed Berry curvature, tensor/symmetry identities and channel-resolved numerical benchmarks with explicit assumptions.
- **Not established:** real material, ferroic stability, actual spinor Wannier operator provenance, generic signs or quantitative driven-DC intrinsic NMEE.
- **Next material gate:** first audit a realistic structure for symmetry and a *same-Hilbert-space* spinor Wannier package supplying H, Sxyz, physical basis overlaps/connection and a full BZ gap. Only after the operator and gauge gate passes should a DFT response be claimed.
- **Reproducibility:** every new stage adds its own `docs/STAGE*`, module, regression tests and machine-readable evidence. The paper's frozen sources and tests remain unchanged. GitHub CI success is a software/model gate, not acceptance of all physical claims.

Updated 2026-10-20.
