# Stage 8c — external 3R-WSe2 relaxed bilayer coordinate acquisition audit

**Status: BLOCKED_NO_VERIFIED_RELAXED_BILAYER_FILE (2026-10-20).** This is an original-research material-source/provenance audit, NOT a Jia et al. PRL reproduction, NOT an actual-material DFT run, and NOT a measured/calculated real-material NMEE. Prior Stage8b P1 CIFs remain explicitly IDEAL_UNRELAXED.

## Audited source candidates

| Source | Verified publicly described content | Exact desired original relaxed AB/BA CIF or POSCAR inspected? |
|---|---|---|
| Pakdel et al., [Nat. Commun. 15, 932 (2024)](https://doi.org/10.1038/s41467-024-45003-w), BiDB [DTU record](https://doi.org/10.11583/DTU.24849819) | 8,451 bilayer stackings studied; published database link claims low-energy structures/properties | **NO**. The target record and six-atom file could not be obtained and independently checked in this pass. |
| Barik & Woods, [Sci. Data 10, 232 (2023)](https://doi.org/10.1038/s41597-023-02146-7), [BMDS](https://doi.org/10.6084/m9.figshare.21799475) and [public GitHub](https://github.com/rkb12/BMDB-database) | 760 relaxed bilayers; SOC, gap, structural parameter statistics, structure-generator code and figures | **NO**. A property spreadsheet or generation example is NOT proof of a particular published relaxed WSe2 AB atomic motif. |
| Li et al., [Nat. Commun. 17, 2457 (2026)](https://doi.org/10.1038/s41467-026-68854-x), [public source data](https://doi.org/10.17861/2e4207fe-a6c0-402c-8468-97839c74c36b) | Experimental AB/BA ferroelectric signatures in R-WSe2 and experimental electronic structure; accompanying source-data XLSX | **NO**. This source-data workbook is not yet verified to contain a DFT-relaxed six-atom CIF. |
| [MoS2 AA / AB-2H / AB-3R Wannier-based SHG dataset](https://zenodo.org/records/21480888) | Metadata expressly names AB-3R/POSCAR and Wannier90 input/output packages | **NO**, currently listed as **Restricted**. Respect access permissions. |
| [MoS2 microscopy structure repository](https://github.com/HyunPhys/DislocationNetworksInMarginallyTwistedBilayerMoS2) | README declares bilayer and domain-wall atomistic structure input files for TEM simulation | **NO** proof of an *ab initio* relaxed six-atom 3R-AB geometry, despite availability of simulation structures. |

A CIF from **bulk** 3R-WSe2 (periodic R3m) is not a relaxed isolated two-layer slab. A cleavage model would need re-relaxation and its own provenance. Preserve explicit M-on-X versus X-on-M registry: 3R naming alone is ambiguous for a bilayer. Compare [Di Berardino et al., PRB 111, 125403 (2025)](https://doi.org/10.1103/PhysRevB.111.125403) atomic registry notation.

## A physically important broad-BZ caution

Experimental R-stacked WSe2 spectroscopies do not imply that both the conduction and valence extrema, let alone the static NMEE integral, are always restricted to K/K'. For example, Li et al. ([2026](https://doi.org/10.1038/s41467-026-68854-x)) associate conduction/valence edges with Lambda/K in their setting; a different quasi-free-standing epitaxial 3R-WSe2 [ARPES/DFT study](https://doi.org/10.1103/xs47-7hky) finds a valence maximum at Gamma. These are **different experimental environments, not contradictory estimates of one identical geometry**. They jointly motivate a future **full-BZ Gamma/Lambda/K hotspot survey**, not extrapolation of the synthetic eight-band K/K' valley dominance.

## Source authenticity gate (fail closed)

A candidate must provide: original uploaded bytes and SHA256; authoritative archive URL and publication DOI; direct evidence that these are relaxed **bilayer** rather than bulk/monolayer or a model; full atomic sites/lattice/vacuum; explicit W-on-Se / Se-on-W stacking convention; and coordinates suitable for the pinned group-theory toolkit. A program-readable geometry file can at most be `READY_FOR_MANUAL_SOURCE_REVIEW_NOT_VALIDATED`; manually entered DOI/URL metadata does NOT authenticate source provenance. Full space-group symmetry determination and relaxation claims require separate checks.

**Local standalone Stage8c artifact**: `BSF_NMEE_STAGE8C_SOURCE_AUDIT_20261020.zip` contains a conservative six-site P1-CIF provenance preflight (`structure_preflight.py`), three passing tests, a Chinese report, JSON source ledger and SHA256 manifest. Six prior *ideal* CIFs are all correctly labelled `MODEL_ONLY_NO_RELAXED_STRUCTURE_CLAIM`. Its strict-format parser is **not a general CIF symmetry-expansion library** and has not validated a real author source. No fabricated CIF or fictitious structure URL was committed.

## Next priority

Obtain the real WSe2 target coordinates from **BiDB** (then BMDS or an explicit author-provided POSCAR) and perform actual complete-motif symmetry audit with the pinned group-theory toolkit before any band/Wannier/response DFT. Current scientific status remains **BLOCKED**, not a passed material gate. Prior frozen Jia paper reproduction stays unchanged.
