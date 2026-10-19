# Provenance and archival policy

- Independent analysis and code accumulated over rounds 1–32 on 2026-10-19. Original 33 archived ZIPs and SHA256 manifest are kept in the separate offline full research archive `NMEE_round32_complete_audit_covariant_20261019.zip`.
- This source tree is a **curated software release**, not the complete internal worklog, and intentionally omits original paper PDFs, publisher figures, intermediate scratch data and private model files.
- The release preserves the mature static, finite-frequency, PT, BGK, Wannier, and Berry-connection modules and their existing tests, with corrected relative imports, symmetry-safe CuMnAs k sampling in the new CLI, explicit protocol metadata and a Stage-1 handshake validator.
- Curated results under `results/` are regenerated using this release's `scripts/reproduce.py`, not copied from the paper's authors.
- Scientific limitations are not patched with empirical per-figure scaling.

The official work and all claims in PRL 136, 226302 belong to the named authors. Any third-party PDF or private Wannier model should remain outside Git history.
