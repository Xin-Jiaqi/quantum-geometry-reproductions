# Stage 6 work log

- Starting point: frozen stage 5, with weak-coupling interference `A s²+B c²+C sc` and fixed-slide sign reversal. Original PRL reproduction files unchanged.
- Implemented a full 2×2 in-plane static NMEE tensor (`spin x/y`, `field Bz`, `electric x/y`), evaluated using independent finite-field occupied-to-empty mixed curvature.
- Found **exact model mirror** at zero slide independently for same- and cross-sublattice hopping. Mirror residual < 4e-18 eV in sampled generic k, proving the xz;x coefficient forbidden there.
- Derived C3/m_z selection rule: trace of 2×2 tensor must be odd and C3-invariant, hence begins cubic for analytic responses; traceless symmetric sector can begin linearly.
- Tested an independently rotated sliding Hamiltonian: C3 covariance at ~3e-11 absolute; angle +60 gives trace sign reversal.
- Validated weak-slide scaling, performed dense valley and magnetic step convergence, and implemented character-weighted symmetry projection for two physically defined hopping perturbations.
- No claim of orbital-character assignment for a real material, no sliding energy landscape, no response-protocol resolution for PRL driven DC.
