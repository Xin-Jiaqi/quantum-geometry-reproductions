# v76 → v77 Handoff | 2026-10-23

**Original Xiang et al. PRB 113 L201406 (2026) status: PARTIAL.** This handoff extends previous v49–v75 source audits and fully tested finite-model methods, without rewriting or purporting to correct the authors' printed formulas.

## Achievements in v76

1. True Bx-shifted equal-energy contours for **both** + and − bands and multiple E shells in tilted 2D Dirac with exact band-dependent DOS. Conduction/valence elastic collision Markov operators satisfy detailed balance, conservation and positive spectra; independent direct/spectral AC and DC tests. The population-only response is invariant under Bx translation if the impurity kernel depends on relative transfer, but the original interband PHE is NOT thereby disproven.
2. A legitimate **two-band interband-coherent** completely positive GKSL with short-range spinor-overlap angular jumps and explicit pure-dephasing sector, restricted to β=0 same-k ±E shells. Both optical E→B and B→E vertex orders, full dense vs spectral Liouville and local band U1 frame covariance. Homogeneous integrated EB signal happens to vanish for this coherent proxy; nonzero angular harmonic quantifies the local response. The proxy's E vertex lacks the source length-gauge diagonal drift, so does not reproduce original SM Eq20 tau pole.
3. **Constructive no-go:** gamma_phi=0 and .14 share precisely the same population transport collision spectrum and tau_tr but give different interband spin and local two-color amplitudes. Scalar transport tau alone CANNOT uniquely identify the full mixed response.
4. **Necessary Bx equilibrium collision Ward:** d_B D_B F0 + D0 d_B F_B=0 is required by magnetic equilibrium. Frozen Bx collision has model stationarity defect ≈0.3467493. Genuine original Bx scattering requires translating/moving equal-energy shells, not a fixed-ring jump approximation.
5. **Sufficient controlled Bz GKSL implementation** with Bz-dependent mass, spinor overlap jump operators and dissipator, exact field-equilibrium stationary F_B. For chosen fixed orbital E vertex and current, independent 5-point full-B derivative agrees with analytic Liouvillian+collision-derivative to 3.7e-11; wrong no-magnetic-collision derivative differs by 0.0796876. Bz and orbital E proxy are source-distinct.

## Scientific limits and v77 required steps

- **Priority 1:** Work in the original source's Bx Zeeman/tilted Dirac k-space and correct length-gauge D_y. Construct **field-translated full two-band elastic GKSL** including the k-moving basis connection, diagonal Boltzmann drift and interband coherence. All Bx response and bath vertices must obey physical continuity and equilibrium Ward. Explicitly distinguish frozen-H0 versus equilibrated-B baths.
- **Priority 2:** Recover SM Eq15 selective electric 1/eta coefficient within this *same* conserving bath, or give a reproducible counterexample to its existence; preserve both EB field orderings and the true epsilon/v/Z printed definitions.
- **Priority 3:** Treat printed Eq20/31 relative sign and main/SM helicities as **two unresolved literal source branches**; do not correct/choose a sign by fitting or through the Bz toy model.
- **Priority 4:** Physical transverse charge transport and magnetization/boundary subtraction, continuum k/q/omega/gamma limit and SI Fig2 remain unproven. Orbital MNHE requires full Wannier magnetic position plus original Eq49/52 factor and m!=l interferences.
- **Release gates:** 52 physics tests + frozen evidence regeneration locally 53/53 PASS; final PR head CI, merge SHA, independent main CI and verified portable zip hashes must be checked and added only after the real results.

**Planning:** paper reproduction not elevated to COMPLETE; numerical tests certify **specific physics constructions**, not paper truth.
