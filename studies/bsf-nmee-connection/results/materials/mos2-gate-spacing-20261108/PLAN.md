# Fixed-registry MoS2: gate versus separation

Before-result plan (2026-11-08). Keep the BiDB 1930 3R A registry, in-plane cell and each monolayer geometry fixed. Zero-probe-field self-consistent references at (spacing increment Å, physical vertical field V/nm): (0,-0.50), (0.15,-0.50), (0.30,-0.50), (0.30,0). Reuse the existing (0,0) reference. Negative field is chosen because the earlier isolated monolayer response is positive for this sign, opposing the negative lower-layer interface response. The -0.50 field is an exploratory bound, not a claim about experimental device tolerance.

The existing hopping multiplier is not calibrated to distance; do not convert it to Å or pressure. 0.15 and 0.30 Å are modest geometric probes relative to the 3.038 Å interface gap, not demonstrated accessible equilibrium displacements. Add only targeted points if a sign bracket or field trend warrants them. No NEB, relaxation, strain scan, new material or twelve-state fit.

For each reference use the existing full-KS frozen-potential spin-only NMEE model at 96 and 112 states and b=0.002 eV; check b=0.001 at sign-sensitive points. Require a positive indirect gap throughout. Record lower, upper and total separately. A total zero is not evidence of a lower-layer reversal. The gate is self-consistent, whereas the infinitesimal probe E/B response is frozen-potential. All VASP and wavefunction processing execute in PBS compute allocations only.

Interpretation limits: PBE fixed reference geometries, coarse 6x6 integration; no converged absolute signal, experimental pressure mapping, or proof that a gate leaves the stacking pinned. If no lower-layer reversal near d0 is found, state that negative result instead of extrapolating a critical field.
