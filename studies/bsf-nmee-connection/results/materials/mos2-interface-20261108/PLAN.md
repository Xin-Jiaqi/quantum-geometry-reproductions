# Fixed-3R MoS2: electrostatic environment versus electronic mixing

Unchanged BiDB 1930 and constituent monolayer 1929; spin-only, fixed-ion insulating EB response. No new materials, strain, NEB or twelve-state model.

1. Generate genuinely zero-B/E, dipole-corrected references with full spinor overlaps and ionic+Hartree potential. Reference jobs 558139 (bilayer), 558140 (monolayer); PBS DP_002, 4 CPU each, 24 GB, 2 h limit. Existing normal zero-B charge is an initial guess only. Self-consistency retained.
2. Reconstruct full-space response and ensure it has the observed signs before interpreting perturbations. Report empty-state sensitivity, forbidden monolayer response and B-step dependence. Frozen-potential model is not identical to self-consistent finite-field DFT.
3. Hold the actual interface environment fixed and scale its projected Mo d/S p interlayer Hamiltonian block. Preserve full-space links and spin matrices. This only removes the selected mixing: alpha=0 is NOT an exact electrostatic-only double layer.
4. Independently extract the plane-averaged interface electrostatic environment and apply it to the isolated monolayer reference. This is a controlled test of the smooth electrostatic-field explanation, not all microscopic electrostatics, exchange-correlation or self-consistent screening. Compare two smooth cutoffs outside the constituent layer; do not fit the field to NMEE.
5. If local electrostatic embedding reproduces the response, describe electrostatic activation. If hopping removal suppresses it while electrostatic embedding is insufficient, identify hybridization-assisted activation within the tested approximation. If baseline/window controls fail or the residual dominates, leave mechanism unresolved. Never label all intralayer changes as electrostatics.

Possible unified picture to test (not a result): interface-induced change adds to the constituent-layer response: activation from a symmetry-forbidden baseline in MoS2, competition/reversal of an existing response in MoSSe. Additive decomposition is descriptive; a common causal electronic channel still needs evidence.

LVHAR exports ionic+Hartree, not exchange-correlation: https://vasp.at/wiki/LVHAR . Exact discrete Bloch links are from VASP/Wannier90 interface: https://vasp.at/wiki/LWANNIER90 .
