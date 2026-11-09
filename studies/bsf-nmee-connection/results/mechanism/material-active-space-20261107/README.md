# Testing the material content of the twelve-state mechanism

This checkpoint does NOT close the material minimal-model gap. It tests and narrows it, using the existing true-zero-field MoSSe S/S reference. No new DFT, strain, material or stacking scan was run. Scientific work ran only in PBS compute jobs; local work prepared inputs, inspected small outputs and tested algebra.

## Response-independent selection

Each isolated-layer d+p block has 14 occupied and 8 empty spin states. Retain the top two valence states, the lowest two empty states, and two combinations of the remaining six empty states. The latter come from the SVD of their coupling to the OTHER layer's top valence states. The amplitude construction retains that direct coupling exactly. A second, independently specified construction weights the coupling by the high-empty minus mean-valence denominator. Neither uses the sign or magnitude of NMEE to select states.

All Hamiltonians, real-space layer spin, Zeeman operators and Bloch-overlap electric coupling are projected consistently from the existing KS reference. Truncated 96-state basis metrics are explicitly orthonormalized. The new tests check invariance of the selected subspace under arbitrary within-layer basis rotations and preservation of the source coupling. The active Hamiltonian projected through the KS lift agrees with the local construction to about 6.4 micro-eV in the main case.

## Main outcomes

For the window8/96-state reference, bottom C in 1e-9 (nm V T)^-1:

| Model/control | Selected interlayer mixing off | Full strength |
|---|---:|---:|
| Complete material reference | +13.498 | -6.602 |
| Active-state hoppings retained in complete background | +13.498 | -30.460 |
| Direct twelve-state truncation | +23.368 | +28.394 |

Removing the top-valence to selected high-empty hopping from the active-in-background control restores +2.532. Thus the identified hopping really drives reversal in the material background, but the twelve-state Hamiltonian plus directly projected operators does not retain the entire response. Both selection methods and all three orbital windows show this failure. An independently determined constant background also fails: it predicts +18.524 at full coupling.

Actual root-mean-square r-to-low-empty/r-to-valence hopping norms are about 0.61–0.76; r-to-r hopping is appreciable. They do not reproduce the toy model's strong selectivity. These norms are not the toy rho parameter or response weights. A matching uniform 6x6 BZ norm check of the toy gives about 0.017–0.019 at rho=0.1; even rho=1 gives about 0.17–0.19 and does not reverse. Norm ratios alone cannot predict the response sign.

Adding deeper occupied states does not provide a stable truncation shortcut. Keeping the top valence pair plus ALL eight local empty states per layer (20 states total) restores +21.456 to -7.234 in the main representation. But window4 gives positive full-strength C (+3.618/+5.303 for 96/112); windows8/12 give negative C. Three twelve-state and three sixteen-state energy-pair selections also fail to reverse at window8. Therefore neither twelve nor twenty is established as a reliable, minimal material response model. The 20-state apparent magnitude match at window8 must not be marketed as quantitative agreement.

## What the reduced model omits

The twenty-state closure controls complete the selected twelve-space with eight orthogonal high-empty states. In this limited representation, removing active/complement zero-field Hamiltonian coupling still leaves C=-2.325; removing the electric connection gives +36.513. Blocking H, Zeeman and electric connections recovers the bare twelve-state result +28.394. These are diagnostic interventions, not experimental knobs.

More importantly, repeat the interventions WITHOUT deleting any of the 96/112 material states. Fix the zero-field active projector P; Q=1-P. H_block, Z_block and S_block replace the corresponding matrix A by PAP+QAQ. E_block replaces each physical overlap link L(k,k') by P(k)L P(k')+Q(k)L Q(k'), retaining the full Hamiltonian, readout and Zeeman vertex. The electric coupling remains the same covariant finite-difference construction. This is a deliberately changed electric vertex, not a gauge transformation and not a self-consistent material operation.

For window8/96, C is -6.602 (all), -14.006 (H_block), +23.829 (Z_block), -28.237 (S_block), +94.224 (E_block). Across ALL six window/cutoff choices, the original and H_block remain negative and E_block becomes positive. Z_block changes sign between windows; S_block stays negative. Thus the robust finding is that the electric connection to the excluded sector is essential in these controls. Do not upgrade the Zeeman control to an equally robust claim. In the complete space Q includes occupied and empty states: the result is not a unique attribution to one empty orbital, one three-band loop, or a new tensor.

The physically useful distinction is between states selected by their direct zero-field interlayer hopping and states needed for the field-induced response. Eliminating the latter generally changes electric/Zeeman/readout vertices and can generate mixed-field effective terms. A bare truncation retains neither all virtual transitions nor their induced effective operators.

## Implications and next test

The earlier twelve-state KM example remains a valid existence demonstration. Its selectivity and g_c^2/Delta rule are NOT established MoSSe conditions. A proposed next step is a field-consistent elimination of the omitted states, deriving the effective electric and mixed EB terms from the material Hamiltonian and physical operators, with no response-sign fitting. Check whether it restores reversal in a compact model before adding further bands or claiming universality. Symmetry classification should then constrain the actual field-linked orbital paths, not only count allowed tensor components. No new PRL-novelty or experimentally realizable switching claim is made here.

## Files and replay

- material12.json: two response-independent active-space constructions, three orbital windows, full/background/bare controls, norms and orbital weights.
- ladder.json: occupied/empty retention ladder, full-band group diagnostics and readout partition. The large band-group cancellations are not a unique microscopic source attribution.
- closure.json: predefined high-empty energy pairs and twenty-state vertex interventions.
- validate20.json: 96/112 cutoffs and three windows; includes the failed window4 case.
- full_vertices.json: the robust, complete-space electric-connection test. Spin readout linearity is checked to about 1e-10 in reported units; Zeeman finite-difference residuals are below 1.6e-4.
- toy_coupling_norms.json: contextual norm comparison, not a fitted material rho.
- scripts/: exact PBS analyses and the local subspace helper; scripts reuse the existing parent model.py and ../zero-field-refined-20261107/full-controls-cache.npz. No transport helper, secret, PAW data or large wavefunction cache is included. This is not a standalone licensed DFT distribution.

Inherited limits: 6x6 material mesh, approximate d+p representation, and imperfect high empty-state convergence in the 112-state cache. This is a qualitative mechanism audit, not a final response magnitude or gate/strain threshold. All claims distinguish failures and passing controls.
