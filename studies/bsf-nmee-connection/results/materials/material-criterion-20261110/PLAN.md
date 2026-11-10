# Bounded KM + three-material criterion test, then discussion

2026-11-10. User requests a candidate material criterion and a pause for physics/manuscript discussion after the KM + MoSSe/InSe/MoS2 iteration. Do not add materials, DFT, strain, fitting states, or automatic further exploration.

Test whether the low-order signed channel dependence exposed in eight-state KM predicts material gate trends. Use the existing fixed separated-layer lower-valence/upper-empty channel. Train a local quadratic in its signed amplitude alpha at -0.25,0,+0.25. Save coefficient predictions before computing new alpha=-0.5,+0.5 states; alpha=1 is a physical-bilayer consistency/extrapolation check, already partly known from past work, NOT a blind prediction. Do not move the training window after seeing failures.

Common three-material diagnostic: frozen zero-gate Hamiltonian, add internal layer-step U*(I-2P_lower)/2, U=0,+/-1meV. This common field is NOT the screened external gate. Remove/scale the same selected channel consistently in H and in the step vertex. Use both existing KS windows (MoSSe/MoS2 96/112, InSe144/160) and frozen-empty4/8 projection definitions.

Bridge to actual screened references: for MoSSe and InSe additionally perform the same channel test on existing independent SCFs at zero,+/-0.10V/nm, with the fixed zero-field selected space transported by existing PAW overlaps. MoS2 does not have matching symmetric small-gate SCFs; do not label the layer-step calculation as an equivalent external-field validation.

All other Hamiltonian blocks and physical vertices retained. No rescreening of artificial controls. This is a local alpha Taylor surrogate around a strongly mixed background; unlike isolated eight-state KM, the selected embedded subspace is incomplete and no global t^2 or pure-interference interpretation is assumed. If the quadratic or transfer of direction fails, report its range of failure; do not turn a fitted coefficient into a universal energy-only rule.

Primary metrics: G signed and |G|, actual C direction, insulating gap, all-window stability, held-out sign and amplitude errors, difference between internal step and screened external response. Log any invalid/gapless point without deleting it silently. One bounded refinement only for any materially ambiguous result.

Remote processing uses already-present caches/overlaps in Japan, avoiding wavefile migration. Scientific Python, hashes and postprocessing only inside PBS compute allocation. Request4CPU/16GB/30min per stage; new output<0.1GB. Group/project balance unknown; user says Japan free. Shared unrelated jobs unchanged.

After fixed-prediction validation: MoSSe/MoS2 half-amplitude tests retain layer signs; InSe has opposite predicted/actual signs in the internal-step +half case and fails full-interface enhancement/suppression. Use the one bounded refinement: highest-band frozen-8 settings, halve B and internal U step, keep all validation alphas. Also record an exploratory energy-weighted mixing diagnostic at zero gate (bare-layer and selected-channel-removed full-space occupied/empty backgrounds) and the signed screened layer energy offset. These descriptors were requested AFTER seeing the failure and cannot be advertised as a blind universal classifier. Do not refit the polynomial or add further materials/jobs after this check.
