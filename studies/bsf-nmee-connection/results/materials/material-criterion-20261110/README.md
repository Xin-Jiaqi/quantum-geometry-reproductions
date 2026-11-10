# KM → MoSSe / MoS₂ / InSe: a tested criterion with a clear boundary

2026-11-10. The bounded three-material iteration is complete. **No universal material criterion has been established.** A low-order signed-channel description works well locally in MoSSe and MoS₂, but fails in InSe, including the direction of some layer changes and the enhancement/suppression of the full interface. This failure is retained rather than repaired by moving the fit window. Stop numerical expansion here for discussion, as requested.

## Physical result

The preceding eight-state KM analysis identified cross terms between interlayer response paths. Here the same selected *class* of channel—lower-layer valence / upper-layer empty mixing—is tested in the existing materials. No new DFT, material, model states, structure or strain was introduced.

| Existing material and control | Low-order prediction for restoring the selected channel | Actual full-interface result | Boundary |
|---|---|---|---|
| MoSSe, internal layer step | upper-layer gate sensitivity increases | increases | local new half-amplitude layer errors <6% |
| MoSSe, screened external gate | upper-layer gate sensitivity increases | increases | local new half-amplitude errors <2.3%; small full lower-layer residual is not accurately predicted |
| MoS₂, internal layer step | upper-layer gate sensitivity increases | increases | local new half-amplitude errors <6.6%; no matching symmetric small-gate SCFs exist here |
| InSe, internal layer step | upper-layer gate sensitivity increases | decreases | at positive half-amplitude even the upper/lower slope signs are wrong |
| InSe, screened external gate | upper-layer gate sensitivity increases | decreases | local slope signs survive, but amplitudes and full enhancement/suppression do not |

“Gate sensitivity” in this table is the magnitude of G=ΔC/Δcontrol, not the magnitude or sign of C itself, and not a percentage divided by a nearly vanishing C. At full restoration, the upper |G| relative to the same-window selected-channel-removed reference is:

- MoSSe internal step: actual2.86–3.34, predicted3.12–3.68.
- MoSSe screened gate: actual6.90–12.11, predicted7.48–13.20.
- MoS₂ internal step: actual11.92–26.11, predicted11.67–27.91.
- InSe internal step: actual0.143–0.159, predicted1.58–1.73.
- InSe screened gate: actual0.315–0.362, predicted2.56–2.76.

Ranges are the four retained band/projection settings, not statistical intervals. Absolute responses and calibrated experimental gains remain unconverged. The known physical α=1 endpoint is a retrospective extrapolation test; it is **not a new blind prediction**. Newly computed ±1/2-channel controls were predicted before their job was submitted.

## Same definition where possible, distinct fields where necessary

Three-material common diagnostic: from each zero-gate SCF cache, apply an internal step

H(U)=H(0)+U W,   W=(I−2P_lower)/2,

where P_lower is the cached real-space half-space operator in the full KS space. U is an energy in eV, evaluated at0 and ±1meV. It is **not** a screened external field. This is the closest common material counterpart of KM's internal layer offset, while retaining the complete material response space and physical vertices.

For MoSSe and InSe, an additional bridge uses the independent fixed-ion SCFs at external fields −0.10,0,+0.10V/nm, with the fixed zero-field selected space transported by the existing PAW overlaps. These slopes include the original self-consistent field-induced change of the reference states and vertices. MoS₂ has no equivalent symmetric small-field data in the archive; its internal-step results cannot be promoted to external-field validation.

C denotes the static Ex Bz-induced layer Pauli-spin density coefficient, in10⁻⁹(nm V T)⁻¹. Internal-step G is in10⁻⁹(nm V T eV)⁻¹. Screened-external G is in10⁻⁹(V² T)⁻¹. Neither absolute slopes nor expansion coefficients are pooled across these two field definitions.

## Fixed candidate and held-out test

Use zero-field separated-layer occupied/empty projectors in the existing selected material orbital space (Mo d/chalcogen p or In s,p/Se p). Scale only the Hermitian lower-valence/upper-empty block by a real α. All other states, Hamiltonian blocks and response vertices remain. Under the internal step, the selected block's U dependence is scaled consistently too. Under screened external fields, scale the transported selected block independently at each existing SCF state. Modified Hamiltonians are **not rescreened**.

At each material/window/field definition, calculate α=−1/4,0,+1/4 and freeze the local expression

G_l(α) ≈ A_l+B_l α+D_l α².

Then predict α=−1/2,+1/2,1 without refitting. There are20 independent window/field groups,40 newly held-out half-amplitude cases and20 previously known or consistency-check physical endpoints. Two band windows and frozen-empty4/8 definitions are retained:96/112 for MoSSe/MoS₂,144/160 for InSe.

This is a **local Taylor hypothesis around a material with the selected channel removed**, not KM's exact weak-total-hopping expansion. Other material channels and unselected states remain strongly coupled. Odd/even terms under changing α therefore cannot automatically be called pure second-order interference/direct contributions. The test asks whether KM's economical response description transfers; it does not assume that it must.

The fixed prediction file was saved and uploaded before the validation job. No training window, polynomial order, channel definition or material was changed after seeing failures. All360 training/validation field-state records are retained, with45 further probe-refinement records. All sampled points remain insulating (minimum gap about0.795eV in the main stages). No unfavorable control was excluded.

## Why InSe is a useful counterexample

Representative largest-window, frozen-8 values:

| InSe control | Upper G predicted | Upper G calculated | Consequence |
|---|---:|---:|---|
| internal step, α=+1/2 | −40.63 | +52.79 | opposite gate-slope signs |
| internal step, α=1 | −854.25 | −70.61 | wrong enhancement/suppression relative to removed-channel G=+492.70 |
| screened gate, α=+1/2 | −21.58 | −10.24 | correct slope sign, poor amplitude |
| screened gate, α=1 | −98.44 | −11.31 | wrong enhancement/suppression relative to removed-channel G=+35.88 |

Use the field-specific units above; these numbers cannot be compared as if all were slopes per physical V/nm. Across the four windows the same physical failure persists. Halving both the magnetic probe and internal-step interval preserves every tested slope sign; the maximum layer slope change is0.0114%. The counterexample is not removed by this numerical refinement.

Thus the material response over this channel-restoration interval contains important behavior beyond the second-order description. This observation does not identify a unique microscopic origin, prove an inevitable strong-coupling law, or invalidate the existing common-channel intervention evidence.

## Exploratory physical applicability diagnostic

After observing the failure, evaluate an energy-weighted mixing measure in the selected-channel-removed full-space background. For its occupied and empty eigenstates define M_vc=(K_selected)_vc/(ε_c−ε_v); record the largest singular value over k. Exact within-degenerate-sector rotations do not change this value. This estimates the size of first-order state mixing, not the signed spin response.

| Representative material | Background minimum direct gap | Largest energy-weighted mixing η |
|---|---:|---:|
| MoSSe | 1.391eV | 0.154 |
| MoS₂ | 1.538eV | 0.143 |
| InSe | 0.882eV | 0.305 |

The bare separated-layer version also ranks InSe higher (0.254 versus0.145/0.134). This supports stronger state mixing as a plausible reason a short expansion is less reliable in InSe. **It does not establish resonance, a closing gap, a universal η threshold, or a causal explanation of the signed response.** It is a post-hoc descriptor in one representative window per material, not an independently validated classifier. Different vertex cancellations can also make response derivatives much more delicate than state-mixing norms.

The screened change of the projected upper-minus-lower mean layer energy is positive in both MoSSe and InSe: about0.094 and0.101eV per external V/nm, respectively. In the InSe α=+1/2 control, the internal-step upper slope is positive while the screened-external slope is negative. Since the mean layer offset shifts in the same direction, a positive scalar conversion of external gate to a rigid layer offset cannot explain this example. Real gate-induced changes beyond the average offset must be retained, including the reference-state and vertex changes already present in the SCFs. This is a conditional model-control example; it is not a proof of one uniquely identified screening or intralayer-polarization contribution, nor of a sign discrepancy at the physical α=1 state.

## What remains defensible for the manuscript

Retain: interfaces rewrite signed layer responses and redistribute which layer is electrically tunable; channel controls show that baseline direction and control direction are distinct. KM supplies a compact cross-path mechanism and genuinely tested weak-coupling predictions. MoSSe/MoS₂ locally support an economical low-order material description; InSe exposes the boundary and the need to retain more complete electronic/gate response.

Do not yet claim: one universal weak-coupling formula for all three materials; an energy-only criterion; an η threshold; a universal interlayer-hybridization amplifier; or a realistic operation that independently flips a selected hopping phase. No extra states were added to rescue the failed rule. We now stop for the requested discussion of the physical axis, manuscript and next directions.

## Reproducibility / continuation

- `PLAN.md` records the original test and the single explicitly bounded diagnostic refinement. `predictions.json` hashes `train.json` and predates validation submission. `summarize.py` never refits predictions using validation.
- `train.json`, `validate.json`, `refine.json` retain all inputs' SHA256 and job/node identities. `all_controls.csv`, `prediction_tests.csv`, `summary.json` retain all controls, failures and errors. `final_checks.json` includes refinement and the exploratory mixing diagnostics.
- The selected-channel-removed screened baselines reproduce the previous archive exactly. Physical screened gate slopes agree with the independently implemented preceding kernel within1.3×10⁻⁹. Layer sums hold; no skipped or invalid points in this run.
- Japan compute jobs:563995(train),564011(validate),564028(refine), all exit0 on csuper010;2m22s,2m17s,50s. Scientific processing and hashes ran only inside PBS compute nodes. Existing caches and PAW overlaps were reused; no DFT or wavefile migration. Keeper updated. No further jobs are planned before discussion.
- Remote reproduction requires the recorded private scientific caches and the existing licensed environment. Run `run.pbs` with STAGE=train/validate/refine on a compute node; **do not run calculate.py on a login node**. Local prediction/evaluation: `predict.py`, `summarize.py`, `check_final.py`. The archived response helper is unchanged from the preceding material implementation.
- New chats should start at repository `AGENTS.md` and the unified worklog's continuation entry. Files preserve verified state; chat memory, external caches and available compute credentials are not automatically inherited or guaranteed.
- The core model and GitHub workflow/notification settings are unchanged. Research archive uses `[skip ci]`; no remote CI pass is claimed.
