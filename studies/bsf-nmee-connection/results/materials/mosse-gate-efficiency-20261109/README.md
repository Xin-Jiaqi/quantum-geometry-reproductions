# MoSSe S/S: interface-enhanced gate sensitivity of NMEE

Model-only result, 2026-11-09. Reuses the unstrained, relaxed, zero-field S/S bilayer (BiDB 7269) and matched lower-layer reference from `mosse-strain-20261108/e000`. No new DFT, relaxation, strain, material, or fitted small model was run.

## Physical conclusion

At the same imposed gate potential, the full bilayer has about **30% larger net NMEE gate sensitivity** than the independent-layer or smooth electrostatic references. The upper layer accounts for about **78%** of the net change, instead of the approximately equal changes in the references. The full bilayer remains lower-negative, upper-positive and insulating throughout the sampled model interval. This is moderate control enhancement, not a same-layer switch or giant resonance.

The additional mechanism check supports a layer-alignment contribution: changing the relative layer energies strongly changes the response when cross-layer occupied–empty mixing is present. Removing the selected mixing makes this channel much weaker. Most of the total gate response still exists in isolated layers; it must not all be attributed to interface mixing.

## Observable and result

`C_l` denotes the layer Pauli-spin density coefficient in `delta s_x,l = C_l E_x B_z`, with the convention and spin-only Zeeman coupling of the existing material response code. `C_net = C_lower + C_upper`. Units of C: `1e-9 (nm V T)^-1`. Gate gain is the **absolute derivative** `G_l = dC_l/dF_model` near zero, units `1e-9 (V^2 T)^-1`. No percentage uses the almost compensated C_net(0) as a denominator.

| Model | G_net | G_lower | G_upper |
|---|---:|---:|---:|
| Full bilayer, 96–112 KS states | −15.56 to −15.66 | −3.37 to −3.36 | −12.19 to −12.31 |
| Matched isolated pair, 64 states/layer | −11.85 | −5.93 | −5.93 |
| Electrostatic pair, 64 states/layer, both cutoffs | −11.91 to −11.98 | −5.95 to −5.99 | −5.95 to −5.99 |
| Selected cross-layer V–C removed, all six settings | −11.79 to −12.78 | −5.73 to −6.01 | −5.82 to −7.05 |

Ranges reflect basis/projection/reference choices, not statistical errors. Smaller 48/56-state reference results are also retained in the raw tables and plotted ranges. Across all reference settings the enhancement is approximately 30–38%; the higher-cutoff comparison gives approximately 30–32%.

### Why a layer-energy shift matters

As a diagnostic, split the **gate operator**, not the NMEE into additive microscopic mechanisms:

`V_linear = V_step + V_residual`, with `V_step = d_Mo/20 * (I - 2 P_lower)` in eV per V/nm, `d_Mo = 6.038885 Å`. This explicitly chosen split uses the existing physical layer projector. Its tangent response is additive; the microscopic channel interventions themselves are not additive or unique.

At gate increments ±0.01 V/nm:

| Model | G_net from whole gate | From layer offset | From residual profile |
|---|---:|---:|---:|
| Full 96 | −15.5592 | −5.2326 | −10.3267 |
| Full 112 | −15.6628 | −5.2143 | −10.4485 |
| No selected V–C, six settings | −11.7860 to −12.7799 | −0.9090 to −1.3342 | −10.8644 to −11.6318 |

Thus the interface adds sensitivity to layer alignment; it does not simply multiply the isolated-layer Stark response. The selected occupied–empty interlayer gaps move in opposite directions under gate bias, and the corresponding virtual admixtures change. This is a change of electronic mixing, not a claim of gate-created free carriers. Nor is it necessary for the bare hopping itself to increase. The decomposition is not a pair of independently realizable external fields.

## Field, references, and intervention boundary

- `F_model` adds **one-electron energy +F_model*z/10 eV** with z in Å, centered in the vacuum cell. It is a model potential slope in V/nm, **not a self-consistently screened applied device voltage**. Do not directly equate its sign/value to the legacy VASP `EFIELD` convention or to earlier gated DFT without conversion.
- The same real-space linear profile acts inside layers in every comparison; independent layers have their own direct Stark response. Replacing the whole field by a rigid layer offset would miss most of the response.
- Lower reference matches the actual bilayer intralayer geometry. Upper reference is its horizontal mirror: actual upper/lower intralayer mismatch is only **0.000691113 Å**, but no separately exact upper-layer DFT is claimed. The upper response follows `C_upper(F) = -C_lower(-F)` for this mirrored reference.
- Electrostatic references use the lower and reflected-upper planar ionic+Hartree environment, subtract the isolated potential, smooth by 0.20 Å and taper at 3.2/3.8 Å. This is an explicit, nonunique static reference, not a self-consistent nonhybridizing material.
- Selected 44 d+p states are embedded in full 96/112-state KS spaces. Each layer has 22 selected states, 14 selected occupied states. Fixed zero-gate valence projectors define VV, CC, and cross-layer VC channels. Three frozen-empty choices (4, 8, 12) check projection dependence. Both H0 and the gate operator undergo the same channel removal. Full physical response vertices/links and remaining KS states are retained.
- “No selected hopping” is not literal isolated layers or pure electrostatics, because other states and vertices remain. No unique percentage decomposition of microscopic NMEE is asserted.
- Same meshes and frozen electronic-response assumptions as the preceding strain materials model; these calculations do not establish complete k-mesh, screening or device convergence.

## Evidence and reproducibility

`scan.json`: 407 accepted curve rows, no invalid rows. All model/reference points remain insulating; minimum gap 0.627551 eV (full bilayer minimum 0.627980 eV). The half-Zeeman-step maximum C change is 6.16e-5 in the tabulated units. Constant-potential PAW operator errors are below 2.7e-6; isolated constant-energy-shift invariance is at ~1e-9. Gate differences ±0.025/0.05 give consistent near-zero slopes.

`gate_direction.json`: 16 diagnostic rows, complete; tangent additivity residual below 2.52e-4 over ±0.01/0.02 gate increments. These residuals are numerical checks, not additional physics.

`curves.csv`, `slopes.csv`, `summary.json`: unrounded tabulations and summary. `prior.json`/`PLAN.md`: hypothesis and scope sealed before outcomes. `operators.json`: operator provenance and reference checks. `scheduler.txt`/`jobs.json`: compute-node completion evidence. `gate_direction.py`, `prepare_operators.py`, `scan_gate.py`, `analyze.py`, `mixing_channels.py`: scientific sources.

Existing source caches and licensed wavefunctions/potentials stay at `/work/sunq/xjq/mosse-strain-20261108/e000`; generated operators are in `/work/sunq/xjq/mosse-gate-efficiency-20261109/operators.npz`. PAW helper imports follow the existing `/work/sunq/xjq/mosse-7269-mechanism-20261104` implementation. A fresh scientific run requires those source data and the existing numpy/scipy/ASE/vaspwfc/PAW environment. Submit `run.pbs` and then `direction.pbs` through PBS; no scientific processing on login nodes. Local summaries/figures can be rebuilt from the archived JSON/CSV using `summarize.py` and `plot_gate.py`.

Both jobs completed on csuper001, DP_002: main 559753 (3 min 19 s), diagnostic 559757 (50 s), exit status 0. No new DFT was launched.

## Figure caption

**Interface-enhanced electrical control in MoSSe S/S.** (a) Calculated net-response change from zero gate, comparing the full bilayer, matched isolated-layer pair and electrostatic pair. The higher-state-cutoff representatives are emphasized; thin lines show other tested settings. Each marker is a calculated point, and joining lines guide the eye. (b) Signed near-zero gate slopes of lower, upper and net responses; removal of selected cross-layer valence–empty mixing weakens the enhancement and upper-layer dominance. Marker centers and spans show midpoints and ranges across tested settings, not statistical errors. Duplicate full results across projection windows are omitted because the full Hamiltonian has no projection intervention. `F` is an imposed one-electron slope, not calibrated external gate voltage. Figure sizes 183×73 mm, SVG/PDF editable text and PNG 600 dpi. Visual inspection passed; source validator reports only that no TIFF submission raster was made; vector PDF/SVG and PNG preview suffice for this worklog figure.

## Consequence for the PRL candidate story

This supports a concrete link between “interface rewrites layer response” and “gate controls the net response”: cross-layer mixing makes layer alignment an additional control channel. It goes beyond a change of cancellation alone, but a ~30% frozen-potential gain is insufficient by itself to establish a PRL-level claim. No dramatic tunable resonance was found over the sampled range.

The next economical independent check, to discuss before launching, is a small symmetric positive/negative gate pair at the **same fixed S/S geometry**, with matched isolated-layer controls and self-consistent screening. Test whether the bilayer retains stronger net sensitivity and upper-layer dominance under the same *applied* field. If screening erases the gain, retain the interface-rewritten susceptibility as established physics and do not claim improved external gate efficiency. No new DFT has been queued.
