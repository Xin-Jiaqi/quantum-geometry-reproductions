# MoSSe S/S biaxial strain

Status: prior material-model hypothesis complete; independent relaxed DFT at −2%, 0, +2% submitted. No strained NMEE result is established yet. Follow the unified manuscript worklog for scientific conclusions.

- `prior_prediction.json`, `prediction.json`: prediction recorded before strained DFT; geometric cross-layer scaling only, not calibrated deformation potentials.
- `strain_point.py`, `point.pbs`: constrained same-registry ionic relaxation, independent zero-B bilayer and geometry-matched isolated-layer SCF, complete retained KS response and smooth electrostatic reference.
- `mechanism.py`, `mixing_channels.py`: selected occupied/empty layer mixing interventions; nonadditive and projector dependent.
- `PLAN.md`: decisions and failure branches fixed before the results.

Run only inside the existing authorized PBS project allocation. VASP and its PAW potentials are not distributed. The exporter reuses the project's existing `spin_transitions.py`, `vaspwfc.py`, and PAW library at the explicit parent path; the prediction requires the previously archived true-zero-field cache. `CASE` and `STRAIN` select each independent case. Original large wave functions remain on the compute resource. This snapshot does not introduce or validate a new response formula.
