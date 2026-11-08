# MoSSe S/S: relaxed biaxial strain and layer NMEE

Completed: independent fixed-cell, all-ion relaxed DFT at −2%, 0, +2%, with matched isolated-layer and smooth screened electrostatic references. All retain the same C3v stacking and a positive computed gap. See the unified manuscript worklog for the physics discussion.

## Physical result

For fixed positive Ex Bz, lower-layer Cxx is −3.36, −6.92, −7.32 at −2%, 0, +2% (96-state representative; units 1e-9 (nm V T)^-1, Pauli spin density). The matched isolated layer is +19.68, +11.43, +5.31 (64 states). Every electrostatic reference is positive. The upper layer remains positive; no sampled point exhibits a local reversal or total zero crossing. Compression roughly halves the lower negative response, while the total negative response grows because the upper positive response also decreases.

Selected interlayer occupied–empty mixing remains a sufficient route to local reversal at all three strains and all six retained-space settings. A stronger hopping norm does not imply a stronger negative response. Compression also changes the positive isolated-layer response and the relative transition energies. Energy/state exchange diagnostics support opposing tendencies for 0 versus −2%; all cross-reference exchanges involving +2% fail the degeneracy-compatibility guard and are explicitly excluded. No forced band rematching was used.

This supports robust interface rewriting of a layer response, not an achieved strain switch or a universal reversal criterion. The prior hopping-only strain hypothesis was uncalibrated, largely a metric-scaling effect, and did not predict the observed tensile trend. No new gate, larger strain, material, or twelve-state search was launched.

## Reading the data

- `responses.csv`: all 33 actual full/isolated/electrostatic results, including retained-state/reference variants. For the isolated layer, `lower` means the entire matched single layer.
- `mechanism.csv`, each case's `mechanism.json`: all 72 selected-channel controls, plus mixing scans. They are Hamiltonian interventions, not additive or unique mechanism contributions.
- `spectrum.json`: complete attempted cross-reference diagnostics, with invalid pairings preserved. `spectrum_partial.json` records the earlier 0/−2 analysis.
- `verification.json`, `checks.json`, `relaxed_geometry.json`: geometry, symmetry, matched layer structure, insulation, links/norms, and probe-step checks.
- `MoSSe_strain_layer_response.{png,svg,pdf}`, `plot.py`: all three physical points; vertical ranges are setting sensitivity, not statistical error. Lines guide the eye, without claiming intermediate results.
- `prior_prediction.json`, `prediction.json`, `prediction_metric_audit.json`, `PLAN.md`: sealed prior hypothesis, its pre-DFT metric audit, and planned decision branches.
- `jobs.json`, `scheduler_status.json`: actual job lineage and terminal status. A historical failed post-relaxation guard and two obsolete dependency jobs are retained, not concealed as successful runs.

Important raw-output convention: for `mono/response.json`, only index 0 of `C_total_lower_upper` is the physical whole isolated layer. Indices 1 and 2 are auxiliary partitions at the single Mo plane, not two material layers; do not interpret their large opposite numbers as a bilayer signal. The CSV, figure, and worklog use index 0.

## Actual method and recovery

PBE-D3(0), SOC, 450 eV, 6×6 mesh, zero B and zero vertical bias for each self-consistent reference. Biaxial strain is relative to the source BiDB in-plane a=3.250850 A, not a newly optimized zero-stress lattice. Cell fixed; all ionic positions relaxed to the original 0.01 eV/A criterion. The matched monolayer retains the bilayer's intralayer geometry. NMEE is the existing frozen-potential static spin EB response in the retained full KS space. 96/112-state amplitudes are not final absolute convergence; no new self-consistent finite probe-field, phonon, or experimental feasibility validation is claimed.

The initial script intended z-only relaxation, but this ASE exporter wrote T T T in every actual relaxation POSCAR. All three structures therefore underwent all-ion relaxation and retained the same stacking (maximum lateral drift 3.73e-7 A). The old `result.json`/`progress.json` method strings for 0/−2 are unmodified raw provenance; this audit and actual inputs supersede their inaccurate z-only description.

Job 558256 completed tensile relaxation, then failed a dimensionless 1e-8 registry assertion. `continue_relaxed.py` reused that converged structure with the Cartesian 1e-5 A registry check; no force threshold was relaxed and VASP ionic relaxation was not repeated. Its introductory legacy phrase “constrained ionic relaxation” refers to the fixed-cell workflow; actual free-ion flags and the detailed metadata govern the method. `strain_point_v1.py`/`v2.py` preserve executed versions; current `strain_point.py` documents a corrected future entry point. `check_registry.py` verifies real saved roundoff/periodic images are accepted and a 0.01-A displacement is rejected. `continuation_audit.json` and `registry_validation.json` retain the evidence. `relax_inspection*.json` are interim diagnostics, not final structures.

All Japan-side scientific runs and postprocessing used PBS compute nodes. Local scripts `summarize.py`, `verify.py`, and `plot.py` operate on collected small files. Licensed VASP/POTCAR and large wavefunctions are not distributed. Export/PAW helpers and the prior material cache are dependencies at the explicit original project paths; this snapshot does not introduce a new response formula. Future reproductions require the licensed inputs and existing parent project. Raw compute directory: `/work/sunq/xjq/mosse-strain-20261108`.
