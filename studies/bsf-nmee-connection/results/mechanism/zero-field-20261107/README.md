# MoSSe S/S A: zero-field reference and response controls

This is a research checkpoint, not a fully converged gate-voltage prediction.

- Observable: bottom-layer static spin response to Ex Bz; values C in 10^-9 (nm V T)^-1. The Zeeman perturbation is b sigma_z with b in eV.
- `initial`: self-consistent BEXT=0 reference from scratch, 112 bands, 6×6 k points. Finite-B wavefunctions and density were not reused.
- `refined`: additional fixed-charge iterations using that zero-field density; new overlaps and operators were exported from the new wavefunctions. `controls96-analysis.json` projects this reference and all interventions onto the first 96 states before evaluating response.
- Both 96/112-state controls and all three d+p windows retain +,+,-,- in the weak/strong-state versus weak/strong-energy comparison. This supports electronic-state reorganization in this material; it is not a unique additive decomposition into energy and geometry.
- U raises the upper layer relative to the lower: H_U=U(1/2-P_bottom). U is an internal energy offset, not an external gate field or voltage. No gate screening or field-induced intralayer changes are included.
- Bottom C stays negative across the scanned U range [-0.6,0.6] eV. Net C crosses zero between U=-0.15 and 0 while the layers remain oppositely directed there.
- Important limitation: the 112-state representation still has a spurious B-even electric-spin response. The 96-state bottom response has a much smaller residual and retains the mechanism. Absolute layer amplitudes are not converged; the 64-state projection is particularly insufficient. Do not use these data to fit a critical experimental voltage.

Scientific execution was entirely on PBS compute nodes. Job identities and checks are in `evidence.json`. The authoritative interpretation is in the unified manuscript worklog.

## Replay / dependencies

The scripts run in the existing remote material project `/work/sunq/xjq/mosse-7269-mechanism-20261104`, in the `zero-field-20261107` and `zero-field-refined-20261107` subdirectories. They require the existing parent `model.py`, `spin_transitions.py`, `vaspwfc.py`, and `pawlib`, plus NumPy, SciPy, ASE and the licensed VASP installation. `refine.pbs` depends on the completed initial seed. PAW datasets and large wavefunctions/caches are retained there, not redistributed here. This snapshot is therefore not a standalone VASP distribution. `scripts/` preserves the actual analysis implementations; no shared response core was changed.
