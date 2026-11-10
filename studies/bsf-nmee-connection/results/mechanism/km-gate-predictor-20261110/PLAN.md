# Gate direction from the existing eight-state KM model

Before new calculations, 2026-11-10.
Question: can zero-bias energies, layer character, and physical spin/electric
vertices predict which layer is strengthened by a small positive layer offset?
No new material, DFT, strain, or twelve-state fit.

Use the existing mirror bilayer, path fraction 0.25, mass45 meV, t10/15/20 meV,
spin and electric probe along 60 degrees, Bz. Gate coordinate U raises upper
on-site energies by U/2 and lowers bottom by U/2. It is an internal potential
energy difference, NOT an unscreened physical field in V/nm.

Differentiate the existing insulating spin-electric kernel analytically in U,
then take its central B derivative. Split the gate derivative into changes of
transition denominators and changes of spin/velocity matrix elements, using
the full layer vertex W=(Ptop-Pbottom)/2. Eigenstate rotation is controlled by
W_mn/(E_n-E_m), with the sign fixed by full products of physical vertices.
This is a derivative diagnostic, not a pair of separately physical materials.
Reject unresolved same-sector degeneracies; retain the full response and audit
phase invariance, derivative closure, magnetic step and integration changes.

First compute predictions at U=0 and save them. Then independently evaluate
U=+/-0.25 and +/-0.50 meV, preserving an insulating common gap. Validate signs
and the local linear prediction. Do not fit to these held-out states. Repeat
representative states at fraction0.24/0.26 only if needed to test a claimed
mechanism, not to chase a universal sign law. Compare qualitatively with the
already-completed MoSSe/MoS2/InSe controls; do not call an instantaneous KM
energy-sector decomposition the same object as a material's fixed separated-
layer channel intervention. Report failures and limits explicitly.
