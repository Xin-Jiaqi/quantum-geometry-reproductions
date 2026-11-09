# Independent gate-SCF test, 2026-11-09

Question: does gate screening preserve the approximately 30% net NMEE gain and
upper-layer dominance predicted by the frozen-potential MoSSe S/S model?

Use the SAME relaxed e000 geometry (BiDB 7269). Fixed ions, insulating charge
state, no gate doping. External physical field F = ±0.10 V/nm. VASP EFIELD=-F/10
eV/Å: https://vasp.at/wiki/EFIELD defines its sign opposite to conventional field.
Use LDIPOL, IDIPOL=3 and the existing vacuum cell. Check actual potential slopes.

For each sign compute full bilayer, exact isolated lower, exact isolated upper.
Do not assume the upper reference is the exact mirror of the lower. Existing
zero-field bilayer/lower data are retained; the central derivative does not
require a new upper zero-field calculation. All calculations use the prior SOC,
450 eV cutoff, 6×6 full mesh, EDIFF1e-9 and 112/64 bands. This is a fast mechanism
test, not a claim of final mesh convergence. Reuse zero-field WAVECAR/CHGCAR for
bilayer/lower initialization with ICHARG=1, then fully converge the new field;
upper begins from atomic densities. Do not alter source calculations.

Measure (C(+F)-C(-F))/(2F) for lower, upper, and net EB response. Compare full
net with the sum of BOTH independently calculated isolated responses. Compare
96/112 and 48/56/64 KS response spaces, halve the magnetic probe step, check
layer sum, insulating gaps, normalization and wave/link gauge consistency.
The gate density is self-consistent; the probe Ex,Bz NMEE definition remains the
existing fixed-potential spin-only interband susceptibility (no new B-probe SCF).

Expected if hypothesis survives: |G_net| exceeds isolated sum across retained
settings, with upper |G| larger than lower. If screening suppresses the gain,
report that negative result and do not increase the field to seek enhancement.
No same-layer reversal target. A second smaller field pair is warranted only
if central-slope nonlinearity or numerical uncertainty obscures this decision.

Two scoped PBS jobs, each 4 CPU,24GB,2h maximum, one field sign per job, three
SCFs sequentially. Scientific input preparation, VASP, export, wavefunction
reading and response all on compute nodes. Frontend only short queue/quota
queries, small-file transfers and submission. No new material/strain/relaxation.
Estimated new storage a few GB, no migration to another cluster. Keep licensed
POTCAR/WAVECAR out of public archive. Update unified worklog with result or
honest pending status and resume metadata if queue/runtime extends this turn.
