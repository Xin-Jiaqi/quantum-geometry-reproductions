# InSe homobilayer NMEE: pre-result plan

2026-11-09. Authorized after selecting BiDB 864, with matched isolated layer
(parent 862), and its m001-related pure-sliding partner. No response-dependent
material selection, relaxation, strain, twelve-state fitting or wider survey.

Observable: fixed-ion insulating Pauli spin density delta s_x,l=C_l Ex Bz.
Perpendicular physical gate F is independent of the two probe fields;
VASP EFIELD=-F/10 eV/Angstrom. Gate is self-consistent; Ex/Bz response is the
existing full-KS frozen-potential spin-only derivative, NOT full three-field
self-consistency or orbital magnetization.

Hypotheses before outcomes:
1. Matched isolated D3h layer has zero C at F=0 and C(F)=-C(-F).
2. A/B partners obey C_lower(B,F)=-C_upper(A,-F), C_net(B,F)=-C_net(A,-F)
   for a free slab. These are symmetry checks, not novelty claims. A allowed
   nonzero C is not guaranteed to have large magnitude or a specified sign.
3. Interface-induced electrical sensitivity may involve Gamma-region In s/p,
   Se pz cross-layer occupied-empty mixing, with px/py needed for spin and
   in-plane electric vertices. This is a falsifiable transfer hypothesis, NOT
   inferred from optical Rashba evidence or from a small gap alone.
4. No prediction of absolute layer sign, favored layer or numerical enhancement
   is made. A failure of enhancement is scientifically useful and retained.

Structures: A identity stacking plus (1/3,2/3), C3v, BiDB sliding stable and
PBE gap .892 eV. It is 1.695 meV/cell above sampled inversion structure 863,
not asserted a unique ground state. Database dipz=0.0 is not a verified finite
polarization. B is the horizontal mirror plus common in-plane translation;
same lower-layer chemical orientation. Cell height 31.734612 A, all centered;
mono uses exact constituent geometry and cell, with In/Se atomic ordering.

Economical staged sequence: independent SOC SCF and full-KS response at mono0,
A0, B0; then A +/-0.10 V/nm and matched mono +/-0.10. Stop a batch if a gap
closes or SCF fails; no signal inferred from an unconverged run. Mirror/single
layer symmetry and response-band cutoffs define the numerical floor before
interpreting sign/gain. Distinct finite-field baselines are retained, not
interchanged with older finite-B density protocols.

Use PBE, In_d/Se PAW PBE.64, SOC, ENCUT450, full Gamma-centered6x6 mesh,
EDIFF1e-9, neutral fixed ions. KPAR4 on4 MPI ranks, NCORE1; bilayer160 spinor
bands and mono80. Report full-KS response at128/144/160 and64/72/80 respectively,
B-probe energy +/-.002 and +/-.001 eV. All occupied bands are retained.
Selected In s,p and Se p orbitals exclude deep In4d/Se4s only for future
Hamiltonian channel interventions; the layer response itself is not downfolded.

Follow up only after baseline: matched electrostatic reference, full versus
VC/VV/CC interventions, and Gamma contribution; use screened gate reference
when claiming practical electrical advantage. No assumption that Rashba spin
splitting, Mexican-hat dispersion or band-edge optics already proves static NMEE.

Resource: Japan IMR tested SOC/W90/PAW stack; first-tier SOC+link workflow not
validated and data/software reuse favors IMR. Two DP_002 allocations at most,
4 CPU/24GB/2h each; mono batch3SCFs, bilayer batch4SCFs. Each step writes its
completion separately. Expected new output <15GB. User work quota caps zero/
unspecified; group/project cap unknown. All VASP, wave processing, profiles,
hashes and numerical work on PBS compute nodes, never the frontend. No other
account jobs changed. Local only prepares inputs and summarizes small results.

## Execution refinement after the first A0 result (2026-11-09)

The zero-field layer values are unexpectedly large and nearly compensate. One
independent 8x8 A0 SCF, initialized from the converged 6x6 charge density, is
added to test integration sensitivity before interpreting a large magnitude.
The primary seven-state gate comparison is unchanged. After the mono allocation
finished, Ap was delegated to a separate allocation (560116); the original
bilayer runner skips its completed record. No duplicate gate SCF is intended.

## Conditional prediction after C(A,0), before inspecting bilayer gate responses

2026-11-09T18:20:27.691364+08:00

The large opposite layer C(0), together with the much smaller isolated-layer
gate sensitivity, motivates a separate possibility: if paired cross-layer
virtual transitions dominate and their gaps are detuned oppositely by F, the
CHANGES of the two layer responses can have the same sign and comparable
magnitudes. Net electrical enhancement could then occur without strong layer
selection. This is conditional on the electronic mechanism, not required by
C3v and not yet a result. It does not predict an absolute sign or gain factor.
Opposite layer slopes, little net gain, or strong layer asymmetry are retained
as alternative outcomes. Gamma concentration and directional VC controls will
test the proposed interpretation. Large zero-field layer C alone proves none
of these gate properties.

## Outcome of the gate prediction (all seven primary points complete)

The symmetry predictions survive the initial calculation. The conditional
same-sign layer-increment branch does NOT occur: positive F makes the negative
lower C less negative and the positive upper C less positive. Both layer
magnitudes decrease, with the upper changing more. Their increments oppose
one another; net gain is about 1.4-1.6 relative to the matched isolated pair
over the tested band cutoffs, rather than a giant gain inferred from C(0).
Neither layer nor the net response reverses sign within +/-0.10 V/nm.

Selected cross-layer occupied-empty mixing matters, but is not the unique
origin of the negative lower C: it stays negative when that selected mixing
is removed. Occupied-occupied/empty-empty mixing also matters. Removing one
directional occupied-empty block reverses the direction in which the gate
tunes each layer, while leaving the static layer-spin signs unchanged. These
are frozen-potential Hamiltonian interventions, not realizable independent
experimental knobs or additive contributions. The real gate-SCF gain and
these interventions must not be conflated.

Execution record: the original 560011 batch completed A0/B0, then was explicitly
stopped during the newly started Am to stay within its allocation. Its Am
data are not used. Fresh Am/Ap completed in 560139/560116. The reference and
channel controls completed in 560207. The independent 8x8 A0 check and a small
Gamma orbital-character diagnostic have separate allocations in jobs.json.
