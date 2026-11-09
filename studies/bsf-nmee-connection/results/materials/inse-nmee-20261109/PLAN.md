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
