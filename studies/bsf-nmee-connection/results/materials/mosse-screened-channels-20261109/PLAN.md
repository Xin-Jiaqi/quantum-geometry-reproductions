# Screened MoSSe channel test — sealed before results

2026-11-09. User authorizes the proposed causal test and raises a concern about
restricting discovery to MoSSe/MoS2. No new material or DFT is launched here.

Question: does the selected cross-layer occupied–empty mixing that changes
zero-gate layer NMEE also support the screened upper-layer gate sensitivity?

Reuse the actual zero and ±0.10 V/nm gate-SCF wavefunctions, potentials, spin
vertices and Bloch links. For each applied field, freeze its original converged
electrostatics while modifying selected Hamiltonian channels. This is a
conditional electronic intervention, NOT a self-consistent material with
hybridization switched off. Screening is neither refitted nor recomputed.

Define 44 d+p orbitals and layer valence projectors ONCE at zero gate. Transport
these physical orbitals into each field-dependent KS space using PAW overlaps
<u_F|u_0> and symmetric orthonormalization; report representation loss. Use the
same zero-gate layer/valence labels for all fields and both magnetic probe signs.
At each field scale that field's selected VC Hamiltonian block consistently;
retain all remaining KS states and physical spin/electric response vertices.
Thus the field-dependent part of the selected channel is scaled as well.

Compare alpha=0,0.5,0.75,1 for VC; remove VV+CC and either direction of VC as
controls. Repeat 96/112 KS states and three zero-gate projection windows (4/8/12
frozen empty states). Compute C at zero gate and the same ±0.10 secant G used
in independent SCF. Reject gap closure or poor subspace representation; report
finite-window and magnetic-step sensitivity. Alpha=1 must exactly recover the
existing SCF results. Do not identify no-VC with isolated layers or electrostatics.

Support: removing the same selected channel weakens both local sign rewriting
and upper-selective gate gain, robustly across projection choices. Opposite
result: if selectivity persists or grows, revise the common-channel story;
neither force coincident alpha thresholds nor infer that C reversal causes G.
Directional controls may identify a sharper hypothesis, not a unique additive
decomposition. A prospective material prediction still requires a separate test.

Broaden after/alongside this diagnostic with a small physically varied set:
same-family stronger SOC, heterobilayer with changed layer alignment, and later
a different symmetry. Material inclusion must not depend on obtaining a large
NMEE. No candidate NMEE is asserted from optical hybridization evidence.

Resource: existing data and tested Python/PAW stack on Japan IMR; CPU DP_002,
4 cores, 8 GiB, 20 min limit; expected additional output below 200 MB. All wave
reading and numerical processing inside PBS, never login nodes. Preflight:
queue enabled, no DP_002 jobs; numeric user quota caps zero/unspecified and
group/project cap unknown. User's free-resource statement retained. No other
account jobs touched.
