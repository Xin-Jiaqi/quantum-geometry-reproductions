# Non-Janus control: 3R MoS2

2026-11-07 Asia/Shanghai. Authorized scope: compare isolated layer, stable 3R bilayer, mirror partner and a few perpendicular biases. Suspend twelve-state work. Do not select a material or stacking by the computed response sign.

## Selection and question

MoS2 is selected over WSe2 because it shares the Mo–S backbone with the established MoSSe case, permits the same PAW protocol, and BiDB 1930 is sliding-stable and the most strongly bound of the database's sampled bilayers (26.5136 meV/A²; nearest 2H 26.5050). This small energy difference is not a claim of an exact global ground-state order. WSe2 remains a future SOC-strength comparison, not part of this batch. Local database copied read-only from the user's WSL E:/databases/bidb.db. All coordinates originate in BiDB. No NEB or phonon calculation is implied by database slide stability.

The physical question is whether an interface creates the fixed-ion, insulating s_x = Cxx Ex Bz response absent in a free, unbiased D3h monolayer, and whether a stable stacking or gate can control it. A nonzero response is symmetry-allowed, not symmetry-guaranteed. Its amplitude and sign are unknown in advance.

## Structures and independent predictions

- mono: BiDB 1929, point group D3h, horizontal mirror. With zero gate, total Cxx must vanish. Real-space halves of a monolayer are NOT separate layers and will not be reported as layer responses.
- A: BiDB 1930, parallel 3R registry, C3v, stacking operator identity plus translation; upper Mo relative to lower Mo at (1/3,2/3) in the exported cell.
- B: horizontal m001 image of A, followed by a common in-plane translation to put bottom Mo at the same origin; upper Mo registry (2/3,1/3). Both constituent monolayers have horizontal mirror, so B retains the same chemical orientation of the bottom layer. It inherits A's stability by symmetry in a free slab.
- At conventional perpendicular field F, C_bottom(B,F)=-C_top(A,-F); C_net(B,F)=-C_net(A,-F). The fields Ex and Bz and laboratory axes are held fixed. This relation permits, but does not require, a measurable nonzero signal.
- For the isolated monolayer, C_mono(F)=-C_mono(-F). Gate-induced monolayer response is an essential control: interface activation alone does not prove a unique hybridization mechanism, since ordinary electrostatic symmetry breaking could also activate it.

## Bounded calculation sequence

First independent mono F=0 and A F=0, k6. After checking the response and numerical floor: A F=+/-0.25 V/nm, B F=0, mono F=+/-0.25 V/nm. Add one response-resolution/mesh check where needed, without a broad field scan. All slabs are neutral, ions fixed, common height 27.25422529 A. Density partition uses the interface midplane and vacuum midpoint, with symmetric boundary weights.

Use self-consistent SOC finite electric enthalpy, all four (Ex,Bz) signs, Ex=0.004 V/A and BEXT=0.016 eV initially. These are numerical differentiation probes, not proposed laboratory fields. 450 eV, Gamma-centered 6x6x1 full mesh, EDIFF=1e-9, Mo_pv/S, VASP6.4.3. Gate is the conventional physical F, implemented as EFIELD=-F/10; this reverses the native VASP force-field sign used by some older project scans. Record both signs explicitly. Convert Pauli spin per area to C in 1e-9(nm V T)^-1 using each structure's actual area and g=2.0023193. This is a spin-only Zeeman response, not orbital magnetism, transport spin accumulation or ion-relaxed response.

## Interpretation decided before results

1. A reliable interface signal absent in unbiased mono supports interface activation. Gate mono control distinguishes this basic symmetry effect from a stronger interlayer-hybridization claim.
2. Opposite measurable total signals in A/B offer stacking control. Mirror relation alone is not a novelty claim or proof of experimental detectability.
3. Large opposite layer signals with weak net signal point toward hidden layer response and electrical imbalance; report the cancellation rather than select only the large local values.
4. A signal below the numerical floor or impractically small is a negative material outcome. Do not turn numerical sign changes into switching.

No twelve-state or universality claim is required for this test. The role of KM will be reconsidered after the material comparison.

Execution refinement: unbiased mono/A/B retain all four probe signs as the explicit time-reversal check. Gate points use +/-Ex at positive B, exploiting time reversal of the nonmagnetic insulating slab; this is the same reduced-sign strategy already used for MoSSe. Their B-even residual is recorded as unmeasured (null), never zero. The original four-sign source is archived separately.

Runtime refinement decided before inspecting any response: the bilayer four-sign runs are too slow for the short DP allocations. A/B will stop only AFTER both positive-B electric-field runs have converged and written their completion records. A separate PBS collector then reads those full densities using the same time-reversal-reduced protocol as the gate calculations and prior MoSSe work. Their negative-B residual remains unmeasured. The zero-gate monolayer retains all four signs. This avoids spending the exploration budget on redundant time-reversed SCFs; it does not replace unconverged densities with results. Both original and collector job IDs are retained.
