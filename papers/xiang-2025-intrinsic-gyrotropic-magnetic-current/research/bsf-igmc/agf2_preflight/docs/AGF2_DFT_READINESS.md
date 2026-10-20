# Stage 4: AgF₂ DFT / Wannier scientific-input gate

**Scope**: make the first candidate calculation reproducible *before* running a heavy first-principles job. Literature validation is separate from independent reproduction of its numerical inputs.

## 1. Sources and already-established (literature) facts

Zhu et al., *Sliding Ferroelectric Control of Unconventional Magnetism in Stacked Bilayers*, *Phys. Rev. Lett.* **135**, 056801 (2025), DOI `10.1103/dmzg-ck2t`.

| Item | Source-reported value | Consequence for this project |
| --- | --- | --- |
| Monolayer AgF₂ LG | `p2_1/b11`, LG17 | Verify atomistic monolayer is in the intended polymorph, not a bulk AgF₂ phase |
| Bilayer LG / PG | `pb11` / `m_100` | Non-centrosymmetric polar `P_z` possible |
| Polar pair | `O1={m_001 |(0.5,-0.16)}`, `O2={m_001 |(0.5,+0.16)}` | These are **stacking operations**, not complete POSCARs |
| Source Pz | `±1.66 pC/m` | Independent DFT polarization reference, not a fitted input |
| Source barrier | `17.1 meV/f.u.` | Literature BSF gate passed **for their setup** |
| Magnetic Type-I vs II, 0% strain | Type-I lower by `0.26 meV/f.u.` | Energetic ranking is extraordinarily sensitive to DFT protocol/convergence |
| SOC anisotropy Type-I | easy axis `[100]`; MPG `m'_100` | Single-state `C_xy, D_xx` **allowed**, not guaranteed finite |
| Type-II | studied at `4%` biaxial tensile strain, connector `m'_001` | Separate structural/magnetic campaign, not interchangeable with zero-strain Type-I |

Those are sourced to the author's published main paper and the user-provided detailed notes; we have not executed an independent AgF₂ calculation. Prior eight-band model's lack of a double well is a **different Hamiltonian** and does not contradict this literature BSF result.

A public bulk AgF₂ crystallographic record (e.g. JARVIS/JVASP) is **not the same as the LG17 monolayer**, so it must not be reused as an unvalidated bilayer POSCAR. The present investigation has *not obtained complete author-relaxed 2D atomic coordinates* or an authenticated pseudopotential/U/vdW recipe.

## 2. Atomic-data blocker and input contract

Before VASP, obtain an actual, provenance-preserved LG17 **monolayer** structure or the exact relaxed B₁/B₂ pair and document its source/license. Explicitly check species, magnetic sublattices, lattice vectors, vacuum axis, fractions and numbers of formula units. For each B state record `POSCAR` SHA256; reject fake coordinates and reject stacking translations as substitutes for atomic positions.

Inspect the author's actual Supplemental Material / available coordinate repository to freeze `PBE/PBEsol` or other functional, Ag 4d Hubbard U / J treatment, vdW correction, PAW POTCAR variants, k density, plane-wave cutoff, vacuum/dipole correction and moment order. Do **not** infer the needed U or cutoff from a different AgF₂ paper. The JSON template explicitly contains `null` for unresolved physical settings and coordinates; the validator returns `BLOCKED_INPUTS` until each required value is fixed.

The preflight script verifies presence/hashes/nominal choices only. It cannot prove the provided structures are equivalent to the paper's atomic configurations, so a separate independent 2D layer-group and structure-mapping check is still mandatory.

## 3. Minimal execution order once real input exists

**A. Structural baseline (zero strain, Type-I)**: same lattice and identical functional/PAW/U/vdW protocol for B₁ and B₂; compare total energies and moment-preserving constrained relaxation; establish B₁/B₂ local minima and ion-relaxed pathway / barrier. Include a check of layer-group compatibility over different `symprec`. Importantly, compare energies normalized to the **same formula unit**.

**B. Magnetic-order gate**: for *both* B states test Type-I/II before SOC, then SOC-magnetic `[100], [010], [001]`. Confirm state-to-state domain connector and single-state antiunitary group **from actual magnetic moments**, not the nonmagnetic structure's point group. Because `0.26 meV/f.u.` is sub-meV, force/SCF precision must be established by **explicit convergence data at a substantially smaller energy scale** and functional/U sensitivity, not a single calculation.

**C. Finite-strain Type-II campaign**: explicitly apply 4% biaxial strain, re-relax both stacks, independently verify two local minima and magnetic order, MAE, energy ordering and corresponding `m'_z` magnetic connector. Do not compare numerical Type-I and Type-II results as if performed at the same strain.

**D. Band and operator gate**: with self-consistent SOC, test whether each magnetic phase is an insulator and whether the relevant excitation frequency lies below the direct optical gap. Use Wannier stage-1 quality gates without selecting Wannier ranks by the desired response. Extend the interface to export **projected spin and full position/Berry operator** in exactly the Hamiltonian's Wannier gauge (and their hashes), validate spinor matrix elements and gauge covariance. Existing `hr.dat + centres` is not a sufficient IGMC input.

**E. IGMC response**: compute phase-resolved clean AC susceptibilities with explicit Zeeman convention (`g`, `μ_B`, units, B amplitude, rad/s vs THz), Fermi occupation, frequency, degeneracy treatment and k-mesh/broadening convergence; compare B₁/B₂ in the **same lab axes**, and include spin-Zeeman vs orbital contributions plus induced electric-field background as distinct physics. Preserve the independent Xiang Fig. 2(d) quantitative issue; never fit a multiplicative factor to paper artwork.

## 4. Current evidence/state

- `PASS_SYMMETRY`: two connection parities and Type-I internal tensor allowedness derive algebraically.
- `PASS_LITERATURE_BSF`: the *published* AgF₂ structural calculations report distinct stable stacking minima and a finite barrier.
- `BLOCKED_ATOMIC_INPUT`: no source-verified B₁/B₂ POSCAR pair has been supplied.
- `BLOCKED_DFT_PROTOCOL`: functional/U/vdW/POTCAR choices are not frozen by authenticated source input.
- `BLOCKED_REAL_MATERIAL_IGMC`: no same-gauge multiband spin/position operators and no numerical current have been calculated.

**Go/no-go**: source-backed AgF₂ is a *sensible first DFT target*, but the present files are a **preflight dossier, not DFT results**. A real calculation should begin only after the atomic coordinates and the computational recipe are resolved.
