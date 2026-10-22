> **2026-10-22 direct-data correction: THIS EARLIER LITERATURE-ONLY SCREEN IS SUPERSEDED.** The user supplied the full BiDB ASE/SQLite database; complete record-level filtering and a real three-operation Seitz check changed the first-choice theory demonstrator to **Pd2Te2 id 980**, with **Pd2Se2 id 975** as a close control. Both are nonmagnetic, metallic, eight-atom, Stable p-6m2 (D3h) bilayers with a directly verified *horizontal layer-swapping mirror*, matching BD27 ideal assumptions. VSeTe is the preferred *magnetic* extension; VTeS has only two Stable entries and neither has pure mz. The precise results, UID, caveats and next experiment are in [the direct database report](BD29_BIDB_DB_DIRECT_SCREEN_20261022.md). No extra files or CIF export are requested from the user; the earlier lack-of-record-access paragraph below is **historical only**.

# BD29 | C2DB / BiDB scouting for a source-backed layer-response material

**Date:** 2026-10-22  
**Status:** MATERIAL SHORTLIST ONLY — no structure-level symmetry validation, response calculation, or novelty claim yet.  
**Purpose:** Find a computationally economical real bilayer for the response-field idea while preserving the user's distinction between a concrete example and a general theory.

## 1. Databases and what was actually verified

### C2DB
The official [Computational 2D Materials Database](https://c2db.fysik.dtu.dk/) provides monolayer structure, stability, electronic, magnetic, and optical properties. Its public 2DHub documentation says the full dataset is provided on request; the browser interface remains useful for record-level searches. The 2024 high-throughput magnetic-order study by Sødequist and Olsen screened stable magnetic C2DB monolayers with one magnetic atom per unit cell using spin-spiral calculations and SOC, reporting 58 ferromagnets, 21 collinear antiferromagnets, and 85 non-collinear ground states, including 15 DMI-driven chiral spirals. This is an important warning not to accept a ferromagnetic starting state as the actual ground state. Sources: [C2DB/2DHub](https://www.2dhub.org/c2db/c2db.html), [Sødequist & Olsen, npj Comput. Mater. 10, 170 (2024)](https://doi.org/10.1038/s41524-024-01318-2).

### BiDB
The [van der Waals Bilayer Database](https://bidb.fysik.dtu.dk/) is the more direct first filter: its paper constructed 8,451 homobilayers from 1,052 monolayers and calculated electronic, magnetic and vibrational properties for 2,586 stable candidates. It reports stacking-dependent magnetic states in 37 monolayer families and specifically discusses AgVP2Se6, FeTaSe3/FeTaTe3 (the paper text is inconsistent about this formula; do not use it until resolved), and Janus VTeS. It also reports that a VTeS stacking can become ferrimagnetic because the Janus dipole drives interlayer charge transfer. Source: [Pakdel et al., Nat. Commun. 15, 932 (2024)](https://doi.org/10.1038/s41467-024-45003-w).

**Access limitation:** I verified the database's published documentation and material claims, but did not extract the exact current BiDB UID/CIF or the full bidb.db rows for VTeS and AgVP2Se6. Thus this is a source-backed shortlist, not a completed record-by-record database screen. The official database's noncommercial license/terms must be followed.

## 2. Ranked candidate shortlist

| Rank | Candidate | What is source-supported | Why it may help | Main risk / next gate |
|---|---|---|---|---|
| **A — first to inspect** | **Janus VTeS homobilayer** | BiDB reports metallic bilayers, stacking-dependent magnetic moments and a ferrimagnetic stacking driven by charge transfer from the built-in Janus dipole. | Compact chemical formula; Te gives a plausible SOC channel; stacking already changes layer magnetism, making layer-resolved response a natural question. | Exact BiDB structure, orientation, primitive-cell size, magnetic order, stability and symmetry are not yet retrieved. Do **not** assume horizontal layer exchange or reuse the hexagonal m_z rule: a Janus layer breaks that operation unless the full bilayer mapping genuinely restores it. Confirm the spin-spiral ground state rather than trusting a collinear FM initialization. |
| **B+ — polar magnetic benchmark** | **AgVP₂Se₆** | BiDB identifies stacking-related FM/AFM states; experimental reports establish a van der Waals ferromagnetic semiconductor with \(T_C\) around 19 K. C2DB spin-spiral work predicts that DMI can destabilize the collinear FM state toward a chiral spiral; it is polar and has in-plane ferroelectric order. | Strong magnetoelectric narrative and layer/registry dependence; relevant to the user's interest in magnetic bilayers and intrinsic ferroic minima. | Low symmetry weakens the elegant high-symmetry zero-set story; magnetic ground-state competition (collinear FM vs spiral) must be resolved; likely more complex than the first minimal model. Sources: [C2DB magnetic-order study](https://doi.org/10.1038/s41524-024-01318-2), [experimental ferromagnetic semiconductor](https://doi.org/10.1002/adfm.201910036). |
| **B — benchmark, not novelty target** | **Bilayer CrI₃** | Bilayer AFM CrI₃ has reported layer-opposite Rashba in-plane spin textures; FM CrI₃ has reported interlayer-coupling-related spin canting; current-induced SOT has been studied in TMD/CrI₃ heterostructures. | Clear SOC and magnetic benchmark; useful to validate layer-resolved operator projections and magnetic symmetry handling. | Crowded prior art in spin-layer locking and SOT; a material-specific novelty claim would be weak. Treat existing results as direct/nearest prior art for those exact observables, not as proof that every general response-field theory is pre-empted. Sources: [AFM CrI₃ layer-spin locking](https://www.sciencedirect.com/science/article/pii/S0022369722005170), [2026 TMD/CrI₃ SOT](https://doi.org/10.1103/llyn-hs8b). |
| **B− — symmetry/insulator control** | **MnPS₃** | Experimentally established layered antiferromagnetic insulator, readily exfoliable; it is a reasonable magnetic-symmetry test system. | Useful for testing frequency-dependent spin-density/torque response without metallic Fermi-surface complications. | Mn/P/S has relatively weak SOC compared with Te/I systems, so a sizable current-induced out-of-plane spin signal is not assured; magnetic unit cell and AFM order increase setup cost. Better as a symmetry/response control than first material. Source: [MnPS₃ material summary](https://quantumlab.uark.edu/mnps3/). |
| **Control only** | **MoSSe, MoTe₂, WTe₂ bilayers** | The 2026 LEE paper explicitly uses these as first-principles examples of layer Edelstein response. | Best for reproducing a published baseline and checking operator conventions, layer projection and current protocol. | The generic LEE/current-induced layer-spin claim is already directly covered; do not use them as the sole evidence for a new phenomenon. Source: [LEE, PRL 137, 076401 (2026)](https://doi.org/10.1103/qhvs-pfkq). |
| **Alternative to check** | **2H-VSeTe Janus monolayer / bilayer** | Literature reports a magnetic 2H Janus phase with substantial Te SOC; recent work studies stacking/strain control of polarization and DMI. | Potential way to retain trigonal geometry while introducing magnetism and polar asymmetry. | Exact stable bilayer registry and orientation must be obtained from BiDB; Janus orientation changes the layer-exchange operation, and DMI/polarization already has nearby literature. Source: [Janus VSeTe stacking/DMI study](https://www.sciencedirect.com/science/article/pii/S0925838826028677). |

## 3. Prior-art roles: do not collapse these categories

### Direct / potentially direct competition
- **LEE, PRL 137, 076401 (2026):** direct for the general nonmagnetic bilayer layer-Edelstein response, static stacking symmetry criterion, and the published MoSSe/MoTe2/WTe2 demonstrations. The exact overlap of its Supplemental Material with a continuous-registry, operator-resolved zero/index classification still needs an equation-level audit.
- **CrI₃ spin-layer locking and TMD/CrI₃ SOT papers:** direct for their own momentum-space spin textures and heterostructure torque calculations, respectively. They make CrI₃ a strong control, not a clean novelty target.
- **BiDB's VTeS / AgVP₂Se₆ results:** direct for the already reported stacking-dependent magnetic state in those material families, **not** for an as-yet-uncomputed current-induced layer-resolved spin-response field.

### Foundational / nearest prior art
- **BiDB (Pakdel et al., 2024):** foundational for screening stacking stability, magnetic-state switching and registry-dependent material properties.
- **C2DB spin-spiral screening (Sødequist & Olsen, 2024):** foundational for deciding whether the monolayer magnetic ground state assumed in a bilayer calculation is plausible.
- **2026 trigonal spin-texture / spin-charge-conversion screening:** close method and response precedent for current-induced spin polarization in monolayer trigonal materials; it does not, from the accessible main text, establish the same bilayer registry-space layer-resolved theorem.

### Introduction context
- General sliding ferroelectricity, magnetic van der Waals materials, and optical shift-current vortices can motivate configuration-space response textures. Their existence does **not** automatically pre-empt a higher-level theory of a different observable. Compare the exact proposition, source/output representation, and assumptions.

## 4. What is the high-value physics question now?

Do not start by running a broad material sweep. First ask:

> **Can a source-validated bilayer's layer/site-resolved current-induced spin response be treated as an equivariant field over registry space, and does the magnetic/Janus layer representation generate an operator-dependent selection rule or response texture that is not already implied by the fixed-registry LEE criteria?**

This is a more precise S2 theory question than “find a new material with spin polarization.” It also allows a material example to serve as a test of a general result rather than being mistaken for the entire novelty claim.

## 5. Minimal, economical workflow

1. **Record-level data gate (no DFT):** obtain exact BiDB UIDs/CIFs for VTeS and AgVP₂Se₆; collect the stable stackings, binding energies, magnetic flags/order, layer group, inversion flag, band gap and primitive cell size. For VTeS, explicitly identify which layer orientation each UID represents.
2. **Magnetic-order gate (no bilayer DFT yet):** cross-check the monolayer ground state against C2DB's spin-spiral data / original literature. If the state is spiral or uncertain, do not silently use a collinear FM state.
3. **Symmetry gate (Python only):** derive the actual unitary and antiunitary magnetic layer group for one selected stacking, including the Seitz translations and layer/site permutation. Decide whether \(S_\pm^z\) are meaningful channels for that magnetic state; do not assume the prior ideal \(m_z\) model applies to Janus VTeS.
4. **Minimal Hamiltonian gate:** only after the actual symmetries and SOC orbitals are source-backed, build a minimal spinful TB/Wannier model and calculate the fixed-current layer-resolved \(S_z\) response over a small set of symmetry-related registries.
5. **Material/novelty gate:** compare the resulting exact observable and symmetry orbit with LEE SM and magnetic spin-texture/SOT papers. Only then decide whether to run a small DFT set or larger registry scan.

**No DFT, WSL, Slurm, or Japan-supercomputer jobs are authorized by this screening note.**

## 6. Where user support would be useful

The main remaining bottleneck is **record-level structure access**, not brainstorming. The web-accessible database documentation and published paper identify the candidate families, but the exact BiDB record UID/CIF and full property row have not yet been extracted in this environment.

If convenient, upload either:
- the BiDB bidb.db export, or
- just the BiDB row/JSON/CIF files for **VTeS** and **AgVP₂Se₆**.

The two candidate records are sufficient; a full C2DB dump is unnecessary. With those files, the next pass can produce a grounded symmetry and computational-cost comparison. Without them, the literature-based shortlist can continue, but I should not pretend the exact database-level screening is complete.
