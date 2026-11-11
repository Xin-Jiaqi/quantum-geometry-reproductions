# Stage 8: literature-based candidate prescreen for sliding-bilayer static NMEE

**Status: LITERATURE PRESCREEN ONLY. No DFT, Wannierization, or material response calculation has been run.** This stage selects realistic structures for a later symmetry/structure preflight; it does not claim a real-material NMEE prediction.

## 1. Freeze the target before screening

The target carried over from Stages 7i–7j is the **gapped, zero-temperature equilibrium/static** Zeeman derivative of the spin-electric response,
\[
\Gamma_{xz;x}=\left.\partial_{B_z}\chi_{S_x,E_x}\right|_{B=0},
\]
not the nonequilibrium driven-DC intrinsic coefficient discussed in Jia et al. (2026). The paper reproduction stays frozen at `papers/jia-2026-nonlin-magnetoelectric-edelstein/`.

We want to test the Stage-6 mechanism in a material family with a high-symmetry nonpolar reference stacking and a sliding path to polar partners:
- Reference: AA-like parallel 1H-TMD homobilayer with layer group `p-6m2`, point group `D3h`, including (C_3) and horizontal mirror (m_z).
- Polar R-stacked AB/BA partners: layer group `p3m1`, point group (C_{3v}); inversion and (m_z) are broken.
- In a fixed (C_{3v}) polar domain, the spatial tensor component `Gamma_xz;x` is allowed; the full tensor must still be computed using the actual relaxed structure and SOC representation.
- The reference (C_3+m_z) symmetry predicts a cubic-or-higher onset of the in-plane trace with displacement, while the traceless symmetric/nematic part can be linear. This is a testable hypothesis for a material, not a claim that the response magnitude or sign is universal.

The stacking symmetry assignments are supported by the 2-layer TMD stacking classification in the supplementary material for [Di Berardino et al. (2025)](https://journals.aps.org/prb/supplemental/10.1103/PhysRevB.111.125403/PRB_DiBerardino2025_SI.pdf) and the explicit (C_{3v}), AB/BA sliding-ferroelectric discussion for bilayer MoS2 in [Scientific Reports (2019)](https://doi.org/10.1038/s41598-019-50293-y).

## 2. Ranked shortlist

| Rank | Candidate | Evidence already in literature | Why it fits our physics | Main caveat / status |
|---|---|---|---|---|
| **A1 — primary** | **3R-WSe2 homobilayer** | Experimentally observed AB/BA ferroelectric domains and field-driven band/domain switching; R-stacked bilayer has an indirect-gap electronic structure. A 2026 ARPES/DFT study reports a (520\pm20\) meV spin-orbit-driven splitting at K in quasi-free-standing rhombohedral WSe2. | Strong W SOC, a semiconductor, real AB/BA polar partners, and a same-material sliding coordinate. Best initial candidate for a strong spin response. | Still must verify the gap over the **entire slide path** with SOC and under all source corners; confirm the exact computational structure is substrate-free and has the intended point group. No NMEE value known. |
| **A2 — benchmark/control** | **3R-MoS2 homobilayer** | AB/BA polar states and switching have extensive theoretical/experimental support. One NEB calculation reports degenerate AB/AC states and a 15 meV switching barrier for its chosen path. | Same clean 1H-TMD symmetry family as WSe2; a well-studied structural/control case for testing (C_3+m_z\) selection rules. | SOC is weaker than in WSe2, so the response may be smaller. A control is still valuable even if its response is weak. |
| **B1 — secondary** | **3R-WS2 homobilayer** | Published DFT/SOC calculations show stacking-dependent band splitting in 3R WS2; it shares the parallel 1H-TMD stacking family. | W supplies strong SOC; provides a chemically related check against WSe2-specific band-edge effects. | This pass has not verified an equally direct experimental sliding-switching/barrier dataset or full-path gap. Keep behind WSe2 until those items are checked. |
| **B2 — exploratory** | **Cd2Se2 bilayer** | A June 2026 preprint screens 280 M2X2 compounds, reports eight stable sliding-ferroelectric candidates, and identifies SOC-induced valence-band spin splitting in bilayer Cd2Se2 that changes with polarization. | Directly connects sliding polarization and spin splitting; potentially high novelty. | Preprint-level evidence; the exact relaxed structure, reference-stack (C_3+m_z) symmetry, full-path insulating gap, material availability, and usable Wannier operators are not established by this prescreen. Do not make it the first DFT target. |
| **C1 — alternate symmetry class** | **1T-PtS2 bilayer** | A 2026 high-throughput device study includes PtS2 among promising sliding-ferroelectric semiconductor channels. The 1T-PtS2 monolayer symmetry is inversion-symmetric but lacks the horizontal mirror characteristic of 1H-TMDs. | Heavy Pt and a potentially useful insulating gap make it attractive for a broader spin-orbit response study. | It does **not automatically test** the specific Stage-6 (C_3+m_z) reference-stack selection rule. Treat it as a separate symmetry-class project unless an explicit stacking analysis establishes the needed operations. |
| **C2 — later extension** | **MoS2/WS2 heterobilayer** | A high-throughput study reports 16 semiconducting TMD heterobilayers with sliding ferroelectricity and notes the experimentally reported MoS2/WS2 system. | Strong W SOC plus a semiconducting heterointerface. | Unlike a homobilayer, exchanging layers is not generally a symmetry of the chemical system; lattice mismatch/commensuration also complicates the clean minimal model. Good follow-up, not the clean first test. |

### Source ledger

1. **R-stacked bilayer WSe2, ferroelectric domains and band alignment:** [Nature Communications (2026)](https://doi.org/10.1038/s41467-026-68854-x).
2. **Rhombohedral WSe2 and large SOC splitting:** [Physical Review Materials 10, 044003 (2026)](https://doi.org/10.1103/xs47-7hky).
3. **MoS2 AB/BA switching and NEB barrier:** [Scientific Reports (2019)](https://doi.org/10.1038/s41598-019-50293-y); additional experiments and switching discussion in [Nature Communications (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11685460/).
4. **High-throughput PtS2 and other sliding-ferroelectric channels:** [Nano Letters (2026)](https://doi.org/10.1021/acs.nanolett.6c02536).
5. **Cd2Se2 sliding-ferroelectric/spin-splitting candidate:** [Hao et al., SSRN preprint (2026)](https://doi.org/10.2139/ssrn.6865959). Treat as preprint, not independent experimental validation.
6. **TMD heterobilayer sliding ferroelectric screening:** [ACS Applied Materials & Interfaces (2025)](https://doi.org/10.1021/acsami.4c19017).
7. **Why symmetry class matters for 1T-PtS2 vs 1H-TMDs:** [npj Computational Materials (2022)](https://doi.org/10.1038/s41524-022-00828-1).

## 3. Screening workflow (no DFT yet)

### Gate 1 — structure and symmetry

For WSe2 first, obtain a published or database-backed relaxed monolayer and generate the parallel-layer AA reference plus the AB/BA translations. Run the group-theory toolkit on the **full atomic motif**, not only the Bravais lattice. Confirm:
- AA reference: (C_3), (m_z), and the predicted tensor constraints.
- AB/BA: expected (C_{3v}) symmetry and `Gamma_xz;x` allowed.
- AB and BA are related by a symmetry mapping that reverses the polar domain; distinguish operations within one fixed stack from operations mapping two different stacks.

If this gate fails, do not proceed under the Stage-6 mechanism label; either select another stacking or register a separate symmetry class.

### Gate 2 — electronic feasibility

From literature/database data only, check the existence of a positive band gap and the known SOC/band-edge character. **A literature band gap at one stacking is not enough**: later DFT must verify the full Brillouin zone for AA, AB, BA and intermediate sliding configurations, and ensure the gap stays open at the (B_z\) and (lambda S_i) source corners.

### Gate 3 — response-operator readiness

Before response calculations, require a spinor Wannier representation with (H(\mathbf{k})), same-frame (S_{x,y,z}(\mathbf{k})), physical basis overlaps/position or Berry connection, covariant velocities, and remote-band-window convergence. `hr.dat` plus a band plot alone is a **BLOCKED** input.

### Gate 4 — physical claims stay separate

- Sliding ferroelectricity: two stable polar minima, polarization reversal, and a verified sliding barrier.
- Our current target: gapped (T=0) equilibrium/static `Gamma_xz;x`.
- Jia et al. driven-DC intrinsic NMEE: a separate nonequilibrium/relaxation-limit gate.

Passing one does not prove the others.

## 4. Decision

**Start with 3R-WSe2; keep 3R-MoS2 as the control.** They share the clean 1H-TMD stacking-symmetry route, while WSe2 offers substantially stronger SOC evidence. WS2 is the secondary same-family check. Cd2Se2 and PtS2 should not displace these until their precise structural symmetry and complete-gap feasibility are established.

Next action is **Gate 1 only**: retrieve the exact structures and perform a symmetry/stacking preflight with the existing group-theory toolkit. This is structure analysis, not a DFT response calculation. No candidate has yet passed the actual-material gate.
