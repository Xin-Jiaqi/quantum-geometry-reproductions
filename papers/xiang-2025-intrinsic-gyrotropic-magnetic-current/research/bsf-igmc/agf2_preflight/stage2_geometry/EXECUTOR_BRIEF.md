# Minimal authorized handoff for an external VASP harness

This directory is already on GitHub. It contains **unrelaxed** C2DB-derived starting coordinates and a verified symmetry relation, **not** author-relaxed structures or DFT results.

**Do not submit computational jobs** until the researcher has independently provided a licensed VASP installation, private Ag/F POTCAR hashes and versions, cluster account/budget and specific execution approval. You may read and validate this folder now.

Once the physical resource gate is explicitly cleared, the ONLY first campaign is 0%-strain Type-I AgF2 **non-SOC, spin-polarized fixed-cell ionic relaxation** from the three **B1** vertical F-F interface-plane separation guesses (2.6, 3.0, 3.4 A); candidate B2 files serve as magnetic mirror partners for independent checks. For Stage-1 Type-I the unitary m_z connector flips in-plane axial spin, hence interlayer Neel vectors are opposite. The geometry is free to slide unless constraints are explicitly reported.

Log all software parameters, computational provenance, converged energies/forces/local magnetic moments, and failed starts. Use a new results folder, do not alter original IGMC reproduction, do not invent absent values, do not start SOC/NEB/Wannier/IGMC yet. Stop if G0 (CALCULATION_GATE.md) cannot be met. Return a lightweight JSON receipt and real output hashes, and wait for independent physics review.
