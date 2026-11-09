# InSe mirror-domain gate test — prediction before new results

2026-11-09. User authorizes a few identical small gates in the mirror partner,
an independent reliability check, and reading field-driven sliding/FQFE work.
No new material, strain, NEB, switching dynamics, or twelve-state fitting.

Keep BiDB864 A and its verified m001-related B at exactly the existing ionic
coordinates and laboratory axes. Observe delta s_x,l=C_l Ex Bz; F is a separate
perpendicular control field. Positive physical F uses VASP EFIELD=-F/10.

Prediction: C_lower(B,F)=-C_upper(A,-F), and conversely. At zero F, both domains
are lower-negative/upper-positive. A small positive F weakens both magnitudes
in A; it should strengthen both in B. The net gate-induced changes should have
the same sign. The exact identity relates opposite fields, not A(+F) to B(+F);
same-field opposite local trends are a near-zero prediction to be tested.
Do not generate B responses by symmetry and present them as independent DFT.

Five new SCFs: B(-0.10), B(+0.10) on 6x6; B(0), B(+0.10), A(+0.10) on 8x8.
Reuse independent completed A(0) on 8x8 and A/B zero + A gates on 6x6.
The 8x8 positive-gate comparison tests the actual tuning directions, rather
than assuming zero-gate mesh stability proves a derivative. Compare each
finite-field state with its own same-mesh baseline. No claim of fully converged
gain factors. All fields are in V/nm. All calculations retain an insulating gap.

PBE+SOC, In_d/Se PBE.64, ENCUT450, full 6x6 or8x8 mesh, EDIFF1e-9, fixed ions,
neutral slab, dipole correction, 160 spinor bands. Full occupied KS space;
128/144/160 response cutoffs and magnetic +/-0.002 and +/-0.001 eV checks.
Record existing SCF forces for context, without mistaking force cancellation
at a symmetric registry for a switching barrier or domain-wall calculation.

Japan PBS compute nodes only. Three resumable batches B6, B8, A8, each at most
4 CPU/24GB/2h; at most two task allocations active. The B8 batch has two SCFs;
completed states are individually checkpointed. Other shared-account jobs are
not altered. Small local inputs/results only on the Mac. New storage expected
below20GB; latest work UID usage6565700748KiB, numeric caps0/unspecified;
group/project cap and balance unknown. User states Japan is free. First-tier
full SOC+PAW/W90 response chain remains unvalidated; reuse verified Japan chain.

Literature question is distinct: a fixed-registry electronic response does not
show that a real gate leaves the ions/domains fixed. Read Born-charge force,
energy barrier, boundary pinning and domain-wall results before deciding what
mechanical check is needed for InSe. FQFE within one monolayer and relative
sliding of two intact layers must not be conflated.
