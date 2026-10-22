# Scientific pivot: from single MnPS3 system to layer-group / spin-layer-group response addressability

**Research memo, 2026-10-22. Status:** mathematical negative control PASSED; useful new paper-level effect UNPROVEN; full 80-LG/528-MLG exhaustive mode screening **NOT STARTED**. **No new DFT or HPC task.** The accepted MnPS3 R1 remains closed with no SCF convergence; the R2-P0 private four-site review is read-only and not yet returned.

## 1. Why this is scientifically better than indefinitely fitting MnPS3

The existing toolchain \`Xin-Jiaqi/group-theory-operations-toolkit\` v0.18.0 already has **80 crystallographic layer groups, 528 magnetic layer groups, 122 magnetic point groups, magnetic operations, layer-stack FE equivalence, layer-resolved Edelstein/IMHE and invariant tensor screens**. Those are useful foundations but not evidence of a new theorem. Current group theory is highly competitive:

- [Zhu et al., PRL 135, 056801 (2025)](https://doi.org/10.1103/dmzg-ck2t): **already survey 80 layer groups** and magnetic FE switching connection operators.
- [Du et al., PRB 113, 214455 (2026)](https://doi.org/10.1103/g471-xrkl): **448 collinear spin layer groups** and general stacking solutions. These are NOT the same mathematical objects as **528 crystallographic magnetic layer groups**.
- [Zhang et al., arXiv:2605.25484v2](https://arxiv.org/html/2605.25484v2) PREPRINT: **40,087 spin layer groups** under \(I_k\le12\), including 448 collinear, 6083 coplanar-noncollinear, 33556 noncoplanar, plus small coreps. No blanket novelty in enumerating SSG layer-group irreps.
- [Chen et al., Nature 640, 349 (2025)](https://doi.org/10.1038/s41586-025-08715-7) already classify unconventional magnon spin-space-group representations.
- [Sun et al., Nat. Commun. 17, 3101 (2026)](https://doi.org/10.1038/s41467-026-69635-2) classify 3 FE altermagnet spin-band switching types.
- [Ni et al., arXiv:2601.16052](https://arxiv.org/abs/2601.16052) PREPRINT screens magnetic layer groups for electric chiral-magnon splitting; [Xie et al., arXiv:2601.15031](https://arxiv.org/abs/2601.15031) PREPRINT reports a chiral-magnon response survey using **464 collinear spin layer groups** (definition differences relative to 448 NOT reconciled here); [Liu et al., Nano Lett. 26, 11017 (2026)](https://doi.org/10.1021/acs.nanolett.6c02164) explicitly predicts FE switching magnon geometry and thermal Hall.

**Consequence:** another tensor-allowed list, or generic FE-controlled magnon chirality/Berry curvature, does NOT pass novelty. The candidate must ask an **unanswered physical response question**, not advertise group-table size.

## 2. New discriminating concept: FE-domain *magnetic-mode addressability*

Fix a real physical magnetic site orbit and local-moment normal-mode vector space \(V_Q\), and a driving field space \(V_E\). For an individual FE/magnetic state with stabilizer \(H\), its symmetry-allowed static force-to-field matrix satisfies
\[
D_Q(h)K=K D_E(h),\quad h\in H,\qquad K\in\mathrm{Hom}_H(V_E,V_Q).
\]
The input **must include physical magnetic-site moment orientation and mode representation**; ordinary layer/point-group names alone do not determine these.

**First HARD NO-GO (verified):** if the two **full SOC magnetic states** are exact unitary symmetry partners \(c\), their normalized responses are related by
\[
K_B=D_Q(c)K_A D_E(c)^{-1}.
\]
Thus **rank, singular values and the dimension of the allowed intertwiner space** cannot change. Fixed-lab \(E_x\) may look OFF→ON under domain rotation, but the corresponding mapped field/mode has unchanged physical strength. This is an **elementary representation-theory theorem**, not a proposed new effect.

The [independent gate](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/domain_covariance_gate.py) actually checks one mirror/C4 rotated partner, Hom dimensions, invariant rank/singular values and an apparent laboratory-axis OFF→ON that is purely reorientation. It also verifies a representative nonsymmorphic glide \(g(k)^2=e^{-ika}\), giving eigenvalues \(\pm i\) at a zone edge. **NO full 80/528 catalog was run.**

This cleanly partitions the most promising classes:
- **0. Trivial conjugate domains:** exact full magnetic symmetry, response only reorients; **do not call an intrinsic new channel**.
- **1. Structural but not SOC-magnetic conjugacy:** a stacking operation pairs ionic structures or nonrelativistic spin-group states, but not *full SOC* magnetic Hamiltonians; now different mode couplings are permitted but must be derived/verified (no magnitude follows from symmetry).
- **2. Different magnetic stabilizers/induced site representations:** switching connects inequivalent magnetic isotropy types or physical moment structures; \(\dim\mathrm{Hom}\) can differ. This may offer robust mode-selective zero vs allowed contrasts, but does not guarantee FE degeneracy or dynamics.

**Candidate physical question:** Can a sliding FE path change the **number or symmetry label of independently electrically addressable, physically real AFM modes**, beyond a change in measurement axes and beyond prior spin/FE switching classification? The answer is currently **unknown**, and most easily tested as a group/no-go exercise.

## 3. Why all-layer groups may actually matter at nonzero wavevector

At \(q=0\), many uniformly driven response tensors collapse to magnetic point-group invariants; a blind 528-MLG scan likely adds no new *tensor*-level result. **Full LG/MLG** are needed for true Wyckoff-site orbit permutations, layer exchange, induced magnetic-mode representations, and at **finite \(q\)** for nonsymmorphic fractional translations.

A glide with \(g^2=T_a\) obeys \(\rho_q(g)^2=e^{-iqa}\). At \(q=\pi/a\), eigenvalues \(\pm i\) demand *little-group/projective* treatment, not naive mirror parity. There may be **momentum-selective layer/sublattice magnetic-mode addressability** missing from ordinary uniform-field lists. But ordinary light at \(q\approx0\) cannot directly inject zone-edge momentum into an ideal translationally invariant sample. A useful story needs a real drive: structured/near-field electric field, piezo/phonon mode, supercell folding or another explicitly justified momentum source. Existing projective SSG and magnon band theory is prior art, not proof this response is new.

## 4. Concrete high-return plan with kill tests (no VASP)

**Gate G0 (completed model proof):** response domain conjugacy/equivariance null test and nonsymmorphic \(g^2=T_a\) phase sanity test. Local NumPy model pass; GitHub group-toolkit [Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32) submitted for mandatory CI/review, NOT merged into default branch. No full enumeration and no novelty.

**Gate G1 (next, cheapest and decisive):** choose **two simple, independently specified magnetic site-orbit toy lattices** with a genuine magnetic order representation, not arbitrary K-matrices: (a) exact unitary FE-domain partner (negative control), (b) a pair where the structural FE operation is **not** a true SOC magnetic mapping or where stabilizers differ. Compute stabilizers, layer-site permutations, induced transverse magnetic-mode representations, and allowed field-to-mode \(\mathrm{Hom}\) spaces. Show **a strict invariant distinction** beyond rotated coordinates, or kill the idea. No need to enumerate 528 types before this.
 
**Gate G2 (conditional):** if G1 survives, calculate mode-addressability equivalence classes across **selected representative layer groups** and relevant Wyckoff magnetic sites; handle antiunitary coreps rigorously before expanding scope. Do not pretend a 528 table can be populated without ordered magnetic sites and mode bases.

**Gate G3 (conditional):** explore 1–2 nonsymmorphic finite-\(q\) exemplars with physically motivated external momentum, rather than another generic 2D magnon band-nodal catalog.

**Gate G4 (only if physically distinct and feasible):** *then* a concrete material with a small validated response computation. Do not reactivate MnPS3 R1 merely because this group-theory idea exists.

## 5. Scope / danger statements

A magnetic point group classification tells only **allowed**, not **nonzero**, **magnitudes**, **energies**, **stable FE**, or **experimental detectability**. The number 528 refers to *magnetic layer groups* and cannot be substituted for the spin layer-group lists (448 collinear; 40,087 under a specified supercell cutoff). Antiunitary operations in causal finite-frequency Kubo generally require complex-conjugation/time-reversal/frequency care; do **not** apply the unitary static singular-value statement unmodified to arbitrary antiunitary AC responses.

**Scientific success criterion:** a demonstrably different domain-/site-/mode-addressability selection principle not reducible to established group-covariant tensor rotation, spin-layer-group irrep classification or published electric-magnon switching. A useful **negative theorem and precise no-go catalog** is a legitimate outcome but not, without additional insight, a stand-alone high-impact physics paper.

## 6. Github location / responsibilities

- Dedicated exploratory group theory code/memo: [group-theory-operations-toolkit Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32) (protected main; no forced commits or CI bypass).
- This memo: independent physics focus under \`quantum-geometry-reproductions/.../research/brainstorm-20261021/\`, not original PRL reproduction.
- Current R2-P0 is a **separate approved READ-ONLY private MnPS3 site/SCF review**; it can supply physical site data later but does not make MnPS3 a required material or allow any job submission.
