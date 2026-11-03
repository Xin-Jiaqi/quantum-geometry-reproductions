# BSF × IGMC in AgF₂: novelty and story gate after the first real-material screening

**Research note, 2026-10-22. Status: PROMISING_HYPOTHESIS / MATERIAL_RESPONSE_UNVERIFIED / NOT_NOVELTY_CERTIFIED.**

## Executive verdict

The project has a coherent physical story worth testing, but **the combination of “sliding ferroelectricity” and “IGMC” alone is not enough to claim high novelty**. Both building blocks are already established in recent papers. The potentially publishable result is a quantitative, source-grounded real-material demonstration that the *magnetic stacking connector* controls the sign parity of a **specific IGMC channel** differently from the accompanying Hall response, and that the prediction survives SOC, magnetic-state, structure, and numerical checks.

The strongest current candidate story is:

> **Connector-selective switching of insulating displacement IGMC in a sliding-ferroelectric antiferromagnetic bilayer.** Sliding reverses \(P_z\), but it does not reverse every magnetic response. The unitary Type-I and antiunitary Type-II interlayer connectors produce different sign rules for the displacement IGMC and the anomalous Hall response; strain changes which connector class is relevant.

This is a symmetry hypothesis grounded in published ingredients, not yet a calculated AgF₂ IGMC signal. The next cheap material gate is a SOC-on static B1 calculation, followed by a response-theory audit before committing to Wannier transport work.

## 1. Evidence already obtained

### Real-material B1 screening

From PR #38 commit 247af15662c160b32e142cb2ecc0ae982aa7af3f:

- 11/11 ionic steps electronically converged; VASP stopped normally at the screening force threshold, max force \(0.0782\ \mathrm{eV/\AA}\).
- Type-I Ag moments \([-0.567,+0.567,+0.566,-0.566]\ \mu_B\) and zero net moment retained.
- Minimum sampled interlayer F–F distance \(3.002\ \AA\); occupation-defined gap \(1.3875\ \mathrm{eV}\) on the nonconverged \(12\times8\times1\) screening mesh.
- One local WSL run, 8 ranks × 1.62 h. No new PBS work.

This is enough to move from “can the initial structure hold a plausible AFM state?” to “does the SOC electronic state support the intended response channel?” It is **not** a converged local minimum, easy-axis determination, magnetic ground-state comparison, Berry-phase polarization, barrier, or material response.

### Existing toy-model evidence in this repository

The independent eight-band model already demonstrates that explicitly sliding-dependent real-space interlayer hopping can produce a nonzero insulating displacement-IGMC geometric kernel and different B1/B2 parity for unitary \(m_z\) versus antiunitary \(m'_z\) connectors, without inserting a phenomenological layer-odd onsite potential. It explicitly lacks a verified total-energy double well and physical units. Therefore, another toy-only parity plot would add little; the real-material channel amplitude and stability are the missing information.

## 2. Nearest literature and what remains distinct

| Work | What is already established | Consequence for our novelty claim |
| --- | --- | --- |
| Xiang et al., PRL 134, 116301 (2025), [DOI](https://doi.org/10.1103/PhysRevLett.134.116301) | Spin-Zeeman quantum geometry (Zeeman QGT) produces conduction and displacement IGMC in SOC-coupled materials; symmetry classification and model demonstrations. | Do not claim to discover Zeeman-QGT IGMC or its general symmetry classification. A new result must be the switchable material/channel prediction and quantitative material evidence. |
| Zhu et al., PRL 135, 056801 (2025), [DOI](https://doi.org/10.1103/dmzg-ck2t) | AgF₂ and Fe₂MoSe₄ bilayers exhibit sliding-ferroelectric control of unconventional magnetism, spin polarization and AHE through a stacking connection operator. AgF₂ B1/B2 polarization and sliding barrier are source-reported. | Do not claim the first sliding-ferroelectric magnetic switching in AgF₂. Their work does not report Zeeman-QGT IGMC. |
| Chakraborti, Ghosh & Nandy, arXiv:2508.14745 (2025), [preprint](https://arxiv.org/abs/2508.14745) | ZQGT-driven IGMC in 2D unconventional magnets, including \(d\)-wave and mixed magnetic orders, with conduction/displacement channels and symmetry signatures. | Theoretical IGMC in 2D unconventional magnets is already an active direction. We need a connector-resolved switching prediction and a material calculation, not another symmetry-only classification. |
| Ghorai, Sarkar & Agarwal, PRB accepted 2 Sep 2026, [DOI](https://doi.org/10.1103/xr5v-gkv5), [arXiv](https://arxiv.org/abs/2601.04787) | An orbital-origin contribution to intrinsic gyrotropic magnetic current, including an orbital displacement mechanism, demonstrated for CuMnAs; response can track Néel reversal. | A spin-Zeeman-only calculation may miss a competing or additive orbital channel. Keep spin and orbital magnetic driving explicitly separated. |
| Chakraborti et al., *Gyrotropic Fingerprints of Magnetic Topological Insulator–Unconventional Magnet Interfaces*, Advanced Quantum Technologies (first published 30 Aug 2026), [DOI](https://doi.org/10.1002/qute.70427) | Zeeman-QGT IGMC at interfaces as a symmetry fingerprint of unconventional magnetic order. | Symmetry fingerprints and magnetic-order diagnosis are also crowded. Our differentiator must be ferroelectric connector switching in a real bilayer, not “IGMC is a probe of magnetism.” |
| Xiang, Jin & Wang, PRL 135, 146303 (2025), [DOI](https://doi.org/10.1103/14mp-263q) | A broader spin-quantum-geometry framework linking spin QGT and Zeeman QGT to spin transport and Edelstein responses. | Do not conflate spin QGT, Zeeman QGT, conventional Edelstein response, and IGMC. Define the measured current and magnetic-drive operator explicitly. |

## 3. The actual predicted fingerprint

The relevant observable is the AC current response to an oscillating magnetic field,
\[
J_i(\omega)=\left[\sigma^C_{ij}(\omega)+\sigma^D_{ij}(\omega)\right]B_j(\omega),
\]
not just the static magnetoelectric polarizability \(\alpha^{PB}_{ij}=\partial P_i/\partial B_j\). The static-like polarization derivative is related most directly to the displacement channel; it cannot replace the full channel-resolved IGMC tensor. See the companion correction note, B1_IGMC_TENSOR_ERRATUM_20261022.md.

Under the project's channel-resolved symmetry convention, conditional B1 magnetic point group \(m'_x\) allows

\[
\sigma^C=
\begin{pmatrix}
0&c_{xy}&c_{xz}\\
c_{yx}&0&0\\
c_{zx}&0&0
\end{pmatrix},\qquad
\sigma^D=
\begin{pmatrix}
d_{xx}&0&0\\
0&d_{yy}&d_{yz}\\
0&d_{zy}&d_{zz}
\end{pmatrix}.
\]

The most experimentally legible **insulating displacement-channel** prediction is the \(D_{xx}\) sign. The source-defined partner connectors give:

| Domain comparison | Type-I: unitary \(m_z\) | Type-II: antiunitary \(m'_z\) |
| --- | --- | --- |
| \(P_z\) | reverses | reverses |
| \(D_{xx}\) | **reverses** | **same sign** |
| AHE \(\sigma_{xy}\), if nonzero at the chosen chemical potential | **same sign** | **reverses** |

This crossed parity is the strongest candidate story: it separates the switching of the insulating displacement IGMC from the AHE response, and distinguishes two magnetic-connector classes. Type-II in the source corresponds to a strained branch, reported around 4% biaxial tensile strain; it must not be presented as a freely accessible second zero-strain magnetic state.

**Important response-theory gate:** the repository intentionally retains both the channel-resolved C/D transformation rule and the printed PRL Eq. (14) alternative in src/igmc_lab/tensor_symmetry.py. The source-level convention must be reconciled against the main text and Supplemental Material before a final symmetry claim. No final paper should silently select whichever convention yields the desired sign.

## 4. What would make the story publishable?

### Minimum material evidence

1. **SOC state:** SOC-on B1 is self-consistent, remains insulating or has a clearly identified Fermi-level regime, and yields a defensible magnetic point group. The current non-SOC collinear run does not determine a physical easy axis.
2. **Real domain pair:** B1/B2 geometries are source-grounded and mapped in the same laboratory axes. The published paper reports a BSF pair, but our current B1 geometry is C2DB-derived and screened only; it is not the authors' relaxed POSCAR.
3. **Channel-resolved response:** compute \(\sigma^D_{xx}(\omega)\) and, where physically relevant, \(\sigma^C_{xy}(\mu,\omega)\) from a validated SOC Wannier Hamiltonian with same-gauge position and spin operators. Keep spin-Zeeman and orbital magnetic couplings separate.
4. **Actual switching:** the response parity is confirmed from both partner Hamiltonians, not inferred solely from a tensor table. If Type-II is used, include the source-relevant strain and verify the magnetic state at that strain.
5. **Numerical trust:** converge k mesh, frequency/broadening, occupied-band summation and operator gauge covariance; report physical units and a response magnitude. Distinguish an insulator's in-gap displacement current from a conduction channel that may vanish in the gap.
6. **Useful comparison:** compare \(D_{xx}\) against AHE or another independent magnetic readout; establish whether their different switching parities are large enough and robust enough to be experimentally useful.

### Novelty judgment

- **Symmetry-only story:** low-to-moderate novelty; likely too incremental by itself.
- **Validated toy model plus correct channel-resolved group theory:** useful theory note, but this is already largely in the repository.
- **SOC/Wannier material prediction with converged finite-frequency amplitude, cross-parity switching, and spin/orbital decomposition:** a credible publishable materials-physics story, subject to the result size and literature overlap.
- **High-impact claim:** not justified until a large/robust response or a genuinely new physical mechanism is demonstrated.

## 5. Next steps, ordered by information gain per compute hour

**Gate 0 — now, no VASP:** preserve the tensor erratum and settle the Eq. (14) convention against the source; keep the toy model as a control, not material evidence.

**Gate 1 — one SOC-on static B1 screening point:** use the actual screened B1 coordinates and locked private POTCAR, run the noncollinear VASP executable with a declared Type-I x-axis initial spin pattern, no ionic relaxation, no automatic retry. Check SCF, site-resolved spin/orbital moments, occupation gap, and whether the prepared state remains physically meaningful. This does not establish the easy axis or MAE; one orientation is only a screening state. Full task is in B1_SOC_STATIC_SCREENING_TASK_20261022.md. **Not executed or authorized by this report.**

**Gate 2 — zero-cost domain map:** generate a symmetry-image B2 POSCAR from B1 only for group/operator mapping, clearly labelled NOT DFT-relaxed; check the actual connector and channel parity. Do not call this a second stable minimum.

**Gate 3 — only if Gates 0–2 pass:** SOC/Wannier Hamiltonian and same-gauge spin/position/B operators, then low-cost coarse \(D_{xx}\) response. If the response is symmetry-forbidden, numerically negligible, gap-unstable, or not convergent, stop before large k meshes and NEB.

**Gate 4 — later, separately authorized:** verify structural double-well/barrier and Type-II strained branch if the scientific claim requires them; quantify finite-frequency response and orbital contributions.

## Decision

**Continue the brainstorm, but do not yet promote it to a materials discovery.** The best candidate story is connector-selective, channel-resolved IGMC switching. The immediate bottleneck is no longer the non-SOC geometry; it is whether SOC and the correct IGMC response theory produce a real, converged, measurable signal in the material.

No VASP, PBS, or Wannier job was launched in preparing this report.
