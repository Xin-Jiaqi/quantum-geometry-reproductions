# BD27 | A registry-torus response field: six symmetry-forced zeros and a conditional index budget

**2026-10-22 · INDEPENDENT THEORY/MODEL RESULT · NOT a material result or a new theorem claim.**

## 1. Research motivation and direct competitors, correctly scoped

The 2023 [Ji–Yu–Xu–Xiang *General Theory for Bilayer Stacking Ferroelectricity*, PRL 130, 146801](https://doi.org/10.1103/PhysRevLett.130.146801) uses the six special hexagonal-registry points **G(0,0), A(1/2,0), B(0,1/2), C(1/2,1/2), M(1/3,2/3), N(2/3,1/3)** (Fig. 2(b)), enumerates their allowed stacking symmetry/polarizations, and is direct/foundational prior art for **registry-space layer symmetry**. The 2026 [Zhou et al. *Layer Edelstein Effect*, PRL 137, 076401](https://doi.org/10.1103/qhvs-pfkq) is **direct** for general current-driven layer-odd spin and the 80-LG bilayer response classification. Whether its supplementary tables/equations already imply the six-location and winding-index rules for the exact observable is still **OPEN**. [Agarwal et al., *Shift photocurrent vortices from topological polarization textures*, npj 2D Mater. Appl. 9, 43 (2025)](https://doi.org/10.1038/s41699-025-00552-6) studies another observable and physical drive; it is **highly relevant Introduction/nearest texture precedent**, **NOT evidence by itself that the current-induced layer-odd spin-response field result is directly published**. Likewise [2025 experimentally observed polar merons](https://doi.org/10.1038/s41467-025-60647-y) are important motivation but a distinct observable.

The intended possible S2 result is an **operator-specific global registry-space constraint** that cannot be formulated merely as an isolated stack's point-group tensor. This test asks what such a claim could look like. The baseline Neumann-covariance result and Poincaré–Hopf theorem are classic mathematics; **their combination is not a new general symmetry principle**.

## 2. Precise physical assumptions and covariance

Consider **two identical, equally oriented 1H-type hexagonal monolayers with compatible vertical-mirror/threefold symmetries and a horizontal layer-exchange transformation**. The interlayer translation \(\boldsymbol\delta=x\mathbf a_1+y\mathbf a_2\) is periodic modulo \(\Lambda\), so the configuration space is the torus \(\mathcal T^2=\mathbb R^2/\Lambda\). Choose
\[
\mathbf a_1=(1,0),\qquad
\mathbf a_2=(-1/2,\sqrt3/2).
\]
The **measured operator** is the layer-odd out-of-plane axial spin \(S_-^z=S_T^z-S_B^z\), and the **source** is uniform in-plane charge current \(\mathbf j\), with an explicitly fixed *current-driven nonequilibrium response protocol* (not the equilibrium DC response). Represent
\[
\delta\langle S_-^z\rangle=\mathbf F(\boldsymbol\delta)\cdot\mathbf j.
\]
The response is a real, continuous, periodic, nonsingular 2D covector, identified with a planar vector using the physical Euclidean metric. For finite-\(\omega\) complex response, independent real/imaginary zeros, causal conditions, and index classification require new analysis; this note applies only to a real DC coefficient.

Let \(C_{3z}\) be an ordinary proper rotation, \(m_v\) an in-plane line mirror, and \(m_z\) the **horizontal layer-swapping mirror**. Their induced actions are:
\[
\mathbf F(C_3\boldsymbol\delta)=C_3\mathbf F(\boldsymbol\delta),\quad
\mathbf F(m_v\boldsymbol\delta)=-m_v\mathbf F(\boldsymbol\delta),\quad
\mathbf F(-\boldsymbol\delta)=-\mathbf F(\boldsymbol\delta).
\]
The second negative sign is due to out-of-plane **axial** spin under a vertical mirror; the third comes from the **layer permutation** under \(m_z\), while in-plane current remains even. These are full-family *mapping* rules; \(m_z\) is not itself a symmetry of a generic individual registry. At a fixed registry, the operation may require a **Seitz in-plane translation** to map the two infinitely periodic layers. A concrete material must pass its actual Seitz/site/orbital test.

**Key limitation:** conductivity, fixed-\(j\) versus fixed-\(E\), the Fermi surface/relaxation protocol, and spin operator embeddings must remain physically consistent across all configurations. No real-material response has been calculated.

## 3. Forced six-zero set, using established BSF labels

There are four solutions to \(2\boldsymbol\delta\in\Lambda\) modulo \(\Lambda\) (the layer-swap involution's fixed points): G, A, B, C. At each, \(\mathbf F=-\mathbf F\), so \(\mathbf F=0\).

There are three \(C_3\)-fixed points, obtained by \(C_3\boldsymbol\delta\equiv\boldsymbol\delta\pmod\Lambda\): G, M, N. At each, \(\mathbf F=C_3\mathbf F\) implies \(\mathbf F=0\) because a 120° rotation has no nonzero real fixed in-plane vector.

The overlap is G, so the **union contains exactly six symmetry-enforced zero locations: G/A/B/C/M/N**. This does **not** say the full physical response has exactly six zeros; additional zeros are possible, and zeros might be degenerate. It also does not claim those locations are low-energy ferroic states. These are the same registry coordinates given in BSF Fig. 2(b), now with a separate *response-specific* constraint.

## 4. A fully explicit *allowed response field* and independently checked index assignment

Use fractional registry \(\mathbf x=(x,y)\), the above Cartesian basis matrix \(A\), and \(J=\begin{pmatrix}0&-1\\1&0\end{pmatrix}\). Define a simple *symmetry-invariant scalar auxiliary function*
\[
\Phi(x,y)=\cos(2\pi x)+\cos(2\pi y)+\cos(2\pi(x-y)),
\quad \boxed{\mathbf F(\boldsymbol\delta)=J\nabla_{\rm cart}\Phi.}
\]

This is only a **Fourier witness of symmetry compatibility**; no physical response is known to be a gradient of this scalar, no \(\Phi\) is claimed to be the stacking energy, and its coefficients do not come from DFT.

The direct script [\`stacking_response_torus_gate.py\`](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/stacking_response_torus_gate.py) (NumPy only) was run and verifies:
- \(C_3\) ordinary vector covariance, \(m_v\) **axial mirror** covariance, \(m_z\) layer-odd antisymmetry, and Bravais translation periodicity on 50 random points;
- exactly the six *forced* G/A/B/C/M/N zeros with numerical field norms below \(10^{-12}\);
- nondegenerate local Jacobian determinant signs (winding/Poincaré indices)

| BSF point | Fractional registry | Why it vanishes | Index **in this witness** |
| --- | --- | --- | ---: |
| G | (0,0) | \(m_z\) and \(C_3\) | +1 |
| A | (1/2,0) | layer swap (with lattice translation) | -1 |
| B | (0,1/2) | layer swap (with lattice translation) | -1 |
| C | (1/2,1/2) | layer swap (with lattice translation) | -1 |
| M | (1/3,2/3) | \(C_3\) | +1 |
| N | (2/3,1/3) | \(C_3\) | +1 |
| **Total** | | | **0** |

A small symmetry-preserving second Fourier harmonic \(\Phi\mapsto\Phi+0.05\sum_{\bf g}\cos(4\pi{\bf g}\cdot{\bf x})\) passes the same tests and preserves the simple-zero indices. This shows robustness to **one specific small perturbation**, not arbitrary microscopic perturbations. The full first-harmonic model has precisely six stationary points (one can derive this from \(\sin(2\pi x)=-\sin(2\pi y)= -\sin(2\pi(x-y))\)); a numerical root scan gives the same six roots, but no experimental response is modeled.

## 5. Stronger general mathematical consequence and its limits

On the torus, a smooth planar response vector \(\mathbf F\) is a section of the tangent bundle under the Euclidean-metric identification. If its zeros are isolated, the classical **Poincaré–Hopf theorem** requires
\[
\boxed{\sum_{p:\mathbf F(p)=0}\mathrm{ind}_p(\mathbf F)=\chi(\mathcal T^2)=0}.
\]
Further, if a zero at a \(C_3\)-fixed registry is *nondegenerate*, \(D\mathbf F\) commutes with a 120° rotation and has the form \(aI+bJ\), whose determinant is \(a^2+b^2>0\). Thus G, M, N necessarily each have **index \(+1\)**, under the stated regularity/nondegeneracy conditions. A, B, C are permuted by \(C_3\), and have the **same** index if simple. If there are **no additional zeros**, the index sum forces A/B/C to be \(-1\) saddles. If additional zeros occur, their total index must compensate the six forced points, and no **universal exactly-six count** or universal half-translation saddle index follows.

These deductions are **standard index theory plus group covariance**, not yet a publishable novelty. If a physical response is nonanalytic because Fermi surfaces reconstruct as registry changes, the smooth-torus premise may fail; finite temperature and controlled lifetime may be needed. If the response is complex at AC frequency, simply treating \(\mathbf F\) as a real tangent vector is invalid.

## 6. Why this might grow beyond an isolated example

The natural proposed generalization is to treat a physically defined response tensor/covector as an **equivariant section over the registry torus** and classify, for each real layer/site output representation and drive representation:
1. symmetry-fixed forced zeros and zero-set orbits;
2. leading allowed Taylor order at each orbit;
3. index/charge-budget constraints for smooth **real two-component** responses;
4. when actual Seitz layer swaps, antiunitary Onsager, and nonsymmorphic phases modify those statements.

This is S2-level methodology if it yields a reproducible classification that is not already explicitly provided by published stacking/LEE theory. It does **not** become a new physics theorem just by combining textbook facts. To reach a strong general story, it needs a response-specific, previously unarticulated selection law, a physically nonzero model \(\mathbf F(\delta)\), and a practical measurement in a plausible bilayer or small-angle moiré setting.

**Priority of the next tests:** (i) Compare with the 2026 LEE Supplemental Material's complete 80-LG response tables and the 2023 BSF Fig. 2(b)/Table S7; (ii) derive the general fixed-locus and index classification for **one** alternative layer group with a *different* symmetry orbit, not a blind 80-LG scan; (iii) separately complete the source-backed MoS\(_2\) orbital SK + SOC model to test if \(S_-^z\) can be actually nonzero in a fixed-\(j\) Boltzmann/Kubo calculation. Reject apparent “new vortices” if they merely restate the published polarization or shift-current textures.

**Science status:** G3 symmetry/test **PASS**; six forced zero locations and the conditional index budget **DERIVED**; microscopic nonzero current-induced spin response and general novelty **UNPROVEN**; **NO new DFT/HPC jobs**.
