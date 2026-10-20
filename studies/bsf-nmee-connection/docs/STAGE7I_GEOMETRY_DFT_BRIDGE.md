# Stage 7i: Static NMEE as mixed Berry curvature; first-principles input gate

**Status:** independently validated eight-band synthetic model; zero-temperature, insulating **equilibrium/static** spin–electric Zeeman derivative only. NOT published driven-DC intrinsic NMEE, ferroelectric bistability, or a DFT material prediction. The prior Stage 7a–h research and frozen paper reproduction are preserved.

## 1. A useful general physics bridge

Introduce an auxiliary spin source \(\lambda_i S_i\) and spin-only vertical Zeeman perturbation \(-g\mu_B B_z S_z\) into the original eight-band Hamiltonian. In its *existing* Pauli spin convention, for the gapped occupied projector \(P(k,\lambda_i,B_z)\),

\[
\mathcal F_{\lambda_i k_j}
=i\,\mathrm{Tr}\!\left(P[\partial_{\lambda_i}P,\partial_{k_j}P]\right),\qquad
\boxed{\Gamma_{iz;j}=\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}
\left.\partial_{B_z}\mathcal F_{\lambda_i k_j}\right|_{\lambda=B=0}.}
\]

This is the same occupied-to-empty static Kubo kernel
\(-2\operatorname{Im}\sum_{v,c}(S_i)_{vc}(v_j)_{cv}/(\epsilon_v-\epsilon_c)^2\), but now a **mixed parameter-space Berry curvature of the entire occupied vector bundle**. It is insensitive to arbitrary unitary mixing within occupied bands. It also directly connects to the modern theory of electric polarization differentiated with respect to a spin source (subject to the physical charge/spin/g-factor convention).

This standard geometry is not a new universal transport effect. In particular, a gapped **equilibrium/static** spin–electric susceptibility derivative is distinct from current-driven Edelstein spin accumulation and from the unresolved nonequilibrium **driven-DC** coefficient in Jia et al.

## 2. Independent Wilson-plaquette verification (no Kubo denominators)

Construct four occupied frames at corners \((k_j\mp\delta k/2,\lambda_i\mp\delta\lambda/2)\), at fixed \(B_z\). The oriented determinant Wilson loop \(W=\prod_{\square}\det(U_a^\dagger U_{a+1})/|\det|\) gives

\[
\mathcal F_{\lambda_i k_j}\simeq\mathrm{Arg} W/(\delta k\,\delta\lambda),\quad
\Gamma_{iz;j}(k)\simeq
\frac{\mathrm{Arg}[W(+\Delta B)W(-\Delta B)^*]}
{2\Delta B\,\delta k\,\delta\lambda}.
\]

A gauge change \(U_a\to U_a Q_a\) with independent \(Q_a\in U(4)\) at each corner leaves the closed Wilson determinant unchanged. There is no individual-band gauge fixing or occupied-subspace near-degeneracy division.

For two K/K' disks of radius 0.5 nm^-1, 22×60 quadrature nodes/valley and \(\Delta B=0.25\) T:

| Slide u | t_perp | Kubo \(\Gamma_{xz;x}\) | Wilson \(\Gamma_{xz;x}\) | Relative error |
|---:|---:|---:|---:|---:|
| .1 | 5 meV | -6.916903698e-7 | -6.916308508e-7 | 0.008605% |
| .4 | 5 meV | +2.164836066e-6 | +2.164701420e-6 | 0.006220% |
| .4 | 20 meV | +9.843337225e-5 | +9.842167207e-5 | 0.011886% |

Halving both source-loop steps cuts the error about fourfold, consistent with a centered quadratic finite-loop error. Every in-plane component \((S_{x/y},k_{x/y})\), the slide-reversal sign and occupied \(U(4)\) gauge invariance were checked independently. **These checks are two local valley disks, not a rigorous full-BZ error bound.** An exploratory direct *uniform whole-BZ* Wilson-loop evaluation was insufficiently converged near sharp valley hotspots and is explicitly **not** used as supporting evidence.

## 3. Portable spatial symmetry filter

A spin response index \(i\) and Zeeman index \(a\) are both **axial**, whereas electric-field index \(j\) is **polar**. Their two determinant factors cancel, leaving an ordinary rank-three polar spatial law:

\[
\Gamma'_{ia;j}=R_{ii'}R_{aa'}R_{jj'}\Gamma_{i'a';j'}.
\]

For \(\Gamma_{xz;x}\) in a *fixed stacking*, inversion and horizontal \(m_z\) each forbid the coefficient; C3 alone does not. Time reversal alone does **not** forbid \(\partial_{B_z}\chi_{x;x}\), because both the spin and B field are odd under time reversal. Do not put operations mapping **different sliding domains** into the little group of a single displaced stack.

If the original reference stack has C3 plus a domain-exchanging mirror, the **scalar trace** of the in-plane sliding tensor is allowed first at cubic displacement order, whereas its traceless symmetric (nematic) part can be linear; sign and magnitude are microscopic, not symmetry-determined.

## 4. What is truly needed from DFT/Wannier

For a candidate **gapped spinor-SOC bilayer**, before any quantitative response claim require:

1. Structural stacking geometry, little-group point operations of *each fixed stack*, possible slide-partner map, full-BZ insulating occupied count, chemical potential and numerical gap guards under all B/λ source steps.
2. Spinor Wannier \(H(k)\) **and full nonconserved SOC spin matrices** \(S_{x,y,z}(k)\); lattice, Wannier/orbital positions and *basis overlap/Berry connection*, or physically covariant velocities. **A bare hr.dat or band energies alone is insufficient.**
3. In a k-dependent Wannier basis, a consistent connection convention: with \(\mathcal A=i\langle e|\partial_k e\rangle\), \(D_kH=\partial_k H-i[\mathcal A,H]\) rather than differentiating a gauge-changed Hamiltonian naively. For a Wilson loop, the link must include the basis overlap between k and k+dk, not just ordinary eigenvector dot products.
4. Independent Wilson plaquette versus static Kubo/projector response under the *same* spin and Zeeman normalization, and convergence under source steps, adaptive k-mesh, Wannier remote-band window, and spatial symmetry operations. K/K' hotspots are model-specific: a material study must scan **the whole BZ** for other hotspots.
5. Careful 2D **sheet-unit** conventions and spin source (\(\sigma\) versus \(\hbar\sigma/2\)); no arbitrary use of supercell vacuum thickness. The modeled B is **spin-only Zeeman**, not the orbital electromagnetic field.

**Separate claims must stay separate:** sliding *ferroelectricity* requires stable pair of minima, polarization reversal and a barrier; *driven DC NMEE* requires its own Keldysh/density-matrix zero-frequency and relaxation protocol. Neither follows automatically from the present static calculation.

## 5. Research implication and next gate

A testable, still model-specific hypothesis: spin-orbit coupling plus nontrivial same-/cross-sublattice interlayer hybridization and virtual remote-band dressing can yield a sliding-linear **nematic** and a sliding-cubic **scalar** contribution that compete, even while the insulator stays gapped. C3 plus layer-exchange symmetry restricts powers; Rashba helicity and complex matrix elements determine their signs.

The immediate DFT preparation is **not** to launch VASP. First, use point-group screening and a spinor-Wannier operator completeness preflight on an actual structural candidate; only after this gate can the same Wilson/Kubo *static* benchmark be reproduced reliably with Wannier data.

Source and reproduction: bsf_nmee/mixed_wilson_geometry.py, bsf_nmee/symmetry_dft_gate.py, scripts/run_stage7i_wilson.py, tests/test_stage7i_wilson_static.py, tests/test_stage7i_dft_symmetry.py. The independent offline Stage 7i archive additionally contains source-step sequences, all 2x2 tensor components, standalone 8-band reference, U(4) tests and a machine-readable DFT input gate.
