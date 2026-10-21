# B+D Gate 4 | Polar registry and symmetry-allowed magnetoelectric readout

**2026-10-21. Independent original-project brainstorming.** Not Xiang/Jin/Wang PRL 135 146303 (2025) reproduction; no changes to its benchmark code, results or scientific acceptance. No DFT, physical material prediction, proven ferroelectric switching, or priority assertion.

**Main result:** Replacing identical one-site layers with inequivalent A/B honeycomb motifs removes Gate 3's structural no-go. One explicit *parameterized* stacking-energy model supports opposite-electronic-dipole AB/BA local minima, with a connection symmetry; a separate **32-band atomic-SOC model** exhibits FE-even total versus FE-odd layer-staggered **in-plane** static spin magnetoelectric response. However the originally proposed \(E_x\to S_z\) linear response vanishes by \(C_{3z}\) **at the polar minima**. The spinful model's nonpolar AA registry is **lower in energy**, so the AB/BA extrema are only *metastable* there. This is a **screening result**, not an original physics discovery.

## Three-part research logic

| Part | Candidate question | Actual result / scope |
| --- | --- | --- |
| Idea | Can a magnetic bilayer have genuine \(\pm P_z\) stable stackings connected by \(m_zT\), with distinct total vs layer-opposite spin responses? | Symmetry-compatible physical atomic registry construction; **toy local extrema only** |
| Test | Derive electronic dipole and energy from A/B bonds, test symmetry, full orbital SOC and insulating spin Kubo, and null controls | Numerical PASS in specified toy parameters; axial \(E\to S_z\) **NO-GO** by \(C_3\); in-plane spin response survives |
| Peer work | Is this a new FE-spin theory, or just symmetry + standard magnetoelectricity? | **Strong direct overlap** with 2025 PRL stacking connection, 2026 *Nature Communications* unified FE-spin symmetry and 2025/26 magneto-optic / charge–spin bilayer works. Novelty **NOT established** |

## A. A physically explicit two-site honeycomb stacking geometry

Let primitive vectors be \(\mathbf a_1=(1,0),\mathbf a_2=(1/2,\sqrt3/2)\) in arbitrary units. Each layer has **two different chemical species**, A and B, at \(\mathbf r_A=0,\mathbf r_B=(\mathbf a_1+\mathbf a_2)/3\), separated from the other identical layer by vertical distance \(h\). The top and bottom lateral origins are shifted by \(-\boldsymbol\delta/2\) and \(+\boldsymbol\delta/2\). Two candidate honeycomb registries:
\[
\delta_{AB}=+\tfrac13(\mathbf a_1+\mathbf a_2),\qquad
\delta_{BA}=-\tfrac13(\mathbf a_1+\mathbf a_2).
\]
A horizontal reflection **exchanges the two layers** and maps AB to BA; for A-type AFM moments \(+\hat z,-\hat z\), the combined antiunitary \(m_zT\) maps the magnetic stackings while preserving the chosen Néel orientation. **Neither daughter must have PT symmetry.**

The first check is a **spinless four-orbital bilayer** (one local orbital per A/B per layer) with in-plane A–B nearest-neighbor bonds and *genuine* \(\delta\)-dependent interlayer hopping
\[
t_{ab}(\mathbf k,\delta)=t_\perp\sum_{\mathbf R}
e^{-|\mathbf R+\delta+\mathbf r_b-\mathbf r_a|^2/(2\ell^2)}
e^{i\mathbf k\cdot(\mathbf R+\delta+\mathbf r_b-\mathbf r_a)}.
\]
The sublattice chemical contrast is \(\varepsilon_A=-\varepsilon_B=1.4\). Band energy is the sum of two occupied spinless bands. A separate empirical pair contribution is explicitly added:
\[
E_{\mathrm{pair}}(\delta)=\sum_{\mathbf R,ab}
\left[q_aq_b \frac{e^{-\kappa r_{ab}}}{r_{ab}}+
A_{\rm core}\exp\left(-\frac{(r_{ab}-h)^2}{2w^2}\right)\right],
\quad
r_{ab}=\sqrt{|\mathbf R+\delta+\mathbf r_b-\mathbf r_a|^2+h^2}.
\]
Here \(q_A=-q_B=0.9,\;\kappa=2,\; A_{\rm core}=0.12,\;w=0.25,\;h=0.75\). The full *toy* energy is \(E_{\mathrm{band}}+E_{\mathrm{pair}}\). **It is not a self-consistent Hartree/DFT functional**, may double-count electrostatic terms, and uses uncalibrated parameters. It is an existence test of a simple symmetry-respecting energy landscape.

### Spinless structural/electronic results

19×19 k mesh: at AB/BA the *orbital-center layer-electron-number difference* \(\Delta n=n_T-n_B\) is \(-0.01510082323\) and \(+0.01510082323\) electrons/unit cell; AA has zero. The out-of-plane electronic dipole proxy is
\[
p_z^{\mathrm{el}}=-e(h/2)\Delta n
\]
provided the basis orbitals' centers are fixed on their layers and their ionic \(z\)-dipoles cancel. **This is not a calibrated full material polarization; no SI number is claimed.** The sign reverses under the physical layer exchange.

- AB and BA energy \(E=-2.54521862748\), identical to numerical precision; AA \(E=-2.50921435111\) (synthetic units).
- Hessian eigenvalues at either AB or BA: \(0.17846,0.53538\), so both are strict local minima **for these chosen parameters**.
- At a candidate stacking saddle, fractional registry \((u,v)=(0,1/2)\), Hessian eigenvalues \(-0.17225,+0.70057\). The saddle energy is approximately \(0.00495643\) above AB. This is a **candidate path barrier/upper bound**, not a converged minimum-energy switching barrier.
- \(\Delta n\) becomes zero in the disconnected-layer limit \(t_\perp=0\), in the identical-species limit \(\varepsilon_A=\varepsilon_B\), and at the symmetric AA registry.
- **Parameter sensitivity matters:** reducing the empirical ionic charge from \(q=0.9\) to \(q=0.45\) or 0 removes AB/BA's positive-definite Hessian. We have not predicted universal structural stability.

These geometric motifs capture symmetry already familiar from hBN; **the toy does not claim magnetic hBN**. For comparison, hBN AB/BA dipoles and stacking landscapes have been experimentally/theoretically established in [Yasuda et al., *Science* 372, 1458 (2021)](https://doi.org/10.1126/science.abd3230), and analyzed in [Bennett et al., *Nat. Commun.* 14, 1629 (2023)](https://doi.org/10.1038/s41467-023-37337-8). The toy's numerical energies and charge transfers are **not** those papers' DFT values.

## B. Atom-local SOC + magnetic response on the same geometry (32-band)

A separate calculation resolves A/B species, \(s,p_x,p_y,p_z\) orbitals, physical spin and T/B layer, giving \(2\times2\times4\times2=32\) states. **All interlayer and intralayer hopping blocks are spin-independent real-space bond-direction Slater–Koster matrices**. Spin mixing arises from onsite
\[
H_{\rm SOC}=\lambda_{\rm at}\,\mathbf L\cdot\mathbf S,\quad\lambda_{\rm at}=0.30,
\]
and an assumed fixed A-type-AFM exchange \(H_m=m\tau_z\sigma_z\), \(m=0.34\) unless specified. Only eight lowest bands are filled. The model has a true insulating direct *and indirect* gap of \(0.12356\) at the polar stacking for the sampled grids. Here magnetic ordering is a hypothetical ingredient, **not a prediction that a B/N honeycomb magnet exists**.

For physical orbital mirror parity \(D_{m_z}=\mathrm{diag}(+,+,+,-)\), define
\[
\mathcal C=m_zT=\tau_x D_{m_z}\sigma_x K.
\]
We directly checked at generic momentum and displacement
\[
\mathcal C H(\mathbf k,+\delta)\mathcal C^{-1}
=H(-\mathbf k,-\delta),\qquad
\mathcal C v_x(\mathbf k,+\delta)\mathcal C^{-1}
=-v_x(-\mathbf k,-\delta).
\]
The mirror exchanges layers; time reversal flips physical spin; together they leave the layer-staggered exchange order \(m\tau_z\sigma_z\) fixed. \(\mathcal C\) is a **connection between two structural states**, not a symmetry of either polar daughter.

**Crucial structural limitation:** With this *32-band* band energy plus the same independent ionic pair term, AB/BA are still **local** minima (coarse-grid Hessian eigenvalues \(\sim0.347,1.041\)), but AA is *lower* in energy: \(E_{AB}\approx-11.71683\) versus \(E_{AA}\approx-11.77679\). Therefore these two polar states are **metastable, not the toy global ground state**. The spinless model's \(0.00496\) saddle difference **must not** be recycled as the 32-band switching barrier.

### Correct symmetry filter before selecting spin observable

At exactly AB and BA, threefold rotation \(C_{3z}\) survives. A \(z\)-component of a spin density (axial \(z\)) is invariant under \(C_3\), whereas \(E_x,E_y\) transform as a 2D vector. Thus **the entire linear response row must vanish**:
\[
\chi_{S^z,E_x}=\chi_{S^z,E_y}=0.
\]
The same conclusion holds for a layer-staggered \(S_-^z\). This is a conventional Neumann/principle-of-symmetry result, not a new theorem. It applies to both static and finite-frequency **linear** spin-density response. Previous Gate 3's nonzero linear \(E_x\to S_z\) at a *generic* registry cannot simply be claimed in an AB/BA ground state.

For the *allowed* in-plane components, compute **static, insulating, occupied-to-empty length-gauge Kubo** using the full 32-band eigenstates and physical velocity \(v_x=\partial_{k_x}H\). No finite-\(\tau\) Drude terms or metallic driven-DC physics enter. Define
\(S_{\rm tot}^{x}=\sigma_x/2\), \(S_-^{x}=\tau_z\sigma_x/2\).
With \(\mathcal C=m_zT\), \(S_{\rm tot}^{x}\) is **connection-even** and \(S_-^{x}\) **connection-odd**, as verified by direct 32-band calculations.

| mesh | \(\chi_{S_{\rm tot}^x,E_x}^{AB}\) | \(\chi_{S_-^x,E_x}^{AB}\) |
| ---: | ---: | ---: |
| \(7^2\) | \(+7.09085484\times10^{-4}\) | \(-8.324333\times10^{-6}\) |
| \(11^2\) | \(+6.54500549\times10^{-4}\) | \(-4.211775\times10^{-6}\) |
| \(15^2\) | \(+6.54176002\times10^{-4}\) | \(-3.501043\times10^{-6}\) |
| \(17^2\) | \(+6.54199346\times10^{-4}\) | \(-3.484312\times10^{-6}\) |
| \(19^2\) | \(+6.54204901\times10^{-4}\) | \(-3.482833\times10^{-6}\) |
| \(23^2\) | \(+6.54206049\times10^{-4}\) | \(-3.482975\times10^{-6}\) |

For the BA partner at every mesh: \(\chi_{\rm total,x}^{BA}=\chi_{\rm total,x}^{AB}\) and \(\chi_{\rm staggered,x}^{BA}=-\chi_{\rm staggered,x}^{AB}\) to numerical precision. Both \(S_{\rm tot}^z\) and \(S_-^z\) linear coefficients are zero to ~\(10^{-17}\) at AB/BA. The calculated occupied-layer charge difference is \(\pm0.076205\) at \(N=19\) (different 32-band Hamiltonian from the spinless \(\pm0.015101\) model); it reverses as expected under the connection.

**True negative controls:** \(\chi_{\rm tot,x}\) at an AB registry falls to numerical zero when setting **either** AFM exchange \(m=0\), atomic SOC \(\lambda_{\rm at}=0\), or all interlayer bonds to zero. Under magnetic-order reversal \(m\to-m\), the *static* in-plane spin response reverses as expected for conventional magnetoelectricity. The response is **not** a proof of a new QGT, nor of a large material magnetoelectric coefficient.

**The fact that a FE-even total response coexists with a FE-odd staggered response is a symmetry-derived diagnostic**, potentially useful for layer-selective readout, but it is a normal magnetic response-tensor consequence. At the model parameters the staggered channel is \(\sim0.5\%\) of the total one. Real-space layer spin accumulation, spin continuity/torques and material sensitivity have not yet been established.

## C. Literature novelty gate: heavy overlap, not yet a new story

Most important direct source overlooked by a superficial review:

- **[A unified symmetry framework for spin–ferroelectric coupling in altermagnetic multiferroics](https://doi.org/10.1038/s41467-026-69635-2)**, *Nature Communications* (2026), full text [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13039994/): classifies **three FE spin coupling regimes** based on FE switching and spin Laue symmetry; has actual DFT six FE minima and multiple FE switching paths in bilayer MnPS\(_3\), plus transport fingerprints. Our elementary \(\mathcal C=m_zT\) parity relation cannot be sold as a comparable new universal classification. Their focus is band/spin reconstruction; an *observable-specific full nonlinear or finite-frequency spin-density response* might be distinct, but needs explicit equation-level demonstration.
- **[Sliding Ferroelectric Control of Unconventional Magnetism in Stacked Bilayers](https://doi.org/10.1103/dmzg-ck2t)**, *PRL* **135**, 056801 (2025): pre-existing connection-operator structural/magnetic rule. Our mathematical symmetry idea must be attributed to this prior art where appropriate.
- **[Unconventional Magnetism, Sliding Ferroelectricity, and Magneto-Optical Kerr Effect in Multiferroic Bilayers](https://doi.org/10.1021/acsami.5c19598)**, *ACS AMI* (2025): FE sliding reverses spin/valley/Kerr signatures. Broad FE spin-optical readout already occupied.
- **[Multistate manipulation of charge-spin conversion in two-dimensional ferroelectric bilayers](https://doi.org/10.1103/np2g-86h3)**, *PRB* (2026): standard Edelstein/Spin Hall spin conversion already tuned by FE.
- **[Photo-magnetization in two-dimensional sliding ferroelectrics](https://doi.org/10.1038/s41699-022-00297-6)**, *npj 2D Materials and Applications* (2022): helicity/illumination-dependent spin+orbital photomagnetization and FE-switching sign control have already been proposed for ZrI\(_2\), WTe\(_2\), MoS\(_2\) bilayers. **Simply pivoting from forbidden linear \(S_z\) to circular-light second-order \(S_z\) is NOT a new mechanism.**
- **[Mixed Weyl semimetals and low-dissipation magnetization control in insulators by spin-orbit torques](https://doi.org/10.1038/s41467-017-01138-7)**, *Nat. Commun.* (2017): momentum/Zeeman Berry curvature and spin magnetoelectric susceptibility standard.
- **[Néel spin-orbit torque in antiferromagnetic quantum spin and anomalous Hall insulators](https://doi.org/10.1038/s41467-025-63171-1)**, *Nat. Commun.* (2025): insulating AFM field-induced staggered torque already exists.

Citations identify **topic/concept overlap**. We have **not** matched the exact 32-band Kubo equation/operator to full supplementals of every competitor; absence of a verbatim matching formula does **not** establish priority.

## D. Gate 4 scientific verdict and next discriminating work

- **PASS (structural):** unlike the one-site model, a two-species atomic motif can have nonzero stacking-induced electric layer dipole with opposite sign, local energy minima and a symmetry connecting them.
- **PASS (symmetry and microscopic method):** atom-local SOC with strictly spin-independent spatial bonds, \(m_zT\) daughter mapping, zero magnetic/SOC/disconnected limits, and k-mesh convergence of allowed response.
- **NO-GO (naive S_z idea):** \(C_3\) at polar local minima forbids *linear* \(E_{x,y}\to S_z\) irrespective of quantum geometry.
- **LIMITATION:** FE local minima are parameter contingent, 32-band polar pair is **metastable** with lower nonpolar AA; neither electronic energy model contains calibrated electron-ion self-consistency. No claimed switching barrier, Curie/Néel temperature, or real-material spin moment.
- **NOT NEW:** layer-dependent static *linear* magnetoelectric coefficient and its \(m_zT\) switching parity are already within conventional response and magnetic-connection frameworks.

**Next Gate 5 (if pursuing B+D):** on one *actual* polar magnetic bilayer candidate where a real FE path and magnetic order have been demonstrated, derive full layer-resolved causal susceptibility \(\chi_{S_\ell^\alpha,E_b}(\omega)\) and explicit exchange/SOC torque response; separate conventional spin-magnetoelectric/equilibrium mixed Berry terms from an **interlayer-irreducible**, experimentally distinguishable spectral/dynamical contribution. Use magnetic connection to establish specific **frequency-resolved** constraints that differ from existing static FE-spin and Kerr classifications. Compare exact formulas with above papers **before** committing to DFT. If the distinctions collapse to symmetry or ordinary Berry curvature, archive this as a negative result.

**Run scripts:** [honeycomb polar structure and energy](../scripts/bd_gate4_polar_honeycomb.py) (NumPy) and [32-band atomic SOC AFM Kubo](../scripts/bd_gate4b_atomic_spin.py). Dense mesh optional: \`python scripts/bd_gate4b_atomic_spin.py --dense\`; may take noticeably longer.

**Strict separation:** these tests are independent paper-inspired research under \`research/brainstorm-20261021/\`, not a reproduced original PRL claim, and do not modify the frozen paper CI.
