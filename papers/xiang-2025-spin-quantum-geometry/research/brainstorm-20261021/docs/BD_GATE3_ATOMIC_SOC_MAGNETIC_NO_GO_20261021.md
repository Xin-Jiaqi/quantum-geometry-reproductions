# B+D Gate 3: microscopic atomic-SOC pathway, magnetic no-go and polar-connection hypothesis

**Date:** 2026-10-21.  
**Scope:** independent physics brainstorming in the Xiang/Jin/Wang 2025 Spin QG vicinity, **NOT the original paper reproduction**. No real-material prediction, claimed novel effect, Wannier fit, DFT, or verified sliding ferroelectricity.  
**Software:** [standalone reproducible NumPy script](../scripts/bd_gate3_atomic_soc.py). The frozen PRL reproduction files, tests, CI and \`PARTIAL\` status are untouched.

## 1. Problem formulation: three independent gates

**Gate M (microscopic origin):** Can a truly spin-independent Slater–Koster two-center interlayer bond, together with local atomic \(\lambda\mathbf L\cdot\mathbf S\), generate the PT-compatible spinful interlayer orbital hopping postulated in the prior eight-band Gate 2?

**Gate N (magnetic indispensability):** Does the finite-frequency layer-opposite spin response actually require AFM exchange, or does it persist when the magnetic order \(m\to0\)? Which genuinely magnetic-odd response survives?

**Gate C (structure and switching):** Can a pair of *polar-perturbed* daughter states be related by a physically motivated **\(m_z T\)** connection, in contrast to our previous nonpolar, \(C_{2z}\)-related, individually-PT states? Does this establish FE? **No**, unless polar asymmetry is derived from the microscopic registry and a stable double well is verified.

## 2. The full 16-band atomic model (not guessed spinful tunneling)

Layers: \(T/B\), orbitals \(s,p_x,p_y,p_z\) on each layer, physical spin \(\uparrow/\downarrow\). The total Hilbert-space dimension is \(2\times4\times2=16\).

- Onsite orbital energies \((-1.15,4,4,+1.15)\). The high-energy \(p_x,p_y\) orbitals are explicitly retained before downfolding.
- Collinear layer-staggered exchange \(m\tau_z\sigma_z\), \(m=0.36\) by default.
- **Only onsite spin-orbit coupling**, \(\lambda_{\rm SOC}(L_x\sigma_x+L_y\sigma_y+L_z\sigma_z)/2\), \(\lambda_{\rm SOC}=0.36\).
- The **entire interlayer hopping is spin independent**. For a bond \(\mathbf d=(R_x+\delta_x,R_y+\delta_y,h)\), \(h=0.95\), \(\hat{\mathbf d}=\mathbf d/|\mathbf d|\), Slater–Koster block in \((s,\mathbf p)\):

$$
T_{\rm bond}(\mathbf d)=
e^{-|\mathbf d_\parallel|^2/(2\ell^2)}
\begin{pmatrix}
V_{ss} & V_{sp}\hat{\mathbf d}^{\,T}\\
-V_{sp}\hat{\mathbf d}&V_{pp\pi}I_3+
(V_{pp\sigma}-V_{pp\pi})\hat{\mathbf d}\hat{\mathbf d}^{\,T}
\end{pmatrix}\otimes I_\sigma,
$$

with \(V_{ss}=0.28\), \(V_{sp}=0.23\), \(V_{pp\pi}=0.04\), \(V_{pp\sigma}-V_{pp\pi}=0.17\), \(\ell=0.68\). The Bloch bond phase is \(e^{i\mathbf k\cdot\mathbf d_\parallel}\); \(v_x=\partial_{k_x}H\) includes its complete derivative. Synthetic orbital dispersions are \(\mp0.12(\cos k_x+\cos k_y)\) for \(s,p_z\).

This is a transparent **model** of virtual atomic SOC paths, not a universal accurate Slater–Koster parameterization of a specific real material. Orbital-dependent vertical radial hopping and hopping-distance relaxation have not been calibrated to first principles.

The natural antiunitary magnetic operator, fixed at each in-plane \(\mathbf k\) in this orbital-embedded convention, is

$$
PT=\tau_x\,{\rm diag}(1,-1,-1,-1)\,(i\sigma_y)K,\qquad (PT)^2=-1.
$$

The \(C_{2z}\) partner operator uses \({\rm diag}(1,-1,-1,+1)\) in orbital space and \(i\sigma_z\) in spin space.

## 3. Gate M: physical virtual path actually generates spinful hopping (PASS within model)

Eliminate high \((p_x,p_y)\) orbitals with a **fixed-energy \(E=0\) Schur complement**:

$$
H_{\rm eff}(E=0)=H_{LL}-H_{LH}H_{HH}^{-1}H_{HL}.
$$

A representative virtual path is

$$
s_T \xrightarrow{\ t_{sp_x}(\mathbf d)\ }p_{x,B}
\xrightarrow{\ \lambda_{\rm SOC} L_y\sigma_y/2\ }p_{z,B},
$$

and its \(p_y\) counterpart. Therefore spin dependence is generated from *spin-independent interlayer hopping plus on-atom \(p\)-SOC* rather than included as a bare hopping assumption.

For \(\mathbf k=(0.5,0.4)\), \(\boldsymbol\delta=(0.17,0.23)\) and the parameters above, the resulting \(s_T\to p_{z,B}\) spin-matrix \(\mathsf T_{s,p_z}\) has Pauli decomposition \(\sum_{\mu=0,x,y,z}t_\mu\sigma_\mu\):

| Component | Complex amplitude (synthetic units) |
| --- | --- |
| \(t_0\) | \(+0.4693082+0.0041220i\) |
| \(t_x\) | \(-0.0035485-0.0001551i\) |
| \(t_y\) | \(+0.0043171-0.0005678i\) |
| \(t_z\) | \(-1.5\times10^{-7}-0.8\times10^{-7}i\) |

**Both necessary counterfactuals:** with \(\lambda_{\rm SOC}=0\), all spin-dependent components disappear; with \(V_{sp}=0\), this \(s\to p_z\) interlayer channel disappears. Full 16-band Hamiltonian and its \(E=0\) Schur complement satisfy their respective \(PT\) checks numerically.

**Limit:** an \(E=0\) Schur complement is NOT an exact frequency-dependent, normalized low-energy theory. To use \(H_{\rm eff}\) for response coefficients one must also downfold the physical spin operator and electric-field vertices, accounting for eliminated-state weight and basis/position-operator connections. Here **all quoted measurable-response proxies below were instead computed in the original full 16-band space**, so this omission cannot silently contaminate those numbers.

## 4. Gate N: the initial AC staggered-spin response fails the magnetic gate

Direct full-band length-gauge retarded Kubo, with **four occupied bands**, uses the physical Hermitian spin-density operators \(S_-^z=\tau_z\sigma_z/2\), \(S_{\rm tot}^z=\sigma_z/2\). All occupied-to-unoccupied transitions of 16 bands are included. For occupied \(n\), unoccupied \(m\),

$$
\chi_{O,E_x}(\omega)=\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}
\sum_{n\in{\rm occ},m\in{\rm emp}}
\left[\frac{-i\,O_{nm}v^x_{mn}/\Delta_{mn}}
{\omega-\Delta_{mn}+i\eta}-
\frac{(+i\,O_{mn}v^x_{nm}/\Delta_{mn})}
{\omega+\Delta_{mn}+i\eta}\right],
$$

up to electron-charge / overall sign convention. Parameters: \(\omega=0.30,\eta=0.02,\delta=(0.17,0.23)\); \(13\times13\) BZ mesh. The indirect gaps remain open (e.g. \(1.08\) at \(m=\pm0.36\), \(1.57\) at \(m=0\)); this is **subgap interband response**, not a metallic RTA current.

| \(m\) | \(\operatorname{Im}\chi_{S_-^z,E_x}\) | \(\operatorname{Re}\chi_{S_{\rm tot}^z,E_x}\) |
| ---: | ---: | ---: |
| \(0\) | \(+1.31922\times10^{-7}\) | \(\sim0\) |
| \(+0.36\) | \(+2.39792\times10^{-8}\) | \(+1.21205\times10^{-5}\) |
| \(-0.36\) | \(+2.39792\times10^{-8}\) | \(-1.21205\times10^{-5}\) |

This destroys a tempting claim: the original **layer-opposite AC signal survives without magnetism**, and is even under Néel reversal. It cannot be called a magnetic-specific mechanism, regardless of how elegantly it is built from microscopic SOC.

The total-spin signal *is magnetic odd*: \(\chi_{S_{\rm tot}^z,E_x}(m)=-\chi_{S_{\rm tot}^z,E_x}(-m)\) and it vanishes at \(m=0\); it also vanishes when interlayer hopping, \(s\!-\!p\) bonding, or local SOC is removed.

### But the magnetic-odd static result is conventional mixed Berry geometry, not novelty

In the **same complete 16-band model** compute the static \(\chi_{S_{\rm tot}^z,E_x}(0)\) both from the exact occupied–empty Kubo sum and an **independent finite-difference occupied-subspace projector**, adding a uniform auxiliary Zeeman \(b_z S_{\rm tot}^z\). Define

$$
F_{k_x b_z}=i\,\mathrm{Tr}\,P\,[\partial_{k_x}P,\partial_{b_z}P].
$$

At \(m=+0.36\), \(N_k=11\times11\):

$$
\chi_{S_{\rm tot}^z,E_x}^{\rm static}=+1.1439529597\times10^{-5},
\qquad
\int_{\rm BZ}F_{k_x b_z}=-1.1439530857\times10^{-5}.
$$

The sum is only \(-1.26\times10^{-12}\) (finite-difference \(h=10^{-3}\)). The \(m<0\) signs invert together. Hence this attractive-looking magnetic-odd insulator signal is **ordinary mixed \((\mathbf k,b)\) Berry curvature**, i.e. the familiar intrinsic **linear spin magnetoelectric polarizability** in a gapped magnetic system. Its exact spectral/geometric representation, sign switch, and dependence on SOC do **not** establish a new effect. The standalone script implements both numerical routes and asserts their agreement.

Established theoretical context includes [Mixed Weyl semimetals and low-dissipation magnetization control in insulators by spin-orbit torques, *Nature Communications* (2017)](https://doi.org/10.1038/s41467-017-01138-7) and the Xiang spin-QG and Jia et al. NMEE literature; their distinct response orders and operators must not be conflated.

## 5. Gate C: PT forbids ferroelectric polarization; use a connection between PT-broken daughters

Because \(T\) leaves electric polarization invariant but inversion sends it to its negative, **an exact \(PT\) symmetry of a state prohibits any ordinary net polarization**. In a 2D finite-thickness sample the out-of-plane polarization \(P_z\) is particularly unambiguous (unlike possible in-plane polarization modulo a quantum). Our earlier exactly-PT bilayer states connected by \(C_{2z}\) are *not* demonstrated sliding-ferroelectric daughters. Moreover, \(C_{2z}\) **does not reverse the polar \(z\) axis**, so its registry pairing cannot alone be invoked as FE-\(P_z\) switching.

A better **symmetry template**, not a self-consistent material, uses the antiunitary structural connection

$$
\mathcal C=m_zT
=\tau_x\,{\rm diag}(1,1,1,-1)\,\sigma_x K
$$

(up to a global phase). It exchanges layers and reverses the physical \(P_z\) orientation while keeping the chosen AFM layer Néel vector unchanged under the combined mirror/time reversal. Explicit numerical checks show

$$
\mathcal C H(\mathbf k,\boldsymbol\delta,+u;m)\mathcal C^{-1}
=H(-\mathbf k,-\boldsymbol\delta,-u;m)
$$

with a **manually imposed** polar layer potential \(u\tau_z\). Either \(u\ne0\) daughter breaks \(PT\) and the old \(C_{2z}\) mapping, as it must. This is **not** an output of the Slater–Koster model's relative shift: the polar potential is an external hypothesis. The approximate occupied-layer charge imbalance is a *proxy* for a dipole, not a full ionic/electronic Berry-phase polarization calculation, and no double well has been found.

Static \(13^2\) mesh, \(m=.36,\ |u|=.08,\ \delta=\pm(.17,.23)\):

| observable response to \(E_x\) | \((+\delta,+u)\) | \((-\delta,-u)\) |
| --- | ---: | ---: |
| \(\chi_{S_{\rm tot}^z,E_x}(0)\) | \(+1.25692149\times10^{-5}\) | \(-1.25692149\times10^{-5}\) |
| \(\chi_{S_-^z,E_x}(0)\) | \(+0.00124707674\) | \(+0.00124707674\) |
| electron layer occupancy \(\langle\tau_z\rangle\) | \(-0.0170841\) | \(+0.0170841\) |

The distinct **odd/even response parity** across the connected states agrees with the antiunitary symmetry of the observable and the electric perturbation. This selection rule is standard magnetic representation theory. A response reaching \(10^{-3}\) in synthetic units is **not** evidence of a material-scale signal or enhanced novel effect. The static staggered response is odd under magnetic-order reversal, and is close in spirit to a Néel torque/magnetoelectric response, which already has strong peer precedent.

## 6. Tight near-neighbor literature, updated this round

1. **Xiang, Jin & Wang**, [*Spin Transport Revealed by Spin Quantum Geometry*, PRL 135, 146303 (2025)](https://doi.org/10.1103/14mp-263q): core spin QG reference; none of these toy results reproduce its metallic or nonlinear results.
2. **Zhou et al.**, [*Layer Edelstein Effect*, PRL 137, 076401 (2026)](https://doi.org/10.1103/qhvs-pfkq): nonmagnetic bilayer layer-odd current-induced spin magnetization, symmetry and bilayer stacking. **Direct warning against calling the \(m=0\) AC layer-spin phenomenon inherently new.** This paper emphasizes current-driven and electrically gated configurations; our insulating finite-frequency observable is not automatically equation-by-equation identical.
3. **Pan, Jiang & Fabian**, [*Multistate manipulation of charge-spin conversion in two-dimensional ferroelectric bilayers*, PRB 114, 115135 (2026)](https://doi.org/10.1103/np2g-86h3): explicit FE bilayer EE/SHE multistate switching; direct FE-spin competition.
4. **Peng et al.**, [*Sliding Ferroelectricity Driven Spin-Layertronics in Altermagnetic Multilayers*, Advanced Science, e76050 (2026)](https://doi.org/10.1002/advs.76050): DFT CuF2 FE-paired bilayers, spin/layer control, and four-state multilayers. **Published full-text close neighbor**, not just a preprint.
5. **Zhu et al.**, [*Sliding Ferroelectric Control of Unconventional Magnetism in Stacked Bilayers*, PRL 135, 056801 (2025)](https://doi.org/10.1103/dmzg-ck2t): connecting-operator framework and magnetic stacking relationships.
6. **Tang et al.**, [*Néel spin-orbit torque in antiferromagnetic quantum spin and anomalous Hall insulators*, Nat. Commun. 16, 7790 (2025)](https://doi.org/10.1038/s41467-025-63171-1): insulating AFM staggered E-driven torque and AC response, direct barrier to rebranding it.
7. **Hanke et al.**, [*Mixed Weyl semimetals and low-dissipation magnetization control in insulators by spin-orbit torques*, Nat. Commun. (2017)](https://doi.org/10.1038/s41467-017-01138-7): mixed momentum/magnetization Berry geometry and insulator torque have long prior art.
8. **2026 antiferromagnetic moiré piezospin response**, [*Quantized piezospintronic effect in antiferromagnetic Moiré systems*, npj Spintronics (2026)](https://www.nature.com/articles/s44306-026-00135-1): mixed \((k,\mathrm{strain})\) Berry curvature and spin-dependent pumping; direct warning for any future "mixed stacking-spin geometry" pump claim.

These papers establish *close conceptual overlap*; a final novelty decision would require matching the actual Kubo observable, symmetry, local operator, and frequency protocol at equation level, not just citing abstracts.

## 7. Provisional verdict and next high-value experiments

- **Gate M: PASS as a toy mathematical mechanism.** Spin-dependent interlayer low-energy tunneling is generated from onsite atomic SOC and spinless Slater–Koster bonds; no artificial eight-band spinful hopping needed.
- **Gate N1: NO-GO for "magnetism-required staggered AC spin response."** Signal survives in the nonmagnetic parent.
- **Gate N2: OVERLAP with ordinary insulator spin magnetoelectricity.** Magnetic-odd static total-spin response is exactly the standard mixed Berry-curvature coefficient numerically.
- **Gate C: SYMMETRY TEMPLATE ONLY.** \(\mathcal C=m_zT\) pairs manually polar-perturbed daughters with distinct spin-response parities. **No physical sliding-ferroelectric pair yet.**

**The next genuinely scientific bottleneck is not another response formula.** Construct two **microscopically inequivalent, noncentrosymmetric monolayer sublattices** with explicit atomic registry, electronic/ionic charges and a total-energy surface \(E(\boldsymbol\delta)\); establish paired minima with opposite calculated \(P_z\) and a finite switching barrier. Only then run the full band spin and torque response, test the disconnected-layer/SOC-off controls, and compare equation-by-equation with 2025 NSOT and 2026 bilayer EE/SHE papers. If symmetry/continuity says the remaining signal is just a known response, archive and redirect without material DFT.

**No modifications to original-paper reproduction code, frozen figures or scientific acceptance.**


### Gate S: a structural reason why tuning this one-site model will never create genuine sliding FE

There is a simple *ionic-geometry* no-go, separate from the electronic \(PT\) constraint. Take identical one-site centrosymmetric Bravais layers at
\(\mathbf r_{T,\mathbf R}=(\mathbf R-\boldsymbol\delta/2,+h/2)\) and
\(\mathbf r_{B,\mathbf R}=(\mathbf R+\boldsymbol\delta/2,-h/2)\).
Ordinary spatial inversion about the midpoint sends
\(\mathbf r_{T,\mathbf R}\mapsto\mathbf r_{B,-\mathbf R}\) **for every \(\boldsymbol\delta\)**, and vice versa. Therefore no amount of changing *only the real-space interlayer hopping magnitudes* between these identical, single-site layers can manufacture a **structurally switchable net \(P_z\)**. In the AFM state, inversion may be broken by staggered spins but its product \(PT\) remains, again forbidding electric polarization. The opposite-parity \(s,p_z\) orbitals on the **same atomic center** do not change this geometric conclusion.

The minimum credible replacement is an explicitly decorated/noncentrosymmetric atomic motif (at least an inequivalent two-site arrangement with correct site symmetry and chemically meaningful positions), possibly with an appropriate relative layer orientation. Only then can different registries possess distinct top/bottom electrostatic environments and an actual calculated electronic+ionic \(P_z\). Even with such a motif, *paired opposite dipoles and distinct minima* must be checked rather than assumed. This is a **model-scoped geometric no-go**, not a universal no-go for single-orbital models with inequivalent chemical layer structures.
