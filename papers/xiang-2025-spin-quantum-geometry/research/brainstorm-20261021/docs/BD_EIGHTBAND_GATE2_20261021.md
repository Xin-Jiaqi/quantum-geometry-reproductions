# B+D Gate 2: two-orbital PT bilayer with registry-generated spin-dependent interlayer hopping

**2026-10-21. Independent exploratory theory, NOT Xiang/Jin/Wang PRL 135, 146303 (2025) reproduction.**  
**Verdict:** a **conditional model-level positive control**, not novelty verification and not a material prediction.

## Research question

Can a spin-density response in a magnetic bilayer be nonzero only when **both spin-independent and spin-dependent physical interlayer bonds are present**, while being symmetry related across two distinct structural registries, and disappearing when either channel is disconnected?

A true result has to pass: explicit microscopic bond Fourier transform; PT and stacking-partner checks; gauge/embedding sewing and U(2) invariance; a full gapped **finite-frequency** susceptibility rather than equating a geometric kernel with DC transport; and elimination of the uncoupled-layer signal.

## Eight-band model and orbital representations

Hilbert space: layer $\tau$ (top/bottom) × two opposite inversion-parity orbitals $\rho$ (even/odd, conceptually $s/p_z$) × physical spin $\sigma$. Represent the antiunitary magnetic operation, acting at fixed in-plane momentum, by

$$
\mathcal A=PT=\tau_x\rho_z(i\sigma_y)K,\qquad \mathcal A^2=-1.
$$

For a pair of layer positions separated by the in-plane shift $\boldsymbol\delta=(\delta_x,\delta_y)$, the interlayer hopping block is assembled by Fourier transforming genuine displaced bonds $\mathbf d=\mathbf R+\boldsymbol\delta$, with $\mathbf R\in\mathbb Z^2$:

$$
t_e(\mathbf k,\boldsymbol\delta)=t_0\sum_\mathbf R e^{-|\mathbf d|^2/(2\ell^2)}e^{i\mathbf k\cdot\mathbf d},\qquad
t_{o,a}(\mathbf k,\boldsymbol\delta)=\lambda\sum_\mathbf R d_a e^{-|\mathbf d|^2/(2\ell^2)}e^{i\mathbf k\cdot\mathbf d}.
$$

The *assumed*, symmetry-compatible spin-dependent interlayer hopping uses $t_{o,x}\rho_x\sigma_y+t_{o,y}\rho_x\sigma_x$, while the spin-independent hopping uses $t_e(\rho_0+0.18\rho_z)$. If $t=a+ib$, the corresponding Hermitian layer block is $a\tau_x-b\tau_y$. This is **not** a claim that a particular $s/p_z$ material realizes these spinful amplitudes; $\lambda$ is an assumed microscopic bond coefficient that must be obtained from an orbital/Slater–Koster/SOC calculation in a real system.

The intralayer Hamiltonian is independent of $\boldsymbol\delta$ and contains the orbital splitting, staggered exchange, and staggered Rashba:
$$
H_{\rm intra}=\left[-0.12(\cos k_x+\cos k_y)\right]I+
\left[\Delta_\rho+0.17(\cos k_x-\cos k_y)\right]\rho_z
+0.30\tau_z(\sin k_y\sigma_x-\sin k_x\sigma_y)
+m\tau_z\sigma_z,
$$
with $\Delta_\rho=1.5,\ m=0.48$. The main example uses $t_0=0.30,\ \lambda=0.22,\ \ell=0.72$. All values are synthetic dimensionless energy/length units.

**Why eight bands matter:** the strict one-orbital-per-layer four-band representation allows only $I,\tau_x,\tau_y,\tau_z\sigma_{x,y,z}$ as PT-even Hermitian products; in particular direct $\tau_{x,y}\sigma_i$ are forbidden. Enlarging to opposite-parity orbitals permits $\tau_{x,y}\rho_x\sigma_{x,y}$, which are explicitly PT even. This is not a universal no-go theorem outside the chosen representation.

## Symmetry and physical connection

For every tested $\mathbf k$ and registry,
$$
\mathcal A H(\mathbf k,\boldsymbol\delta)\mathcal A^{-1}=H(\mathbf k,\boldsymbol\delta),\quad
H(-\mathbf k,-\boldsymbol\delta)=C_{2z}H(\mathbf k,\boldsymbol\delta)C_{2z}^{-1},
$$
with $C_{2z}=i\sigma_z$. Also $v_x=\partial_{k_x}H$ obeys the corresponding negative-sign velocity sewing. $k_x\to k_x+2\pi$ obeys **orbital-embedding unitary sewing**, not strict same-basis periodicity; numerical checks explicitly include the embedded layer displacement.

**Important limit:** $C_{2z}$ relates two **structural registries**, but has not been shown to reverse a polarization $P_z$, yield an energetically stable double well, or constitute a sliding-ferroelectric switch.

## Observable and Kubo protocol

Use a *physically local* layer-staggered spin-density operator rather than naming a conventional spin flow a measurable current:
$$
S_-^z=\tfrac12\tau_z\sigma_z,\qquad \delta S_-^z(\omega)=\chi_{S_-^z,E_x}(\omega)E_x(\omega).
$$
The response is evaluated at half filling, **inside a true band gap**, using all occupied–unoccupied transitions of the eight-band independent-electron, length-gauge retarded formula. For valence $v$ and conduction $c$, define
$$
\Delta_{cv}=\varepsilon_c-\varepsilon_v,\quad
r^x_{cv}=-i(v_x)_{cv}/\Delta_{cv},\quad
Q_{vc}=(S_-^z)_{vc}r^x_{cv}.
$$
Then, up to an overall chosen sign/elementary-charge convention,
$$
\chi(\omega)=\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}
\sum_{v,c}\left[
\frac{Q_{vc}}{\omega-\Delta_{cv}+i\eta}-
\frac{Q_{vc}^\ast}{\omega+\Delta_{cv}+i\eta}
\right].
$$
For exact Kramers doublets, the same result is independently checked using rank-two spectral projectors $P_a$ and traces $\mathrm{Tr}(P_a S_-^z P_b v_x)$. The spin-density observable does not use an ambiguous conventional spin-current definition; **its induced local exchange/SOC torques still require a complete dynamical continuity analysis**.

**Crucial distinction:** this is a finite-$\omega$ *interband insulator* susceptibility, not an intrinsic $\tau^0$ metallic DC transport coefficient or the full original Xiang et al. nonlinear spin-QGT response. A nonzero real interband kernel by itself is not a positive quantum metric or a measured conductivity.

## Numerical acceptance and nontrivial negative controls

At $\delta=(0.17,0.23),\ \omega=0.30,\ \eta=0.02$, using the gap parameters above:

| 2D k mesh | $\operatorname{Re}\chi_{S_-^z,E_x}$ | $\operatorname{Im}\chi_{S_-^z,E_x}$ |
| --- | ---: | ---: |
| $11\times11$ | $-2.560967\times10^{-8}$ | $3.260802\times10^{-7}$ |
| $17\times17$ | $-2.583495\times10^{-8}$ | $3.282797\times10^{-7}$ |
| $25\times25$ | $-2.583936\times10^{-8}$ | $3.283163\times10^{-7}$ |

The $25^2$ calculation has an indirect gap $\simeq0.65746$ and a minimum direct gap $\simeq0.83726$: $\omega=0.30$ lies below the interband threshold. At $\eta\to0$ the small real contribution tends toward zero while the phase-quadrature imaginary part remains finite: for $\eta=0.04,0.02,0.01,0.005$, $\operatorname{Im}\chi$ approaches $3.288\times10^{-7}$. No subgap absorption claim follows from this calculation.

At each mesh, *three negative controls* vanish within numerical precision:

- $t_0=0,\ \lambda=0$: all physical interlayer hopping off.
- $t_0=0.30,\ \lambda=0$: only spin-independent interlayer hopping.
- $t_0=0,\ \lambda=0.22$: only spin-dependent interlayer hopping.

The response survives **only when both independent physical bond channels coexist** in this particular model. Comparing only the spin-independent hopping strength while leaving a disconnected-system magnetic term was the *failure* of Gate 1; that artifact is absent here.

Under the actual partner map $\boldsymbol\delta\mapsto-\boldsymbol\delta$, 
$$
\chi_{S_-^z,E_x}(-\boldsymbol\delta)=-\chi_{S_-^z,E_x}(+\boldsymbol\delta).
$$
At the symmetric registry $\delta=(0,0)$ the same component vanishes. **These are C2z tensor consequences**, not newly discovered universal selection rules.

Checked numerically: Hermiticity; $(PT)^2=-1$; exact PT and C2z Hamiltonian/velocity relations; momentum-boundary orbital embedding; Kramers degeneracy; independent $U(2)$ rotations in each doublet; eigenbasis Kubo versus rank-two projector traces (maximum local difference $\approx1.2\times10^{-16}$ at $11^2$); actual indirect gap; $\pm\delta$ parity; the full interlayer/local-torque *operator commutator identity*; and real-space-bond registry reindexing at $R_{\max}=6$. Code truncates at $R_{\max}=4$ near the small-displacement baseline; the larger truncation was explicitly used to check integer-lattice registry invariance.

**This is a test of the model equations and numerical implementation, NOT a proof that the spin density is experimentally detectable.** Magnitudes lack physical SI calibration and material-specific lifetimes.

## Known peer work and immediate novelty blockers

Our earlier [17-item peer-work audit](PEER_WORK_AUDIT_20261021.md) already identifies at least four direct competition classes:

1. [Zheng et al., *Nano Lett.* 2024](https://doi.org/10.1021/acs.nanolett.4c01657): extended-space geometry and stacking-dependent electric multipoles.
2. [Sliding FE / layer spin Hall, *J. Phys. Chem. Lett.* 2024](https://doi.org/10.1021/acs.jpclett.4c01479): sliding and layer-spin Hall locking.
3. [Sun et al., *Advanced Science* 2025](https://doi.org/10.1002/advs.202505860): magnetic interlayer coupling controls Berry- and quantum-metric-driven Hall mechanisms in CoPSe3 bilayers.
4. [Tang et al., *Nature Communications* 2025](https://doi.org/10.1038/s41467-025-63171-1): insulating AFM staggered current-induced spin accumulation and Néel torque.
5. [Zhu et al., *PRL* 2025](https://doi.org/10.1103/dmzg-ck2t) and [Du et al., *PRB* 2026](https://doi.org/10.1103/g471-xrkl): stacking connection and magnetic spin-layer-group theories.
6. [Spin geometric tunneling in CrPS4, *Nature Physics* 2024](https://doi.org/10.1038/s41567-024-02675-x): spin geometric phases in interlayer transport (different experiment but a warning against claiming spinful interlayer tunneling is new).

These sources have not been exhaustively matched to the exact new *eight-band* spin-density Kubo equation and orbital-bond representation; **originality is OPEN**, not established by the presence of a nonzero coefficient.

## Next falsification gates (do not leap directly to DFT)

- **G3, independent physical microscopic mechanism:** derive the odd-parity spinful two-center / SOC tunneling from specified atomic orbitals, local magnetic order, hopping paths and Wannier symmetry; do not merely fit $\lambda$ to obtain nonzero response.
- **G4, symmetry beyond C2z:** ask if genuine BSF two-state switching (inverting *measured* polarization) permits the same pair, and identify its magnetic connection operator, not just an arbitrary $\delta\to-\delta$ pair.
- **G5, full response and observability:** compare the local spin density with layer-resolved exchange/SOC torque and finite-frequency dynamics, contact terms, and spin relaxation; evaluate an experimentally addressable probe.
- **G6, peer equation-level comparison:** especially insulating Néel torque, layer Edelstein, magnetic-coupling quantum metric, and proper spin-current literature; identify a genuinely different measurable prediction.

**Decision:** keep B+D alive as a mathematically consistent, physically explicit *toy hypothesis*. No paper-style discovery language, no claim of ferroelectricity, and no changes to the original PRL reproduction tests or scientific acceptance.


## Supplementary phase-gauge and amplitude counterfactual (2026-10-21)

On the same $17^2$ mesh, keeping $\lambda=0.22$ and varying the ordinary $t_0=0,0.10,0.15,0.20,0.25,0.30,0.35$ gives $\operatorname{Im}\chi\approx0,\ 9.2299\times10^{-9},\ 2.55084\times10^{-8},\ 6.03346\times10^{-8},\ 1.373751\times10^{-7},\ 3.282797\times10^{-7},\ 9.645773\times10^{-7}$. Thus the signal does not require a registry-dependent **intralayer** odd-dispersion term and goes continuously to zero when ordinary interlayer hopping is removed. This is **only model interference**, not an identified quantum-geometric mechanism or originality proof.

For fixed $t_0=0.3$, the response is an **even function of $\lambda$** within numerical noise. This is expected: the global opposite-orbital phase redefinition $\rho_z$ sends $\rho_x\to-\rho_x$ and hence $\lambda\to-\lambda$, while leaving the other Hamiltonian terms and observable invariant. Similarly, the sign of $t_0$ alone is unobservable in this model by an allowed combination of orbital/layer phase conventions. **Do not interpret the sign or evenness of an individual hopping coefficient as a stacking-switchable physical observable.** The physically distinct $\pm\boldsymbol\delta$ registries have instead been checked with the actual $C_{2z}$ connection operator.
