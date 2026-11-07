# Stage 9: complete stacking path and the physical meaning of bilayer NMEE

## 1. The Hamiltonian, before the response

Basis: (top, bottom) × (A, B) × (spin up, spin down). Spin Pauli matrices are σ; the A/B Pauli matrices are τ. The layer displacement **d** is Cartesian; the code stores its components in primitive lattice vectors.

$$
H(\mathbf k,\mathbf d;B_z)=
\begin{pmatrix}h_+(\mathbf k)&T(\mathbf k,\mathbf d)\\T^\dagger(\mathbf k,\mathbf d)&h_-(\mathbf k)\end{pmatrix}
-g\mu_B B_z\,\mathbb 1_{\rm layer}\otimes\mathbb 1_{\rm sub}\otimes\sigma_z.
$$

Each layer contains ordinary nearest-neighbour hopping, staggered onsite energy, intrinsic SOC and Rashba SOC:

$$
h_\ell=\begin{pmatrix}m\mathbb 1+d_I(\mathbf k)\sigma_z&F_\ell(\mathbf k)\\F_\ell^\dagger(\mathbf k)&-m\mathbb 1-d_I(\mathbf k)\sigma_z\end{pmatrix},\qquad \ell=\pm1,
$$

$$
d_I=-2\lambda_I\sum_{j=1}^3\sin(\mathbf k\cdot\mathbf b_j),\quad
F_\ell=-t\sum_j e^{i\mathbf k\cdot\boldsymbol\delta_j}\mathbb 1
-i\ell\lambda_R\sum_j e^{i\mathbf k\cdot\boldsymbol\delta_j}
\frac{\sigma_x\delta_{jy}-\sigma_y\delta_{jx}}{a}.
$$

Here δ/a = (0,1), (√3/2,−1/2), (−√3/2,−1/2); b = (δ₂−δ₃, δ₃−δ₁, δ₁−δ₂). Thus h−=σz h+ σz in spin space: the layers are horizontal-mirror partners with opposite Rashba SOC. The onsite term and Rashba terms lower the symmetry of the ideal graphene monolayer. Calling this simply layer-group-80 graphene with identity stacking would be incorrect.

The spin-independent interlayer hopping is

$$
T_{ab}=t_\perp\sum_{\mathbf R}
e^{-|\mathbf R+\mathbf r_b-\mathbf r_a+\mathbf d|^2/\rho^2}
e^{i\mathbf k\cdot(\mathbf R+\mathbf r_b-\mathbf r_a+\mathbf d)}\mathbb 1_{\rm spin},
\quad \mathbf r_A=0,\quad\mathbf r_B=\boldsymbol\delta_1.
$$

Both the hopping magnitude and Bloch phase change with sliding. The electric vertex is ∂H/∂k in this same orbital embedding. An in-plane uniform electric field is treated in linear response, not inserted as a periodic onsite scalar potential. No Coulomb self-consistency, structural force field, ionic relaxation or spin-dependent interlayer hopping is included.

Numerical values: t=0.85 eV, m=10 meV, λI=10 meV, λR=20 meV, a=0.142 nm, ρ=0.16 nm. The effective g=10 follows the existing program convention −gμB Bσ; it does not include a further 1/2. Tunneling amplitudes are model parameters, not fits to a material.

## 2. What is actually measured/calculated

The in-plane response is δsᵢ=Bz Σⱼ Cᵢⱼ Eⱼ, with i,j=x,y. The output is spin expectation per area in the existing Pauli convention, not automatically SI magnetic moment. Multiplication by ℏ/2 gives angular momentum density; converting to physical magnetization requires the chosen moment convention. Program units of C are (nm V T)⁻¹.

Operational definition: prepare a weak static Bz, then apply a weak low-frequency E. Differentiate the complete insulating spin–electric susceptibility with respect to B:

$$
C_{ij}=\left.\partial_{B_z}\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}
\left[-2\operatorname{Im}\sum_{v\in\mathrm{occ},c\in\mathrm{empty}}
\frac{(\sigma_i)_{vc}(\partial_{k_j}H)_{cv}}{(\epsilon_v-\epsilon_c)^2}\right]\right|_{B_z=0}.
$$

The electron charge is accounted for by using eV and V/nm. All calculations here are at T=0 with four occupied bands and an open indirect gap. μ is placed inside the gap for each parameter set; this does not change carrier filling. Electric and magnetic fields are infinitesimal response probes. The finite difference B step is numerical, not an experimental switching field.

`prepared_B_response` supplies the corresponding finite-frequency susceptibility. Its zero-frequency limit agrees with the static kernel, and layer-resolved results sum to the total. For AB, going from zero to ℏω=3 meV changes Cxx by about 3.7% while the interband gap is about 33 meV. This supports an explicit low-frequency readout protocol.

An independent mixed-Liouville EB solver is also checked at finite adiabatic rates (`bilayer_stage9_limit_check.json`). Equal-rate extrapolation from η=0.1 and 0.03 meV agrees with the static Cxx within 0.00031% for both AB and the generic small displacement. The initial much-smaller-rate scan is retained in `bilayer_stage9_protocol.json`: some values acquire numerical errors rather than converging monotonically. Those entries must not be interpreted as evidence of a physical order-of-limits effect. Neither scan establishes equivalence term by term with the published intrinsic/extrinsic split or repairs the author's omitted SM term.

The electronic layer dipole is Dᶻ/(e d_layer)=−(Ntop−Nbottom)/2 per cell. Its areal value is a 2D polarization (dipole per area). A physical layer spacing has not been invented. Within fixed-height point orbitals and rigid in-plane sliding, the ionic z dipole does not change, so the electronic change is meaningful. Absolute material polarization including relaxation is outside the model.

## 3. Coordinates and the complete AB→BA path

The code's a₁=δ₁−δ₂ and a₂=δ₁−δ₃ have a 60° angle. Define conventional 120° vectors A₁=a₁ and A₂=a₂−a₁. Then fractional coordinates transform as (u,v)₁₂₀→(u−v,v)₆₀.

Therefore conventional AB=(1/3,2/3) and BA=(2/3,1/3) are code coordinates (−1/3,2/3) and (1/3,1/3). The straight segment between these representatives has length one nearest-neighbour bond (0.142 nm). AB/BA labels depend on layer and orbital convention; using the coordinates without declaring the lattice basis is unsafe.

We sampled 17 points on this full segment at each of t⊥=5 and 20 meV. Both paths remain gapped (minimum sampled indirect gaps 31.83 and 29.40 meV). Both the dipole and C change sign under λ→1−λ. At the midpoint both vanish. At the endpoints Cxx=Cyy and off-diagonal terms vanish within numerical precision. In this model the endpoint tensor is isotropic; C3 alone would also allow an antisymmetric component, which additional point symmetries exclude here.

For 5 meV, the finer-grid AB values are Cxx=−2.47130×10⁻⁵ (nm V T)⁻¹ and Dᶻ/(e d_layer)=+2.33732×10⁻⁷ per cell; BA reverses both signs. Thus ordinary AB/BA switching is not evidence for a polarization-independent spin switch.

The occupied-electron energy is recorded for diagnostics only. It is not a total structural energy and must not be used to label AB/BA stable minima, compute barriers or infer switching fields.

## 4. Tiny sliding: what changes at AA versus AB

For a fixed generic direction through the nonpolar AA reference, C3 plus layer-exchanging mirror gives the previously verified minimum allowed powers Dᶻ∼u³, Tr C∼u³ and Cxx−Cyy∼u. Special directions or vanishing microscopic coefficients can suppress those leading terms further.

For a displacement q about polar AB, the correct statements are instead

$$
D^z(\mathrm{AB}+q\hat{\mathbf n})-D^z(\mathrm{AB})=p_2q^2+O(q^3),\qquad
\Delta\operatorname{Tr}C=c_2q^2+O(q^3),\qquad
C_{xx}-C_{yy}=n_1q+O(q^2).
$$

C3 forbids a scalar's linear gradient at AB. The layer mirror sends AB to BA and therefore does NOT make the local AB expansion odd in q. The traceless symmetric pair (Cxx−Cyy, Cxy+Cyx) can vary linearly. These sum/difference operations concern response directions, not the even/odd decomposition between magnetic domains.

Two independent fractional directions were tested. For direction (0.23,0.41), even-in-q dipole change divided by q² is −9.50, −9.50, −9.49 ×10⁻⁷ at q=.01,.02,.04; the odd directional difference divided by q is 5.635, 5.632, 5.623 ×10⁻⁶. Refining the grid changes the leading dipole coefficient by 0.066% and the spin coefficients by less than 0.0003%.

This can describe a small imposed shear or a displacement along a shear phonon near a polar structure. It does not prove that AB is stable or that the signal is large. At q=.04 in this direction, the anisotropy relative to |Tr C(AB)| is only about 0.45%; the polarization change is about 0.65%. An asymptotically lower power does not guarantee a larger effect at a finite displacement. Nor does u³ imply an exclusively third-harmonic response under oscillatory u.

## 5. Microscopic mechanism and its limits

### Opposite layer responses leave a small residual

At AA and t⊥=5 meV, Cxx(top)=−0.00972153 and Cxx(bottom)=+0.00972153; the total vanishes. At AB they become −0.00950962 and +0.00948491, leaving −2.47128×10⁻⁵. The net is approximately 0.26% of the magnitude of one layer contribution. BA exchanges the magnitudes. These are responses to EB, not zero-field magnetic moments or spontaneous magnetic order.

### Selection rules and exact comparisons

| Controlled change | Electronic dipole | In-plane EB spin response | Interpretation |
|---|---|---|---|
| Decouple layers | zero | total zero | opposite mirror partners cancel |
| Set Rashba to zero | need not vanish | zero | spin-z conservation blocks the in-plane spin vertex |
| Set staggered onsite m to zero | zero | zero | layer/sublattice inversion is restored even at generic displacement |
| Reverse both layers' Rashba signs | unchanged | changes sign | global σz spin rotation preserves energies/charge, reverses measured in-plane spin |
| Reverse m | changes sign | changes sign | inversion relates these Hamiltonians |
| Set intrinsic SOC to zero | may remain finite | may remain finite | intrinsic SOC tunes the response but is not necessary in the tested gapped model |

The Rashba-reversal comparison is exact in this model and proves that a spectrum and a layer-charge imbalance alone do not determine C. It is not a demonstrated experimental knob for reversing Rashba without changing other parameters.

### Why two hopping channels cannot simply be added

In the fixed physical A/B basis, write T=t_s V_s+t_c V_c (same- and cross-sublattice channels). The exact projected resolvent has self-energy Σtop(z)=T(z−hbottom)⁻¹T†. It contains terms proportional to t_s², t_c² and t_s t_c. The mixed term represents alternate virtual interlayer paths that connect the same states. Their relative phases change with stacking.

The response also involves spin and velocity vertices, which must be transformed together with the states; the self-energy alone is not a complete response formula. In the analytic weak-coupling regime the total response is even under simultaneous (t_s,t_c) reversal, and begins with a quadratic form in these amplitudes. Earlier Stage 7 independently recovered competing even and mixed contributions in a valley k·p calculation. Stage 8 confirmed that, near AA, adding the two isolated channels misses over 99% of the leading directional response. This is a controlled orbital-basis mechanism, not a uniquely observable decomposition into two geometric tensors.

With same-sublattice tunneling alone, T=f(k,d) times the sublattice/spin identity. A bottom-layer σz rotation makes the diagonal layer Hamiltonians identical, and a scalar layer phase makes f real at each k. The transformed Hamiltonian commutes with layer exchange; the occupied projector has equal layer charges pointwise. Thus the layer dipole is exactly zero in this limit. Do not use this pointwise transformation to drop the k-dependent connection from the electric vertex: it does not prove the spin-electric response is zero.

### A distinct effect at a fixed stacking

At λ=1/4 along AB→BA, changing t⊥ from 5 to 20 meV gives

| t⊥ | Cxx | Cyx | Dᶻ/(e d_layer) per cell | indirect gap |
|---|---:|---:|---:|---:|
| 5 meV | −1.77239×10⁻⁵ | −2.88609×10⁻⁷ | +1.68612×10⁻⁷ | 32.23 meV |
| 20 meV | −1.59427×10⁻⁴ | +1.14187×10⁻⁵ | +1.89445×10⁻⁶ | 31.12 meV |

All C values are (nm V T)⁻¹. The transverse component changes sign between 10 and 12.5 meV, while P keeps its sign and the gap stays open. The spin tilts to opposite sides of the −x direction (about +0.93° and −4.10° using arctan(Cyx/Cxx)); it does not fully reverse. This is a hybridization-controlled response anisotropy, not a topological transition. The local quadratic approximation in t is insufficient to describe the reversal. Because both state matrix elements and energy denominators change, this scan alone is not a quantitative separation of their roles. Earlier fixed-spectrum controls show that geometry matters, while the full sign change also depends on energy weights.

Increasing t⊥ is only a hypothetical proxy for changing layer separation/pressure. Real compression also changes other Hamiltonian terms; no pressure scale or feasible experimental amplitude is claimed.

## 6. What drives the structure, and what reads it

The model is H(k,d), not an equation of motion for d. The readout fields Eparallel and Bz are not assumed to cause sliding. A possible separate write operation uses perpendicular Ewrite or mechanical shear. With the 2D polarization pz and structural energy per area U₀, the gate contribution is U=U₀−pz Ewrite and the shear force per area is −∇d U₀+Ewrite ∇d pz.

Neither U₀ nor its barriers are provided by the current electronic model. At an exactly C3-symmetric polar endpoint, ∇d pz=0: a perpendicular field can favor one domain without directly providing a nonzero uniform first-order shear force at that ideal point. Fluctuations, a symmetry-breaking perturbation, and domain-wall motion can therefore matter. Wang and Dong's hBN first-principles study explicitly treats this issue through off-diagonal Born effective charges and assisted switching, rather than assuming field-driven sliding [PRB 111, L201406 (2025)](https://doi.org/10.1103/PhysRevB.111.L201406); [accessible manuscript](https://arxiv.org/html/2505.09084v1).

A sequential write/read proposal would apply a structural control pulse, remove it if a metastable state exists, and then read with weak Eparallel and Bz. If the gate stays on, its direct effects on layer potential and Rashba SOC must be added to H. Current calculations do not establish metastability, deterministic switching, or stability of intermediate stackings.

## 7. What is and is not a paper claim

Established: complete Hamiltonian and coordinate convention; clean insulating EB response with a readout protocol; AA and AB expansions; full AB/BA sign relations; microscopic null controls; layer cancellation; an interlayer-hybridization-dependent transverse sign change at fixed stacking.

Not established: a new quantum-geometric tensor, a universal law for all layer groups, a fully analytic few-band closed coefficient, real material amplitudes/stability, or quantitative write-field requirements. The published [NMEE work](https://arxiv.org/html/2507.23415v1) already establishes EB-induced spin response in time-reversal-invariant inversion-broken insulators. Merely adding a bilayer and polarization reversal would not by itself establish a strong new novelty claim.

The promising question is narrower: can a symmetry-resolved EB spin response diagnose or control changes of interlayer hybridization that polarization alone does not specify, with a measurable amplitude in stable structures? The current model supplies candidate signatures and falsifiable limits, not a finished affirmative material answer.

## 8. Reproduction and evidence

Run from repository root with PYTHONPATH including `papers/jia-2026-nonlin-magnetoelectric-edelstein` and `studies/bsf-nmee-connection`. Use one BLAS thread. Execute `scripts/run_stage9_path.py`, `scripts/run_stage9_protocol.py`, `scripts/run_stage9_coupling.py`, and `scripts/run_stage9_limit_check.py` under the study directory. These create the corresponding JSON files. Run `scripts/plot_stage9_physics.py` for the figure and source CSV.

The six new direct physics tests are `tests/test_stage9_protocol.py`; saved-evidence convergence and interpretation checks are in `tests/test_stage9_evidence.py`. Earlier stage manifests remain untouched. Main grid: C3-closed exterior 61, radial 60, angular 96; finer grid: 101/100/144. Finite B step 0.0625 T, with coupling fine controls also halving it. Layer/frequency and independent-solver checks use 41/45/72 and are numerical protocol checks, not the source of the main high-accuracy figure.
