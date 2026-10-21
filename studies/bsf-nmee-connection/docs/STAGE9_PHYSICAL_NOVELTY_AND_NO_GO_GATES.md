# Stage 9 concept gate: what physics is genuinely new?

**Status:** Critical research planning and exact limiting-case analysis, not a new real-material calculation. The original Jia et al. PRL NMEE reproduction remains frozen and its driven-DC/static discrepancy unresolved. Source of Stage8g energies: user-reported Japan IMR VASP OUTCAR audit, NOT independently opened original OUTCAR by this repository note. No new HPC jobs.

## I. Three physically different quantities must not be conflated

**LEE, transport:** current-induced NONEQUILIBRIUM layer-resolved spin accumulation in a metallic or carrier-doped state; depends on occupations, relaxation/broadening, and the current/electric-field normalization. Refer to Zhou et al., *Layer Edelstein Effect*, PRL 137, 076401 (2026), https://doi.org/10.1103/qhvs-pfkq .

**Published NMEE:** cooperative electric and magnetic fields producing spin magnetization; Jia et al., PRL 136, 226302 (2026), https://doi.org/10.1103/r3x2-kb93 . Intrinsic and extrinsic terms have their own dynamical/relaxation definitions. The existing frozen reproduction explicitly identifies an unanswered diagonal density-matrix / static-vs-driven issue. Do not silently identify either of them with the next observable.

**Our validated model observable:** gapped, 0-K, **equilibrium/static** spin-response variation with a uniform electric field and a spin-only perpendicular Zeeman perturbation, via mixed Berry curvature on the *entire occupied manifold*:
\[
\Gamma^{\rm stat}_{iz;j}=\partial_{B_z}\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}\,
i\,\mathrm{Tr}P[\partial_{\lambda_i}P,\partial_{k_j}P].
\]
This is closer to a derivative of modern electric polarization with respect to a spin source and magnetic field, or equivalently an equilibrium electric-field modulation of spin magnetic susceptibility (after correct electron-charge, spin, g-factor and orbital magnetic operator conventions), than to a demonstrated DC current-induced torque. In a gapped equilibrium insulator there is no bulk steady conduction current. The spin-only source is not a complete physical magnetic-field treatment, which can also include orbital coupling.

A **new mandatory theory gate** is to derive the frequency/order-of-limits relation between thermodynamic static, clean driven DC and disorder/relaxation-limited responses. Test any proposed identity in the full existing Kane–Mele model, not by importing a formula whose diagonal term was previously unresolved. If protocols are inequivalent, rename and develop the static nonlinear magnetoelectric response honestly; do not claim a direct reproduction of driven-DC NMEE.

## II. A strict spin-conservation null test overlooked by SOC-based candidate ranking

Consider \([H(\mathbf k),S_z]=0\) for all k, where \(S_i\) are physical Pauli spin operators. Since \(v_j=\partial_{k_j}H\), \([v_j,S_z]=0\), and with a perpendicular spin Zeeman \(H_B=-g\mu_BB_z S_z\) the same blocks are preserved. For transverse \(S_x\), only *opposite* spin sectors have nonzero spin-flip matrix elements, while \(v_j\) connects states with the *same* spin. Hence for each occupied-empty Kubo pair
\[
\langle v|S_x|c\rangle\langle c|v_j|v\rangle=0
\quad\Longrightarrow\quad
\mathcal F_{\lambda_xk_j}(B_z)=0,\quad
\Gamma_{xz;j}^{\rm stat}=0.
\]
This remains true for layer-odd transverse spin \(S_x^{(-)}=\Pi_T S_x-\Pi_B S_x\) **provided physical layer projectors commute with \(S_z\)**.

This is a mathematical selection rule **independent of whether the stack has C3v and independent of an Ising SOC splitting magnitude**. A realistic candidate needs actual nonconservation of \(S_z\), e.g. Rashba/non-Ising SOC and/or spin-mixed orbitals. W's strong atomic SOC and a large valley spin splitting alone DO NOT certify a sizable \(S_x\) response. The 2019 paper by Lin et al., PRB 100, 155408, https://doi.org/10.1103/PhysRevB.100.155408, establishes that electronic ferroelectricity can couple to Rashba splitting in bilayer TMDs, so non-Ising physics is plausible, but quantitative spin mixing must be *measured* from genuine SOC spinors rather than presumed.

Cheap **pre-Wannier go/no-go**: on actual relaxed AB/mid/BA WSe2 with SOC, inspect convergence of full-BZ gap, \(\langle S_x\rangle,\langle S_y\rangle,\langle S_z\rangle\) spin textures and off-diagonal spin/velocity transition-matrix norms (especially Gamma and Lambda as well as K/Kprime). A spin expectation texture alone is not sufficient to compute the Kubo response; the off-diagonal operators and remote-band window are indispensable.

## III. Layer spin is a physical spatially resolved observable, not an arbitrary Wannier label

Layer-even/odd spin operators require a defensible physical projection into top/bottom slabs; changing a Wannier orbital naming convention is NOT a physical change, but choosing different spatial layer regions can affect reported magnitudes when wavefunctions hybridize. Require (i) the same physical spinor Hilbert space, (ii) a fixed real-space layer window or PAW atom-projected spin-density definition, (iii) a sensitivity check to the layer partition near the interlayer region, (iv) the correct exchange of layer projectors by the real sliding domain symmetry, and (v) total = top+bottom after accounting for interstitial contribution. Refuse quantitative 'hidden layer spin' if only diagonal layer weights from HR are available.

Stage8f glide rule is valid for a truly layer-exchange-symmetric *actual* structure and a symmetry-covariant spin-density definition:
\[
M_{ij}^{(+)}(s=1/2)=0,\qquad M_{ij}^{(-)}(s=1/2)\ {\rm symmetry\ allowed}.
\]
Allowed is NOT necessarily nonzero; one must not infer material magnitude from a toy model.

Stage8h must independently verify post-relaxation midpoint glide, not just the ideal POSCAR. High midpoint energy means it is likely a transient sliding intermediate or a domain-wall-like local stacking, NOT an available stable memory state. Domain walls and devices may break local glide through strain/substrate/electric field; an isolated commensurate-cell symmetry theorem is not automatically a device-scale domain-wall theorem.

## IV. Novelty check versus existing 2019–2026 literature

- Lin et al. (PRB 100, 155408, 2019) already predict FE-controlled Rashba spin splitting/texture in MX2 bilayers: https://doi.org/10.1103/PhysRevB.100.155408 .
- Zhou et al. (PRL 137, 076401, 2026) introduce and calculate nonequilibrium bilayer hidden layer Edelstein responses: https://doi.org/10.1103/qhvs-pfkq .
- Pan et al. (PRB 114, 115135, 2026) address multistate ferroelectric control of Edelstein and spin Hall charge-spin conversion: https://doi.org/10.1103/np2g-86h3 .
- Jia et al. (PRL 136, 226302, 2026) propose intrinsic/extrinsic NMEE: https://doi.org/10.1103/r3x2-kb93 .
- Li et al. (Nat. Commun. 17, 2457, 2026) report experimental AB/BA 3R WSe2 domain polarization and valley band edges: https://doi.org/10.1038/s41467-026-68854-x .

Therefore the novelty **cannot** be merely '3R stacking breaks inversion', 'AB/BA flip a response', 'strong W SOC', or 'hidden layer polarization allowed by glide'. Stronger potential story: *nonlinear (electric × perpendicular magnetic) equilibrium spin susceptibility with domain-odd TOTAL and domain-even LAYER-STAGGERED sectors, a symmetry-protected total-spin node but an actual finite staggered signal, and an independently validated microscopic quantum-geometric mechanism with nontrivial remote-band cancellation*. Distinguish this from LEE current-driven transport and from known FE-Rashba spin textures. Necessarily calculate an actual finite staggered quantity, not only a symmetry allowance.

## V. Proposed minimal stop/go decision

**Gate 0 (no new DFT):** independent raw OUTCAR/CONTCAR symmetry and force audit from Stage8h. Confirm midpoint glide and net upper-layer in-plane force; do not reinterpret 24.117 meV/cell three-point increase as NEB barrier.

**Gate 1 (small numerical task):** convergence of basic constrained PES differences to a few meV per 2-formula-unit cell, including midpoint residual Fz and an explicit one-arm mirror relation; for any truly quantifiable model, full ionic endpoints and the glide stability after relaxation.

**Gate 2 (most important electronic-structure falsification):** single properly relaxed 3R-WSe2 AB spinor-SOC electronic structure and a symmetry-compatible midpoint, with **spin-mixing and full-BZ insulating gap** examined BEFORE attempting costly spinor Wannier interpolations. If transverse source response has negligible matrix elements, switch material/observable rather than enlarge WSe2 batch.

**Gate 3 (defining manuscript):** compare total and physical layer-resolved \(\Gamma^{\rm stat}\) across AB → glide midpoint → BA with full Brillouin zone, remote-band and gauge-covariant Wilson/Kubo agreement. Check domain parity, scalar/nematic irreps, robustness to layer projectors and sensible experimental sheet units. An honestly meaningful nonzero result plus a mechanism is needed for material novelty.

**Gate 4 (optional higher-ambition):** establish rigorously whether/which dynamic clean/disordered \(\omega\to0\) limit coincides with equilibrium static susceptibility. Only then use 'NMEE/Edelstein' as a transport interpretation. Experimental protocols may instead target \(F_{\mathrm{mix}}\sim-\Lambda (\mathbf E_\parallel\cdot\mathbf B_\parallel)B_z\), where the spin-only vs full orbital response and device boundary conditions are distinguished.

### Decision rule

Do **not** automatically prioritize 12 more z-only slide jobs, CI-NEB, broad material screening, or deep Wannier convergence merely because the pipeline runs. We have structure/symmetry feasibility and a first real energy profile, but not a response magnitude or direct link to measured Edelstein spin accumulation. **The decisive question is whether real WSe2 generates a finite, converged, layer-resolved equilibrium static cross response with nontrivial origin that was not already established by FE-controlled Rashba studies.**

**Scientific status:** mechanism hypotheses physically consistent, WSe2 structural calculations reported but raw symmetry check pending, novelty conditional on future quantitative evidence, and no actual real-material NMEE yet.
