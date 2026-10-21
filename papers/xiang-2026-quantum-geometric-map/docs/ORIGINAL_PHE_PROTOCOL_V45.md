# v45 | Original spin-induced PHE: explicit static-B / RTA protocol audit

**Source:** Xiang, Jia, Xu & Wang, *Quantum geometric map of magnetotransport*, PRB **113**, L201406 (2026), DOI [10.1103/kd12-4ycj](https://doi.org/10.1103/kd12-4ycj). Independently checked main Eqs. (2), (7), (8), (10), (11) and the final RTA discussion against [arXiv:2510.02661v2](https://arxiv.org/pdf/2510.02661). **Scientific status remains PARTIAL.** Do not infer paper error or a unique experimental protocol without controlling all terms.

## 1. Correct scientific target

The paper's linear-in-B bilinear conductivity is \`j_a=σ^(i)_{ab,c} E_b B_c\`. The **MNHE** (\`τ^0\`) is antisymmetric in indices \`a,b\`; the **PHE** (\`τ^1\`) is **symmetric**. The specific main Eq.(8)/Eq.(11) tilted Dirac PHE target is \`σ^(1)_{xy,x}=σ^(1)_{yx,x}\`. Previous v41-v44 QWZ **antisymmetric linear Hall** constitutes an independent mathematical methodology, NOT an acceptance of the original paper's symmetric mixed-field PHE. Crucially, the authors explicitly state their map uses **phenomenological relaxation-time approximation** and that microscopic bath-coupled descriptions lie beyond its scope. A conserving SCBA-Středa solver is thus a *separate extension*, not a prerequisite for reproducing the printed RTA theory.

## 2. Exact uniform spin-field momentum shift

For main Eq.(10) at **fixed** tilt, including the source's spin perturbation \`H_B=-gμ_B B_x σ_x\`:

\[
H(\mathbf k,B_x)=t_xk_x I+v_F(k_x\sigma_y-k_y\sigma_x)-g\mu_B B_x\sigma_x
=H_0(k_x,k_y+\delta),\quad \delta=g\mu_B B_x/v_F.
\]

This identity is exact as a two-by-two matrix; the alternate SI Eq.(62) helicity gives the opposite sign shift. Therefore a translation-invariant equilibrium conductivity integrated over the **entire** continuum k-plane and computed with its **actual B-dependent velocities and equilibrated occupations** cannot depend on a constant Bx in this fixed-tilt ideal model. In particular the finite-T constant-τ Boltzmann Drude cross-channel integral

\[
\sigma_{xy}^{D}/(e^2\tau)
=\sum_s \int \frac{d^2k}{(2\pi)^2}v^x_sv^y_s[-f'(\varepsilon_s)]
\]

satisfies \`∂σ_xy^D/∂B_x=0\` after \`q_y=k_y+\delta\`. Direct finite-field quadrature tests this, and an independent integrand derivative **separates a nonzero velocity change and equilibrated-population change that cancel as a total q_y derivative**. At β=.125, μ=.30, T=.055, δ=.10, a 301×301 grid on [-2.5,2.5] has centered B response about \`-1.87e-19\` in normalized model units. At a 400×400 grid the derivative pieces are velocity \`+0.00496902533\`, population \`-0.00497012860\`, sum \`-1.10327e-6\`; its discretization residual decreases with a 600×600 grid. A tight hard cutoff produces spurious nonzero finite-field response that decays as the domain extends.

This is a **conditional no-go for a fully equilibrated static-B Drude protocol**, not an automatic refutation of the original *partial two-color Liouville* Eq.(8), whose electric/magnetic perturbation orderings and population preparation may define a different quantity. The physical distinction, including field-dependence of tilt, UV completion, field preparation, collision model and possible contacts, requires additional derivation. The full linear response at finite-B under any protocol cannot be inferred from only one selected geometric term.

## 3. Independent main Eq.(8) angular integration

Using main Eq.(10) direct eigenstates, \`s=±1\`, the Zeeman curvatures exactly agree with those printed below Fig. 2:

\[
Z_{s,-s}^{xx}=s k_xk_y/k^3,\quad Z_{s,-s}^{yx}=-s k_x^2/k^3,
\quad v^x_s=t_x+s v_F k_x/k,\quad v^y_s=s v_F k_y/k.
\]

Substituting the literal curvature and band velocity into **main Eq.(8)**, with the negative T=0 Fermi derivative and \`v_F=1\`, yields the independent angular integral (normalized by \`τ g μ_B\`):

\[
C_8(\beta,\mu>0)=\frac{1}{2(2\pi)^2}\int_0^{2\pi}
\frac{\beta\cos^2\theta-\cos\theta+2\cos^3\theta}
{1+\beta\cos\theta}\,d\theta
=\frac{\beta}{4\pi(1+\sqrt{1-\beta^2})^2}.
\]

The **printed main Eq.(11)** coefficient is the opposite sign,
\[
C_{11}=-\frac{\beta}{4\pi(1+\sqrt{1-\beta^2})^2}.
\]
At β=1/8, μ>0, \`C8=+0.0025064159066917\` while \`C11=-0.0025064159066917\`; both flip sign for μ<0. Independent eigenbasis matrix elements and angular quadrature match to 1e-12 across multiple beta. This supplements the earlier v37 *literal SI Eq.(63)* sign audit. The precise signed **physical** coupling \`-gμ_B B\`, manuscript helicity convention, and transport preparation may affect comparison: preserve both branches rather than silently “correcting” the authors.

## 4. Acceptance and next gate

New tests: 10. They cover main Eq.(10) + SM helicities, eigenvector Zxx/Zyx, Eq.(8) angular integral vs Eq.(11), beta/mu odd parity, full-plane finite-B Drude, UV cutoff artifacts, independent velocity/population cancellation, and invalid inputs. Run:

~~~bash
python -m unittest tests.test_dirac_phe_protocol_v45 -v
python scripts/reproduce_dirac_phe_protocol_v45.py --beta .125 --mu .3 --temperature .055 --zeeman-shift .1 --grid 301
~~~

**PASS:** explicit original-PHE formula and model shift identities under stated conventions. **OPEN:** a fully prepared stationary RTA mixed EB Liouville response with the magnetic equilibrium source, order of limits and gauge/current contacts; Fig. 2(b) literal Eq.(11) SI ordinate remains about \`2π\` above publisher curve (v37) and is **not** fixed by v45. Source paper itself does not prove a microscopic bath theory. **Next v46** should construct the *same-protocol* complete length-gauge mixed-field current, compare equilibrated-B and adiabatic-B protocols, test k-shift Ward and the published Eq.(8) branch without fitted prefactors. Do not contact authors or modify other paper modules.

## Follow-up v47 signed tensor reconciliation

v47 corrected a later (v46) electric-source derivative error and demonstrated that the **symmetric** cross-conductivity in the explicitly declared frozen-H0 RTA tends to the signed **main Eq.(11)** coefficient throughout beta=.025–.6. This does not itself resolve the independently computed literal main Eq.(8) signed branch, which may involve further field/connection conventions; do not erase the prior v45 calculation. See [v47 correction](ORIGINAL_PHE_SYMMETRIC_SOURCE_CORRECTION_V47.md). Source SM collision formula remains independently unavailable.
