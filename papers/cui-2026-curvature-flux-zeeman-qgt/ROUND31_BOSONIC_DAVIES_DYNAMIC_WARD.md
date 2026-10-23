# Round31 | Bosonic thermal correlation spectrum to interacting Fock collision and finite-frequency Ward

**Original Cui paper PARTIAL.** This is a **specified weak-coupling harmonic-boson model**, not an ab-initio phonon calculation or a diagrammatic proof of microscopic material KME.

## Narrow physical advance after Round30

Round30's within-number-sector transition rates were *postulated* as \(W_{mn}=\gamma|A_{mn}|^2e^{-\beta(E_m-E_n)/2}\). This guaranteed a KMS ratio but did not provide a bath correlation function or frequency-dependent microscopic absorption/emission weights.

In Round31 the **environment Hamiltonian and coupling are explicitly specified**, albeit as a simplified continuum model:
\[
H_{\rm tot}=H_{\rm Fock}+\sum_q\omega_qb_q^\dagger b_q+
A\otimes \sum_q g_q(b_q+b_q^\dagger),\qquad
A=n_{\rm L}-n_{\rm R}+0.29(n_{0\uparrow}-n_{0\downarrow}).
\]
Here \(H_{\rm Fock}\) is the **same 16-state interacting fermion model** (four physical spin modes with SOC, Hubbard \(U=0.48\), intersite \(V=0.19\)) accepted in Round30. \(A\) is a number-preserving local electron density coupling and \([A,N]=[H_{\rm Fock},N]=0\).

### Thermal bath spectral derivation

The explicit spectral-density **model assumption** is a smooth Ohmic continuum
\[
\mathcal J(\nu)=g_2\nu e^{-\nu/\omega_c},\quad \nu>0,\qquad
g_2=0.045,\quad\omega_c=2.4.
\]
For a thermal harmonic oscillator, the unsymmetrized bath correlator is
\[
C(t)=\int_0^\infty d\nu\,\mathcal J(\nu)
\left[(n_B(\nu)+1)e^{-i\nu t}
+n_B(\nu)e^{i\nu t}\right],
\quad n_B(\nu)=\frac1{e^{\beta\nu}-1}.
\]
An independent **finite 48-level oscillator Hamiltonian and thermal trace** checks these two Fourier weights, their KMS ratio and the exact matrix-time two-point function. This is an independent boson-space oracle, not a fit to the fermion response.

The two-sided Fourier transform of the continuum correlation, within the **Born, Markov, weak-coupling and fully secular** assumptions, gives a transition \(|n\rangle\to|m\rangle\), \(\Delta=E_m-E_n\),
\[
\Gamma(\Delta)=2\pi\,\mathcal J(|\Delta|)
\begin{cases}
 n_B(|\Delta|), & \Delta>0\quad\text{(boson absorption)},\\
 n_B(|\Delta|)+1, & \Delta<0\quad\text{(boson emission)}.
\end{cases}
\qquad W_{mn}=|A_{mn}|^2\Gamma(\Delta).
\]
Thus \(\Gamma(+\nu)/\Gamma(-\nu)=e^{-\beta\nu}\) follows from a **declared oscillator thermal correlator**, rather than being imposed by hand on an otherwise arbitrary transition rate. The Ohmic zero-frequency limit is \(\Gamma(0)=2\pi g_2/\beta\), included as a *single collective* dephasing jump from the energy-diagonal component of \(A\). A small per-number-sector nondegeneracy gate avoids arbitrary eigenvector gauges at exact degeneracies; clustered/near-degenerate nonsecular physics is unresolved.

The pairwise **total number \(N\)-conserving** jump Lindbladian is Gibbs stationary for the interacting \(e^{-\beta(H-\mu N)}\), but **does not select a unique chemical potential** or supply a particle reservoir. Electron bath coupling is given explicitly, but the material-dependent \(g_q\), \(\omega_q\), phonon eigenvectors, scattering self-energy, Lamb shift and nonsecular memory are **not** computed from first principles. Naming this model "ab initio electron-phonon collision" would be false.

Standard weak-coupling Davies and thermal-bath correlator context: [Albash et al., *New Journal of Physics* 14, 123016 (2012)](https://doi.org/10.1088/1367-2630/14/12/123016); [OpenQuantumTools bath-correlator notes](https://uscqserver.github.io/OpenQuantumTools.jl/dev/basics/bath.html). Recent nonsecular caveats: [*Physical Review X* 16, 011040 (2026)](https://doi.org/10.1103/sfp3-3sqf). These are **background references**, not source-Cui claims.

## Strongly independent finite-Fock and charge checks

- Exact 16-state four-mode fermionic Fock basis and 0–4 particle-number sectors (1,4,6,4,1) are retained, with *physical spin* SOC/Hubbard operators.
- Pairwise \(\Gamma(+\nu)/\Gamma(-\nu)\) KMS, rate matrix detailed balance in **grand** Gibbs weights, explicit collective \(\omega=0\) jump, stationary interacting Gibbs density and exact conservation of total \(N\), including \(N^2\).
- The bath transfers fermions **between sites** in the reduced subsystem. Therefore the sitewise charge Ward reads
\[
\partial_t\langle n_L\rangle+\langle J_{L\to R}\rangle
=\mathrm{Tr}(n_L\mathcal D[\rho]),
\]
and similarly for \(R\). The explicit nonzero bath source must be retained. We do **not** invent a separately derived bath bond-current operator to hide it.
- A cross-check against Round30's \(e^{-\beta\Delta/2}\)-weighted *postulated* rates verifies actual energy dependence of the boson spectrum beyond merely satisfying detailed balance.

## Dynamic finite-(q,omega) local electromagnetic Ward and bath gauge vertex

On the two-site torus, let \(G=0.43n_L-0.32n_R\), \(U=e^{i\phi G}\). For a time-dependent **pure longitudinal gauge** with \(e^{-izt}\), \(z=\omega_{\rm drive}+i\eta,\ \eta>0\), the electronic Hamiltonian vertex is
\[
 V_G=i[G,H]+izG,\qquad
 (\partial_\phi\mathcal D_\phi)[\rho_0]=
 i[G,\mathcal D(\rho_0)]-\mathcal D(i[G,\rho_0]).
\]
The **same** full bath jump operators are co-rotated by \(U\). With this explicit bath-vertex source the full retarded Fock Liouvillian response
\[
 [-iz-\mathcal L_0]\,\delta\rho=
 -i[V_G,\rho_0]+(\partial_\phi\mathcal D_\phi)[\rho_0]
\]
must reproduce the exact gauge rotation \(\delta\rho=i[G,\rho_0]\). **Both** site densities and the physical coherent bond current including \(i[G,J]\) contact must vanish for pure gauge. The model's nonzero bath source is included in local continuity, and both omission of the bath-gauge vertex and freezing the jump basis are negative controls.

This is a genuine test of **finite nonzero \(z\)** for a toy Fock model, rather than claiming that adding \(i\eta\) to a collisionless response constitutes a microscopic scattering theory.

## Magnetic finite-frequency source and independently differentiated thermal collision vertex

For a physical Pauli-spin Zeeman \(V_B=\partial_B H=-0.7\,c^\dagger\sigma_y c\), construct **independent \(B=\pm\epsilon\)** interacting Fock Hamiltonians and *rederive* the Ohmic-boson jump matrices from their field-dependent energy gaps. The resulting
\[
 \mathcal D_0[\partial_B\rho_{\beta,\mu}]
 +(\partial_B\mathcal D_B)[\rho_{\beta,\mu}]=0
\]
is separately checked, as is the **nonzero-frequency** response
\[
[-iz-\mathcal L_0]\delta\rho_B
=-i[V_B,\rho_0]+(\partial_B\mathcal D_B)[\rho_0].
\]
The charge Ward **with the bath source** is tested at finite \(z\). Omitting the B-dependent bath contact must change the response significantly. This is an **equilibrium and driven synthetic Zeeman response**, NOT the author's KME current tensor or an SI units result.

## Independent failure modes and acceptance scope

The gate rejects invalid inverse temperature, negative bath spectral density, invalid boson cutoff and non-Hermitian/number-nonconserving Fock inputs. It compares finite oscillator exact correlation vs spectral weights, KMS vs system rates, source-complete continuity, co-rotated jump covariance, finite gauge vertex vs exact rotation, finite-\(z\) gauge density/current and magnetic contact vs omission.

**Still missing:** actual phonon dispersions and electron–phonon matrix elements, a demonstrated self-energy/vertex Ward-Takahashi identity for a translational metallic continuum, Lamb shifts, nonsecular dynamics near degeneracies, unique \(\mu\) relaxation, experimentally meaningful SI KME and author's Eq.(19), Table I or material figures. The results are **finite synthetic many-electron and model-bath calculations**. The full Cui source remains **PARTIAL**.

Reproduction:
~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/audit_bosonic_davies_dynamic_ward.py
~~~
