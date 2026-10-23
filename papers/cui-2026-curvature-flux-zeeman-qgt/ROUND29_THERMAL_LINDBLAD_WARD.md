# Round29 | Finite GKSL: thermal detailed balance versus local charge conservation

**Whole original Cui paper: PARTIAL.** This is an independently constructed *finite one-particle* spinful SOC open-system model, not a microscopic electron-phonon KME.

## Motivation and precise scope

Round28's affine linearized \(-\gamma Q(\delta\rho-i[G,\rho_0])\) is phenomenological; it is **not** an autonomous, trace-preserving completely-positive thermal bath. Here we explicitly compare **two different genuine GKSL dissipators**, proving which mathematical constraints each satisfies and locating the missing source and gauge vertex. This contrast is **not a universal no-go theorem** for all collision kernels. It cannot settle the author's Eq.(19), materials or SI response.

A two-site single-particle \(4\times4\) Hamiltonian explicitly acts on *site* \( \otimes\) *physical spin-1/2* and includes spin-mixing SOC hopping and onsite spin fields. Site projectors \(P_0,P_1\) commute with an arbitrary local charge-gauge generator \(G=g_0P_0+g_1P_1\). All quantities use synthetic natural units \(\hbar=1\). The reference thermal state is the **normalized canonical one-particle** Gibbs state
\[
 \rho_\beta=e^{-\beta H}/\mathrm{Tr}[e^{-\beta H}],
\]
**NOT** the Fermi-Dirac many-electron covariance \(f(H)\). In particular these GKSL jumps are not being asserted as fermionic bath scattering operators for the actual metallic KME.

## GKSL and independently constrained baths

For jumps \(L_\ell\),
\[
 \mathcal D[\rho]=\sum_\ell\left(L_\ell\rho L_\ell^\dagger
 -\tfrac12\{L_\ell^\dagger L_\ell,\rho\}\right).
\]
This is a completely-positive, trace-preserving Markovian generator with coherent \(\dot\rho=-i[H,\rho]+\mathcal D[\rho]\).

**A. Local onsite dephasing:** \(L_i=\sqrt{\gamma_{\rm loc}}P_i\). It conserves charge on each site *in the dissipator*:
\[
 \mathrm{Tr}[P_i\mathcal D_{\rm loc}(X)]=0 \quad\text{for all }X.
\]
For generic hopping \([P_i,H]\ne0\), however, \(\mathcal D_{\rm loc}(\rho_\beta)\ne0\); it is unital and does not establish finite-temperature thermal detailed balance with respect to \(H\). Indeed the model's finite-\(\beta\) equilibrium energy drifts under dephasing.

**B. Canonical detailed-balance transitions:** let \(H|n\rangle=E_n|n\rangle\), and choose a declared Hermitian synthetic bath-coupling matrix \(A\). Introduce 12 rank-one jumps for \(m\ne n\),
\[
 L_{m\leftarrow n}=\sqrt{W_{mn}}|m\rangle\langle n|,\qquad
 W_{mn}=\gamma |\langle m|A|n\rangle|^2
 e^{-\beta(E_m-E_n)/2}.
\]
The rates obey the explicit *pairwise population detailed balance*
\[
 W_{mn}p_n=W_{nm}p_m,\qquad p_n\propto e^{-\beta E_n}.
\]
Therefore \(\mathcal D_{\rm th}(\rho_\beta)=0\), with exact trace preservation and complete positivity. These **postulated** KMS-type rates are **not** derived from a phonon spectral density, microscopic bath correlation functions, Fermi factors, Pauli blocking or disorder diagrams. Calling them an experimentally established Davies bath or material relaxation time would be unjustified.

Because rank-one transitions are nonlocal in the site basis, \(\mathrm{Tr}[P_i\mathcal D_{\rm th}(X)]\) need not vanish. The correct finite-site continuity law for an arbitrary \(X\) is
\[
 \partial_t n_i+\operatorname{div}J_i
   =\mathrm{Tr}[P_i\mathcal D_{\rm th}(X)].
\]
The bath source sums to zero globally, but must **not** be discarded in a local Ward identity. No invented bath bond current or closed-system Ward cancellation is used to hide this source.

## Gauge covariance and missing collision vertex

For a diagonal real local charge-gauge generator \(G\), \(U_\phi=e^{i\phi G}\), the thermal jumps must jointly transform \(L_\ell^\phi=U_\phi L_\ell U_\phi^\dagger\) with the Hamiltonian and state. Then
\[
 {\cal D}^\phi[X]=U_\phi{\cal D}[U_\phi^\dagger XU_\phi]U_\phi^\dagger,
\]
and its explicit linear *bath vertex* is
\[
 \left.\partial_\phi\mathcal D^\phi[X]\right|_0=
 i[G,\mathcal D(X)]-\mathcal D\big(i[G,X]\big).
\]
At thermal equilibrium \(\mathcal D(\rho_\beta)=0\), this contact cancels \(\mathcal D(i[G,\rho_\beta])\) for a pure gauge. **Freezing the jump matrices while rotating the state is an intentionally wrong-gauge negative control.** Gauge covariance here requires the *assumption* that the synthetic bath coupling also transforms consistently; a physical microscopic electromagnetic bath model remains unbuilt.

## Independent numerical acceptance gates

1. Two-site physical spin SOC and nondegenerate spectra; real onsite projectors and gauge transformations.
2. Separate GKSL action versus independent column-major **Kronecker Liouvillian**, independent classical population master equation, plus RK4 versus finite-dimensional matrix exponential.
3. Positive semidefiniteness, trace and Hermiticity under both model evolutions; pairwise KMS balance and Gibbs stationarity of the thermal jump bath.
4. Positive/negative contrasts: local dephasing has **zero bath site-charge source** but breaks finite-\(\beta\) Gibbs stationarity; thermal jumps preserve Gibbs but have a **nonzero site bath-charge source** for a localized state.
5. Gauge controls: simultaneous rotation of H, state and jumps is covariant; stale jumps are not. Independent finite-amplitude rotation derivative checks the **bath-vertex contact**, and its deliberate omission gives a finite pure-gauge residual.

The executable audit is [audit_thermal_lindblad_ward.py](scripts/audit_thermal_lindblad_ward.py), and the frozen checkpoint, once created from a first verified remote run, is [round29_thermal_lindblad_ward_audit.json](results/frozen/round29_thermal_lindblad_ward_audit.json).

## Hard limits

This is **two sites**, **one normalized canonical particle**, and chosen GKSL matrices, not a translationally converged metal, quantum Boltzmann kinetic theory, Ward-Takahashi self-energy and vertex, nonlinear magnetoelectric Edelstein, finite-frequency Eq.(19) KME or original author's real-material calculation. It demonstrates a practical design tension in these two chosen toy baths, not an impossibility theorem. The earlier 257 Cui tests are preserved unchanged. New tests and actual Python 3.11/3.12 counts must come from GitHub logs, not this document.

## First independent CI measured numbers (2026-10-23)

The first GitHub-hosted CI model audit (Python 3.11) produced the frozen [raw Round29 JSON](results/frozen/round29_thermal_lindblad_ward_audit.json). The frozen-data test is added only *after* these values were observed, and the final isolated test target is **274 = 257 prior + 16 new + 1 frozen regression**, requiring successful Python 3.11 **and** 3.12 actual CI jobs.

At β=0.8, γ_th=0.16, γ_loc=0.18 in the synthetic four-dimensional spin×site toy:

- Canonical pairwise detailed balance residual **3.47e-18** and thermal GKSL Gibbs stationarity **3.69e-17**.
- Onsite CP dephasing Gibbs violation **0.0232482** with positive initial canonical energy drift **0.0569797**.
- A localized input has thermal-jump site bath source **0.0748172** while the global source cancels to **2.78e-17**. Site continuity **including** the source is exact; *omitting* it produces **0.0748172** false residual.
- Jointly rotated jump matrices have U(1) covariance residual **8.18e-18**. Freezing the jumps yields **0.0108446** erroneous local gauge covariance; omitting the collision vertex produces pure-gauge residual **0.0178684**.
- Independent finite-amplitude bath-vertex derivative residual **5.11e-13**; dense GKSL exponential vs RK4 **2.86e-10**, smallest evolved eigenvalue **0.0364169** (finite-time positivity in the tested example).

The observed values certify this exact model protocol, not a material bath or a universal collision no-go. The source-paper status remains **PARTIAL**. The final CI status is recorded separately after both versions complete.
