# Round30 | Interacting fermionic Fock thermal bath, Pauli exclusion and gauge/magnetic vertices

**Accepted scope only after remote CI: exact four-mode interacting finite Fock-space model. Original Cui et al. paper remains PARTIAL.**

## Why this round

Round29 used a normalized **single-particle canonical** Gibbs reference. Such a model cannot by itself test Fermi occupation, Pauli exclusion, a genuine multi-electron interaction, or number-conserving collisions acting on more than one particle. Round30 passes into an explicit **four fermion modes, 16-dimensional many-body Fock space**, not merely four single-particle orbitals relabeled as many-body states.

This closes a concrete upstream theoretical gate, but does **not** claim that a postulated finite-size GKSL bath is a microscopically derived electron-phonon self-energy or the conserving vertex needed for the measured metallic KME.

## Exact CAR and physical interacting Hamiltonian

Let \(c^\dagger_{i\sigma},c_{i\sigma}\) create/annihilate electrons on two sites \(i=0,1\) and physical spin \(\sigma=\uparrow,\downarrow\). Jordan-Wigner occupation-basis signs \((-1)^{\sum_{j<a}n_j}\) are implemented **exactly**, so

\[
 \{c_a,c_b^\dagger\}=\delta_{ab},\quad
 \{c_a,c_b\}=\{c_a^\dagger,c_b^\dagger\}=0,\quad (c_a^\dagger)^2=0.
\]

Use a spin-orbit-coupled spinful hopping block \(h\) from the independently accepted Round29 toy and interaction

\[
 H=\sum_{ab} h_{ab}c_a^\dagger c_b
  +U (n_{0\uparrow}n_{0\downarrow}+0.87n_{1\uparrow}n_{1\downarrow})
  +Vn_{\mathrm L}n_{\mathrm R},
 \qquad N=\sum_a n_a,
\]

with declared dimensionless natural units \(U=0.48,V=0.19\). The Zeeman vertex is the **physical spin** \(\partial_B H=-\kappa\sum_{i} c^\dagger_{i\sigma}(\sigma_y)_{\sigma\sigma'}c_{i\sigma'}\), \(\kappa=0.7\), not an orbital pseudospin.

One-body noninteracting control: for \(U=V=0\), exact Fock grand-canonical \(\langle n_a\rangle\) equals the eigenstate Fermi-Dirac one-body density \(\big[f_\beta(h-\mu)\big]_{aa}\). Turning on U,V changes these occupations and the many-electron equilibrium, which cannot be inferred from a single-particle Fermi matrix alone.

## Explicit grand-canonical thermalized fermion sectors

The true interacting **grand canonical** reference density is

\[
 \rho_{\beta,\mu}=\frac{\exp[-\beta(H-\mu N)]}
 {\operatorname{Tr}\exp[-\beta(H-\mu N)]}.
\]

The number sectors \(N=0,1,2,3,4\) have dimensions **(1,4,6,4,1)**. An independently declared Hermitian number-preserving density coupling \(A=n_{\mathrm L}-n_{\mathrm R}+0.29(n_{0\uparrow}-n_{0\downarrow})\) generates **within-sector** rank-one transitions \(|m,N\rangle\langle n,N|\) of rate

\[
 W_{mn}=\gamma\,|\langle m,N|A|n,N\rangle|^2
 e^{-\frac{\beta}{2}(E_m-E_n)}.
\]

At \(\beta=0.83,\gamma=0.14,\mu=0.16\), these postulated rates satisfy
\(W_{mn}p_n=W_{nm}p_m\), where \(p_n\propto e^{-\beta(E_n-\mu N_n)}\). Hence the corresponding completely positive, trace-preserving GKSL generator obeys \(\mathcal D(\rho_{\beta,\mu})=0\). The exact CAR algebra enforces Pauli blocking for a transition to an already-occupied mode.

**Crucial limitation:** all jumps preserve **total** \(N\), so the bath cannot change the relative probabilities of number sectors. Therefore **every** Gibbs state with another chemical potential is also stationary at fixed \(\beta\); this bath does **not** independently thermalize to a prescribed \(\mu\). No grand-canonical particle reservoir or unique equilibrium attractor is claimed.

## Local continuity and bath-mediated transfer

Let \(n_\mathrm L=n_{0\uparrow}+n_{0\downarrow}\), \(n_\mathrm R=N-n_\mathrm L\), and coherent intersite charge current \(J_{\mathrm L\to\mathrm R}=-i[H,n_\mathrm L]\). For any density X,

\[
 \partial_t\langle n_\mathrm L\rangle+\langle J_{\mathrm L\to\mathrm R}\rangle
    =\operatorname{Tr}[n_\mathrm L\mathcal D(X)],
\]
\[
 \partial_t\langle n_\mathrm R\rangle-\langle J_{\mathrm L\to\mathrm R}\rangle
    =\operatorname{Tr}[n_\mathrm R\mathcal D(X)].
\]

The many-electron KMS bath conserves **total** charge \(\operatorname{Tr}[N\mathcal D(X)]=0\), and preserves the whole N distribution, but it can redistribute **local** site occupation. An electron-only continuity law without the explicit bath source is wrong. No unspecified bath current or microscopic lead is invented to cancel the source.

## Electromagnetic and Zeeman bath vertices

Local U(1) generator \(G=0.43 n_\mathrm L-0.32n_\mathrm R\) is diagonal in the Fock occupation basis. The gauge-covariant jump bath must transform its jumps and Hamiltonian **together** under \(U_\phi=e^{i\phi G}\). The exact contact is

\[
 \partial_\phi\mathcal D^\phi[X]\big|_{\phi=0}
  =i[G,\mathcal D(X)]-\mathcal D(i[G,X]).
\]

At equilibrium the contact cancels the unphysical dissipative response to a pure-gauge density-matrix rotation. Independently finite-difference actual gauge-rotated GKSL jumps. The intentionally stale jump set violates covariance.

**New magnetic contact:** The equilibrium \(\rho_{\beta,\mu}(B)\) and thermal jumps constructed from \(H(B)=H+B\,\partial_B H\) change **together**, giving the static differentiated stationarity identity

\[
 \mathcal D_0\!\left(\partial_B\rho_{\beta,\mu}\right)
 +(\partial_B\mathcal D_B)_{B=0}(\rho_{\beta,\mu})=0.
\]

This is independently checked using separate \(B>0\) and \(B<0\) diagonalizations and thermal jump construction, versus differentiation of \(\rho\). Omitting \(\partial_B\mathcal D_B\) gives a finite failure. **This is an equilibrium Zeeman thermal-collision contact, NOT a measured finite-frequency KME or source Eq.(19) response.**

## Independent numerical oracles and negatives

- Fermion anticommutation and explicit **blocked versus allowed** hopping in occupation basis.
- Grand Gibbs Fock partition vs **independent one-body Fermi-Dirac** in the noninteracting limit; interacting orbital-occupation change.
- Fock sector construction and exact zero cross-N jump rates, per-jump [L,N]=0 and general-state \(\mathrm{Tr}(N\mathcal D)=\mathrm{Tr}(N^2\mathcal D)=0\).
- Pairwise KMS master equation vs stationary interacting density; nonunique \(\mu\) negative control.
- Explicit local bath charge source in a two-electron input, full sitewise continuity, and deliberately omitted-source failure.
- Independent 256-dimensional **Kronecker Liouvillian action** vs direct looped GKSL.
- Joint local-gauge unitary covariance, stale jump negative, analytic vs finite-amplitude bath vertex, pure-gauge stationary cancellation.
- Independent finite-Zeeman-field Gibbs + KMS-jump derivative for equilibrium **magnetic** bath vertex, with omission negative.
- Reject non-Hermitian, number-violating and invalid temperature/rate input. A narrow nondegenerate-per-N regime is enforced; non-Abelian thermal jump treatment remains unresolved.

Numerical results and the mandatory regression: [round30_interacting_fock_bath_audit.json](results/frozen/round30_interacting_fock_bath_audit.json), measured from the first remote audit before freezing. Reproduction:

~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/audit_interacting_fock_bath.py
~~~

## Not accepted

No microscopic phonon or disorder bath correlation function, Fermi-golden-rule electronic collision integral derived from a material, Ward-Takahashi self-energy/vertex ladder, conserving finite-q metallic KME, real-source Wannier physical spin+position provenance, thermodynamic BZ limit, SI units, original Table I or Eq.(19). All computed values are **finite, interacting, synthetic many-electron diagnostics**. **Whole source paper PARTIAL.**
