# Round28 | Local-charge-conserving, gauge-covariant *linearized* collision

**Paper status PARTIAL. Scientific stage: a finite synthetic collision sanity check, not a microscopic bath.**

## What was missing after Round27

Round27 proved finite-(q,omega) electromagnetic pure-gauge invariance without collisions. Replacing z by omega+i(gamma+eta) does **not** prove that a collision operator respects local charge continuity, and damping an equilibrium density-matrix *response* directly is not gauge covariant. We isolate these two logically independent failure modes before any attempted dynamic Zeeman-induced current or metallic KME claim.

## The explicit physical and mathematical protocol

Use the accepted Round27 spin-1/2×site 3×3 periodic SOC lattice with n=18 one-particle basis states. Projector P_i is the two-spin-dimensional projector onto site i. Let rho0=f(H), T>0, and the local-block map

\[
 {\cal D}(X)=\sum_i P_i X P_i,\qquad {\cal Q}=1-{\cal D}.
\]
\[
  \mathrm{Tr}[P_i{\cal Q}(X)]=0\qquad\text{for every site }i.
\]

The explicitly declared **affine linearized** phenomenological collision is

\[
 {\cal C}[\delta\rho;G]=-\gamma {\cal Q}\{\delta\rho-i[G,\rho_0]\}.
\]

G is the real diagonal **Peierls longitudinal gauge generator** for the external vector potential. Under a site-local electromagnetic gauge transformation by real diagonal L,

\[
 G\mapsto G+L,\quad V\mapsto V+i(z-\operatorname{ad}_H)L,
 \quad \delta\rho\mapsto\delta\rho+i[L,\rho_0].
\]

The local collision projection commutes with U=e^{iL}; thus both collision covariance and local continuity are preserved. This is **more restrictive than saying gamma is a scattering time**. Q is a site-dephasing projector, but the additional rho0/G-dependent affine forcing is inserted to fix the chosen thermal reference. The resulting full driven equation is **not claimed to be a completely-positive autonomous Lindblad thermal bath, a detailed-balance electron-phonon collision integral, or a Ward-derived microscopic self-energy/vertex pair**.

For drive convention e^{-izt}, z=omega_drive+i eta, Im z>0 and hbar=1, solve

\[
  (z-\mathrm{ad}_H+i\gamma {\cal Q})\,\delta\rho
    =[V,\rho_0]+i\gamma{\cal Q}\bigl(i[G,\rho_0]\bigr).
\]

This has a strict local lattice Ward consequence,

\[
 \operatorname{div}\delta J-i z\,\delta n_i=0,
\]

where delta J includes the explicit Peierls current contact i[G,J] when the vector generator G is applied. The source-paper's omega_source=-omega_drive convention from earlier work is **not** silently changed by using z here.

## Two genuinely independent numerical solvers

1. **Spectral fixed-point:** diagonalize H and use R0(X)=(z+i gamma-ad_H)^(-1)X. Iterate X=R0(b+i gamma D(X)) to the stated convergence threshold, with b=[V,rho0]+i gamma Q i[G,rho0]. For eta>0, the iteration bound gamma/(gamma+eta)<1.
2. **Orbital Liouville:** construct the full n²×n² complex Kronecker commutator and independently assemble Q=I−sum_i P_i^T⊗P_i, then directly solve the dense linear system. This does not reuse the spectral fixed-point kernel.

Exact sitewise charge-collision traces, spectral-vs-Liouville agreement, finite nonzero scalar/Zeeman response, collisionless Round27 limit, variable gamma, local U(1) covariance and the full bond-current contact are tested.

## Critical negative controls

- **Naive RTA**: use Q=I with no covariant reference. Then the collision *directly* changes local charge and violates div(delta J)-iz delta n=0 under a non-pure longitudinal scalar perturbation.
- **Locally conserving but non-gauge-covariant**: use Q=I−D but mistakenly set the reference i[G,rho0]=0. Sitewise continuity still holds, yet the full pure gauge produces an erroneous physical current. Charge conservation **alone** does not imply electromagnetic gauge covariance.
- **Genuine physical-spin magnetic source**: apply -kappa B_i sigma_y at each site. Verify nonzero magnetic density-matrix response and local *charge* Ward. Add a pure electromagnetic gauge to this magnetic drive: the charge/current observables must be unchanged, including Peierls contact. This is a **finite synthetic Zeeman charge-response check**, not source KME Eq.(19) in a thermodynamic material.

## Boundaries and acceptance

- *Positive claim*: a conserving, gauge-covariant **linearized model relaxation** and corresponding finite-dimensional microscopic Ward on a synthetic SOC torus, valid for Im z>0 and the indicated source/Peierls protocol.
- *No claim*: completely positive thermal Markovian bath; self-energy/vertex diagrams satisfying microscopic Ward-Takahashi identity; detailed balance, electron-phonon/impurity kinetics, dynamical full-metal intraband Drude, thermodynamic BZ convergence, real H/S/r Wannier operators, SI GHE or experimentally measurable KME.
- This also does **not** promote the entire original Cui paper above **PARTIAL**.

Reproduce:
~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/audit_conserving_collision.py
~~~

Current new tests and frozen benchmark acceptance must be verified against **actual remote Python 3.11/3.12 jobs**, in addition to all existing Round17–27 audits. No old Round21 occupation_partition changes or rewrites.
