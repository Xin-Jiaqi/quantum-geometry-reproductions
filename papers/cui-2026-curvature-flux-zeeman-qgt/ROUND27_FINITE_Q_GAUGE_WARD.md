# Round27 | Finite-(q,ω) microscopic longitudinal U(1) Ward

**Scientific acceptance: synthetic finite-periodic-lattice charge-gauge method only. The Cui paper remains PARTIAL.**

## Actual gap

Rounds 23–25 verified the complete *static equilibrium* Zeeman-current no-bulk-current Ward. Round26 checked a declared same-Wannier-gauge H/physical Pauli spin/position input interface on a toy. **None** verified longitudinal finite-(q,ω) charge conservation and electromagnetic contact terms. This is a necessary guard before claiming a response in a measurable electromagnetic protocol, but it does **not** supply a conserving scattering kernel or turn the old clean interband gyrotropic susceptibility into macroscopic KME.

## Explicit physical protocol and sign convention

Use a **finite 3×3 periodic two-dimensional torus** with one physical spin-1/2 at each site, nearest-neighbor hopping T_x=t_x I+iα σ_y and T_y=t_y I−iβ σ_x, onsite mσ_z and weak nonuniform scalar potential. SOC spin is not conserved, but U(1) particle number is. Operators act on the 18-dimensional *single-particle* site×spin Hilbert space. Work in natural units ħ=1; charge is absorbed into the potential/Peierls gauge generator; energies and complex frequency have the same model units. No physical conductivity volume or SI units.

For a nonuniform real site-phase L_i=cos[2π(q_x x/N_x+q_y y/N_y)+φ], q_x,q_y integers and q≠0, let U(t)=exp[iλ(t)L]. From the undriven H obtain an **exact pure-gauge Hamiltonian**

\[
H_g(t)=U(t)HU^\dagger(t)-\dot\lambda(t)L.
\]

At first order in the drive λ(t)∝e^{-izt}, z=ω_drive+iη with η>0, its *two necessary terms* are

\[
K=i[L,H]\quad (\text{Peierls/vector}),\qquad
V_{\phi}=izL\quad(\text{scalar}),\qquad
V_g=K+V_{\phi}=i(z-\operatorname{ad}_H)L.
\]

The source-paper Fourier variable in earlier Cui response comparisons is **ω_source=−ω_drive**. This Round27 check uses **drive** z throughout; it does not identify z with the source-paper positive-frequency retarded quantity.

With f(H) the thermal equilibrium one-particle density matrix and \(\operatorname{ad}_HX=[H,X]\), the exact collisionless Liouville result is

\[
\delta \rho=(z-\operatorname{ad}_H)^{-1}[V_g,f(H)]
           =i[L,f(H)].
\]

As the local projectors n_i commute with L, their total pure-gauge density response vanishes. The physical oriented bond-current is J_ij=i(H_ij |i⟩⟨j|−h.c.), with an *explicit Peierls current-operator contact*

\[
\delta J_{ij}^{\rm contact}=i[L,J_{ij}],
\qquad
\delta\langle J_{ij}\rangle=
\mathrm{Tr}(J_{ij}\delta\rho)+
\mathrm{Tr}(\delta J_{ij}^{\rm contact}f(H))=0.
\]

The exact local operator continuity equation holds before integrating or imposing symmetry:

\[
i[H,n_i]+\sum_{j\in{\rm outgoing}}J_{ij}
 -\sum_{j\in{\rm incoming}}J_{ji}=0,
\quad
\operatorname{div}\delta J-i z\,\delta n=0.
\]

The second relation includes the appropriate explicit contact term when the perturbation enters hopping. The complex z is an analytic retarded probe; **η is not a conserving scattering bath** and cannot supply Drude relaxation. No real-measurement conductivity is extracted from a pure gauge.

## Independent numerical and failure oracles

- Kubo spectral pair differences vs independent **orbital-basis N² Liouville matrix inversion**, no shared diagonal response formula.
- Full local operator continuity at each site and complex-frequency continuity separately for scalar, vector and their sum.
- Independently finite-difference **exact finite-amplitude Peierls Hamiltonian and bond-current operators**; separately test time-dependent gauge-rotated density matrix against the full Liouville equation.
- Several nonzero lattice q and complex drive frequencies, SOC torque nonconservation, error gates on invalid q and retarded frequency.
- **Negative controls:** omit scalar potential and obtain finite spurious charge density; omit Peierls current contact and obtain finite spurious current and a local continuity violation. Merely projecting a current onto a symmetry-allowed tensor would conceal these errors.
- Frozen [raw audit](results/frozen/round27_finite_q_gauge_ward.json), executable [script](scripts/audit_finite_q_gauge_ward.py), **16 additional Cui-specific tests**.

For the 3×3 torus at (q_x,q_y)=(1,1), z=.42+.09i:
- spectral vs full orbital Liouville residual 3.73e-15;
- pure-gauge density response 2.22e-15; current response 3.43e-15;
- scalar and vector *separate* density amplitudes ~.31981, mutually canceled by gauge completion;
- missing current contact: spurious bond current max **1.34659**, local finite-q continuity residual **3.59468**;
- contact-complete vector continuity residual 2.39e-15.

These are dimensionless/natural-unit **finite system** diagnostics, not material GHE/KME.

## Non-claims

No genuine sample Wannier physical spin or position, magnetic-field Zeeman KME at finite q, transverse Maxwell/diamagnetic tensors outside the longitudinal pure-gauge channel, metal transport/Drude limit, collision self-energies with Ward-consistent vertices, thermodynamic BZ convergence, actual source Fig./Table I/Eq.(19) SI observable. A gauge Ward necessary condition is not a sufficient test of the full conductivity. **Original Cui reproduction: PARTIAL.**

## Reproduction

~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/audit_finite_q_gauge_ward.py
~~~

Two-version GitHub CI acceptance applies **only after both actual run logs show 243/243 tests and Round27 audit success**. Earlier 227 Cui and repo-wide tests are not conflated.
