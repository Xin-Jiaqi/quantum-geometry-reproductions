# Hamiltonian and validation contract

## 1. Eight-state basis and microscopic motif

Orbital basis `(top A↑, top A↓, top B↑, top B↓, bottom A↑, bottom A↓, bottom B↑, bottom B↓)`. The spatial motif has **two distinct orbital centers in every layer**:

\[
\mathbf r_A=(0,0),\qquad \mathbf r_B=(0,0.37)a.
\]

Bottom / top are separated by an out-of-plane distance not parameterized for the in-plane hopping; the in-plane relative translation is \(u\hat y\). In a cell-periodic orbital gauge,

\[
T_{ab}(k,u) = t_0\sum_{R_x,R_y=-2}^{2}
 e^{-\kappa[(R_x+x_b-x_a)^2+(R_y+y_b-y_a+u)^2]}
 e^{i(k_xR_x+k_yR_y)}.
\]

All interlayer orbital, interlayer spin-flip SOC-like, and spin-dependent tunneling blocks use the same **position-derived** \(T_{ab}(k,u)\), with separately frozen scalar couplings. This is a two-center overlap **ansatz**, not a Slater–Koster fit to an actual element or DFT hopping calculation.

The Hamiltonian has the block structure

\[
H=\begin{pmatrix}h_{\rm top}&V(k,u)\\V^\dagger(k,u)&h_{\rm bottom}\end{pmatrix},\quad
V_{ab}=T_{ab}(k,u)\,[s_0-i\lambda_{\rm inter} s_x+\xi_{\rm AF}s_z],
\]

where \(\xi_{\rm AF}=0\) in the FM-like branch. On each layer the orbital block is \(m\rho_z + 2t_{AB}\cos k_y\,\rho_x\), identical top/bottom. (This *intralayer* mixing is a tight-binding ansatz, not recalculated from the atomic motif.) The remaining blocks are spin-orbit and magnetic terms obtained by lifting the previous four-band allowed-operator model to the orbital identity:

\[
\begin{aligned}
H_{\rm so}&=\lambda\,\tau_z(\sin k_y s_x-\sin k_x s_y)+\lambda_z\sin k_x s_z,\\
H_\mathrm{AF}&=J\tau_z s_y+c s_z,\\
H_\mathrm{FM}&=J s_y + c\tau_z s_z.
\end{aligned}
\]

The orbital, layer and spin identity factors have been suppressed in these expressions. **There is no \(\delta(u)\tau_z\) or any other separately imposed odd-in-u onsite potential.** This is enforced in the public API. Orbital energy splitting \(m\rho_z\) is layer independent and therefore does not represent a polar field.

## 2. Symmetry representation (verified numerically)

Take \(U_{m_z}=\tau_xs_z\), \(U_{m'_z}=U_{m_z}(is_y)K\), and internal \(U_{m'_x}=s_zK\), with orbital identity. We verify

\[
H(k,-u)=U_{m_z}H(k,u)U_{m_z}^\dagger \quad (\mathrm{AF}),
\]
\[
H(-k,-u)=U_{m'_z}H(k,u)^*U_{m'_z}^\dagger \quad (\mathrm{FM}),
\]
\[
H(k_x,-k_y,u)=s_z H(k_x,k_y,u)^*s_z\quad (\mathrm{both}).
\]

The `K` is implicit when equations explicitly show the complex conjugation. Expected combined time-reversal phase convention of the low-frequency displacement channel is the independently tested convention from the IGMC paper reproduction; do **not** interpret it as a statement about an arbitrary complex causal response at one frequency. Connection alone does not establish that two candidate structures are energetically stable.

## 3. Geometric coefficient and polarization checks

Define zero-T occupied/empty pairs \(n<m\), velocity \(v_x=\partial_{k_x}H\), physical Pauli spin \(s_x\). The declared reduced coefficient is

\[
\frac{D_{xx}}\omega = \Bigl\langle\sum_{n<m}(f_n-f_m)
\frac{-2\operatorname{Im}[v^x_{nm}s^x_{mn}]}{(\epsilon_m-\epsilon_n)^2}\Bigr\rangle_{\rm BZ}.
\]

The independence check diagonalizes \(H\mp b s_x\) at each k point, forms the **occupied four-band non-Abelian** Wilson determinant along x and differentiates the Berry polarization wrt b. The reference mesh convergence verifies that the value approaches the direct eigenspace-pair expression. There is no reused call to the direct formula inside the Wilson routine. Orbital x positions are equal, so the stated x-direction cell-periodic convention does not omit a nonzero intracell x-coordinate correction in this specific model; this simplification is **not** transferable to arbitrary Wannier centers.

The layer-charge proxy \(\langle\tau_z\rangle\) is the total occupied density difference between layers. It is odd under the two connectors, even though displacement current is odd in the unitary branch and even in the antiunitary branch. It is **not** the complete Berry-phase switchable polarization of a real ionic crystal.

We require globally positive **indirect band gaps**, and four occupied bands everywhere at each chosen fixed chemical potential. In the present parameter set, \(\mu_\mathrm{AF}=0\) and \(\mu_\mathrm{FM}=-0.03\) eV (same chemical potential in each branch's ±u pair). The values are taken inside the corresponding half-filling gap, not optimized to fit response curves.

## 4. Limitations / falsification gates

- The electronic **band-energy** sum is not the crystal's total energy; neither \(u=+0.19\) nor \(u=-0.19\) has been shown to be an energy minimum. No genuine BSF has yet been predicted.
- Intralayer \(2t_{AB}\cos k_y\rho_x\) is phenomenological, as are SOC/exchange/tunnel prefactors. Only **the u-dependence of all interlayer hoppings** derives from specified orbital-center displacements.
- Spin-dependent tunneling in the AF-like branch is crucial to the magnitude; physical origin and relation to actual magnetic orbitals require DFT.
- The full orbital magnetic moment and Maxwell–Faraday-induced E fields are not included. This work probes the **spin-Zeeman** IGMC model, not the complete electromagnetic transport current.
- Bands are insulating over the sampled u interval but k meshes are finite; no topological gap guarantee away from mesh samples is claimed. The energy and charge proxies are not experimental observables.
- Original Xiang PRL Fig.2(d) discrepancy and antiunitary convention are still unresolved independently; this study does not retrospectively alter that paper's `PARTIAL` status.
