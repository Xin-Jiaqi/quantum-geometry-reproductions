# Stage 7 continuation: a controlled four-band downfolding gate

**Status:** new analytic identities + independently checked spectrum/projector diagnostics, **not** a completed two-/four-band NMEE response formula. The existing Stage 7 eight-band K/K' expansion, near-complete *signed* valley capture, and orbital interference tests remain frozen in STAGE7_PHYSICS.md.

**Physics scope:** synthetic eight-band bilayer Kane–Mele model; zero-temperature, gapped equilibrium/static finite-Zeeman NMEE only. Neither the driven-DC formula nor ferroelectric double-well nor any real material is established.

## 1. Why an accurate eight-band valley expansion does not validate a naive four-band model

At each exact microscopic K/K' center, choose P0 as the middle four eigenstates of the **decoupled** (t_perp=0) bilayer, and Q0 as the other four. At fixed u and t_perp, partition the original eight-band Hamiltonian:

\[
H(\mathbf k)=\begin{pmatrix}A&B^\dagger\\ B&D\end{pmatrix},\quad
A=P_0HP_0,\ B=Q_0HP_0,\ D=Q_0HQ_0.
\]

The diagnostic ratio is r=||B||_2/Δ, Δ=min{|a-d|: a in spec(A), d in spec(D)}. It is **not a rigorous SW error bound**; no claim is made when the block spectra interleave. We also compare the exact *middle four-band* projector to the original P0 through the smallest principal-angle cosine squared.

Independent numerical calculation using the repository's original microscopic formulas (bond=0.142 nm, t=0.85 eV, intrinsic SOC 10 meV, Rashba 20 meV, staggered mass 10 meV, interlayer range 0.16 nm, lattice shell 5, chemical potential 0.02541381265 eV):

| u | t_perp | qx (nm^-1) | ||QHP|| (meV) | Δ (meV) | r | smallest principal cos² |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 5 meV | 0 | 3.536 | 12.875 | 0.275 | 0.9991 |
| 0.1 | 5 meV | 0.10 | 13.823 | 9.309 | 1.485 | 0.8927 |
| 0.1 | 5 meV | 0.75 | 93.319 | 16.001 | 5.832 | 0.2536 |
| 0.1 | 20 meV | 0 | 14.143 | 8.514 | 1.661 | 0.0000057 |

At 20 meV, even at q=0, the four lowest/nearest eigenstates have changed their orbital-space content enough that the **fixed t_perp=0 four-state projector** is no longer a valid small-mixing starting point. At 5 meV the same fixed projector also becomes unreliable on a finite 0.75 nm^-1 disk. **This does not rule out other k-dependent, symmetry-adapted four-band projections.** The sample is only q=(qx,0), not an upper bound throughout the K/K' disks. Mirror u↔−u and TR matrix checks of the reconstructed model reach ≤1.3e-17 eV and 0, respectively.

The exact Feshbach operator is H_eff(E)=A+B†(E−D)^(-1)B when E lies outside spec(D). It is energy dependent; the transformed spin, velocity, and Berry/position operators must also be kept. Merely diagonalizing the 4×4 block A or differentiating a truncated H_eff is **not** a validated response calculation.

## 2. Basis-consistent interference and projector formula

At fixed declared A/B orbital embedding, let

\[
H(\mathbf k;s,c,B_z,\lambda)=H_0+sV_s+cV_c-g\mu_B B_z S_z+\lambda S_i.
\]

The current frozen code computes

\[
\chi_{ij}(B_z)=-2\mathrm{Im}\sum_{v\in occ,c\in emp}
\frac{(S_i)_{vc}\,(\partial_{k_j}H)_{cv}}{(\epsilon_v-\epsilon_c)^2},\qquad
\Gamma_{iz;j}=\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}\partial_{B_z}\chi_{ij}.
\]

For the gapped occupied projector P, the *same convention* may be written without denominators within nearly degenerate occupied or empty spaces:

\[
\chi_{ij}=i\operatorname{Tr}\{P[\partial_\lambda P,\partial_{k_j}P]\}_{\lambda=0},\quad
C_{sc}^{ij}=\left.\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}
\partial_{B_z}\partial_s\partial_c
i\operatorname{Tr}\{P[\partial_\lambda P,\partial_{k_j}P]\}
\right|_{s=c=B_z=\lambda=0}.
\]

With R(z)=(z−H)^(-1), P=(2πi)^(-1)∮R(z)dz and
∂s∂c P=(2πi)^(-1)∮(RVsRVcR+RVcRVsR)dz, which is an exact, gauge-subspace-invariant differentiation algorithm as long as the occupied/empty gap stays open. The k derivative also acts on V_s,V_c through the original Bloch phases: these velocity terms must **not** be discarded. This is an operator expression, not yet a compact analytic coefficient for a 2-/4-band model.

**Exact total hopping parity:** the layer sign matrix L=τ_z⊗I4 maps H(s,c) to H(−s,−c), also mapping velocities consistently and leaving the spin operators invariant. Thus Γ(s,c)=Γ(−s,−c). Since Γ(0,0)=0 (two mirror-related disconnected layers cancel), analytic weak hopping gives

\[
\Gamma(s,c)=A s^2+B c^2+C sc+O((|s|+|c|)^4).
\]

All odd **total** hopping orders are absent. This imposes no universal sign on C, nor does C require [V_s,V_c] ≠ 0. The published Stage 7 cross-sign-odd finite-t response includes higher-order odd-in-c contributions; its division by t² approaches C only in the weak-t limit.

## 3. Reproduction and next decision

New study-only diagnostic: bsf_nmee/sw_projection_diagnostics.py (returns both valleys × u={0.1,0.4} × t_perp={5,20} meV × qx={0,0.1,0.25,0.5,0.75} nm^-1 = 40 sampled rows). Unit tests: tests/test_stage7_sw_projection_gate.py. To run from studies/bsf-nmee-connection: python -m unittest discover -s tests -v.

**Next step:** search for a k-dependent symmetry-adapted low-energy projector, explicitly bound/eliminate remote-band mixing, and separately benchmark the *transformed observables* against the eight-band static response. If no controlled four-band window covers the response-dominant valley region, retain the accurate eight-band theory. No DFT and no claims of a universal PRL-level mechanism.
