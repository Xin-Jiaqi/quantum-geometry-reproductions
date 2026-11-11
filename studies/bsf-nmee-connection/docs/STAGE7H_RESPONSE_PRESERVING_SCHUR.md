# Stage 7h: response-preserving Feshbach elimination and remote-sector contacts

**Scientific status:** An **exact eight-band algebraic reformulation** plus stringent local numerical verification; **not** yet a closed or error-controlled four-band NMEE theory. Scope is the established synthetic gapped T=0 **equilibrium static** Zeeman-derivative NMEE, not the unresolved driven-DC PRL, DFT, or ferroelectric energetics.

## Physics motivation

Stage7f–g showed that the positive scalar AB/BA hopping coefficient depends on large compensating low/remote band contributions, and a bare middle-four-band Kubo sum can have the wrong sign. Now test what is truly required after eliminating the complementary Q space.

Fix at each k a reference P0 of four middle bands of the **decoupled** bilayer; Q0 denotes the complementary four. This splitting need **not** be perturbative; it is an exact coordinate partition. In P/Q blocks write

\[
H=\begin{pmatrix}A&B^\dagger\\B&D\end{pmatrix},\
r(z)=(z-D)^{-1},\
\Sigma(z)=B^\dagger r(z)B,\
g(z)=[z-A-\Sigma(z)]^{-1}.
\]
With \(\mathcal U=(I,rB)^T\), \(\mathcal V=(I,B^\dagger r)\), and \(\mathcal R=\mathrm{diag}(0,r)\), the **full eight-band** resolvent is exactly

\[
\boxed{G(z)=(z-H)^{-1}=\mathcal U g\mathcal V+\mathcal R.}\tag{1}
\]

For any physical source vertex \(O\), especially \(S_x,v_x=\partial_{k_x}H,S_z\), define the energy-dependent *four-dimensional* dressed vertex

\[
\widetilde O(z)=\mathcal V O\mathcal U
=O_{PP}+O_{PQ}rB+B^\dagger rO_{QP}+B^\dagger rO_{QQ}rB.
\]
It equals \(\partial_\lambda[A+\Sigma(z)]\) for a source \(H+\lambda O\) at fixed k and reference P/Q frame. The **source derivative** is NOT just \(\mathcal Ug\widetilde O g\mathcal V\). There is a rigorously required decomposition

\[
\boxed{\begin{aligned}
GOG ={}& \mathcal Ug\widetilde O g\mathcal V\\
&+\mathcal Ug\mathcal V O\mathcal R
 +\mathcal R O\mathcal Ug\mathcal V
 +\mathcal R O\mathcal R.
\end{aligned}}\tag{2}
\]

The last three terms are called *contact/remote-source pieces* in this study, **purely as an algebraic grouping**. Omitting them is an uncontrolled source truncation, *not* an independently well-defined or physical four-band theory.

The zero-T insulating occupied projector and its derivatives are recovered by contours enclosing all occupied poles (the numerical validation uses two separate isolated occupied doublet contours):
\[
P=\oint_{\mathcal C}\frac{dz}{2\pi i}G,\quad
\partial_O P=\oint_{\mathcal C}\frac{dz}{2\pi i}GOG,\quad
\chi_{x;x}=i\operatorname{Tr}\left(P[\partial_{S_x}P,\partial_{v_x}P]\right),
\quad \Gamma_{xz;x}=\partial_{B_z}\chi_{x;x}.
\]
The physical \(v_x\) and spin matrices are calculated in the **original Bloch embedding** first, then transformed at fixed k; the k-dependent P/Q basis is **not** used to differentiate a naive \(4\times4\) Hamiltonian.

## Local sign-reversal counterexample

At the first K center, slide u=0.4, t_perp=20 meV, finite Zeeman step ±0.25 T, values of the **local** \(\partial_B\chi\) (NOT BZ integration, original program units):

| Evaluation | Local response |
|---|---:|
| Original full eight-band Kubo | **-0.028550806875** |
| Exact Feshbach all four source terms | **-0.028550806874** |
| Only dressed-vertex source terms | **+0.036199325084** |
| Omitted remote/contact corrections | **-0.064750131958** |

At u=.4, t=5 meV the exact result is +0.005652317063, dressed only +0.008359236405, and omitted correction -0.002706919340.

**Hence, the four-dimensional effective propagation \(g\) and dressed source vertex are not sufficient by themselves to reconstruct the full static NMEE; the omitted three terms can even reverse the sign of a local diagnostic.**

For two **small** K/Kprime radius-0.15 nm^-1 disks, u=.4, t=20 meV, radial×angular=12×36/valley, the signed response is full +1.02868162114e-4, dressed only +1.95265130628e-4 and omitted correction -9.23969685144e-5. The three omitted terms reduce the approximate magnitude by ~90% of the actual result. The disk is **not** the entire valley or BZ; at u=.1/t=20 meV the analogous correction is positive +1.70378e-6. No universal contact sign is claimed.

## Exact source metric and energy-separated vertices

A physical eigenvector \(\psi=(p,q)\) at an energy \(E\) off the poles of D obeys \(q=r(E)Bp\). Full normalization induces the **energy-dependent P-space metric**
\[
\boxed{p^\dagger Z(E)p=1,\quad Z(E)=I+B^\dagger(E-D)^{-2}B=I-\partial_E\Sigma(E).}\tag{3}
\]
For two exact eigenstates with generally distinct energies,
\[
\widetilde O(E_m,E_n)=O_{PP}+O_{PQ}r(E_n)B+
B^\dagger r(E_m)O_{QP}+B^\dagger r(E_m)O_{QQ}r(E_n)B,
\quad \langle\psi_m|O|\psi_n\rangle=p_m^\dagger\widetilde O(E_m,E_n)p_n.
\]

At K+(0.05,0.02) nm^-1, u=.4, t=20 meV, comparison of the middle four **exact** eigenstates yields bare P-block maximum source matrix errors 0.557 for Sx, 0.891 for Sz and 0.0859 eV·nm for vx, versus exact dressed errors of order 1e-13. With t=5 meV the analogous bare errors are Sx 0.0311, Sz 0.0454 and vx 0.00234 eV·nm. Z-normalization and full Q-component reconstruction are checked numerically. This **does not** restore excluded Q→Q or Q→P Kubo virtual transitions if one only sums four low eigenstates.

## Reproducibility and limits

- Independent original-equation eight-band implementation and repository-native module give exact resolvent/source identities, tested at several complex energies, with direct matrix-inverse errors <~1e-12.
- Two occupied doublet contour residues reproduce the full occupied projector at ~1e-15 and the local eight-band Kubo static derivative to ~1e-11 absolute program units.
- Contour 24/48/96 points and B steps 0.125/0.25/0.5 T preserve the contact sign effect; standalone two-valley small-radius integration independently reproduces full eight-band direct integration.
- **Crucial:** This method currently uses original eight-band eigenvalues to position numerical contours **for validation**, and represents the contacts in 8D embeddings. It is NOT a computationally independent efficient four-band approach, not an error-controlled SW projection, and not a full BZ calculation.
- Fixed reference P/Q source partition and its dressed/contact components are **not independently observable or k-dependent-gauge invariant**; only the sum and original response are physically meaningful.

Next single physics step: exploit valley C3/spin/orbital symmetry to isolate the smallest **matrix-element phase combinations** inside the energy-dependent spin, Zeeman and velocity vertices, then test whether those govern scalar-cubic versus nematic-linear terms and their signs. Do not claim new universal transport effect from the standard Schur identity alone.

**Added study-only source:** bsf_nmee/feshbach_response.py, tests/test_stage7h_feshbach_response.py, scripts/run_stage7h_feshbach.py, results/sliding_stage7h_feshbach.json (benchmark numbers). Complete independent executable archive and the Chinese physics report are separate. Frozen Jia paper and earlier Stage7 research files have not been modified.
