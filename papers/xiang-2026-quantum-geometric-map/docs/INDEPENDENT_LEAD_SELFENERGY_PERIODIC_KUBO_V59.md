# v59 | Finite-band causal electrode self-energies and independent periodic-bulk Kubo benchmark

**Physics-first aim:** replace v57–58's ideal flat wide-band Γ by an *exactly solved microscopic lead model*, and independently validate the *periodic-bulk charge Hall* geometric/Kubo interface. These are TWO DISTINCT independently tested capabilities. **Do not interpret a four-terminal finite-device voltage as bulk orbital MNHE/PHE**, or infer the original Xiang paper's contested formulas from a QWZ Chern insulator.

## A. Exact semi-infinite 1D electrode

An actual semi-infinite uniform one-dimensional tight-binding electrode has on-site energy ε_l, hopping magnitude t_l>0 and device-to-lead hopping w_l. Its endpoint retarded surface Green function satisfies

\[
g_{11}^r(E)=\big[E-\varepsilon_l-t_l^2g_{11}^r(E)\big]^{-1},\quad
\Sigma_l^r(E)=|w_l|^2g_{11}^r(E),\quad
\Gamma_l(E)=-2\,\mathrm{Im}\,\Sigma_l^r(E)\ge0.
\]

The physical *retarded* branch, x=E−ε_l, is

\[
g_{11}^r=\frac{x-i\sqrt{4t_l^2-x^2}}{2t_l^2},\quad |x|<2t_l;
\]
outside the band it is the real root behaving as 1/x, and Γ is exactly zero. This explicitly retains **both** the energy-dependent real energy shift ReΣ and the finite-band broadening Γ(E). It satisfies causal Im g≤0 and the exact surface Dyson equation, including right at both band edges (where the energy derivative diverges and is separately rejected).

For an explicitly user-specified *field-dependent lead on-site energy* ε_l(B)=ε_l(0)+s_l B (optional modeling input), **both** Σ'_B and Γ'_B must appear in the B derivative:

\[
\partial_Bg=-s_l\,\partial_Eg,\qquad
\partial_Eg=-g^2/(1-t_l^2 g^2).
\]

This optional linear *scalar magnetogating* is a mathematical sensitivity model. A spinless ordinary scalar on-site energy linear in B generally does NOT obey microscopic time-reversal/Onsager by itself. Spinless Onsager is therefore asserted/tested only for B-independent electrode parameters. **It is not claimed to describe a fully gauge-covariant microscopic orbital lead in a real perpendicular field.**

## B. Complete Caroli and magnetic-contact derivative

With diagonal local, site-resolved Σ_l(E,B), an actual finite v55 Peierls device H(B) has

\[
G^r(E,B)=[EI-H(B)-\sum_l\Sigma_l^r(E,B)]^{-1}
\]
and true lead-resolved Caroli transmission

\[
T_{ab}(E,B)=\mathrm{Tr}[\Gamma_a G^r \Gamma_b G^a]
=\Gamma_a(E,B)\Gamma_b(E,B)|G^r_{s_a s_b}|^2\quad(a\ne b).
\]

The full derivative must retain all reservoir contributions:

\[
\partial_BG^r=G^r\left[\partial_BH+\partial_B\Sigma^r\right]G^r
\]
\[
\partial_BT_{ab}=
(\Gamma'_a\Gamma_b+\Gamma_a\Gamma'_b)|G_{ab}|^2+
2\Gamma_a\Gamma_b\,\mathrm{Re}(G_{ab}^*\partial_BG_{ab}).
\]

Thermal integration is adaptive over -f'(E), **split at all actual lead band edges and finite Hamiltonian eigenenergies**. The routine fails closed if its numerical quadrature does not converge or if lead current conservation is violated. It uses v57's exact four-terminal floating transverse probes, and an *independent five-point finite-B rerun* of the entire Green function/lead self-energy/voltage-solve to verify the derivative. Both magnetic gauges are compared at the level of **individual T_ab and d_B T_ab**. Real energy-dependent leads obey the exact retarded/advanced optical theorem \(i(G^r-G^a)=G^r\Gamma G^a\), and lead currents sum to zero.

### Original finite-device reference

For the same v57 four-site ring at B=.31, q/ħ=-.87, μ=.15, T=.16, with *four identical physical semi-infinite chain leads* t_l=2, w_l=sqrt(.55):

| Hall-device observable | Exact 1D finite-band leads | Prior v57 wide-band Γ=.55 |
| --- | ---: | ---: |
| Floating Hall voltage V_H/ΔV | **−0.0661265705180** | about −0.0664226587 |
| Analytic d_B(V_H/ΔV) | **−0.1555717994672** | about −0.1619849894 |
| Independent full 5-point B derivative error | **3.89×10^-13** | previously tested |
| Symmetric vs Landau T and dT disagreement | numerical zero | previously tested |

For t_l=2,5,20,100 with **w_l²/t_l=.55/2 fixed**, V_H/ΔV runs −.06612657, −.06638331, −.06642234, −.06642310, approaching the original wide-band result; the corresponding B slopes approach their wide-band values as well. This is a **nontrivial independent continuity control**: physical finite-band reservoirs reduce to the established ideal lead model in the explicitly defined t→∞ scaling.

The optional s_l=.12 magnetic on-site shift gives a different physical *prescribed* reservoir model and a different B slope. The full Σ'_B and Γ'_B analytic calculation agrees with full five-point B differences; we explicitly DO NOT impose Onsager for this scalar B-odd reservoir shift.

## C. An independent periodic bulk *charge* Hall Kubo oracle

The independent 2D Qi–Wu–Zhang two-band periodic lattice Hamiltonian

\[
H(\mathbf k)=\sin k_x\sigma_x+\sin k_y\sigma_y+
[m+\cos k_x+\cos k_y]\sigma_z
\]

has analytic velocity operators \(d_iH\). Its band Berry curvature integrated over the Brillouin zone by the v54 spectral geometry code gives a **periodic intrinsic charge anomalous Hall (AHE)** Chern number:

\[
C=\frac{1}{2\pi}\int_{\mathrm{BZ}}\Omega_{xy}^{\mathrm{occupied}}(\mathbf k)\,d^2k.
\]

An INDEPENDENT Fukui–Hatsugai–Suzuki **discrete occupied eigenstate link/Wilson plaquette method** yields the same topological invariant. The FHS link phase has the opposite sign to A=i⟨u|du⟩ by construction, and the code documents this sign convention explicitly; it is not an empirical sign choice to make a plot agree.

| Mass m | Spectral integrated C (31×31) | Independent Wilson C | Insulating phase |
| ---: | ---: | ---: | --- |
| −3 | +6.88×10^-10 | 0 | trivial |
| −1 | +0.9999999981 | +1 | nontrivial |
| +1 | −1.0000000019 | −1 | nontrivial |
| +3 | +6.88×10^-10 | 0 | trivial |

The true gap-closing transition masses −2,0,+2 fail closed: a Chern integer of a **gapped** occupied bundle is not assigned there. The periodic analytic Bloch derivatives agree with independent central differences. This is a **2D zero-external-field intrinsic CHARGE anomalous Hall bulk benchmark**, NOT d_B σ(E) and NOT the target orbital intrinsic magnetoelectric Hall of the original paper. The device Hall voltage and periodic Chern conductivity have different boundary/observable definitions and are **never compared by numerical equality**.

## D. Physics acceptance and remaining work

The 21 new independent tests cover lead causal branch/Dyson equation, exact spectral function/band edges and optical theorem, lead-E and lead-B derivatives, full Caroli trace, full five-point B differences, gauge and conservation, physically compatible Onsager, recovery of v57 wide-band contacts; independently the full periodic QWZ Kubo Berry integral, discrete FHS/Wilson invariant, k-grid refinement, exact analytic Bloch velocities and gapless-phase rejection.

| Capability | Status |
| --- | --- |
| Exact retarded 1D semi-infinite chain Green function, Σ(E), Γ(E), d_BΣ and d_BΓ | **VERIFIED for the specified lead Hamiltonian** |
| Finite open-device coherent Caroli four-probe response with specified physical dispersive chains | **VERIFIED numerical/model algorithm; CONDITIONAL finite-device prediction** |
| Independent periodic QWZ intrinsic **charge** AHE Chern benchmark | **VERIFIED bounded bulk oracle** |
| Arbitrary realistic Wannier principal-layer electrode, explicit disorder/dephasing, conserving microscopic bath | **UNRESOLVED** |
| Intrinsic bulk orbital mixed E/B MNHE/PHE of the Xiang paper | **UNRESOLVED**; original source audit **PARTIAL** |

### Reproduce (isolated PRB subtree)

    python -m unittest tests.test_lead_selfenergy_v59 -v
    python scripts/reproduce_lead_bulk_v59.py --output results/frozen/v59_lead_and_bulk_benchmark.json

Do not claim that a 1D local endpoint lead self-energy is a general material electrode. A true Wannier principal-layer surface matrix self-energy, electrode orbital embedding, scattering channels and disorder/dephasing still need separate input and verification.

**Suggested v60:** independently test general principal-layer matrix surface Green functions and field-dependent \(\Sigma(E,B)\) on genuinely 2D devices, then develop a **same-Hamiltonian** finite-open vs periodic Kubo comparison under controlled fixed carrier density, before targeting full quantum-geometric orbital magnetic Hall.
