# v57 | Independent open-lead four-terminal Hall device (NEGF and Landauer-Büttiker)

**Scientific goal:** A real finite-device current with explicitly connected electrodes, not the scalar RTA bath of v56 and not fitted to disputed original-paper SM expressions. **This is not a verified intrinsic bulk orbital MNHE/PHE.**

## Physical assumptions and response

Reuse v55's Peierls finite open graph H(B) and its analytical magnetic vertex M=d_B H. Specify 2 or 4 **distinct local ideal wide-band terminals** at graph sites with B-independent energy-independent diagonal broadenings Γ_a>0. The retarded lead self-energy is Σ^r=-iΓ/2, Γ=sum Γ_a; the retarded Green function is

\[
G^r(\varepsilon,B)=[\varepsilon I-H(B)+i\Gamma/2]^{-1},\quad
T_{\alpha\beta}=\mathrm{Tr}[\Gamma_\alpha G^r\Gamma_\beta G^a].
\]

At finite T integrate each T_ab with the actual Fermi-window derivative −f'(ε−μ). Results are in **q²/h conductance units**; model magnetic field and graph area define the dimensionless flux qB area/ħ. Local B-independent Γ permits **analytic magnetic differentiation**:

\[
\partial_B G^r=G^r M G^r,\quad
\partial_B T_{\alpha\beta}
=2\,\mathrm{Re}\,\mathrm{Tr}[\Gamma_\alpha G^r M G^r\Gamma_\beta G^a].
\]

No spurious scalar scattering rate is added to the central-device Green function: coupling to ideal leads supplies its linewidth. This remains a **wideband coherent** lead approximation, not a full microscopic lead or inelastic bath. B-dependent or nonlocal general lead self-energies require extra B derivatives and contact/gauge terms, and are **not supported**.

For a four-terminal measurement, ordered Left,Top,Right,Bottom, fix V_L=+.5,V_R=−.5 and solve I_T=I_B=0 using the Landauer–Büttiker current Laplacian

\[
I_\alpha/(q^2/h)=\sum_{\beta\ne\alpha}T_{\alpha\beta}(V_\alpha-V_\beta).
\]

Differentiate **the complete floating-probe linear system** in B to obtain d_B(V_T−V_B), not merely a selected transmission derivative. Independently recompute *all* transmissions and *all* floating voltages at five nearby finite magnetic fields.

## Exact numerical benchmark

Symmetric four-site diamond ring, q/ħ=−.87,B=.31, μ=.15,T=.16, each Γ_a=.55; all numbers model-unit normalized:

| Physical 4-lead output | Value |
| --- | ---: |
| Hall ratio V_H/ΔV | **−0.06642265873** |
| Magnetic Hall slope d_B(V_H/ΔV) | **−0.16198498942** |
| Longitudinal effective device conductance in q²/h | **0.51928029946** |
| Floating and total current conservation error | order 10⁻¹⁶ |
| Independent analytic-vs-five-point field derivative | verified to numerical differentiation accuracy |

## Twenty independent physical gates

Original Caroli matrix trace independently matches the local lead matrix-element formula γ_aγ_b|G_ab|². The retarded/advanced **optical theorem** i(G^r−G^a)=G^r Γ G^a implies exact equality of incoming/outgoing transmission totals; all lead currents and their B derivatives sum to zero. Every T_ab and d_BT_ab agrees individually in symmetric vs Landau magnetic gauges, as do the floating-probe outputs. Complete analytic magnetic derivatives agree with five-point finite-field differences of the *full* Landauer–Büttiker calculation.

Onsager–Casimir T_ab(B)=T_ba(−B) passes. **Two-terminal** conductance is even in B and has zero B-linear slope at B=0. The mirror-symmetric **four-terminal** ring has an odd-in-B Hall voltage, zero at B=0. Temperature integration, origin translations, hermiticity and contact validity are tested and fail closed where appropriate. **20/20 new physical tests passed locally.**

## Lead and geometry dependence are not errors to fit away

The Hall slope for equal Γ=.2,.55,1.1 is respectively −.08773,−.16198,−.17398. The Hall voltage ratio at T=.04,.16,.35 is respectively −.11990,−.06642,−.02605. These are physical device/contact-dependent values, not a single intrinsic material coefficient.

A small **sparse concentric-ring** sequence with 1/2/3/4 rings gives B slopes approximately −.162,−.682,−6.365,+6.418, showing strong mesoscopic dependence. **THIS IS NOT A THERMODYNAMIC-SCALING STUDY:** area, effective site density, flux and lead geometry change together, and this family is not a fixed-density planar lattice. No claim of bulk convergence or material size scaling is supported.

## Formal acceptance

- **VERIFIED physics algorithm:** two-/four-terminal coherent wideband local-lead Caroli transmission, analytic B derivative, Landauer–Büttiker floating voltage solution, optical theorem/current conservation, magnetic gauge invariance, Onsager and finite-field oracle.
- **CONDITIONAL device prediction:** finite coherent Hall voltage and its B derivative for particular specified leads and finite graph.
- **UNRESOLVED:** controlled thermodynamic bulk limit, real energy-/field-dependent lead self-energies, realistic materials embedding, conserving inelastic scattering, intrinsic bulk MNHE/PHE and full two-frequency nonlinear EB.

Unlike v56, v57 has explicit source/drain leads and voltage probes. v56 AC scalar-relaxation *closed-graph* numbers and v57 DC Landauer *open-device* numbers refer to **different observables**, so no artificial numerical matching is attempted.

Run in the isolated paper package:

    python -m unittest tests.test_lead_transport_v57 -v
    python scripts/reproduce_lead_transport_v57.py
    python scripts/reproduce_lead_transport_v57.py --graph-npz examples/finite_four_probe_ring.npz
    python scripts/reproduce_lead_transport_v57.py --rings 2 --gamma-lead .2

External input NPZ keys: hoppings (N,N Hermitian), positions (N,2 real), lead_sites (2 or 4 distinct integer indices); loading disables pickle. Model-unit results are not experimental SI predictions. Detailed full-fidelity local handoff includes NPZ and frozen original solver output; paper repo freezes the scientific benchmark and parameter dependence. Whole original Xiang paper source-audit state remains PARTIAL.
