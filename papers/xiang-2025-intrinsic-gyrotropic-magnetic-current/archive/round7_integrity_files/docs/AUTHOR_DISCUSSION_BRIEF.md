# Potential author clarification — draft for researcher review ONLY

Paper: Xiang et al., *Phys. Rev. Lett.* 134, 116301 (2025). This is a **draft scientific discussion brief**, not a sent letter. No contact has been initiated.

## Priority 1: normalization used to produce published Fig. 2(d)

The article defines `Q_nm^{αβ}=(r_nm^α σ_mn^β+r_mn^α σ_nm^β)/2=Re(r_nm^α σ_mn^β)` in main Eq. (5). In SM §III, the intermediate line after Eq. (34) instead uses `Q_nm^{αβ}=2 Re(r_nm^α σ_mn^β)` while retaining the same prefactor of its preceding oscillatory-wavepacket expression. Under `H1=-g μ_B B σ`, main Eq. (12) gives approximately `0.017754 A/(m T)` for the in-gap Fig. 2(d) case (`g=2`, `f=1 THz`, all other model parameters from the caption). The *literal intermediate* would double this to approximately `0.035508 A/(m T)`, closer to a **rough visual read** of Fig. 2(d), but it is inconsistent with independently computed finite-field Wilson polarization, orbital-basis Liouville and an eigenvector-free Pauli-trace derivation.

**Suggested focused question:** Which normalization of `Q_nm` and the Zeeman term was used in the original code generating Fig. 2(d), and would the authors be willing to share its raw plotted values (or data/parameter script) to distinguish a source convention from a plot normalization?

## Priority 2: antiunitary sign in Eq. (14)

Main Table I and SM §V distinguish T-even Zeeman curvature / C and T-odd Zeeman metric / D. In our independent fixed real-phase computation, H and its TR magnetic partner give `K_T(ω)=K_H(ω)*`: the real conduction-like quadrature is even and the imaginary displacement-like quadrature is odd, for both 1 and 3 THz. However Eq. (14) prints one factor `η_T=-1` for both C/D.

**Suggested focused question:** Is Eq. (14)'s antiunitary sign intended for the displacement sector only, or does it implicitly act on a complex frequency-dependent tensor with an additional conjugation or ω-sign convention?

## Lower priority / do not mix with the two main questions

SM Table II row `-42'm'`: its displayed unconstrained `σ_xy`, `σ_yx` do not exhibit the equality required by `\bar4`, though the allowed/nonallowed classification is unchanged. See `docs/ROUND5_SM_TABLE_SOURCE_AUDIT.md`.

## What is not asserted

No claim is made that the authors used the factor-two expression in their code, that Fig. 2(d) contains a proved error, or that a real frequency-domain antiunitary response can be treated without causality/conjugation. No email is sent and no publisher PDF is redistributed.
