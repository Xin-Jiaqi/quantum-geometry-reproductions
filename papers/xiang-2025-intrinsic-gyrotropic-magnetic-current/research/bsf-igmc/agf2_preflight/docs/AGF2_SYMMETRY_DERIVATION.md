# AgF₂ Type-I / Type-II crossed-parity IGMC proposal

## Conventions

Cartesian laboratory axes as used for the magnetic order of the Zhu et al. AgF₂ example, `x=[100]`, `y=[010]`, `z=[001]`. The response is a *polar* current to an *axial* magnetic drive, `J_i(ω)=χ_{ij}(ω)B_j(ω)`.

For a connector `(R, θ)` connecting **two distinct ordered domains/stackings** (not a fixed-frequency complex conjugation identity), the **declared real quadrature** convention of the existing Xiang model reproduction is

`χ_X(B2) = η_X^θ det(R) R χ_X(B1) R^T`, `η_C=+1, η_D=−1`.

This channel parity is part of the present response convention, not a resolution of the published Eq. (14) antiunitary-convention question. Separate physical tensor allowedness from the partner-to-partner sign.

## Internal single-state magnetic group `m'_100`

The published AgF₂ Type-I SOC magnetic symmetry contains the antiunitary mirror `m'_x`, with `R=diag(−1,1,1)`, `det R=−1`, and `θ=1`. For C, invariance requires `χ^C=−R χ^C R^T`, giving **four allowed** basis entries `xy, xz, yx, zx`. For D it requires `χ^D=+R χ^D R^T`, giving **five allowed** entries `xx, yy, yz, zy, zz`.

Thus `C_xy` and `D_xx` are both **symmetry-allowed in this magnetic group**. At zero temperature, C is a Fermi-surface response and vanishes if the material is insulating; D is an insulating-eligible displacement response. Neither allowedness statement establishes a nonzero coefficient for AgF₂.

Do not assume without checking that the *strained Type-II* configuration has the identical internal magnetic group after relaxing atoms and selecting a SOC-consistent magnetic domain. Its **inter-partner mapping** follows the reported `m'_z` connector; internal allowedness is an independent gate.

## Structural and magnetic connection

- Structural pairing: `t1=(0.5,−0.16)` and `t2=(0.5,+0.16)`; a swap of layer roles plus the in-plane translation vector `(1,0)` implements the `m_z` structure connection: `t2=−t1+(1,0)` modulo lattice.
- Type-I magnetic connector: `m_z` (unitary; `R=diag(1,1,−1)`, `det R=−1`).
- Type-II (4%-strain) magnetic connector: `m'_z` (antiunitary; same R, `θ=1`).
- Both reverse polar `P_z` because `R_zz=−1`.

With no mixed in-plane component under this diagonal operator, `C_xy` switches odd for both connectors; `D_xx` switches odd under unitary `m_z`, even under `m'_z`. The `z`-axial anomalous Hall pseudo-vector is **T-odd**, transforms as `(-1)^θ det(R) Rzz`, and switches **even / odd**, respectively. This makes AHE and D_xx a *crossed-parity* candidate; unlike both C and D, AHE is **not** transformed as the polar-current←axial-magnetic-field susceptibility.

| Connection | `P_z` | `C_xy` | `D_xx` | AHE `σ_xy` |
| --- | --- | --- | --- | --- |
| Type-I `m_z` | − | − | − | + |
| Type-II `m'_z` | − | − | + | − |

Prediction is symmetry-exact **if** the two states really are related by the specified connector, share matched magnetic/lattice settings, and the responses use the declared real-phase convention. Unrelated doping or different magnetic chemical potentials invalidate comparing numerical magnitudes as if they were the same physical switch. In an insulator C=0 and the `C_xy` sign is then not an experimentally useful toggle.

## Why this is potentially more than a symmetry exercise

A quantitative DFT/Wannier result could demonstrate **magnetic-order-selective ferroelectric readout**: the same BSF two-state structure reverses the quadrature IGMC in Type-I but not in Type-II, while AHE has the opposite switching parity. A successful claim requires *finite, converged* response from self-consistent electronic structure, a measured or justified magnetic order, and a verified full projected spin/Berry-position operator chain. This document alone establishes **no new real-material current**.
