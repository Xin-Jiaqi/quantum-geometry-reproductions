# BSF × IGMC Stage 5: physical readout interpretation

**Status: symmetry prediction, not a material-response calculation.** This is a separate extension under `research/bsf-igmc/agf2_preflight/`, not part of the frozen Xiang reproduction. The low-frequency spin-Zeeman, real-phase convention is explicit; the full-frequency antiunitary ambiguity in the Xiang paper remains open.

## What is driven?

In the insulating, clean, low-frequency regime, the relevant piece of the Zeeman quantum metric response is a *field-induced electronic polarization*:

\[
P_i^{(B)}(t)=\alpha^{PB}_{ij} B_j(t)+O(\omega),\qquad
j_i^{D}(t)=\partial_t P_i^{(B)}(t).
\]

For `B_j(t)=B_0 cos(ωt)`, the leading displacement response is in the **quadrature** of the magnetic drive, `j_i^D(t)=-ω α^{PB}_{ij}B_0 sin(ωt)`, with overall conventions for the spin Zeeman coupling absorbed in `α^{PB}`. In particular it can exist without mobile carriers at zero temperature. The net charge pumped over a full closed sinusoidal cycle is zero: `∫_0^{2π/ω}j_i^D dt=0`. This is an *AC sensing mechanism*, not a DC nonvolatile current by itself.

The polarization `P_i` is polar and even under time reversal, whereas the magnetic field `B_j` is axial and odd. Consequently, the *real adiabatic magnetoelectric polarization coefficient* transforms between properly symmetry-related magnetic domains as

\[
\alpha^{PB}_{(2)}=(-1)^\theta\det(R)\, R\alpha^{PB}_{(1)}R^T.
\]

This statement follows from transformation of the physical fields; it does not rest on the contested blanket Eq. (14) rule for both complex finite-frequency C/D sectors. For the measured reactive IGMC `j^D=\dot P^B`, this implies the same inter-domain switching parity in the stated low-frequency convention.

## AgF₂: crossed switching of two physical probes

For the AgF₂ bilayer considered by Zhu et al., the ferroelectric domains `B₁` and `B₂` are *source-reported* as stable, with `P_z` opposite. They are related by the spatial connector `m_z` for the Type-I ordering, while the Type-II spin-domain connector after **4% biaxial strain** is `m'_z`. Work in the same laboratory axes for both domains and take `R(m_z)=diag(1,1,-1)`.

- In-plane polarization `P_x` and `P_y` are unchanged under `m_z` while in-plane **axial** `B_x` and `B_y` reverse. Thus the longitudinal magnetoelectric polarization coefficients `α^{PB}_{xx}` and `α^{PB}_{yy}`, and hence the reactive `D_xx` and `D_yy`, **reverse** between Type-I ferroelectric domains.
- With `m'_z=m_z T`, the additional time reversal flips `B` again. Therefore the same in-plane `α^{PB}_{xx,yy}` and leading reactive `D_xx,yy` are **even** under the Type-II inter-domain switch.
- The anomalous Hall vector `σ_xy^{AHE}` is **axial and time odd**, so it is even under `m_z` but odd under `m'_z`: a complementary switching character to the longitudinal insulating-eligible IGMC signal.

| Measured quantity | Type-I, `m_z` | Type-II under 4% strain, `m'_z` | Regime |
| --- | ---: | ---: | --- |
| `P_z` | odd | odd | structural FE domains |
| `D_xx` (`j_x` from AC `B_x`) | **odd** | **even** | insulating eligible, **not calculated for AgF₂** |
| `D_yy` (`j_y` from AC `B_y`) | **odd** | **even** | insulating eligible, **not calculated for AgF₂** |
| `D_yz` (`j_y` from AC `B_z`) | even | odd | symmetry-allowed in Type-I; out-of-plane drive has orbital backgrounds |
| AHE `σ_xy` | **even** | **odd** | only if AHE itself is nonzero at the selected Fermi level |

Type-I's **single-domain** magnetic point group `m'_x` permits `D_xx, D_yy, D_yz, D_zy, D_zz` and `C_xy, C_xz, C_yx, C_zx` under the separately audited real-phase convention. **Type-II single-domain allowedness is not yet independently established** from a relaxed magnetic crystal. Symmetry does not guarantee nonzero responses or a useful experimental magnitude. For a conventional topologically trivial insulator, in-gap **DC AHE can vanish** even if symmetry allows it; do not promise two nonzero readout signals at one and the same chemical potential.

## Physics and computational priorities

1. **Prioritize `D_xx` and `D_yy` together.** AgF₂ Type-I easy axis is reported as `[100]`; the spin Zeeman coupling for `B_y` can probe different interband spin matrix elements from `B_x`. Symmetry makes both available, but only SOC-inclusive same-gauge matrix elements can establish which is larger. Do **not** assert `D_yy > D_xx` before computation.
2. **One controllable material comparison, two separate campaigns.** Type-I source phase is zero strain; Type-II source phase is 4% biaxially strained. Comparing their parity tables is a *mechanism comparison*, not evidence of a magnetic configuration switch in an unchanged sample.
3. **Verify a true material response.** Source-verified relaxed AgF₂ `B₁/B₂` structures, spin magnetic order, SOC, the band gap, and the Wannier-gauge spin/position operators are needed for a multiband `α^{PB}` or IGMC; `wannier90_hr.dat` plus centers alone cannot do it.
4. **Experimental background is essential.** Time-varying B creates an electric field through Faraday's law, and orbital magnetic moments can contribute to the response. Estimate their size before declaring a purely spin-Zeeman measurement.

**Claim allowed today:** an AgF₂-based *crossed-parity magnetoelectric AC readout proposal* with exact low-frequency inter-domain transformation laws **conditional on the reported symmetry/domain mappings**.

**Claims not allowed today:** DFT-predicted nonzero AgF₂ IGMC, a measured signal, a universal magnitude enhancement, a simultaneous AHE/IGMC readout in the same insulating window, or independent reproduction of the source's FE double-well.

### Sources and dependencies

- Xiang et al., *Phys. Rev. Lett.* **134**, 116301 (2025), DOI 10.1103/PhysRevLett.134.116301: current decomposition into conduction and magnetic-field-induced-polarization displacement terms, especially Eqs. (9)–(12).
- Zhu et al., *Phys. Rev. Lett.* **135**, 056801 (2025), DOI 10.1103/dmzg-ck2t: AgF₂ layer group, the B₁/B₂ ferroelectric pairing, Type-I/Type-II magnetic configurations and their strain dependence.
- Xin-Jiaqi/group-theory-operations-toolkit v0.18.0: magnetic operations and tensor-allowedness infrastructure; present analytic mirror derivation is independently coded and does not claim an installed upstream API comparison.
