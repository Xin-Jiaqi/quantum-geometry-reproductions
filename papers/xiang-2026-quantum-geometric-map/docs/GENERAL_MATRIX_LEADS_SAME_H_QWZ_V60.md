# v60 | Independent principal-layer matrix electrodes and same-H finite/periodic topology

**Physics-first:** the original Xiang paper's disputed supplemental formulas are not a truth criterion. This release proves two bounded abilities, NOT a completed intrinsic bulk orbital EB/MNHE response.

## I. True semi-infinite matrix principal-layer lead

With Hermitian orbital-layer block H00 and **general noncommuting** interlayer hopping H01, retarded z=E+iη, η>0:

\[
g_s(z)=[zI-H_{00}-H_{01}g_s(z)H_{01}^\dagger]^{-1}.
\]

Accelerated López-Sancho decimation is tested against the exact one-orbital retarded root, a **separate direct sparse inversion of a 160-layer two-orbital chain**, explicit Dyson residuals and nonnegative spectral density, and arbitrary lead-orbital unitary frame changes.

The full device-contact matrix C_a (nDevice × nLead) gives \(\Sigma_a=C_ag_sC_a^\dagger,\ \Gamma_a=i(\Sigma_a-\Sigma_a^\dagger)\). Nonlocal couplings may yield **off-diagonal Γ matrices**. All Caroli transmission and magnetic derivative calculations use full Γ matrices, not silently diagonal linewidths.

For central gauge transformation G(B), contact C'_a(B)=G(B)C_a, so \(\dot C'_a=(\dot G)C_a\), and

\[
\dot\Sigma_a=\dot C_ag_sC_a^\dagger+C_ag_s\dot C_a^\dagger,\quad
\dot\Gamma_a=i(\dot\Sigma_a-\dot\Sigma_a^\dagger),\quad
\dot G^r=G^r(\dot H+\sum_a\dot\Sigma_a)G^r.
\]

Tests include a **true multi-site contact, non-diagonal Γ**, and a deliberately translated four-site graph (where G(B) is not accidentally the identity). **Omitting d_BC fails the magnetic gauge derivative**, while the full calculation passes independent finite-B Caroli differences and full floating-probe Hall differences. The one-orbital matrix lead converges toward the earlier v59 reference at decreasing positive η.

Reference for two-orbital layers, η=.003 and declared coherent four-terminal parameters: Hall V_H/ΔV **-0.0680077538749482**, derivative **-0.15990871963649705**, longitudinal current per bias **0.5285663187319317 q²/h**. Two magnetic gauges on shifted positions agree at T matrix ~2.4e-13 and d_BT ~8.6e-13.

**Strict assumptions:** η is explicitly finite, all internal lead H00/H01 are B-independent and the d_BC interface supports central-frame gauge changes only. This is not a full magnetic Peierls model of arbitrary three-dimensional electrodes or an implementation of real-material Wannier lead dispersions.

## II. SAME QWZ Hamiltonian finite-open versus periodic geometry

Periodic QWZ \(H(k)=\sin k_x\sigma_x+\sin k_y\sigma_y+(m+\cos k_x+\cos k_y)\sigma_z\). Its *identical Fourier hoppings* produce the finite open lattice: on-site mσz, \(T_x=(\sigma_z-i\sigma_x)/2\), \(T_y=(\sigma_z-i\sigma_y)/2\). Tests independently reconstruct the Bloch matrix from these hoppings at multiple k points.

For occupied projector P (Q=I-P), define Bianco–Resta local marker

\[
C(r)=-4\pi\,\mathrm{Im}\sum_{\alpha\in r}
\langle r\alpha|P X Q Y P|r\alpha\rangle.
\]

The finite-open total marker sums to zero by edge compensation; the **bulk interior** converges toward the v59 separately verified *periodic* Chern:

| size (m=-1) | finite central marker | periodic Chern |
| --- | ---: | ---: |
| 5x5 | +0.985370245 | +1 |
| 7x7 | +0.996590312 | +1 |
| 9x9 | +0.999174829 | +1 |
| 11x11 | +0.999796737 | +1 |

The m=+1 phase is -1, m=±3 trivial. Additional controls confirm finite trace cancellation, complex conjugation flipping the marker, site-orbital unitary covariance and input rejection.

**Scope warning:** this is a genuine same-H **occupied projector topology** comparison, not evidence that finite Hall voltage and bulk Kubo conductivity numerically agree; current-contact/scattering limits are still missing.

## III. Physical acceptance

| Physical capability | Status |
| --- | --- |
| Causal matrix surface Green, Γ spectral, exact decimation and sparse oracle | VERIFIED in specified η>0 semi-infinite block chain |
| Full matrix Γ/Σ Caroli and magnetic gauge contact derivative | VERIFIED specified finite coherent model |
| Same-H QWZ open projector Chern marker versus periodic Chern | VERIFIED topology/geometry |
| Finite four-probe voltage given explicit electrode parameterization | CONDITIONAL |
| Same-H **transport** thermodynamic limit, general B-dependent Wannier lead, intrinsic orbital EB/MNHE/PHE | UNRESOLVED |

21 independent local tests; full remote PRB regression required before merge. The original paper's scientific source audit remains PARTIAL.

Run from the PRB paper subtree:

    python -m unittest tests.test_matrix_lead_marker_v60 -v
    python scripts/reproduce_v60.py

See the physics capability registry and frozen v60 JSON.