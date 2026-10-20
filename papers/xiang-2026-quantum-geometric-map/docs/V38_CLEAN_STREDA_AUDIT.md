# v38: Antisymmetric clean Kubo–Středa sea term as an independent scientific benchmark

**Reference context:** Xiang et al., *Quantum geometric map of magnetotransport*, PRB **113**, L201406 (2026), DOI `10.1103/kd12-4ycj`. This is an **independent extension** motivated by the previously checked v18 symmetric spectral response; it is not a formula quoted from the paper or a claim about its original DC PHE.

**Acceptance:** `PASS`, but only for the **clean, nondegenerate, two-band, zero-temperature, fixed scalar retarded regulator Fermi-sea benchmark**. Full metallic/disordered Kubo–Středa and original paper's double-adiabatic mixed `EB` response remain `PARTIAL / OPEN`.

## 1. Distinguish tensor channels

The previously migrated `spectral_kubo_v18.py` evaluates a **symmetric dissipative spectral kernel** from `RA+AR−RR−AA`, including scalar-disorder SCBA/ladder. Its completion does not calculate the **antisymmetric Hall sea response**.

For an independent two-band periodic test system, define

\[
 H_{\mathrm{QWZ}}(\mathbf k)=\sin k_x\,\sigma_x+\sin k_y\,\sigma_y+[m+\cos k_x+\cos k_y]\sigma_z.
\]

With \(G^{R,A}=(\epsilon-H\pm i\eta)^{-1}\), positive constant \(\eta\), and exact velocities \(v_{x,y}=\partial_{k_{x,y}} H\), calculate the *bare clean sea*:

\[
K_{xy}^{\mathrm{sea}}(\mu,\eta)=\frac{1}{4\pi}\int_{-\infty}^{\mu}\!d\epsilon\,\mathrm{Tr}\{v_xG^Rv_y\partial_\epsilon G^R-v_x\partial_\epsilon G^Rv_yG^R -(R\to A)\}.
\]

This is **a precisely defined normalized Green-function sea kernel**. The program does not pretend that a self-energy, ladder, true Středa surface term, and contact terms have been included.

For the two eigenvalues \(\epsilon_l<\epsilon_u\), \(\Delta=\epsilon_u-\epsilon_l\), \(A=\mathrm{Im}[(v_x)_{lu}(v_y)_{ul}]\), \(z_{l,u}=\mu-\epsilon_{l,u}+i\eta\), the independently derived exact energy primitive is

\[
J=-\frac{1}{\Delta}\left(\frac1{z_l}+\frac1{z_u}\right)+\frac{2}{\Delta^2}(\log z_l-\log z_u),\quad K^{\mathrm{sea}}(\mathbf k)=-\frac{A}{\pi}\mathrm{Im}J.
\]

The zero-integration constant is fixed by the limit \(\epsilon\to-\infty\) in the upper complex half plane. **In a gap** and \(\eta\to0^+\), this tends to the occupied-band Berry curvature \(\Omega_l=2A/\Delta^2\), with the sign orientation fixed by the independently calculated Fukui–Hatsugai–Suzuki (FHS) lattice Chern number.

BZ conversion: \(\sigma_{xy}/(e^2/h)=2\pi\langle K_{xy}(\mathbf k)\rangle_{BZ}\). The clean QWZ model has \(C=-1\) for \(-2<m<0\), \(C=+1\) for \(0<m<2\), and \(C=0\) for \(|m|>2\) in this orientation. We explicitly refuse the degeneracy thresholds.

## 2. Three independent clean-gap checks

1. **Direct trace energy integral:** numerical quadrature of the exact matrix Green functions and \(\partial_\epsilon G=-G^2\), with no eigenvector-based analytic primitive, agrees at representative \(\mathbf k\) points.
2. **Analytic primitive:** integrates the rational two-pole Green functions without a numerical energy mesh.
3. **Topological invariant:** gauge-invariant FHS link-plaquette phase sum, independent of velocity matrices and Green-function response integrals; additionally tests arbitrary Bloch U(1) phase choices and a constant SU(2) basis change.

| QWZ mass \(m\) | FHS \(C\) | occupied-band Berry \(C\) | bare sea at \(\eta=.02\) |
|---:|---:|---:|---:|
| −3 | 0 | 0 | −0.000770 |
| −1 | −1 | −1 | −0.975094 |
| +1 | +1 | +1 | +0.975094 |
| +3 | 0 | 0 | +0.000770 |

For \(m=-1,\mu=0\), bare sea approaches \(-1\) as \(\eta\) decreases: \(-0.803745\) at \(\eta=.16\), \(-0.925426\) at .06, \(-0.975094\) at .02 and \(-0.990035\) at .008. A symmetric Fermi-surface dissipative control vanishes instead approximately as \(\eta^2\) inside the gap.

## 3. What changes when the chemical potential enters the band

For \(m=-1, \mu=1.1,\eta=.02\), with a 101×101 k grid:

| Quantity | Dimensionless \(e^2/h\) value |
| --- | ---: |
| **bare sea only** | +2.447259 |
| compensating rational **boundary-pole diagnostic** | −2.711836 |
| sum, equivalent to the log-occupation part | −0.264577 |
| sharp occupied-band Berry curvature integral | −0.219181 |

This illustrates why a **bare sea kernel** cannot be identified with the full metallic Hall conductivity. The rational `1/z` pole piece has a sharp Fermi-level contribution as \(\eta\to0\), while the logarithmic term approaches the occupied-band Berry sum. The code exposes the exact **algebraic compensating boundary** to diagnose the decomposition. **Crucially, this compensating term is NOT claimed to be the physically complete Kubo–Středa Fermi-surface term**: deriving that requires the full response convention, frequency-dependent scattering self-energy and corresponding vertices.

The finite-\(\eta\) difference between the broadened log-occupation response and the zero-broadening occupied Berry sum is expected. There is no fitted correction, no SI PHE extrapolation and no hidden factor of \(2\pi\).

## 4. Reproduce and reuse

```bash
cd papers/xiang-2026-quantum-geometric-map
python -m pip install -e .
python -m unittest tests.test_streda_clean_sea_v38 -v
python scripts/reproduce_streda_v38.py --mesh 61 --mass -1 --mu 0 --eta .02
python scripts/reproduce_streda_v38.py --mesh 101 --mass -1 --mu 1.1 --eta .02
```

The CLI emits protocol-labelled JSON with the FHS Chern, Berry sum, bare sea and (when requested) algebraic boundary diagnostic. The JSON file `results/frozen/v38_clean_streda_reference.json` captures representative inputs and numbers.

**Next scientific gate:** establish the full metallic Kubo–Bastin/Středa decomposition in a declared microscopic disorder model. When \(\Sigma^{R,A}(\epsilon)\) depends on energy, \(\partial_\epsilon G=-G[1-\partial_\epsilon\Sigma]G\), and consistency demands companion vertices/contact terms. Compare with an exact disordered finite-system Kubo benchmark and the physical \(q\to0,\omega\to0\) ordering before returning to the original PRB mixed \(EB\) PHE. This v38 result **does not** settle the paper's Eq. (11) vs Fig. 2(b) \(2\pi\) difference or orbital three-band physics.
