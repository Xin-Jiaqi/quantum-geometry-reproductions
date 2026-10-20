# Round 6 | Eigenvector-free two-band proof, causal quadratures and an independent BZ integral

Date: 2026-10-20. Primary source: Xiang et al., *PRL* **134**, 116301 (2025), Eqs. (5)–(14), (17), Fig. 2(d); source SM Eq. (34). This is an independent paper-model audit, **not** author-provided plotting code. Paper PDFs and figures are not republished here.

## 1. A new derivation that never computes Bloch eigenstates

Consider **exactly the printed** two-band Hamiltonian, `H(k)=d(k)·σ`, with `E=|d|`, `n=d/E`, and `v_α=∂_αd·σ`. Zeeman perturbation is `H1=-b σ_β`, `b=gμ_B B`, not `-b σ_β/2`. For the occupied `-E` band and empty `+E` band, the gap is `Δ=2E`.

The Pauli trace identity for the band-pair product `A = < - |v_α| + >< + |σ_β| - >` is

\[
A=\partial_\alpha d_\beta-(\mathbf n\cdot\partial_\alpha\mathbf d)n_\beta
-i\,\mathbf n\cdot(\partial_\alpha\mathbf d\times\hat\mathbf e_\beta).
\]

Because the main text uses `r_-+ = i v_-+/Δ`, **without constructing any eigenvectors** we obtain

\[
\boxed{Z_{-+}^{\alpha\beta}=-\partial_{k_\alpha}n_\beta},\qquad
\boxed{\frac{2Q_{-+}^{\alpha\beta}}{\Delta}
=\frac{\mathbf d\cdot(\partial_{k_\alpha}\mathbf d\times\hat\mathbf e_\beta)}{2|\mathbf d|^3}}.
\]

These identities are independent of arbitrary band phases. Importantly, no extra factor of 2 can be inserted into the `Q` expression while keeping the fixed `H1=-b σ` convention. The code path `src/igmc_lab/d_vector_independent.py` implements these **directly from the three real Cartesian components of d**, without importing the original model's eigenvectors, geometric tensor sums or spectral Kubo routines.

Against the older eigenstate implementation, testing all `α=x,y`, `β=x,y,z`, and both `±` magnetization domains at **437 random k-points per domain** yields absolute per-point maxima of `1.78e-15` for Z and `2.66e-15` for `2Q/Δ` (floating point).

## 2. Full nonresonant complex response, not just static tensor parity

For filled lower vs upper occupations `Δf=f_-−f_+`, the band-pair finite-frequency expression is algebraically

\[
\boxed{K_{\alpha\beta}(\Omega)=\int_{\mathrm{BZ}}
\frac{d^2k}{(2\pi)^2}\,\Delta f\;
\frac{Z_{-+}^{\alpha\beta}+i\Omega\,2Q_{-+}^{\alpha\beta}/\Delta}
{1-(\Omega/\Delta)^2}},\quad \Omega=\hbar\omega.
\]

This identity follows by **combining both ±frequency denominators** of the original Eq. (8); for real sub-gap `Ω`, its real quadrature is the conduction-like term and its imaginary quadrature is the displacement-like term. The finite-frequency correction is explicitly `1/[1−(Ω/Δ)^2]`, and the low-frequency limit gives paper Eqs. (11) and (12). The sign is tied to the paper's `e^{-iωt}` convention and `H1=-bσ` perturbation.

The TR partner `H_T(k)=σ_y H^*(-k) σ_y` maps the Dirac Eq. (17) to `β,M -> -β,-M`, i.e. `mass_sign=-1`. Since the drive has spin-pseudovector parity, these **fixed real-phase response coefficients** obey, at nonresonant positive frequencies,

\[
\mathrm{Re}\,K_{H_T}(\Omega)=\mathrm{Re}\,K_H(\Omega),\quad
\mathrm{Im}\,K_{H_T}(\Omega)=-\mathrm{Im}\,K_H(\Omega),\quad
K_{H_T}(\Omega)=K_H(\Omega)^*.
\]

The code verifies both quadratures for **six Cartesian components at 1 and 3 THz**, 12 complete complex cases, with maximum TR residual `1.74e−17 Å⁻¹` and maximum difference from the separate spectral Kubo implementation `5.57e−17 Å⁻¹` (N=81²). This explains why a **single real `η_T=-1`** in printed Eq. (14) does not describe *both* real, phase-resolved channels as defined in Eqs. (11)–(12). It does **not** by itself prove how the author intended their full complex-frequency/antiunitary convention to be applied, and the printed Eq. (14) branch remains distinct.

## 3. BZ numerical integral by a genuinely different quadrature

With the unchanged **Fig. 2** parameters `β=0.25 eV`, `v_F=0.825 eV Å`, `M=0.1 eV Å²`, `a=1 Å`, `g=2`, and **cyclic** `f=1 THz` (`Ω=0.0041356676943898555 eV`), use tensor-product **Gauss–Legendre integration** across the full square Brillouin zone to evaluate the eigenvector-free `D/Ω` formula. Compare with the earlier uniform shifted midpoint grid, which *did* diagonalize H.

| Gauss order per axis | D/Ω [Å⁻¹ eV⁻¹] |
|---:|---:|
| 16 | 0.00267335568406 |
| 24 | 0.01186407559728 |
| 32 | 0.01439655200313 |
| 48 | 0.01518703722847 |
| 64 | 0.01523154610901 |
| 96 | 0.01523407088874 |
| 128 | **0.01523407801563** |
| Older shifted midpoint, N=301² | **0.01523407803570** |

Last two Gauss orders differ by **7.13×10⁻⁹** in the reported units; independent Gauss vs midpoint differs by **2.01×10⁻¹¹**. The marked underconvergence at N=16/24 demonstrates why low-resolution agreement is not an acceptance criterion.

Restoring physical factors without any fitted normalization:

\[
-g\mu_B\sigma^D_{xx}(\mu=0,f=1\,\mathrm{THz})
=\boxed{0.01775389319\;\mathrm{A/(m\,T)}}.
\]

The **literal doubled** middle line after SM Eq. (34) gives `0.03550778639 A/(m T)` but **is not supported** by the independent finite-field Wilson-polarization derivative, orbital Liouville solution, or the new Bloch-sphere derivation. Fig. 2(d) visually appears near `0.04 A/(m T)`, but this **manual figure reading** carries image-resolution uncertainty and is not deposited author data. The physical origin/provenance of the published plotting amplitude remains **OPEN**. No amplitude parameter was changed.

The new Gauss proof is restricted to **T=0 with μ in the gap**. Metallic occupations are discontinuous in k, and this Gauss convergence table cannot be reused as evidence that the doped `μ≈±1 eV` values are converged. Existing midpoint scans in metals are diagnostics only.

## 4. Source consistency and remaining questions

1. Main Eq. (5) uses `Q=Re(rσ)`; the SM §III intermediate after Eq. (34) uses `Q=2Re(rσ)` without compensating a prefactor, causing a doubled cosine term compared to both its preceding line and the main Eq. (12). Both the new eigenvector-free calculation and **two older independent response derivations** support the **main** normalization.
2. The phase-resolved C and D tensors transform differently under magnetic-domain reversal. Main Eq. (14) prints a common real antiunitary sign; ask about the intended frequency/response convention before labeling it a typographical error.
3. Prior SM Table II `\bar4 2'm'` equality omission (round 5) is an **independent lower-priority table issue**, and does not fix Fig. 2(d) amplitude.
4. The **only decisive input to adjudicate what produced the published Fig. 2(d)** would be the authors' original figure data or the exact plotting implementation with spin, Zeeman and current SI conventions. Without those, the validated subset can be frozen while full-paper figure acceptance stays `PARTIAL`.

## 5. Reproduce without changing old checkpoints

```bash
PYTHONPATH=src python scripts/audit_round6.py
PYTHONPATH=src python scripts/audit_round6.py --no-write
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python scripts/verify_release.py
```

Evidence: `results/raw/round6_analytic_dirac_audit.json`, `round6_gauss_convergence.csv`, `round6_complex_TR_quadrature.csv`, `results/figures/round6_independent_gauss_convergence.png`. A new `results/frozen/round6_dvector_reference_v1.json` guards its values and the SHA256 of all prior frozen baselines. None of the preceding frozen files is modified.
