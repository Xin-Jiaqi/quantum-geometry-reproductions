# Round 8 | SM §III–IV: wavepacket and lesser-Green's independent derivation

**Source**: Xiang et al., *Phys. Rev. Lett.* **134**, 116301 (2025), main Eqs. (5)–(12), SM §III pp. 4–6 Eqs. (20)–(34), SM §IV pp. 6–7 Eqs. (35)–(47). Both original PDFs were consulted; publisher pages are not bundled in the public-facing code package. Results refer solely to the published **two-band Dirac model** Eq. (17), 1 THz and the **clean subgap** regime. They do not establish a general multi-band, dissipative, resonant or degenerate proof.

## 1. Direct wavepacket expansion (SM Eqs. 21, 24)

Take an occupied band `n`, an unoccupied band `m`, `Δ=ε_m−ε_n>0`, photon energy `Ω=ℏω`, the perturbation `H1=-b cos(ωt) σ_β`, `b=g μ_B B`. The magnetic field's **two Fourier components are `b/2` each**. SM Eq. (21) yields

\[
M_{mn}(t)=-\frac{b}{2}\sigma_{mn}\left[\frac{e^{-i\omega t}}{\Omega-\Delta}+\frac{e^{i\omega t}}{-\Omega-\Delta}\right].
\]

From Eq. (24), `a_n^(1)=2 Re[M_mn r_nm]`. Define `A=r_nm^α σ_mn^β=x+iy`, as in **main** Eqs. (5)–(6): `Q=x`, `Z=i(A-A*)=-2y`. Evaluating the two harmonics directly gives

\[
\frac{a_n^{(1)}}b=\frac{2\Delta Q}{\Delta^2-\Omega^2}\cos(\omega t)-\frac{\Omega Z}{\Delta^2-\Omega^2}\sin(\omega t).
\]

As `Ω→0`, the **cosine** coefficient approaches `2Q/Δ`, the source of the main Eq. (12) displacement term. Importantly, it uses `Q=Re(A)`, not `2 Re(A)`.

**SM Eq. (34) internal mismatch**: its first line correctly has the wavepacket `2Re[ ... Bbar(t)]`, where each harmonic's `Bbar` contains `−b/2`. The second line extracts `−b Re(A) × [frequency bracket]` but labels the geometric factor `Q=2Re(A)`, thus doubling its cosine contribution. The main Eq. (5) `Q=Re(A)` and independent Wilson-loop, Liouville, Pauli-trace and new wavepacket checks agree. Literal SM ×2 remains an **invalidated diagnostic**, not our result.

**Additional finite-frequency sign discrepancy**: Main Eq. (10) sets `F_n=−∑Ω Z_nm/Δ²` and `P=−b(N cos+F sin)`, implying `P_sin/b≈+ΩZ/Δ²`. Direct SM Eqs. (21),(24) instead give `P_sin/b≈−ΩZ/Δ²`. These are opposite signs in a convention where `Z` is defined by main Eq. (6). However, the sine term is **order Ω in polarization** and contributes first at **order Ω² to the current**. The original main Eq. (12) explicitly drops terms of order `ω²`; thus this discrepancy does not overturn the accepted leading displacement result or explain Fig. 2(d)'s factor-of-two scale. It is an issue with the subleading quoted polarization expression and/or convention, not automatically a claim of erroneous published observables.

## 2. Matrix-resolvent lesser-Green's route (SM Eq. 44)

Instead of rewriting the previously used pair `rσ` band formula, the new module assembles retarded/advanced **orbital-basis** Green resolvents at their unperturbed spectral poles:

\[
G_1^<(\Omega)=i\sum_n f_n\{G^r_0(\varepsilon_n+\Omega)\Sigma P_n+P_n\Sigma G^a_0(\varepsilon_n-\Omega)\},\quad
G^{r/a}_0(E)=(E-H_0\pm i\eta)^{-1}.
\]

Here `P_n=|u_n><u_n|`, `Σ=−b σ_β`; the dimensionless computed current coefficient is `K(Ω)=−J(Ω)/b`. The terms within the same band have opposite denominators `1/(Ω+iη)+1/(−Ω−iη)=0` and are projected out **only after the explicit algebraic cancellation**, avoiding catastrophic floating-point subtractive cancellation. Cross-band matrix inversions are evaluated without geometric band-sum functions. Independent projectors and matrix products are kept in the original orbital basis.

For finite `η` we also compare to the **original spectral Eq. (8)** with the same `+iη` prescription. For `η=0`, reconstructing physical current from SM Eq. (32), i.e. `J(t)/b=−Z cos + ∂_t(a^(1)/b)`, gives the exact same complex K as the lesser Green solution:

\[
K_{\mathrm{GF}}(\Omega)=K_{\mathrm{wavepacket}}(\Omega)
=\frac{\Delta^2 Z+i\,2\Omega\Delta Q}{\Delta^2-\Omega^2}
\]

for a single occupied–unoccupied band pair, with the **fixed** main text `e^{-iωt}`, spin and Zeeman conventions. For all 6 tested tensor components on a `101×101` lattice grid, the maximum per-k complex difference is `1.39×10^-13` in the model units (numerical inversions at high matrix condition numbers). This is an algebraic/numerical equivalence, not independent confirmation of the author's original plotting scale.

At `f=1 THz`, `μ=0`, computed using orbital resolvents at `101²` points: `Im K_xx=6.30109046490295e−5 Å⁻¹`, giving `0.01775609692120911 A/(m·T)`. This differs slightly from the independently converged Gauss baseline `0.01775389319` **solely due to the coarser 101² mesh**. This datum is a cross-method check; it does **not** replace the converged Round 6/7 benchmark.

## 3. SM Eq. (46) Fourier convention

The SM declares `B(ω)=∫dt exp(iωt) B cos(ω1 t)` but prints `2πB[δ(ω−ω1)+δ(ω+ω1)]`. Under this exact Fourier convention, the mathematical transform is instead

\[
\boxed{B(\omega)=\pi B[\delta(\omega-\omega_1)+\delta(\omega+\omega_1)]}.
\]

An independent DFT/FFT of a precisely sampled single cosine cycle returns complex positive- and negative-frequency amplitudes **0.5 B** each (when normalized by sample length), while an incorrectly doubled transformation would return `B` each. **SM Eq. (47) immediately restores the correct prefactor `1/2`** on the current, inconsistent with the preceding Eq. (46) as printed. This is a clear **Fourier normalization typo/local error**, not evidence of what the original Fig. 2 plotting program used. It does **not** explain a factor-two difference in the accepted Eq. (12), which has already been independently checked.

## 4. Scientific decision and remaining blockers

- `PASS`: two-band clean finite-frequency **wavepacket ↔ lesser-Green's ↔ orbital Liouville ↔ spectral Kubo** model kernel in the fixed physical protocol. All tested 6 tensor components agree pointwise; finite positive `η` also tested against spectral expression.
- `PARTIAL`: general multi-band semiclassical treatment, Green's self-energy for actual disorder, and finite-frequency polarization preparation are beyond this model validation.
- `OPEN`: main Eq. (10)'s subleading sine term sign; main Eq. (14)'s antiunitary complex response convention; the source of the main Fig. 2(d) amplitude.
- `SOURCE ERRORS FLAGGED`: SM Eq. (34)'s doubled Q, SM Eq. (46)'s doubled Fourier amplitude. Both are recorded as source-algebra issues; **do not multiply the validated numerical output by two**.
- Source Table II `−42'm'` symmetry relation issue remains independent.

## 5. Reproduce

From the paper directory:

```bash
PYTHONPATH=src python scripts/audit_round8.py --no-write
PYTHONPATH=src python -m unittest tests.test_round8_green_wavepacket -v
PYTHONPATH=src python scripts/verify_release.py
```

Frozen data: `results/frozen/round8_green_wavepacket_reference_v1.json`; audit: `results/raw/round8_GF_wavepacket_audit.json`; visualization `results/figures/round8_GF_wavepacket_audit.png`.
