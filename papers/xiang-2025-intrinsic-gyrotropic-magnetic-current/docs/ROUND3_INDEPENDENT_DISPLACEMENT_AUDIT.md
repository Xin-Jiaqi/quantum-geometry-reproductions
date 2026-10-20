# Round 3 | Independent displacement IGMC prefactor audit

Date: 2026-10-19. Source: Xiang et al., *PRL* **134**, 116301 (2025), main Eqs. (5)–(12), (17), Fig. 2(d); **SM p. 6, Eqs. (33)–(34)**. All source equation transcriptions were inspected directly against the supplied PDFs. The current model methods are **independent reproduction**, not author-provided code/data.

## 1. A source-local twofold discrepancy is real, but source figure provenance remains unknown

Define `z_nm = r_nm^α σ_mn^β`. The main PRL defines

\[
Q_{nm}^{αβ}=(z_{nm}+z^*_{nm})/2=\operatorname{Re}z_{nm},\qquad
Z_{nm}^{αβ}=i(z_{nm}-z^*_{nm})=2\operatorname{Re}(iz_{nm}).
\]

However, the text **immediately after SM Eq. (34) on p. 6** prints `Q_nm=2 Re(z_nm)` and `Z_nm=Re(i z_nm)`. These are respectively `2 Q_main` and `Z_main/2` if the same symbols mean the main article's observables. This is *not* merely a matter of converting Pauli `σ` to a spin operator `σ/2`, because that would change both channels simultaneously.

The cosine term in SM Eq. (34) can be audited **without assuming either definition**. The first line of the SM wave-packet expression contains

\[
\mathrm{Re}\, z\Bigl[\frac{e^{-i\omega t}}{\omega-\Delta}
               +\frac{e^{+i\omega t}}{-\omega-\Delta}\Bigr],
\quad \Delta=\varepsilon_m-\varepsilon_n.
\]

Its exact cosine coefficient is

\[
\operatorname{Re}(z)\left(\frac{1}{\omega-\Delta}
                       -\frac{1}{\omega+\Delta}\right).
\]

The next displayed SM expression replaces `Re(z)` with its printed `Q_SM=2 Re(z)` **but keeps the full same frequency bracket**. Read literally, that is **twice the first line**. In the low-frequency limit the correct coefficient is `−2 Re(z)/Δ`, whereas that intermediate line becomes `−4 Re(z)/Δ`. The main-text `N_n = −∑_m 2Q_main/Δ` and main Eq. (12) agree with the **first line**, not this literal doubled intermediate substitution.

For the sine coefficient, the SM notation `Z_SM=Re(iz)` combined with its different bracket is algebraically consistent with the first line: do **not** double the conduction channel. The localized issue is specifically the `Q` prefactor in the intermediate cosine line. The SM *last* line labels its result as main Eq. (10), but this identification and the preceding literal intermediate cannot both be true without an unstated normalization adjustment.

**Careful status:** the source contains a factor-two tension. Whether the authors' Fig. 2(d) code used the main formula, the displayed SM intermediate or a different normalization **cannot be established from PDFs alone**. We retain the two paths as distinct `paper_main_eq12` (independently verified) and `sm_eq34_literal_intermediate_diagnostic` (**not validated physical theory**). No fitted scaling is used in either.

## 2. Independent verification 1: finite-field Wilson-loop polarization (no QGT formula)

The source uses `H_1(t) = −b(t) σ_x` with `b=g μ_B B`. Freeze a small **static** `b`, diagonalize the original lattice-Dirac Hamiltonian at every `k`, and calculate the Berry-phase polarization of the **occupied lower band**, taking gauge-invariant discrete Wilson loops along `k_x`:

\[
\phi(k_y;b)=-\arg\prod_{j=1}^{N}
\frac{\langle u_{j,k_y}(b)|u_{j+1,k_y}(b)\rangle}
{|\langle u_{j,k_y}(b)|u_{j+1,k_y}(b)\rangle|}.
\]

Then compute the BZ-integrated Berry-connection derivative via a **central finite difference** in `b` (individual phase differences unwrapped modulo `2π`). This calculates `∂_b ∫ A_x d²k/(2π)²` with no `v_nm`, `r_nm`, `Q_nm` or `Z_nm` in the implementation. Independent perturbation theory predicts this derivative to equal the main `∑ f_n 2Q_nm/ε_mn`, or `σ_D/ω` in paper units.

| N × N | Wilson derivative (Å⁻¹/eV) | main Eq.(12) / ω (Å⁻¹/eV) | absolute relative error |
|---|---:|---:|---:|
| 41² | 0.01518277072 | 0.01523664701 | 0.354% |
| 101² | 0.01522542432 | 0.01523407804 | 0.0568% |
| 201² | 0.01523189300 | 0.01523407804 | 0.0143% |
| 301² | 0.01523310377 | 0.01523407804 | 0.00640% |

`b = 10⁻⁴ eV`; `b=10⁻³ eV` changes the result only ~1.8 × 10⁻⁸ in these units. The shrinking difference is dominated by the **finite Wilson-link discretization**, not a magnetic-field fit. This independently supports **main Eq. (12)** with its original prefactor, and does not support using an extra factor two in the physical polarization derivative.

## 3. Independent verification 2: direct orbital-basis Liouville equation

Without computing eigenvector Berry connections or projecting into an eigenband interband sum, construct the equilibrium *orbital-basis* density matrix `ρ₀=(I−H/E)/2` for the gapped half-filled two-band Hamiltonian and solve, at each k,

\[
[\Omega\mathbb 1-\mathcal L_H] \delta\rho_\Omega
=-[\sigma_x,\rho_0],\qquad \mathcal L_H X=HX-XH,
\quad \Omega=\hbar\omega.
\]

Then `K(Ω)=−∫ Tr[v_x δρ_Ω]d²k/(2π)²` (field amplitude normalized to 1 energy unit). This uses a direct 4 × 4 matrix solve at each momentum, not a resubstitution of the main geometric formula. For an AC driving field the **quadrature displacement amplitude** is `Im K(Ω)`.

At `f=1 THz` with `Ω=0.004135667694 eV`, **N=101²**:

- orbital-basis Liouville `Im K = 6.301090464897298 × 10⁻⁵ Å⁻¹`;
- independently coded spectral Kubo `Im K = 6.301090464913824 × 10⁻⁵ Å⁻¹`;
- **absolute difference = 1.65 × 10⁻¹⁶ Å⁻¹**.

The leading-order main Eq. (12) gives `ω × 0.0152340780 ≃ 6.300 × 10⁻⁵ Å⁻¹`, differing from the finite-ω result only by an expected `O(ω³)` correction. It is incompatible with an *additional physical* factor 2 under the stated Zeeman convention.

## 4. Fig. 2(d) SI comparison and limits of source visibility

Parameters exactly as Fig. 2 caption: `β=.25eV, v_F=.825 eVÅ, M=.1 eVÅ², a=1Å`, and explicitly declared `g=2`, `f=1 THz`. The SI multiplier is the independently calibrated `281.794032638 A/(m T)` per Å⁻¹.

| Chemical potential | Main Eq. (12) (A/m/T) | SM Eq. (34) literal middle ×2 (diagnostic only) | Source Fig. 2(d), visual reading (± rough uncertainty) |
|---|---:|---:|---:|
| 0 eV (in gap) | +0.017754 | +0.035508 | about +0.040 ± 0.003 |
| +1 eV | −0.06419 | −0.12837 | about −0.135 ± 0.010 |

The source ordinate explicitly reads `−gμ_B σ_D,xx` in `10⁻² A m⁻¹ T⁻¹`. Readings are **manual approximate** from a PDF figure, not deposited data or verified author parameters. Literal SM ×2 is **numerically close** to the original graph, but not quantitatively identical. **This agreement is circumstantial, not a license to use ×2 as the correct physical prefactor.**

## 5. Scientific verdict / next decisive dependency

- **PASS:** clean gapped lattice-Dirac linear response with paper `H1`, main Eq. (12) prefactor (2 independent *methods*: Wilson finite-field phase and orbital Liouville, plus original spectral Kubo), SI/frequency conventions unchanged, mesh convergence evidenced.
- **SOURCE-INCONSISTENT:** local line following SM Eq. (34), where the symbol `Q` has a literal ×2 mismatch with main Eq. (5) and the immediately preceding first line of SM Eq. (34).
- **OPEN:** exact provenance of the plotted Fig. 2(d) data, including whether its plotted prefactor arose from a different convention or the intermediate SM expression. No original source code/raw plotting data in attached PDFs.
- **NOT IN SCOPE FOR THIS ROUND:** full SM magnetic point group enumeration, clean-metal finite-frequency intraband protocol, real material Wannier/experimental claims.

The next useful step is source-level evidence (author code or raw data, if publicly released). Without it, close the **validated-main-formula subset** and leave paper graph discrepancy open rather than changing the independently validated expression.

## Reproduce

```bash
python scripts/audit_round3.py
python -m unittest discover -s tests -v   # needs PYTHONPATH=src or package install
python scripts/verify_release.py
```

Machine evidence: `parameters_round3.json`, `results/raw/round3_independent_audit.json`, `round3_wilson_convergence.csv`, `round3_displacement_curve.csv`, `results/figures/round3_independent_factor_audit.png`. Prior frozen references preserved unchanged.
