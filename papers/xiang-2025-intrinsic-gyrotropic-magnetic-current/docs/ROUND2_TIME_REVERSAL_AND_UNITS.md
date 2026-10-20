# Round 2: time reversal, SI calibration, unresolved Fig. 2(d)

Date: 2026-10-19. Source: Xiang et al., PRL **134**, 116301 (2025), main Eqs. (5)-(17), Table I, main Figs. 1/2; supplementary §V and §VI. Independent model calculations, not an author-validated correction.

## 1. Source-defined response protocol

At nonresonant low frequency, paper Eqs. (9)-(12) define

\[
J_\alpha(t)=-g\mu_B\,B_\beta\left[\sigma^C_{\alpha\beta}\cos\omega t+\sigma^D_{\alpha\beta}\sin\omega t\right],\qquad
\sigma^C_{\alpha\beta}=\sum_{n m}\int_{\mathbf k} f_n Z_{nm}^{\alpha\beta},\qquad
\sigma^D_{\alpha\beta}=\sum_{n m}\int_{\mathbf k} f_n\frac{2\hbar\omega}{\varepsilon_{mn}}Q_{nm}^{\alpha\beta}.
\]

**Bookkeeping:** our code uses `omega_eV = hbar*omega` as energy, including the paper's `hbar=1` setting. There is no experimental broadening or fitted prefactor. The **published plotted quantity** is `-g mu_B sigma`, not unadorned `sigma`.

## 2. Antiunitary selection rules and the printed Eq. (14)

`j` and `B` are both T-odd. Thus in the clean intrinsic nonresonant response:

- The in-phase `B(t)` channel is T-even as a material response: `sigma_C -> +sigma_C`.
- The quadrature `dB/dt` channel is T-odd as a material response: `sigma_D -> -sigma_D`, since `dB/dt` is T-even.

Under point matrix `R`, the polar-axial rank-two tensor transforms as `det(R) R sigma R^T`. For antiunitary `RT`, apply the channel-specific parity above. This corresponds to SM §V assigning separate Jahn symbols **eV² for C** and **aeV² for D**, as opposed to the **common antiunitary minus sign printed for both in main Eq. (14)**. Both source versions are retained in `tensor_symmetry.transform(protocol=...)`, not silently conflated.

### Independent numerical audit

For paper Dirac `H_+(k)=d_x sigma_x+d_y sigma_y+d_z sigma_z`, define `H_-(k)=d_x sigma_x+d_y sigma_y-d_z sigma_z`. Exactly,

\[
H_-( -\mathbf k)=\sigma_y H_+(\mathbf k)^*\sigma_y.
\]

This is a true time-reversed *Hamiltonian*, not a second implementation of the same geometric formula. In generated paired matrices, `Z_-( -k)=Z_+(k)` and `Q_-( -k)=-Q_+(k)` to ≲10⁻¹⁵. Full finite-frequency Kubo also has even real and odd imaginary parts across the partners. At `mu=+0.85 eV`, `N=101²`: `C_xy(H+) = C_xy(H-) = -0.080552331572...` natural units; at `mu=0`, `N=101²`, `D_xx/omega = ±0.0152340780358` natural units. This is a decisive **channel-resolved protocol check**, but it does not by itself establish author intent behind printed Eq. (14).

### MPG / Neumann checks against the supplied SM tables

A separate nine-component real linear constraint solver yields:

| Generators | C | D | Source comparison |
|---|---|---|---|
| 4mm1′: C₄z, mₓz, pure T | 1 component, `C_xy=-C_yx` | none | SM Tables I/II: consistent |
| 4m′m′: C₄z, mₓz T | 1 component, `C_xy=-C_yx` | 2 independent, `D_xx=D_yy`, `D_zz` | SM Tables I/II: consistent |
| inversion P | all zero | all zero | paper: both inversion-odd |

Applying main Eq. (14) literally to pure T would incorrectly forbid the allowed Rashba C channel in 4mm1′. **Full 122 magnetic point group re-enumeration is still NOT STARTED.**

## 3. Frequency: SM §VI resolves the intended low-THz photon energy

Main Fig. 2 caption prints `omega=10^12 Hz`, which conflates angular omega and ordinary Hz. However SM §VI explicitly estimates `hbar*omega=4×10^-3 eV` for a THz magnetic field. This matches a **cycle frequency** `f=10^12 cycles/s` with `omega=2*pi*f`:

| Input convention | photon energy |
|---|---:|
| `f=10^12 Hz`, `omega=2πf` (supported by SM VI) | 0.004135667694 eV |
| `omega=10^12 rad/s` (counterfactual) | 0.000658211957 eV |

This resolves the **intended physical frequency interpretation**, but does **not prove** what was actually passed to the numerical plotting script.

## 4. Absolute SI, independently cross-checked on Rashba

The source model has `k` in Å⁻¹ and `∫d²k/(2π)²`. Its `sigma_C/D` have units Å⁻¹ once Eq. (12)'s `hbar*omega / delta_epsilon` is accounted for. To recover the ordinate **`-g mu_B sigma` in A/(m·T)**, multiply the raw kernel by

\[
\left[\frac{e\,g\mu_B}{\hbar}\right] \frac{1}{10^{-10}\mathrm m}
=281.794032638\,\frac{\mathrm A}{\mathrm{m\,T}}\,/\mathrm{\AA}^{-1}\qquad(g=2).
\]

From PRL Eq. (16) (`lambda_R=0.6 eV Å`, `m=1.1 m_e`), the Rashba conduction plateau is **3.88457 A/(m·T)**. Independent *direct SI* expression `e*g*mu_B*lambda_R*m/(2π*hbar³)` agrees (3.88457). Main Fig. 1(b) visually approaches about **4 A/(m·T)**. This is a source-to-SI cross-check without fitted factors.

With exactly the source Dirac parameters (`beta=.25 eV`, `v_F=.825 eV Å`, `M=.1 eV Å²`, `a=1Å`, `g=2`, `f=1THz`):

| Dirac case | Numerical predicted `-gmu_B sigma` | Source Fig. 2 visual estimate |
|---|---:|---:|
| Conduction `C_xy`, `mu=1.00 eV`, N=401² | **−33.49 A/(m·T)** | graph peak around −35 A/(m·T) |
| Displacement `D_xx`, `mu=0`, N=401² | **+0.017754 A/(m·T)** | plateau around **+0.04 A/(m·T)** |
| Displacement `D_xx`, `mu=1.00 eV`, N=401² | **−0.064186 A/(m·T)** | trough around **−0.13 to −0.14 A/(m·T)** |

The reference estimates above are **approximate visual readings**, *not* author data or precise image digitization. Differences are large enough to remain notable despite this uncertainty. In particular the in-gap plateau discrepancy is about **2.25×** (0.040/0.01775), while the conduction absolute scale is already consistent. Do not rescale the Q term by two to make the graph look closer; that would violate the presently implemented Eqs. (5), (8), (12).

### Numerical convergence excludes the obvious culprits

At `mu=0`, changing the BZ mesh `N=101,201,301,401` shifts the in-gap SI D only at the ≲10⁻¹³ numerical level. The **full interband Kubo** at 1 THz gives `+0.01775610`, within **0.013%** of the low-frequency `+0.01775389`; thus omitted higher-frequency terms also cannot cause the ≈2× source-curve difference. At `mu=1eV`, D changes by less than 0.2% across the tested grids (not proof of all-point convergence at sharp Fermi edges).

**Open physical/numerical question:** Why does source Fig. 2(d) use roughly twice the amplitude implied by the main-text Eq. (12), source frequency implied by SM VI and the SI prefactor independently anchored by Fig. 1(b) and Fig. 2(b)? Further checks should focus on author's plot-specific frequency/Zeeman convention, matrix-element normalization, and source data rather than fitting scales.

## 5. Runnable evidence

- `scripts/audit_round2.py`, `results/raw/round2_si_audit.json`, `results/raw/round2_si_symmetry.csv`, `results/raw/round2_SI_convergence.csv`, `results/figures/round2_TR_and_SI.png`
- `src/igmc_lab/core.py` (`mass_sign`, explicit frequency and SI), `src/igmc_lab/tensor_symmetry.py`
- `tests/test_round2.py`, `tests/test_symmetry_tensors.py` and first-round tests
- new model operations only; existing NMEE upstream and original first-round folder are not modified.
