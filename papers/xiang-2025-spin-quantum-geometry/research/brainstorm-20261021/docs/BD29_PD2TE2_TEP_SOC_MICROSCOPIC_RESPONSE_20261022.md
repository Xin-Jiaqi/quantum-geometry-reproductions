# BD29 | Pd2Te2-derived 24-band spin-orbit model: first nonzero current-driven layer spin

**2026-10-22. Status: MODEL NONZERO PASS for layer-odd Sz, parity/model controls PASS, layer-even magnitude not converged, MATERIAL & NOVELTY UNVERIFIED. No DFT or external HPC.**

## 1. Source-derived vs chosen quantities

**BiDB facts:** id 980 Pd2Te2 is an 8-atom stable-sliding metallic/nonmagnetic-tagged D3h (p-6m2) structure. Its actual atom coordinates obey C3z, vertical mirror and horizontal layer-exchanging mirror to ~0.0013 Å. BiDB itself is user supplied, read-only, and not uploaded to GitHub. Details in [BD29 direct BiDB record report](BD29_BIDB_DB_DIRECT_SCREEN_20261022.md).

**Model choices, NOT fitted material parameters:** retain only four Te atoms and Cartesian p-orbitals with spin (24 bands). Pd-d orbital contributions are completely OMITTED. Project the slightly sheared BiDB relaxed cell (difference ~0.002997 Å) onto exact hexagonal lattice; symmetrize equivalent Te positions. Geometry-dependent spin-independent Slater–Koster bonds:
\[
t({\bf d})= e^{-(|{\bf d}|-3\mathrm{\AA})/0.8\mathrm{\AA}}
[t_\pi I_3+(t_\sigma-t_\pi)\hat{\bf d}\hat{\bf d}^{\mathsf T}],
\quad t_\sigma=1,\ t_\pi=-0.23\ \mathrm{eV},
\]
with 7.2 Å cutoff and Bravais hopping phase exp(i k·d). Te onsite atomic SOC is \(\lambda{\bf L}\cdot{\bf S}\), \(\lambda=0.48\) eV; inner/outer chalcogen onsite energies 0.25/−0.12 eV. An illustrative chemical potential \(\mu=1.5\) eV and electronic thermal broadening \(T=0.13\) eV are used. These are NOT source-obtained Pd2Te2 hopping/SOC/charge-density fits.

**No ad-hoc spinful hopping vertex** was added to force the desired layer Sz. The SOC is ordinary onsite atomic L·S; registry dependence is solely the real bond geometry in spin-independent p-p SK hopping.

**Transport definition:** for \(S_\pm^z=\sum_{\mathrm{top Te}}S_i^z\pm\sum_{\mathrm{bottom Te}}S_i^z\), a band-diagonal, constant-relaxation-time intraband/Boltzmann response at fixed current is
\[
Q^{(\pm)}_a=\sum_{nk}(-f'_n)\langle S_\pm^z\rangle_{nk}v_{nk,a},\
D_{ab}=\sum_{nk}(-f'_n)v_{nk,a}v_{nk,b},\
\chi_{\pm,b}=\sum_a Q^{(\pm)}_a(D^{-1})_{ab}.
\]
This uses the FULL 2×2 conductivity matrix to distinguish fixed-j from fixed-E response; common relaxation time/charge constants are suppressed. Numbers are in ARBITRARY MODEL UNITS and cannot be converted to experimental spin density without a calibrated model, constants and transport lifetime.

## 2. Validation, including a caught mistake

Full **orbital p + spin-1/2 + site-permutation** Hamiltonian operations were tested for arbitrary in-plane k and relative layer slide u:
\[
H(gk,gu)=U_gH(k,u)U_g^\dagger.
\]
Horizontal mz swaps Te sites and sends u to −u; vertical mirror and C3 rotate u geometrically. The measured maximum relative Hamiltonian covariance errors are:

| Full symmetry | Relative error |
|---|---:|
| C3z | \(5.51\times10^{-16}\) |
| vertical mx | \(1.58\times10^{-16}\) |
| horizontal layer swapping mz | 0 |
| antiunitary spin time reversal on H | 0 |
| time reversal on the two physical velocity matrices | 0 |

A raw relaxed BiDB lattice initially caused C3 error ~0.0023; **after documenting the ~0.003 Å exact-hex projection**, the whole orbital/SOC Hamiltonian passed to machine precision. It would have been incorrect to call the raw-cell model exactly C3 invariant.

Another independent error: a Gamma-centered mesh samples exact band degeneracies. For the **band-diagonal** Boltzmann expression, individual spin/velocity expectations can depend on the diagonalization gauge inside a degenerate eigenspace. It generated spurious tiny or even sizeable forbidden response components. The accepted calculation uses a generic shifted mesh avoiding those exact high-symmetry points and an independent 6-member C3v k-orbit quadrature. It verifies both the u→−u layer parity and SOC-off negative controls. This is not a proof the Boltzmann formula is reliable for all accidental near-degeneracies; fully degenerate/non-Abelian treatment remains an open robustness gate.

## 3. First NONZERO spin-response result

For a physical relative slide \(u_x\) in Å and \(u_y=0\), output is a row covector \((\chi_{\pm,j_x},\chi_{\pm,j_y})\).

| \(u_x\) Å | layer-odd \(\chi_{-,j_y}\) | layer-even \(\chi_{+,j_x}\) |
|---:|---:|---:|
| −0.100 | +0.006511314 | −0.000101975 |
| −0.050 | +0.003245626 | −0.000034302 |
| 0 | <1e−15 | <1e−15 |
| +0.025 | −0.001621678 | −0.000009065 |
| +0.050 | −0.003245626 | −0.000034302 |
| +0.100 | −0.006511314 | −0.000101975 |

These 16×16 shifted-grid results obey **layer-odd χ(−u)=−χ(+u)** and **layer-even χ(−u)=χ(+u)** with errors <1e−14; SOC=0 gives numerically zero both channels. This is a **nonzero microscopic counterexample to the concern that the symmetry-allowed Sz coefficient necessarily vanishes**. Note that the older LEE Eq. (3) simplified Hamiltonian gave Sz=0; this p-orbital SOC model does not have that accidental in-plane-helicity conservation.

The first-order layer-odd coefficient is near \(\gamma\approx-0.065\) in model response units per Å; changing chemical potential and SK parameters will change this magnitude and possibly sign. The layer-even channel is even and small; although data qualitatively support quadratic onset, **its numerical coefficient β is not yet quantitatively converged**.

### k-mesh convergence, u_x = +0.10 Å (same model parameters)

| Nk×Nk | χ layer-odd Sz/jy | χ layer-even Sz/jx |
|---|---:|---:|
| 8×8 | −0.006582671 | −0.000062307 |
| 12×12 | −0.006524787 | −0.000085770 |
| 16×16 | −0.006511314 | −0.000101975 |
| 20×20 | −0.006507386 | −0.000113395 |
| 24×24 | −0.006506043 | −0.000120108 |

Layer-odd changes about 0.08% between 16×16 and 24×24. The small layer-even number changes appreciably and cannot yet be promoted to an actual winding prediction.

**Exact 11-point raw computed data (not simulated by interpolation):** [BD29 p/SOC current response CSV](../data/BD29_PD2TE2_TEP_SOC_BOLTZMANN_20261022.csv).

## 4. What this does and does NOT establish

**Established at the toy-model level:**
1. A realistic BiDB **source structure** with genuine layer-swapping geometry has an explicitly SOC-covariant p-orbital tight-binding realization.
2. Ordinary p–p geometry-dependent hopping PLUS atomic L·S is sufficient to give nonzero current-induced layer-odd out-of-plane spin accumulation on a displaced bilayer, with predicted odd-in-slide symmetry.
3. The leading layer-odd response is numerically stable, survives full fixed-current conductivity conversion and vanishes at SOC=0.

**NOT established:**
1. Actual Pd2Te2 quasiparticle dispersion, Fermi level, orbital character, relaxation time, calibrated spin accumulation or experimental feasibility. Parent monolayer e-hull ~0.134 eV/atom; Te-only model omits important Pd-d bands.
2. A new spin mechanism. An ordinary symmetry-allowed Edelstein response in a shifted SOC bilayer may already follow the 2026 *Layer Edelstein Effect* framework. The specific continuous-registry global zero/index and activation-order classification remains a distinct **THEORETICAL QUESTION**, not a proven original theorem.
3. A measurable layer-even winding −2. The selected weak even coefficient is slower to converge than the robust layer-odd response.
4. A full magnet, Néel torque, domain wall, phonon, SOC Wannier or DFT calculation.

**Next high-value test:** using the same data, improve the layer-even near-G convergence (including careful degeneracy handling) and verify the full hex-registry F±(δ) pattern with a source-validated spinful Pd-d + Te-p Hamiltonian if publicly available. Independently inspect the complete 2026 LEE PRL SM for continuous-registry zero/index predictions; this is the **direct general-theory scope check**. Do not undertake DFT/HPC unless a more specific, high-value test survives these gates.

**Code:** [BiDB local-only geometry/SOC/Boltzmann full symmetry script](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/pd2te2_te_p_soc_boltzmann_gate.py), plus [BD29 database record gate](BD29_BIDB_DB_DIRECT_SCREEN_20261022.md). The raw user database stays private.


## Interlayer hopping off negative control

The local 24-band implementation was re-evaluated after deleting every Te p-p bond between the upper and lower layers, retaining atomic SOC and intralayer hopping. At ux=+0.10 Angstrom, the current-induced layer-odd Sz/jy coefficient falls from -0.006582671 (8x8) to 4.21e-16, and from -0.006524787 (12x12) to 1.07e-16. The layer-even spin signal also falls below 1e-15. Hence the model's nonzero response requires both atomic SOC and interlayer hybridization. This is a model-level control, not an original physical mechanism. The research-branch code now includes an explicit decoupled-layer test.
