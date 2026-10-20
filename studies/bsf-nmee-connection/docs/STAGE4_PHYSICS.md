# BSF × NMEE stage 4: orbital-channel interference and a mixed-curvature slope

## One central result

The stage-3 **gapped sliding-induced sign reversal** is not simply a generic consequence of breaking layer symmetry or turning on interlayer tunneling. In a fixed, explicitly specified eight-band Kane–Mele bilayer, the sign at small positive sliding depends critically on the **relative phase and interference between same-sublattice and cross-sublattice interlayer tunneling**. This is a model mechanism, not a general BSF selection rule.

For the same physical parameters (`t⊥=20 meV`, `λR=20 meV`, `λSO=10 meV`, zero interlayer bias, `T=0`, fixed insulating chemical potential, `Bz` central-difference step `0.25 T`), define four independently evaluated Hamiltonians: `full`, `same-only` (AA/BB), `cross-only` (AB/BA), and `decoupled`. The inclusion-exclusion nonadditivity is

\[
\Gamma_\mathrm{mix}(u)=\Gamma_\mathrm{full}(u)
  -\Gamma_\mathrm{same}(u)-\Gamma_\mathrm{cross}(u)+\Gamma_0(u).
\]

Here \(\Gamma_0\) is zero to numerical precision by the unbroken layer-exchange mirror at zero coupling. Unlike matching only the total hopping norm, the controls **preserve the magnitude of each retained hopping channel**, while removing the others.

| Shift \(u\) | Full | Cross only | Same only | Nonadditivity |
|---:|---:|---:|---:|---:|
| 0.10 | −1.69808×10⁻⁵ | +1.23519×10⁻⁶ | +4.03×10⁻⁹ | **−1.82200×10⁻⁵** |
| 0.20 | +9.32171×10⁻⁷ | +8.97465×10⁻⁶ | +7.83×10⁻⁹ | −8.05031×10⁻⁶ |
| 0.40 | +9.85021×10⁻⁵ | +5.32965×10⁻⁵ | +1.39×10⁻⁸ | **+4.51918×10⁻⁵** |

All values refer to \(\Gamma_{xz;x}\), in the code's static-susceptibility normalization, not a measured NMEE coefficient. At \(u=0.10\), cross-only and same-only responses are positive, whereas the complete Hamiltonian's response is negative: **there is a substantial nonadditive interference contribution**. The same inclusion-exclusion term changes sign before \(u=0.40\). This is a stronger explanatory clue than simply finding a large signal after adding an orbital-hybridization term.

### Relative hopping-sign counterfactual

Changing **only** the sign of all AB/BA interlayer tunneling amplitudes, while keeping AA/BB amplitudes and monolayers fixed, reverses the negative low-sliding response: at \(u=0.10\), the full physical model gives \(-1.69808\times10^{-5}\) whereas this altered Hamiltonian gives \(+2.96566\times10^{-5}\). At \(u=0.40\), both are positive (\(+9.85021\times10^{-5}\) and \(+1.78216\times10^{-4}\)). This counterfactual does **not** by itself represent a physical sliding path, nor is it a mere k-dependent change of eigenvector phase. It supports the statement that the *relative phase* of interlayer orbital couplings, not just their norm, is causally important in this model. Electronic dispersions can change too, so it is not yet a pure gauge-invariant geometric selection rule.

## Underlying geometry: a field derivative of mixed Berry curvature

Write \(H(B_z,h_x)=H_0-g\mu_B B_z S_z-h_x S_x\), with \(v_x=\partial_{k_x}H\) and Pauli spin convention. For a gapped zero-temperature insulator, the **interband static** spin–electric susceptibility kernel used here can be written exactly as

\[
\chi_{x;x}(\mathbf k,B_z)
=-2\sum_{v\in\mathrm{occ},c\in\mathrm{emp}}
  \frac{\operatorname{Im}[(S_x)_{vc}(v_x)_{cv}]}
       {(\epsilon_v-\epsilon_c)^2}.
\]

Up to the adopted sign and normalization of the auxiliary Zeeman parameter, this is the mixed \((k_x,h_x)\) Berry curvature of the occupied projector. The investigated response is its **\(B_z\) derivative**:

\[
\Gamma_{xz;x}^{\mathrm{static}}\ \propto\
\left.\partial_{B_z}\int_{\mathrm{BZ}}
\Omega_{k_x,h_x}(\mathbf k,B_z)\,\frac{d^2k}{(2\pi)^2}\right|_{B_z=0}.
\]

This formulation is independently implemented in `bsf_nmee/response_mechanism.py`; it agrees with the earlier finite-field `chi_k` evaluation to ~\(10^{-18}\) at the checked points. It is **not** a derivation of the authors' unresolved *driven DC NMEE* formula and does not imply that a bulk equilibrium charge-polarization measurement would necessarily produce this response.

For intuition only, the finite-\(B_z\) product rule partitions the numerical derivative of the mixed-curvature summand into **matrix-element variation** and **energy-denominator variation**. At \(u=0.10\), the former contributes \(+1.94\times10^{-5}\) and the latter \(-3.64\times10^{-5}\); at \(u=0.40\) the respective values are \(-5.96\times10^{-5}\) and \(+1.58\times10^{-4}\). Their sum is the physical kernel. **These two pieces individually are not a unique gauge-invariant decomposition**, particularly near degenerate bands: no standalone universal "quantum geometry changed sign" claim follows from them. The exact, gauge-invariant *total* and the orbital-channel counterfactuals are the evidence.

## Symmetry explains the leading power, not the sign

The layer-sign unitary \(W=\mathrm{diag}(+I_4,-I_4)\) maps \(H(t_\perp)\to H(-t_\perp)\), leaving total spin and external Zeeman operators covariant. Therefore \(\Gamma(t_\perp)=\Gamma(-t_\perp)\). In an **analytic, gapped** weak-coupling region,

\[
\Gamma(\tau,t_\perp)=t_\perp^2\,C_2(\tau)+O(t_\perp^4).
\]

This explains why a small interlayer coupling produces a *quadratic* total response, while the physical mirror/layer-exchange symmetry imposes \(\Gamma(-\tau)=-\Gamma(\tau)\). Neither symmetry fixes the sign of \(C_2(\tau)\) throughout a sliding trajectory. A model-specific zero crossing can therefore occur **without closing an insulating gap**. Direct calculations at \(u=0.10\) with \(t_\perp=2.5\) and 5 meV give \(\Gamma/t_\perp^2\approx-0.02698\) and \(-0.02752\) in the same units, consistent with the approach to a quadratic regime. Exactly changing the sign of *all* hopping blocks leaves the computed total response unchanged.

## Robustness beyond a single parameter set

At fixed \(t_\perp=20\) meV, the previously frozen zero is \(u\approx0.1974\). Reducing interlayer tunneling to 10 meV preserves the negative sign at \(u=0.10\), positive sign at \(u=0.40\), and a gapped path, while moving the zero to \(u\approx0.2676\) (coarse root bracket). The sampled minimum distance from the chemical potential rises from about 2.43 meV to 9.44 meV. The response zero therefore does **not require** a tiny gap or a particular 20-meV interlayer coupling.

Exploratory comparisons at weaker Rashba SOC, stronger intrinsic SOC, and longer tunneling decay length also retained endpoint sign reversal under sampled-gap tests, but their complete validated scans are **not** in the stage-4 release and should not be quoted as frozen benchmarks until independently archived.

## Scientific boundaries and next decisive test

This is an intentionally simple synthetic Hamiltonian. The sliding path is not a computed ferroelectric switching trajectory; the model does not determine energy minima, polarization, actual material bands or equilibrium interface geometry. The static diagnostic's electric-field response protocol is not the same as the paper's disputed driven DC coefficient. The exact numeric sign depends on monolayer SOC, hopping phase conventions and chemical potential; the signed channel modification is a controlled alternative Hamiltonian, not a demonstrated real-material operation.

**Next:** obtain a basis-independent perturbative expression for the \(t_s t_c\) interference coefficient in the gapped weak-coupling limit, and test whether a symmetry or orbital-winding selection rule constrains it across other monolayer Hamiltonians. This could turn model-specific interference into a PRL-relevant *general mechanism*; neither claim is yet established.
