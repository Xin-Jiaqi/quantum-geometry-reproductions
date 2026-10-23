# AgF₂ SOC-on first-principles evidence: independent magnetic-group and response-feasibility review

**2026-10-23 | Review of commit fee4e629aed4ac47cc1b234ed23673d503214a55 (PR #38)**  
**Verdict: SOC_SCF_SCREEN_PASS / m'_x_MOMENT_MAPPING_PASS_AT_OUTPUT_PRECISION / RESPONSE_NOT_COMPUTED / SOURCE_EQ14_UNRESOLVED / RESOURCE_CAP_EXCEEDED.**  
This is **read-only**: no VASP, PBS, Wannier, re-SCF, magnetic rotation, MAE, B2 or response run.

## 1. Actual input/output evidence (and what we independently checked)

The GitHub run receipt reports one local WSL VASP 6.6.0 \`vasp_ncl\` SOC-on single point from the **already screened** C2DB-derived Ag₄F₈ B1 snapshot, \`12×8×1\` Gamma, 64 iterations to \`EDIFF=1E-6\`:

| Quantity | SOC-on run | Interpretation |
| --- | --- | --- |
| TOTEN | −35.55001581 eV / 12-atom cell | Source-reported converged SCF; versus non-SOC −35.42794489 eV, −0.12207092 eV at fixed snapshot. **Not MAE and not a stability comparison between different candidate structures; moreover this is not a clean isolated SOC energy because the two input protocols differ in more than LSORBIT (e.g. noncollinear treatment, symmetry and GGA_COMPAT choice).** |
| Occupation gap | 1.3616 eV (HOMO −5.7746, LUMO −4.4130 eV) | Screening diagnostic on limited k grid and finite smearing; not independent band-path / k-converged indirect-gap proof |
| Ag spin moments | dominant \(\pm0.568\,\mu_B\) along x | Type-I seed survives SOC; y about ±0.019, z about −0.005 on Ag |
| Ag local PAW orbital moment | \(|L|\simeq0.143\,\mu_B\), \(L_z\approx+0.069\ldots+0.070\,\mu_B\) | PAW-sphere angular-momentum projection; NOT orbital-B perturbation matrix and not orbital IGMC |
| Reported total spin net | \((0,-0.001,-0.0378)\,\mu_B\) | Small out-of-plane component. Projected 12-ion spin sum is \((0,0,-0.035)\,\mu_B\); the −0.0028 difference is consistent with non-site-projected/interstitial contributions or rounding and is not silently replaced |
| Source/origin | private OUTCAR, OSZICAR, vasprun SHA and sanitized excerpts in receipt | I checked committed CSV/JSON; I did **not** re-open the private WSL original logs, so raw-output cross-validation remains Harness-reported |

**Important VASP distinction:** \`LORBMOM=.TRUE.\` calculates orbital angular momentum **inside PAW augmentation spheres**, not the full gauge-covariant orbital-magnetization or orbital Zeeman operator. Site-projected total \(\Sigma_{12}L_z=+0.246\) and Ag-only \(\Sigma_4 L_z=+0.278\) (values from rounded CSV) must **not** be presented as a bulk physical \(\mathbf M_{\rm orb}\), experimentally measured net magnetization, or orbital IGMC.

## 2. Critical correction: SOC canting and common orbital Lz **do not break m'_x**

The earlier \`RUN_RECEIPT.json\` and \`SOC_STATIC_VERDICT.md\` suggested that the ~2% canting and nonzero \(M_z\) necessarily break exact \(\,m'_x\,\). This is **not supported** by the reported moment pattern.

For \(\,g_x=\{m_x|(1/2,1/2,0)\}\,\), the submitted POSCAR was previously shown to obey exact **geometric** site mapping to ~\(1.28\times10^{-5}\,\mathrm Å\):

- Ag pairs **1↔2, 3↔4**;
- F pairs **5↔8, 6↔7, 9↔12, 10↔11**.

Magnetic **spin and orbital moments** are both time-odd axial vectors. With \(R_x=\operatorname{diag}(-1,1,1)\), \(\det R_x=-1\), the antiunitary operation transforms them by

\[
\mathbf M_j=(-1)\det(R_x)\,R_x\mathbf M_i
=\operatorname{diag}(-1,+1,+1)\,\mathbf M_i.
\]

Thus for \(m'_x\), **\(M_x\) must reverse between mapped partners, but \(M_y,M_z\) are allowed to keep the same sign**. A common \(M_z\) or \(L_z\), or the observed layer-dependent \(M_y\), is **not automatically symmetry breaking**.

**Independent calculation from committed \`MOMENTS.csv\`:**

| Pair | \(S_x\) expected / observed | \(S_y,S_z\) expected / observed | Orbital L vector |
| --- | --- | --- | --- |
| Ag 1→2 | −(−0.568)=+0.568 / +0.568 | +0.019,−0.005 / +0.019,−0.004 | exact at printed precision |
| Ag 3→4 | −(+0.568)=−0.568 / −0.568 | −0.019,−0.005 / −0.019,−0.005 | exact |
| F 5→8, 6→7, 9→12, 10→11 | spin x reverses | spin y,z conserved | exact |

For all 12 ions, max absolute partner deviation in the three **spin** components is **0.001 μB**; for **orbital** moments it is **0.000 μB** at the reported 0.001 μB digit precision. Site forces, magnetization density grid and symmetry of the full spinor density were **not** independently tested.

**Correct verdict:** \`MAGNETIC_mx_PRIME_CONSISTENT_WITH_ALL_REPORTED_SITE_VECTORS\`. Full microscopic MPG classification remains provisional because of finite output precision, \`ISYM=-1\`, k mesh and possible density-sector deviations, **not because this canting contradicts \(m'_x\)**. No forced symmetrization should be performed merely to remove allowed weak canting.

**A useful additional switchable proxy:** out-of-plane magnetization \(M_z\), local orbital \(L_z\), and AHE pseudo-vector \(z\) have the same **axial, time-odd** partner parity: **Type-I unitary \(m_z\): even; Type-II antiunitary \(m'_z\): odd**. Under the project's channel-resolved convention, insulating displacement IGMC \(D_{xx}\) is instead **odd for Type-I and even for Type-II**. This is a potentially measurable crossed-parity fingerprint, but the SOC B1-only run has not supplied any B2 moment, AHE or IGMC value.

## 3. Channel-resolved response-theory gate (not yet passed)

The relevant paper source map labels Xiang PRL (2025) main-text Eq. (11) conduction and Eq. (12) displacement contributions; Eq. (13) shows a Fermi-surface form of C and Eq. (14) gives a magnetic point-group transformation. The **source-mapping file explicitly flags an unresolved** Eq. (14) convention conflict with channel-by-channel P/T classification.

The reproducible \`src/igmc_lab/tensor_symmetry.py\` in the repository implements both named alternatives:
- \`protocol="paper_printed_eq14"\`: antiunitary parity sign − for both C and D;
- \`protocol="channel_resolved_symmetry"\`: antiunitary parity + for C, − for D.

The two predict different C tensor allowedness under \(\,m'_x\). Consequently, **the conduction IGMC selection rules and crossed-parity story must remain conditional** until this source-level issue is resolved from the primary text and Supplemental Material, including frequency conventions and whether one considers conductivity to \(B\), conductivity to \(\dot B\), or real phase-resolved quadratures. Do not silently choose the appealing one.

Under the repository's **declared, conditional channel-resolved convention**, the \(m'_x\) matrices are:

\[
C=\begin{pmatrix}0&C_{xy}&C_{xz}\\C_{yx}&0&0\\C_{zx}&0&0\end{pmatrix},
\qquad
D=\begin{pmatrix}D_{xx}&0&0\\0&D_{yy}&D_{yz}\\0&D_{zy}&D_{zz}\end{pmatrix}.
\]

\(D_{xx}\) is allowed, **not necessarily nonzero**. In a clean zero-temperature band insulator with the chemical potential inside the gap, main-text Eq. (13) identifies C as a Fermi-surface term, so its **bulk conduction IGMC vanishes**, even if symmetry allows C components. AHE need not vanish simply because the state is insulating, but the actual topological/transport regime was not established in this screening.

### Next Hamiltonian/operator contract for an eventual material IGMC

- **Spin-Zeeman**: \(H_B^{\rm spin}=-\mathbf B\cdot\hat{\mathbf M}^{\rm spin}\), with the physical spin g factor, \(μ_B\), and same-gauge off-diagonal spin matrix elements \(S_{nm}^j(\mathbf k)\). Site \`MAGMOM\` initialization and PAW projected moments are not these matrices.
- **Orbital magnetic driving**: a separate field-coupling operator \(\hat M^{\rm orb}_j\), including the periodic-crystal orbital-field subtleties; \`LORBMOM\` site projections **do not** supply it. Keep orbital IGMC separate from Xiang's spin-only response, see Ghorai et al. DOI 10.1103/xr5v-gkv5.
- **Position/velocity**: full same-gauge \(v_i=(1/\hbar)D_{k_i}H\) and position/Berry-connection matrix elements, including Wannier embedding/disentanglement contributions as appropriate. \`wannier90_hr.dat\` alone cannot give the full physical IGMC.
- **Mixed quantum geometry**: with \(m\ne n\), one can independently form interband \(\langle u_m|\partial_{B_j}u_n\rangle=\langle u_m|\partial_{B_j}H|u_n\rangle/(\epsilon_n-\epsilon_m)\) and \(\langle u_m|\partial_{k_i}u_n\rangle\) in a carefully specified smooth/covariant gauge. The real/imaginary parts generate metric/curvature-like mixed terms, but need proper degeneracy handling, occupied-band summation, and clear sign and unit conventions.
- **Independent cross-check**: a finite-Zeeman-field **non-Abelian Wilson-loop polarization derivative** on the same insulating toy/matched Wannier model should agree with the static low-frequency displacement kernel, after frequency factors and field units are fixed.
- **Observable**: for \(B(t)=B_0\cos\omega t\), a pure polarization displacement signal scales as \(J^D\propto\omega\). Account for Faraday-induced E fields and distinguish driven displacement from DC transport.

The current SOC single point has no WAVECAR/CHGCAR saved (\`LWAVE=.FALSE.\`, \`LCHARG=.FALSE.\`), no Wannier spinor Hamiltonian, and no k-resolved position/spin operators. **Therefore one cannot extract a quantitative AgF₂ IGMC coefficient or even infer a lower bound from these outputs.** This is a positive screening, not a response calculation.

## 4. Resource authorization discrepancy: must not recur

The run contract approved **8 MPI ranks × 4h (32 rank-hours)**. The result is **8 ranks × 4.658h = 37.27 rank-hours**, exceeding the cap by **0.658h (5.27 rank-hours)**. The executor disclosed this transparently, but the receipt's \`no_compute_beyond_scope=true\` / \`stayed_within_scope=true\` must not be taken as an assertion that the clock budget was respected. These fields are internally inconsistent with the clock evidence if interpreted as including the resource cap.

**Proposed control for any later job (not authorization):** use a verified supervisor enforcing an absolute wall-clock deadline; reserve sufficient time to request an early graceful stop and escalate to controlled termination at the hard deadline. If monitoring connectivity drops, the hard stop must still be enforced locally. No “only another 40 minutes” discretionary extensions, no automatic PBS fallback or retries. If the physical calculation cannot finish in budget, preserve a PARTIAL receipt and seek renewed authorization.

Avoid inferring total 8-rank RAM from rank0's 0.93 GB peak alone; need actual per-rank/whole-job peak to claim memory safety if that matters to a later run.

## 5. Scientific decision and concrete next deliverables (NO VASP)

**Do not move yet to expensive SOC/Wannier or B2/NEB.** The highest-information free next step is:

1. Run a deterministic **site-wise magnetic symmetry audit** on the existing CSV and POSCAR at several tolerances, recording the exact \(\,m'_x\,\) transformation and discrepancies. The provided 12-site table already supports this strongly.
2. Directly inspect the already-archived, privately held **7-page main text + 11-page SM** for Xiang Eq. (14) and SM Tables I–III, reconcile the C/D antiunitary tensor convention **without postulating an error**. Treat the two named protocols separately if unresolved.
3. Read existing GitHub model and operator code for spin-Zeeman derivative, finite-field Wilson-check availability, gauge covariance, and all missing material Wannier operator files. Return a testable \`RESPONSE_FEASIBILITY_MATRIX.md\` with statuses \`EXISTS_AND_TESTED / EXISTS_UNVERIFIED / MISSING\`.
4. Rank project stories:
   - **primary hypothesis:** Type-I B1/B2 sliding switches finite insulating \(D_{xx}\) while the weak out-of-plane axial response is even;
   - **alternative under strained Type-II:** \(D_{xx}\) even but \(M_z\)/AHE odd;
   - **fail path:** SOC gap/target MPG survives but quantitative D response vanishes or source convention invalidates the proposed fingerprints; demote AgF₂ and do not spend further resources.

**Go/no-go for costly work:** a SOC state preserving the target symmetry is *necessary, not sufficient*. A material response calculation requires a source-validated response protocol and feasible same-gauge operator chain, and a realistic projected signal or clear reason to expect it.

### Sources

- Xiang et al., *Phys. Rev. Lett.* **134**, 116301 (2025), DOI 10.1103/PhysRevLett.134.116301. Internal Eq./SM inventory: \`papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/docs/SOURCE_MAP.md\`.
- Zhu et al., *Phys. Rev. Lett.* **135**, 056801 (2025), DOI 10.1103/dmzg-ck2t.
- Ghorai, Sarkar & Agarwal, orbital IGMC (2026), DOI 10.1103/xr5v-gkv5.
- VASP LORBMOM documentation: https://vasp.at/wiki/LORBMOM — PAW-projector orbital moment only.
- VASP LSORBIT documentation: https://vasp.at/wiki/LSORBIT — Cartesian orbital moments versus SAXIS spinor coordinates.
- Original data: \`RUN_RECEIPT.json\`, \`MOMENTS.csv\`, \`BANDS_GAP.json\` of this same directory. Earlier source-C/D convention: \`src/igmc_lab/tensor_symmetry.py\`.

**No VASP/PBS job was launched; no original scientific output was rewritten.**
