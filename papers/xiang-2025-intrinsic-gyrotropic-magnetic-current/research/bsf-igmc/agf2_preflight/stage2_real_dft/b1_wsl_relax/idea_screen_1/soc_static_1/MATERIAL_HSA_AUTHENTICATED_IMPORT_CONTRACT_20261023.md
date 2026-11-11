# AgF₂ H/S/A materials-operator interface: fail-closed source and gauge contract

**2026-10-23 | ENGINEERING DESIGN ONLY / NO MATERIAL OPERATOR IMPORT, NO DFT OR WANNIER EXECUTED**

## Executive recommendation

The corrected four-band synthetic gate has earned **permission to implement a provenance-checked importer**, **not** permission to run new VASP/Wannier and not an AgF₂ IGMC calculation.

The next low-cost development stage must be **synthetic file serialization → strict reader → reconstruction of identical same-Wannier-gauge H(R), physical spin S(R), full position A(R) → reproduce Dxx/ω≈+0.0342937475** with source-independent and intentionally corrupted fixtures. Do not attempt to declare material readiness until actual Ab-initio matrix-element provenance is proven.

Crucially, the previous AgF₂ B1 SOC static run used LWAVE=.FALSE., LCHARG=.FALSE., and did **not** produce a Wannier interface/WAVECAR/spin matrices; site-projected LORBMOM and MAGMOM are not sufficient. This is a **missing upstream data** issue, not something that can be fixed with a Python parser.

## 1. Official actual file roles: source Bloch gauge ≠ final Wannier gauge

| File / method | Actual role | Required for material chain |
| --- | --- | --- |
| \`seedname_hr.dat\` | Matrix elements of the **Hamiltonian in FINAL Wannier basis** at R, together with an R-list and degeneracy weights (subject to write conventions) | YES |
| \`seedname_r.dat\` (requires \`write_rmn=true\`) | Matrix elements <m,0| r |n,R> in the **FINAL Wannier basis**; contains off-diagonal position terms, not just centers; has no independent degeneracy block | YES (or a separately authenticated mathematically equivalent full position connection) |
| \`seedname.spn\` (formatted preferred) | Physical Pauli spin matrices in the **ab-initio Bloch eigenstate basis**, typically produced by a compatible ab-initio interface, not yet aligned to the Wannier matrices | YES (or an authenticated same-gauge projected physical-spin alternative) |
| \`seedname_u.mat\` and \`seedname_u_dis.mat\` (or \`seedname_v.mat\` if present) | Actual \(\,U(k)\,\) and disentanglement \(\,U^{dis}(k)\,\); map from ab-initio Bloch eigenstates to FINAL Wannier gauge | YES whenever raw Bloch-gauge spn is used |
| \`seedname.mmn\` / \`seedname.amn\` | Neighbor-k overlaps / initial projections **in source Bloch gauge**; useful for Wannier provenance and overlap reconstruction, but **AMN is not the final U(k)** | YES as source provenance/input where relevant; NEVER a substitute for U(k) |
| \`seedname.eig\`, \`.win\`, \`.wout\`, complete k-list, \`wsvec.dat\` where applicable | Energy/k order, input settings, Wannier convergence, centers/spreads, Wigner–Seitz phase corrections | YES for reliable validation and matched physics |
| VASP \`LWRITE_SPN=.TRUE.\` with \`LWANNIER90=.TRUE.\`, \`vasp_ncl\` | VASP 6.6 can **if compiled with compatible Wannier90 interface** write formatted \`wannier90.spn\` directly for noncollinear spinors. VASP documentation restricts it to **default SAXIS**. | Candidate future producer **only**; must verify build and file generation. The existing SOC task used SAXIS=(0,0,1), but did not enable the interface |

**Do not assert that all requested source files already exist.** Read or generate nothing from private licensed VASP inputs until the user separately authorizes that specific future job.

## 2. Exact algebra and compatible spin normalization

Set \(\,V(k)=U^{dis}(k)U(k)\,\) (rectangular \`num_bands×num_wann\` when disentangled; use U alone when no disentanglement). Then

\[
S_j^{W}(k)=V^\dagger(k)\,S_j^{Bloch}(k)\,V(k),\qquad
H^W(k)=V^\dagger(k)\,{\rm diag}(\epsilon_{n k})\,V(k),
\]

with exact source **band order, band exclusions, k order, spinor axis, and matrix conventions**. Reconstruct R-space matrices by the matching discrete inverse Fourier transform and verify Fourier-back values at **every actual mesh k**.

- Source \`spn\` often stores \(\langle \sigma_j\rangle\), **dimensionless physical Pauli**, not angular momentum \(\hbar\sigma/2\), and not a magnetic moment in μB. The synthetic contract uses \(\sigma_j\) directly and \(H_b=H-b\sigma_j\) with b in energy units. For real \(\mathbf B\), the spin Zeeman coupling is \(\,H_B=-\mathbf B\cdot\boldsymbol\mu_s\,\), with a declared sign / \(\mu_B\) / g-factor choice (for conventional electron spin moment, the sign must be derived explicitly). Do **not** silently insert a factor 2 or interpret b as tesla.
- Under noncollinear VASP, document \`SAXIS\`, Cartesian-vs-spinor components and whether the code emitted Cartesian Pauli; the previous B1 static task used default \`SAXIS=(0,0,1)\`. Reject other axes until an explicit tested rotation exists.
- Hermitian \(\,S_j(k)\,\) in Bloch and Wannier gauge, and \(\|S_j^W(k)\|\le1+\epsilon\) for a projected physical Pauli. A mere input manifest asserting \`spin_kind=physical\` is **not** proof of its electronic origin.
- Wannier gauge rotations are only meaningful when H,S,A transform **jointly**; test one independent joint constant rotation and one k-dependent gauge with correct inhomogeneous connection law, or mark the latter **UNTESTED** and refuse broad gauge claims. The existing 4-band synthetic gate tested constant rotations only.

## 3. Position / Wannier-center / WS conventions: reject ambiguous inputs

Write a clear \`X(R)\` carrier convention with H in eV, position A in Å, physical spin dimensionless Pauli, 2D integer R=(Rx,Ry) and a 2×3 Å lattice row matrix. The actual Wannier file has 3-component lattice R, whose Rz must be constrained to 0 for the slab interpolation. Do not flatten Rz without recording the fixed slab setting.

The real-space \`r.dat\` supplies the **full** position matrix; the centers are the diagonal on-site position entries and should be included **once**. Reject adding centers a second time or discarding off-diagonal position elements.

**Subtle but critical**: \`hr.dat\` contains Wigner–Seitz degeneracies. \`r.dat\` has **no separate degeneracy block**. Depending on \`use_ws_distance\`, \`write_ndegen_applied\` and possibly \`wsvec.dat\`, treating \`hr.dat\` and \`r.dat\` as parallel identical weighted R lists can be wrong. The safest first supported real-data route should *explicitly require* \`write_ndegen_applied=true\`, check all degeneracy entries and R list identity, and reject complex WS-folding layouts until a dedicated reference-consistent mapping test is available. If \`write_ndegen_applied=true\` yields all degeneracy values 1, **do not divide them again or apply wsvec twice**.

\[
H(k)=\sum_R e^{i k\cdot R}H_R,\quad S_j(k)=\sum_R e^{ik\cdot R}S_{j,R},\quad
A_i(k)=\sum_R e^{ik\cdot R}A_{i,R}.
\]

The synthetic carrier's \(D_iH=\partial_{k_i}H-i[A_i,H]\) uses fixed basis and physical position conventions; check energy×length dimensions and periodicity. A genuine Wannier gauge connection may have k-dependence; link transporter and velocity must refer to the same gauge.

## 4. Versioned, fail-closed provenance manifest (required fields)

A suggested \`material_operator_provenance.json\` **schema**, not fake material evidence:

- \`schema_version\`, \`producer_code\`, \`producer_version\`, \`producer_build\`, \`producer_method\`, \`spin_source\`, \`spin_quality\` (\`paw_aware\` vs \`pseudo_only\` vs \`unknown\`), \`orbital_B_supported=false\`;
- \`run_kind\`: exactly \`synthetic_fixture\`, \`source_verified_material\` or \`unverified_external\`; **only** a validated producer may promote source_verified_material;
- \`structure_sha256\`, \`incar_sha256\`, \`kpoints_sha256\`, \`potcar_sha256\` (hash ONLY), \`wannier_win_sha256\`, \`files_sha256\` mapping for each actual byte-level input \`eig,mmn,amn,spn,u,u_dis,v,hr,r,wsvec,wout\` as applicable, without private paths;
- \`lattice_angstrom\`, \`num_kpts\`, \`num_bands_raw\`, \`num_wann\`, \`excluded_bands\`, \`num_electrons/filling\`, \`spinors=true\`, \`saxis=[0,0,1]\`, \`kpoints_frac\` + checksum and ordering convention;
- \`projection_origin\`, \`u_matrix_orientation\`, \`disentanglement\`, \`wannier_centers_angstrom\`, \`position_center_once=true\`, \`write_ndegen_applied\`, \`use_ws_distance\`, \`wsvec_applied\`, \`degeneracy_convention\`;
- \`physical_spin_operator_source\`, \`spin_pauli_normalization\`, \`zeeman_sign\`, \`g_factor\`, \`energy_unit=eV\`, \`length_unit=angstrom\`, \`B_unit=tesla\`, \`response_per_B_vs_per_b\`, and \`input_metadata_claims_verified\` per item;
- \`source_blockers\` and \`authentication_status\`: \`SYNTHETIC_ROUNDTRIP_PASS\`, \`MATERIAL_SOURCE_PENDING\`, \`MATERIAL_AUTHENTICATED\`, \`BLOCKED_MISSING_SPIN\`, \`BLOCKED_GAUGE_MISMATCH\`, \`BLOCKED_POSITION_WEIGHTING\`, etc.

Never put POTCAR content, raw charged wavefunctions, auth tokens, personal absolute paths or raw copyrighted proprietary files in GitHub. A checksum gives **integrity** but cannot prove physical correctness by itself.

## 5. Acceptance gates, in increasing evidence strength

**Gate A: synthetic roundtrip, allowed NOW with ordinary Python.**

- Serialize the current **v2_positive four-band synthetic** \`H(R),S_j(R),A_i(R)\` into a small, precisely specified test carrier (JSON or NPZ + checksums) with no false W90 origin. Parse it into \`WannierOperatorBundle\` with all exact units, centers and translations.
- Assert **element-wise** max matrix mismatch ≤1e-10, all \(-R\) partner Hermiticity, b=±1e-4 gap, and re-run direct/Wilson/Liouville in the **imported bundle**. Target direct≈+0.0342937475, Wilson≈+0.0341517445 at 41² (do not force equal finite-N values).
- Negative cases must fail: missing S/A, wrong shape/k order, changed hash, wrong center/extra center, wrong units/Pauli normalization, non-Hermitian matrix, reversed \(-R\), spin gauge changed independently, wrong phase convention. Distinguish mismatched-spin-gauge **numeric wrong response** from structurally invalid metadata.
- Do not mark synthetic roundtrip as an authentic Wannier90 raw-file import.

**Gate B: Wannier90 genuine-file parser with *synthetically authored conforming files*.** Do NOT run actual Wannier90.

- Implement a strict version-specific parser for \`_hr.dat\` and \`_r.dat\`, \`u.mat\`, \`u_dis.mat\`/\`v.mat\`, **formatted** \`spn\` and \`eig\`; treat \`amn/mmn\` as raw-source validation, not transformation substitutes. Include \`write_ndegen_applied\`/WS-weighting guards.
- Construct **small fixtures in actual documented file syntax** with known U, Udis, Bloch S, connection, and known H_R. Obtain independent golden reconstruction from synthetic bundle. Assert S_W(k)=V† S_Bloch V at all k and matching H(R), S(R), A(R) with no fitted gauge rotation.
- Explicitly reject unformatted binary spn in this first parser, partial/missing U, wrong 1-based/0-based k indexing, inconsistent k grid, missing spinor declarations, wrong win flags, bad band exclusion and unfamiliar Wigner–Seitz conventions.

**Gate C: real material source authentication, BLOCKED (no data or compute authorization).**

- Only if new *explicit* compute authorization is granted: generate a consistent VASP/Wannier90 dataset. VASP 6.6 has \`LWRITE_SPN\` with \`LWANNIER90\` for noncollinear spinors **if the VASP binary was compiled with VASP2WANNIER90**, supporting only default SAXIS; verify the local build first.
- Audit the source Bloch eigenstate gauge, actual \`spn\` provenance (including PAW/pseudo-wavefunction approximation status), U/Udis mapping and Wannier center/position matches. Compare interpolated H spin-orbit bands and local spin expectation with **independent** DFT data before any IGMC.
- The prior SOC static run was not a Wannier producer; its local moments cannot fill any missing H/S/A matrix. No new material job is authorized by this document.

## 6. One independent trap: 'embedded Wilson phase winding' detector

The hardened synthetic test is good evidence for Dxx, but its \`wilson_branch_check\` still defines

\[
\mathrm{winding\_wrapped}=\operatorname{wrap}(\gamma_{\rm last}-\gamma_{\rm first}),
\quad \mathrm{winding\_chern\_like}=[|\mathrm{winding\_wrapped}|\ge\pi].
\]

Since \(\operatorname{wrap}\) **always** lies in \([-\pi,\pi)\), this \(\ge\pi\) Chern/winding flag is mathematically unreachable except convention/rounding edge cases. Similarly, \(|\operatorname{wrap}(\Delta\gamma)|\ge\pi\) cannot detect most branch crossings. The genuinely exercised regression guards (**row within 0.05 rad of ±π** and **raw finite-b drift ≥π**) do detect the two injected faults, and the v2_positive phases are far from cut, so the **positive model comparison is not thereby invalidated**.

**Fix before material use:** construct an **unwrapped** \(\gamma(k_y)\) from the wrapped phase sequence, use its **integer net \(2\pi\) winding** across the periodic k_y cycle with the closing link, and test explicitly a positive synthetic model with nonzero Chern/winding as a negative phase-processing control. If the integer winding is nonzero, reject bulk polarization unless the adopted formalism handles Chern sectors. Do not use \(\operatorname{wrap}(\gamma_N-\gamma_1)\) as a winding number. The R4 regression text also says \"2.78 rad >= π\" incorrectly: it injects \`fakeb=3.5 rad\` so the **raw** drift is 3.5>π, while the displayed wrapped value is ~2.78<π. Correct explanatory text only, preserving historical result.

## 7. Research-stage decision

**GO for building/validating the synthetic and genuine-file reader. NOT GO for obtaining new SOC materials data or claiming nonzero AgF₂ response.**

The decisive physics story remains an actual, source-grounded, switchable insulating \(D_{xx}\) in two stacking domains whose axial \(M_z/L_z\) parity is distinguishable. The importer is an **infrastructure prerequisite**, not the scientific result.

### Documentation

- Official VASP \`LWRITE_SPN\` (VASP >6.4.2, noncollinear, compatible build, default SAXIS): https://vasp.at/wiki/LWRITE_SPN
- VASP \`LWANNIER90\` and compile flag: https://vasp.at/wiki/LWANNIER90
- Wannier90 \`.spn\` and \`u.mat/u_dis.mat/v.mat\` file roles: https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/
- Wannier90 \`write_rmn\`, \`write_ndegen_applied\`: https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/parameters/
- Wannier90 \`u_dis*u\` mapping: https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/
- The optional alternative \`wannierberri.utils.vaspspn\` uses pseudo-wavefunctions from WAVECAR and documents that approximation; it is **not** a full PAW-verified physical-spin source: https://docs.wannier-berri.org/en/latest/docs/vaspspn.html

**No simulation was run in writing this interface contract.**
