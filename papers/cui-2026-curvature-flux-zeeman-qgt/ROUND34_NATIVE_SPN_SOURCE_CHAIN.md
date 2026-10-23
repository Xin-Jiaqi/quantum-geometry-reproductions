# Round34 | Native physical-spin Bloch→Wannier chain and one-shot local harness

**Whole Cui source-paper reproduction: PARTIAL.** This round provides tested *native input parsing and mathematical projection*, while an independent real-material spin provenance check remains blocked until actual same-run input is recovered or a documented small calculation creates it.

## 1. Not just another synthetic spin model

Round26 tested supplied same-gauge H/S/A numerics. Round33 parsed actual native Wannier90 **H and r file formats** and warned that \`seedname.spn\` is in the Bloch gauge. Round34 implements source-format read-only parsing for:

| File | Supported format | Native contract |
|---|---|---|
| \`seedname.spn.fmt\` | Wannier90 formatted 3-Pauli input; additionally \`seedname.spn\` if explicitly ASCII | header num_bands num_kpts; m outer, n≤m inner, x/y/z fastest, each complex Re Im; \`S_{n,m}\` upper triangle, reverse by Hermitian conjugation |
| \`seedname_u.mat\` | Wannier90 output | nk,nwann,nwann, per-k 3 fractional coords, **Fortran column-major** rows then columns |
| \`seedname_u_dis.mat\` | Wannier90 output, required for entangled Bloch bands | nk,nwann,nbands, actual matrix size nbands×nwann; columns orthonormal |
| \`seedname_v.mat\` | Wannier90 output, when supplied | independent evidence for \`V=Udis U\` and the *same k order* |
| \`seedname.eig\` | native Bloch band energies | band index, k index, eigenenergy (eV), exact ordering |
| \`seedname_hr.dat\`, \`seedname_r.dat\` | Round33 native parsers | independent Wigner–Seitz-weighted H(k) and position matrix, no position-center double count |

Format authority: [Wannier90 file specifications](https://github.com/wannier-developers/wannier90/blob/develop/docs/docs/user_guide/wannier90/files.md), [Wannier90 U/V code](https://wannier-developers.github.io/wannier90/proc/plot_u_matrices.html), [official spin output of pw2wannier90](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/postproc/). For unformatted \`seedname.spn\` the parser **does not guess machine-specific Fortran record markers**: use the upstream \`w90spn2spn.x\` converter, or a genuinely formatted input.

## 2. Exact spin chain and independent Hamiltonian check

The **spin operator** is Pauli \(\sigma_a\), not a fabricated orbital spin label or an implicit factor of \(\hbar/2\). The file records Bloch eigenstate matrix
\[
S^{B}_{a,nm}(\mathbf k)=
\langle \psi_{n\mathbf k}|\sigma_a|\psi_{m\mathbf k}\rangle.
\]

For disentangled Wannier bands, concatenate **same-run** matrices \(V(\mathbf k)=U^{\mathrm{dis}}(\mathbf k)U(\mathbf k)\), where \(U^{\mathrm{dis}}\) has shape \((N_{\mathrm{Bloch}},N_{\mathrm{Wannier}})\) and \(U\) is \((N_{\mathrm{Wannier}},N_{\mathrm{Wannier}})\). For isolated bands \(N_{\mathrm{Bloch}}=N_{\mathrm{Wannier}}\), \(V=U\).

\[
\boxed{S^{W}_a(\mathbf k)=V^\dagger(\mathbf k)S^B_a(\mathbf k)V(\mathbf k)}
\]
and, independently of the Bloch Pauli data,
\[
H^{W}_{\mathrm{eig}}(\mathbf k)=
V^\dagger(\mathbf k)\,\mathrm{diag}(\epsilon_{n\mathbf k})\,V(\mathbf k).
\]

Independently reconstruct \(H^{W}_{hr}(\mathbf k)=\sum_R e^{i2\pi\mathbf k\cdot R}H_R/d_R\) using **the native Round33 H/r file reader** and check \(H^{W}_{\mathrm{eig}}-H^{W}_{hr}\) at the original input ordered \(\mathbf k\) mesh. The independent file \`_v.mat\`, if available, must equal \(U^{\mathrm{dis}}U\) to numerical tolerance. These are strong *mathematical matching tests*, but matching H cannot prove a separately supplied \`spn\` was computed from the same physical spinor Bloch states.

Inputs are bounded and validated: native length/count, triangular order, index order, complex Hermiticity, physical Pauli operator norm ≤1, k coordinates/order between U/UDIS/V, semior-thonormality, independent V product and energy-H Fourier match. The Round33 unsupported **unexpanded wsvec** regime remains rejected.

## 3. Independent positive and negative physics gates

A *separately constructed* finite spinor fixture uses **four Bloch bands, two Wannier bands and three k points**, a genuine 4-dimensional spin Pauli space, nontrivial Bloch phase gauges, entangled-subspace isometry and a nontrivial final Wannier rotation. Genuine Wannier90 formatted ASCII syntax is emitted to **temporary synthetic files**.

- Native spn triangle and x/y/z order must reconstruct *the independently known input physical Pauli matrices*.
- Native column-major U, Udis, V and eig must match independent arrays; \(V=UdisU\) must hold.
- Projected spin equals the separately computed \(V^\dagger\sigma V\) for all spin directions, all k.
- The physically permitted band-dependent Bloch **phase** transformation must leave **matched** \(S^W\) invariant.
- Deliberately freezing the old Bloch \`spn\` while changing Udis must produce a measurable, **wrong physical spin** despite matching dimensions. A wrong row-major U read must also fail.
- Independent source \`hr.dat\` Bloch-energy Hamiltonian test must reject an altered (but Hermitian) onsite H.
- Isolated-band \(U\)-only mode must pass, and omitting required Udis must fail.
- Bad V, k order, energy index ordering, spin diagonal imaginary contamination, unphysical Pauli norm and truncated input must fail.

The frozen real remote audit is \`results/frozen/round34_native_spin_chain_audit.json\`. CI acceptance demands Python 3.11 and 3.12, all historic Round17–33 audits, and the first-remote-values frozen regression.

## 4. One-shot local DeepSeek harness

The actionable instruction is committed beside this report:

**[ROUND34_LOCAL_DEEPSEEK_HARNESS.md](ROUND34_LOCAL_DEEPSEEK_HARNESS.md)**

It tells the local harness to: check existing actual noncollinear SOC QE/Wannier90 outputs in bounded named research folders; scan without copying or uploading originals; convert machine-dependent \`.spn\` to formatted **only via the official converter on an isolated copy** when available; test the native H/spin/gauge chain; if none exists and an installed official small tutorial can run inside a strict resource bound, generate **one** auditable small source example, not a bulk material calculation; produce one final local JSON receipt. It should not request many separate copy/paste steps from the user.

A real chain requires independently verified original \`.win\`/\`.wout\`, same-run upstream spinor calculation and \`pw2wannier90.x\` parameters, \`.spn.fmt\` plus \`u/u_dis/v/eig\`, and real \`hr/r\`. Even those are not yet real Cui KME: measured SI, material-specific EPW g, and an appropriate conserving collision vertex still remain missing.

## 5. Hard status bounds

There is no authenticated real source in this GitHub repository and no access to the local harness from this conversation. An execution instruction is **not proof the local operation ran**. The parser can verify geometry/operator algebra, but neither a hash nor a parsed \`spn\` proves its generation provenance; actual spin authentication remains \`false\` until the source chain is independently reviewed. No original paper Eq.(19), finite-frequency metallic KME, source Table I or genuine material SI response is accepted. **PAPER PARTIAL**.

## First independent GitHub CI, physical numbers, and strict frozen regression (2026-10-23)

The first Round34 native parser CI surfaced a *genuine test-fixture error*: in the **isolated** band case the fixture kept an extra Bloch gauge rotation outside the actual \`U.mat\`, although \`Udis.mat\` is absent. This was repaired by writing the **full Bloch→Wannier unitary into U** for that fixture. No physical validation threshold was relaxed. The independently rerun [GitHub Actions 39915478996](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39915478996) then passed **370/370** in **Python 3.11 and 3.12**, including all Round17–33 and the new Round34 audit.

The first clean run supplied [the frozen native formatted spin audit](results/frozen/round34_native_spin_chain_audit.json); one subsequent mandatory frozen-data test makes the exact final target **371**, to be accepted only on a freshly CI-tested exact SHA and independently again on merged \`main\`.

Measured **native-format synthetic 4 Bloch-band, 2 Wannier-band, 3-k spinor** data:
- \`spn\` upper triangle, xyz order, and complex Pauli parser round trip residual **0**; independent native \`V\` vs \`Udis U\` **0**; projected spin vs independently computed \(V^\dagger S V\) **0**.
- Physical spin invariant under **matched** Bloch-band phase changes to **2.61845576667e-16**; intentionally stale spin (Udis transformed but spn not) gives **0.5436847695** wrong Pauli matrix; intentionally reading column-major U as row-major gives **1.1991698578** wrong matrix.
- Isolated-band true full native U pathway spin residual **1.2412670766e-16** (no Udis required).
- Native \`eig\` projected Hamiltonian vs Round33 independent **native hr** Wigner–Seitz Fourier error **8.77708367144e-17 eV**.
- Physical projected spin norm max **0.99491408 ≤ 1**.
- Every source status flag remains **\`true_source_spn_verified=false\`**, **\`actual_material_spin_provenance_authenticated=false\`**, **\`real_material_KME_ready=false\`**, original Cui paper **PARTIAL**.

The local DeepSeek harness handoff is a **read-only task instruction** committed into GitHub; this GitHub Actions verification does **not** run the user's local harness and does **not** authenticate any local first-principles calculations.
