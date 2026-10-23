# Round33 | Native Wannier90 dense H/r parsing, Wigner–Seitz weights, and physical-spin stop gate

**Whole source Cui paper: PARTIAL. No authentic material Wannier90 hr.dat, r.dat or spn physical-spin source files were found in the repository.**

## Why this is a genuinely distinct stage

Round26 supplied a **synthetic in-memory H/S/A same-gauge gate**; Round32 supplied an explicit **normalized** EPW output schema with anti-spoofing checks. Neither actually reads original Wannier90 **dense text** data. Round33 implements real **file-format parsers** for the standard Wannier90 \`seedname_hr.dat\` and \`seedname_r.dat\` output formats and performs independent matrix/Fourier/weight tests. The **test files themselves are generated, not actual material calculations**.

Authoritative upstream specification: [Wannier90 file format documentation](https://github.com/wannier-developers/wannier90/blob/develop/docs/docs/user_guide/wannier90/files.md). \`hr.dat\` contains num_wann, nrpts, a block of Wigner–Seitz degeneracy integers, then per-R row/column 1-based complex Hamiltonian elements. **\`r.dat\` has no degeneracy block** and must use the corresponding \`hr.dat\` weights unless the newer \`write_ndegen_applied\` output is already expanded. The official \`wsvec.dat\`/use_ws_distance conventions can require orbital-pair-specific shift vectors and must not be guessed or applied twice. [Wannier90 parameter guide](https://github.com/wannier-developers/wannier90/blob/develop/docs/docs/user_guide/wannier90/parameters.md).

## Supported exact format and restrictions

The native text parsers read, with bounded allocation and file sizes, a two-dimensional ordinary Wannier model with \(R_z=0\). \`hr.dat\` rows \`Rx Ry Rz m n ReH ImH\` and \`r.dat\` rows \`Rx Ry Rz m n ReX ImX ReY ImY ReZ ImZ\` are accepted, including Fortran D exponents. Indexing is **1-based**, R blocks must be complete and unique, and both files must have matching R order, number of Wannier orbitals and number of translation blocks. Refuse missing or extra records, duplicate elements, bad degeneracy integers, mismatched inverses, missing \(R=0\), non-Hermitian \(-R\) partners, imaginary onsite centers and any R with \(R_z\ne0\).

Given per-R degeneracy \(d_R\), the **basic non-expanded source** is interpreted as
\[
  H^{\rm effective}_R=H^{\rm raw}_R/d_R,\qquad
  A^{\rm effective}_{a,R}=r^{\rm raw}_{a,R}/d_R,
\]
and
\[
H(q)=\sum_R e^{2\pi i R\cdot q} H^{\rm effective}_R,\quad
A_a(q)=\sum_R e^{2\pi i R\cdot q} A^{\rm effective}_{a,R}.
\]
The onsite projected position matrix yields Wannier centers \(\mathrm{diag}\,A_a(R=0)\), with **no second addition of orbital centers**.

For \`write_ndegen_applied=true\`, demand \`hr.dat\` degeneracy all ones, and accept the already expanded \`hr\`/\`r\` pair with no extra division. For \`use_ws_distance=true\` with **unexpanded** files, fail closed because \`wsvec\` corrections are missing; the pairwise Wigner–Seitz extra mapping cannot be reconstructed by blindly multiplying each R by a scalar degeneracy. Source flags are **explicit caller declarations**, not independent forensic verification of a raw \`.win\`/\`.wout\`; actual material provenance remains unverified.

**Most important boundary:** these H/r parsers cannot authenticate the **projected physical Pauli spin**. Wannier90 \`*.spn\` is an input Bloch-gauge spin matrix and must be transformed with the actual \`U(k)\`/disentanglement from the matching calculation. A filename, matching band dimension, matching SHA, pseudospin, or made-up \`basis_id\` cannot replace this. The importer returns \`PARSED_H_R_SPIN_NOT_AUTHENTICATED\`, \`physical_spin_available=False\`, \`kme_ready=False\`, original paper \`PARTIAL\`. See [official Wannier90 postprocessing spn notes](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/postproc/) and [WannierBerri Wannier-gauge spin helpers](https://docs.wannier-berri.org/en/master/docs/w90_files.html).

## Independent positive and negative controls

The standalone audit [audit_native_wannier90_hr_r.py](scripts/audit_native_wannier90_hr_r.py) generates native-format **fixture files only**, then:
- verifies R degeneracies \([2,1,2]\), H/r raw values and exact per-R effective reconstruction;
- independently evaluates analytic complex \(\pm R\) Fourier sums at a nontrivial fractional q and Hermiticity, including all Cartesian position components;
- **omits** Wigner–Seitz division on purpose and demonstrates a finite **wrong-H** result;
- separately writes **expanded** \`write_ndegen_applied\` files and verifies no double division;
- deliberately corrupts inverse R Hermiticity, native record count, one-based orbital indices/duplicate entries, degens, raw r pairs, R group/order, wsvec-required flags, and confirms they are rejected;
- returns a still **blocked** material-spin/KME verdict even when synthetic H/r parsing succeeds.

The frozen numerical output is [Round33 raw synthetic file-format audit](results/frozen/round33_native_wannier90_hr_r_audit.json), to be created from the first independent GitHub CI after the tests succeed and then enforced as a regression. **No actual EPW or materials quantum transport rates are computed in Round33.**

## Next source inputs needed

Genuine *same-seed* original \`hr.dat\` and \`r.dat\`, associated \`win/wout\` and any \`wsvec.dat\` plus an independently validated Bloch→Wannier spin provenance chain (\`spn\`, \`u.mat\`/disentanglement) are required **before** the Round26 H/S/A or Round32 complex e–ph gate can be considered on real material. Real matched complex EPW \(g\), DFPT phonons and convergence still remain absent. This is a **native source decoding upgrade**, not a real material KME reproduction.

To reproduce:
~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/audit_native_wannier90_hr_r.py
~~~

The original Cui paper remains **PARTIAL** irrespective of synthetic test count.

## Remote first-run evidence, exact fixture numbers, and frozen acceptance

Initial GitHub PR CI [39911971914](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39911971914) independently passed **344/344** Cui isolated tests and all existing Round17–32 plus new Round33 audits in **Python 3.11 and 3.12**. New frozen JSON [round33_native_wannier90_hr_r_audit.json](results/frozen/round33_native_wannier90_hr_r_audit.json) records actual Python 3.11 remote output before it was added as a regression. The final commit must pass **345/345** and all audits in both interpreters before merging.

Native **generated-format** 2-band, three-R tests show:
- Wigner–Seitz weights exactly **[2,1,2]** and native hr and r weighted reconstructions **0**.
- Independently analytical complex H(q) Fourier error **1.5515838458e-17 eV**, all Cartesian position Fourier error **0**, k-space H Hermiticity error **1.55e-17 eV**, and onsite centers matching raw r **0 Å**.
- **Deliberately omit the hr-only degeneracy weight:** H(k) differs by **0.1059858083 eV**, proving nontrivial scientific importance.
- Already-expanded write_ndegen_applied control: residual **0**; unsupported unexpanded wsvec case rejected.
- Recorded invariants \`authentic_material_input_found=false\`, \`wannier_physical_spin_provenance_authenticated=false\`, \`spn_bloch_to_wannier_gauge_implemented=false\`, \`kme_real_material_ready=false\`, original paper **PARTIAL**.

This stage parses bona fide **text grammar**, not bona fide first-principles materials source content.
