# DeepSeek Harness: bootstrap genuine Wannier90 H/S/A importer, ZERO material compute

**2026-10-23 | GO for CODE AND SYNTHETIC-ONLY FILE INTERFACE | No VASP, PBS, Wannier or new AgF₂ job is authorized.**

## 0. Operating contract

Use this existing [PR #38](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38) and its branch **research/agf2-b1-wsl-relax-20261021** as the only agent-to-agent exchange. All commits are additive or precisely corrected with history preserved. Avoid Git force pushes, rewriting DFT receipts, merging unrelated branches or asking the researcher to relay outputs.

**Read before coding:** MATERIAL_HSA_AUTHENTICATED_IMPORT_CONTRACT_20261023.md, MULTIBAND_GATE_HARDENING_20261023.md, scripts/test_multiband_igmc_gate.py, scripts/wannier_operator_contract_pinned26.py, and the frozen Cui Round26 source commit cdac6418166466b5243e73d8eb291d6dea8d4772.

The last scientific gate was \`PASS\` only for a corrected gapped **synthetic four-band** direct/Wilson/Liouville response, not for AgF₂. Do not reinterpret that as permission to compute new DFT or Wannier. Previous AgF₂ SOC SCF does **not** contain the necessary WAVECAR, W90 H/S/A or source Bloch spin matrix. Do not try to reconstruct them from local moments.

## 1. Scope: first support a synthetic roundtrip, then strict genuine-format parsing

**Phase A — smallest versioned synthetic carrier**:

1. Serialize only the already-tested \`v2_positive\` and \`v2_spinZ_pos\` **complete** H(R), S(R), A(R), lattice, integer R, centers, site/basis identity, and explicit eV/Å/Pauli units into a small JSON/text fixture; no proprietary/magical synthetic \`wannier90.spn\` claims.
2. Implement a strict reader that validates source file hashes, types/shapes, finite complex values, R ↔−R Hermiticity, spin Hermiticity and norm bound, \`centers=diag A(R=0)\`, spin axis and units, then returns the existing \`WannierOperatorBundle\`. Distinguish **source-integrity-check** from actual verified Bloch spin provenance.
3. Re-run the existing **imported** direct Dxx/ω and rank-two A-embedded Wilson slopes, including the fixture where omitting A produces the known strong discrepancy. Require the recovered matrices match bitwise or ≤1e−10 and the 41² values \`+0.0342937475\` and \`+0.0341517445\` within set tolerances.
4. Negative tests: swap S gauge, double centers, corrupt file/hash, remove matrix, wrong basis/source k order, wrong units or Zeeman factor, wrong lattice orientation, inconsistent R list, one missing −R, non-Hermitian S/H/A. A negative control must actually trigger **nonzero exit**; avoid mere \`detected=true\` without enforcement.
5. Label the status \`SYNTHETIC_CARRIER_PASS\`, never \`MATERIAL_AUTHENTICATED\`.

**Phase B — genuine Wannier90 syntax, authored synthetic files only (NO Wannier run):**

1. Verify official Wannier90 version/file format before implementing. Begin with \`seed_hr.dat\`, \`seed_r.dat\`, \`seed_u.mat\`, \`seed_u_dis.mat\` or \`seed_v.mat\`, \`seed.eig\`, **formatted** \`seed.spn\`; inspect \`.amn/.mmn/.win/.wout\` as source-provenance and k/source definitions. \`amn\` is the **initial projection, NOT** final Wannier rotation. Explicitly document whether each is an input or an output of W90.
2. Strictly declare the first support case: 2D slab, 3D orthorhombic lattice, Rz=0, \`write_ndegen_applied=true\`, all degeneracy=1, no ambiguous WS fold, \`SAXIS=(0,0,1)\`, a single consistent spinor Bloch basis/k ordering; reject anything else rather than silently dropping phases, indices or double-dividing R weights.
3. Use actual Bloch→Wannier transform \`V(k)=U_dis(k)U(k)\` (or U alone for isolated band group) to move formatted spn matrices. Check V†V≈I, \`S^W(k)=V†S^Bloch(k)V\`, same-gauge Hermiticity/spin normalization, and reconstruct S(R) by the **matching discrete Fourier convention**. Independently check k-by-k input/output equality and matched H(R) with \`eig\`. The position matrix \`r.dat\` is **already final Wannier gauge**; check its R list/units/centers/weighting, not a second U transformation.
4. Build tiny **synthetic** examples in official file syntax including a *nonidentity complex unitary* U and a disentangled rectangular U_dis, with known Pauli and off-diagonal A. First test the identity gauge, then nonidentity U. Test wrong U orientation / transpose and mismatched k order fail closed, not merely yield different response.
5. Integrate the resulting bundle with the **existing** validated multiband Dxx response tests. Use max 41² and bound ordinary CPU to **≤20 minutes** for the entire bootstrap. If impossible within this resource limit, deliver a PARTIAL narrow parser with precise blockers rather than claiming comprehensive support.
6. Publish no raw private VASP data, wavefunctions, POTCAR, private source paths, secrets, or external licensed PDFs in GitHub.

**Phase C — small scientific correction while nearby:** the hardened synthetic \`wilson_branch_check\` labels \`abs(wrap(last−first))>=pi\` as Chern/winding. Since wrap ∈[−π,π), that predicate cannot reliably detect integer winding. Replace or augment it with genuine unwrapped Wilson-flow \(2\pi\)-winding over the periodic ky cycle including closure. Add **one synthetic nonzero-winding regression** and confirm the existing positive fixture is still safe. The R4 comment \"2.78 rad >= π\" is erroneous, since the injected **raw** drift was 3.5 rad >π and its wrapped equivalent ~2.78 rad <π. Correct only the comment; no need to rerun DFT. Do not silently claim Chern=0 based on a wrapped endpoint alone.

## 2. Absolute scientific acceptance criteria

The material importer is considered **only** \`SYNTHETIC_FORMAT_PARSER_PASS\` if:

- a manifest-hashed synthetic roundtrip exactly preserves four-band H/S/A and units;
- strict official-format source files were parsed with the **correct raw Bloch→final Wannier** V(k) and matched R/weight/center conventions;
- direct/Wilson reference numerical response remains consistent, without fitting;
- each negative case fails in the intended place, including an intentionally bad source spin gauge and a double-counted position center;
- all tests are actually run and pass; package and file-format support limits are explicitly declared;
- a nonzero Wilson winding test correctly rejects an invalid polarization extraction;
- no new VASP/Wannier electronic-structure task was launched.

\`AUTHENTICATED_REAL_AGF2_OPERATORS\` is **impossible** in the current phase because raw AgF₂ W90 Bloch/spin/gauge inputs have not been generated or authenticated. A manifest checksum alone cannot promote them.

## 3. GitHub deliverables

Commit to \`idea_screen_1/soc_static_1/\`:

- \`MATERIAL_HSA_SYNTHETIC_IMPORT_AUDIT.md\` and \`.json\`: supported W90 version/settings, actual success/failure evidence, tensor conventions and remaining source blockers.
- \`scripts/import_hsa_operator_bundle.py\` and \`scripts/test_import_hsa_operator_bundle.py\` with executable assertions, not pseudocode. Optionally separate parsing and manifest helper modules.
- Small **synthetic-only**, clearly labeled reference fixtures with SHA. For larger fixtures, generate deterministically from code and commit generator rather than binary blobs.
- Document schema/source provenance constraints as actual implemented guards; no \`real material PASS\` unless it is genuinely available in a **future, separately authorized job**.
- Short PR Conversation summary with a definite \`PASS_SYNTHETIC_FORMAT / PARTIAL / FAIL\`, list of supported official file formats, one highest-value future material production gate.

Before proposing any later real material run, **check whether local vasp_ncl supports \`LWRITE_SPN=.TRUE.\` + \`LWANNIER90=.TRUE.\`** with the correct compile flag and default SAXIS; this can be a read-only binary/documentation inspection, not an execution of VASP. VASP 6.6 supports the feature in principle, but actual build must be verified. If unverified, keep the source producer marked UNKNOWN.

## 4. Cost and operational guardrails

No VASP/PBS/Wannier jobs, no restart, no convergence/MAE, no automatic installing heavyweight packages, no new material data generation. Strictly no reusing the overspent 4h SOC approval. Use ordinary Python/NumPy within the 20-minute local CPU target. If no valid source parser/gauge proof can be completed, write PARTIAL and stop; do not relax the science to claim PASS.

The goal is an **auditable, minimum-interface gate** that can someday ingest real H/S/A without an untraceable spin/position gauge mistake. Avoid building a large general workflow framework. Results must be committed to PR #38; the user does not transfer logs between AI agents.

### One-line handoff

“Read PR #38 \`HARNESS_HSA_SYNTHETIC_IMPORT_BOOTSTRAP_20261023.md\`; develop and run the restricted synthetic H/S/A importer and tests, prove positive/negative controls, and commit directly to the same PR. No VASP/PBS/Wannier.”
