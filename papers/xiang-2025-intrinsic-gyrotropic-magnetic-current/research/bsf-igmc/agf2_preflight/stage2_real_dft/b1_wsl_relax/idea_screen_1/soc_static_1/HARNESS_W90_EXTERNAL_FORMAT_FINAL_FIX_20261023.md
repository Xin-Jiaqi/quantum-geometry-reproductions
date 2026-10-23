# DeepSeek Harness: final pre-archive correction of synthetic Wannier90 H/S/A importer

**Date 2026-10-23. Task scope: small Python/NumPy-only patch and genuine-file-format tests. No VASP, PBS, Wannier, AGF2 calculations or new licensed inputs.**

## Decision gate

A previous runner proposed closing PR #38 with \`PASS_SYNTHETIC_FORMAT_PARSER\`. The independent reviewer discovered **three definite reasons to postpone archive closure**; see \`POST_1A70A23_EXTERNAL_FORMAT_REVIEW_20261023.md\`. The Phase A synthetic-carrier result is accepted, but the external Wannier90 U-matrix parser is **not compatible with the documented ordering** as written. The goal is one small code correction and rerun, then make the archive decision. Do not reopen older physics projects or consume compute.

**Repository:** https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38  
**Branch:** research/agf2-b1-wsl-relax-20261021  
**Original pre-fix commit:** 1a70a23289eb37277b54171078af8232389e619b

### FIX 1 (mandatory): correct external U-file matrix storage order

Official Wannier90 \`seedname_u.mat\`, \`seedname_u_dis.mat\`, and \`seedname_v.mat\` files are **column-major: rows change fastest, then columns**. Current code is row-major on *both* sides:

- \`scripts/wannier90_synthetic_fixtures.py\`, \`mat_block\` uses \`for r: for c\`;
- \`scripts/import_hsa_operator_bundle.py\`, \`parse_kpoint_matrix_file\` also uses \`for r: for c\`.

Change both to **\`for c: for r\`**. Read and independently verify the source spec: https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/

**Independently authored golden fixture**: create a tiny complex nonsymmetric 4×4 U with nontrivial imaginary phases and a rectangular 6×4 U_dis. Record expected values and physical row/column indices explicitly in test code **before** using your own writer. Confirm a hand-written standard column-major text block parses exactly, and a row-major block fails a strict validation comparison; also re-run the old "wrong U transpose" fault and ensure it is really rejected.

After the fix run Phase A, Phase B identity, k-dependent nonidentity, six-band disentangled, Phase C winding, previous four-band direct/Wilson and regression suite. No fitting of physics response factors.

### FIX 2 (mandatory): real-space spin Hermiticity is only pairwise

\`read_hsa_carrier\` already validates **S_j(-R)=S_j(R)†**, but also incorrectly rejects any **S_j(R) != S_j(R)†** at *all* R. Only R=0 is individually Hermitian. Remove the overstrong off-site requirement; keep the R↔−R condition.

Add a **positive** fixture with complex off-site S(+R) non-Hermitian, S(-R)=S(+R)† that must import; and a **negative** fixture with mismatched ±R spin that must fail. Keep Hermiticity and physical Pauli norm checks for S(k) on a k sample, and explicit source-evidence limitations (numerical norm is not provenance).

### FIX 3 (mandatory): explicit k-list validation and honest spn limitation

- Parse k-coordinate arrays for \`seed_u.mat\`, \`seed_v.mat\`, \`seed_u_dis.mat\`, compare each with \`seed.win\` explicit kpoint list and with one another (orientation, counts and numeric tolerance). Stop discarding returned k lists with \`[1]\`.
- If both \`v.mat\` and \`u_dis.mat\` exist, verify **V = U_dis U** rather than silently choosing v and ignoring u_dis.
- Formatted \`spn\` has k **blocks but no independent per-block k-coordinate labels**. A user or program that permutes valid spn blocks may create a structurally legal file. The old test actually IMPORTED the permuted file and only detected response drift against a known synthetic reference. This is NOT intrinsic fail-closed importer behavior. No parser can unconditionally infer a private spn file's real k mapping from those blocks alone.
- Require a trusted shared source-producer manifest with file hashes, ordered-k checksum and independently validated spin/k-texture oracle before using \`MATERIAL_SOURCE_AUTHENTICATED\`. For unsupported external input mark \`SPN_K_PROVENANCE_UNVERIFIED\`, not successful authentication. In synthetic regression, retain the detectable wrong-response negative control, but **relabel it honestly**.
- Fix gauge law comments in \`test_import_hsa_operator_bundle.py\` and \`wannier90_synthetic_fixtures.py\`: with \`D_k H = partial H - i[A,H]\`, the transformed connection is \`A' = W† A W + i W† partial W\`, **plus i**. The current numeric writer already adds the correct analytic term; test and preserve it.

## Run and scope constraints

Only existing Python/NumPy, at most 10 minutes local CPU. No installing packages, no network copy of private data, no VASP/PBS/Wannier executable, no WSL electronic structure, no new AgF₂ material calculation. Run a handful of synthetic unit assertions, including the imported v2_positive direct +0.0342937475 and Wilson +0.0341517445 reference; the corrected parser must preserve these only when generated files obey the **official** format.

If the standards-verified corrections cannot be completed within this budget, publish **PARTIAL_EXTERNAL_FORMAT_PARSER** with the remaining failure and stop; **do not close PR as PASS**.

## Required deliverables in same PR #38

- Correct \`scripts/import_hsa_operator_bundle.py\`, \`scripts/wannier90_synthetic_fixtures.py\`, and small tests with explicit source-provenance/date annotations. Do not rewrite old SOC/DFT receipts, delete research history or alter unrelated code.
- New \`MATERIAL_HSA_EXTERNAL_FORMAT_FIX_20261023.md\` summarizing corrected U-file order, positive off-site-spin test, k validation/provenance limitation, genuine-format golden values and **actual test exit codes**.
- Short machine-readable QA summary, with new SHA, parser PASS/PARTIAL and material_authenticated **false**. Keep only synthetic fixtures on GitHub.
- PR Conversation final review response and exact next step. Closing PR can then be considered **as an archived, unmerged experiment** only after accounting for branch preservation. No automatic PR merge, delete-branch or costly compute action.

**One-line launch message:**

“Execute PR #38 task \`HARNESS_W90_EXTERNAL_FORMAT_FINAL_FIX_20261023.md\` strictly with Python synthetic files. Fix actual Wannier90 U column-major order, off-site spin Hermiticity, and k-list/provenance guards; run regressions and push to same PR. No VASP/PBS/Wannier, no merge/close until independent final review.”
