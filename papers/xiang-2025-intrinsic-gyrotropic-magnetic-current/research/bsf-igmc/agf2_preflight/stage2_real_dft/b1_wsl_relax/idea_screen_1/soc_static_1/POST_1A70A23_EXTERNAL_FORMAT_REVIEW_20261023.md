# Independent review of 1a70a23: genuine Wannier90 file-format acceptance BLOCKED

**Date 2026-10-23 | PR #38 | Reviewer decision: PHASE_A_ACCEPTED / PHASE_B_NOT_FULLY_WANNIER90_COMPATIBLE / PHASE_C_SUPPORTED_ON_SYNTHETIC / DO_NOT_CLOSE_AS FULL-PARSER-PASS YET.**

**No electronic-structure, VASP, PBS, Wannier or material operator generation was run by this review.** This dated note does not overwrite the original 1a70a23 receipts. The PR remains open for one small synthetic-only format/serialization correction.

## 1. Parts that do survive

- Phase A: manifest-hashed synthetic H(R)/S(R)/A(R) carrier roundtrip reproduced the original four-band v2_positive direct 0.0342937475 and Wilson 0.0341517445 on 41², including synthetic negative tests. The source reports exact 0.0 matrix-recovery error. This warrants **SYNTHETIC_CARRIER_PASS**.
- The synthetic four-band low-ω response remains validated against direct/Wilson/Liouville baselines. This is not any AgF₂ material IGMC signal.
- Phase C: the old mathematically unreachable \`abs(wrap(last-first))>=π\` winding predicate was replaced by an unwrapped Wilson-flow integer test, including a nonzero-winding two-band Chern model. This is a real improvement (subject to the finite-mesh and cut convention documented in the submitted tests).
- Physical gauge algebra: under \(\,H'=W^\dagger HW\,\), \(\,D_kH=\partial_k H-i[A,H]\,\), the compatible position connection is \(\,A'=W^\dagger AW+iW^\dagger\partial_kW\,\). The *numerical writer* adds the positive analytic derivative contribution, consistent with this convention, **but code/docstring comments in test_import_hsa_operator_bundle.py and wannier90_synthetic_fixtures.py incorrectly state minus-i**. Correct the comments/convention in one place; do not change a working numerical sign without algebraic test.

## 2. A definite Phase B file-format error: U-matrix serialization and parse are both row-major

**Wannier90's official file documentation** for \`seedname_u.mat\`, \`seedname_u_dis.mat\`, \`seedname_v.mat\` states that matrix values are written in **column-major order, cycling over ROWS first and then COLUMNS**.

But the submitted source does the opposite:

1. \`scripts/wannier90_synthetic_fixtures.py\`, \`mat_block\` around lines 209–219:

~~~python
for r in range(nrows):
    for c in range(ncols):
        z = M_list[k][r, c]
~~~

2. \`scripts/import_hsa_operator_bundle.py\`, \`parse_kpoint_matrix_file\` around lines 320–362:

~~~python
for r in range(nrows):
    for c in range(ncols):
        m[r, c] = read_next_complex()
~~~

Both implement **row-major** traversal (the column index c changes fastest), while actual Wannier90 file syntax is **column-major** (the row index changes fastest). Thus the two scripts can perfectly reconstruct each other's **synthetically misformatted files** and pass response tests, while silently transposing/misassigning a genuine external Wannier90 U or U_dis. For a generic complex, nonsymmetric unitary, this is a severe **physics/gauge error**, not merely a metadata style difference. The old "wrong U orientation" negative control may also be misleading because the test generator shares the same convention bug.

**Mandatory correction:** Both parser and writer must iterate **c outer, r inner**, and validate against a small *independently hand-authored* official-syntax golden matrix that is complex and non-symmetric (e.g. two blocks with different k points). Never use a serializer written by the same developer as the sole oracle of the external file standard.

Primary source: [Wannier90 official file format](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/), sections \`seedname_u.mat\`, \`seedname_u_dis.mat\`, \`seedname_v.mat\` ("column-major order").

## 3. A second strong limitation: k-order mismatch is NOT structurally rejected

The original audit text says \`wrong U orientation and spn k-order fail closed\`, but actual \`scripts/test_import_hsa_operator_bundle.py\` Phase B negative test **swaps two formatted spn k-blocks, successfully constructs a bundle**, and then merely checks \(\,|D_\mathrm{wrong}-D_\mathrm{golden}|>10^{-6}\,\), with reported difference 8.1×10⁻⁶.

That is a **synthetic response negative control**, not fail-closed **import authentication**. In \`reconstruct_bundle_from_w90\`:

- \`parse_kpoint_matrix_file\` returns k-points, but the caller uses \`[1]\` (matrices only), discarding its k-point list;
- no explicit per-k comparison of the \`u.mat\`, \`v.mat\` and \`u_dis.mat\` k coordinates with each other or \`win\` is performed;
- a formatted \`spn\` supplies only an ordered stack of k blocks, without its own independent k-coordinate labels. Without external provenance or independent spin observables, an arbitrary permutation of spn blocks cannot always be **detected solely from the file contents**.

**Mandatory scope correction:** explicitly validate the k coordinates from each U-family matrix file against \`win\` within documented tolerances and against each other. For spn, require producer-level common-k-order provenance and a source hash or an independently sourced spin-texture oracle. Label \`SPN_K_MAPPING_PROVENANCE_REQUIRED\` until that exists. Do not claim an internally permuted spn is *guaranteed* to be rejected by physics alone.

## 4. Validation path, zero VASP/Wannier

1. Correct parser/writer U-family traversal to \`for c: for r:\` and add a **hand-authored external-golden** asymmetric complex U matrix case whose expected entries are known without this project's writer. Test a 4×4 and rectangular 6×4 matrix. Verify the genuine W90 printed order using official docs or existing fixture files.
2. Re-run Phase A/B/C, the original four-band 41² response/negative tests, disentangled 6-band case, and the winding regression. An honest failure after correct format parsing is grounds for PARTIAL, not for restoring the wrong interpretation.
3. Compare \`u.mat\` and \`v.mat\` file k lists with \`win\` explicitly. If \`v.mat\` and \`u_dis.mat*u.mat\` both exist, verify consistency (rather than choosing one and ignoring the other). Reject any mismatch.
4. Replace "spn permutation fail-closed" with precise evidence and a new **SOURCE_K_MAPPING_NOT_VERIFIABLE_FROM_SPN_ALONE** failure/limitation in material authentication.
5. Correct the documented A transformation sign; only if a separate analytic covariant derivative test fails should implementation be changed.
6. Post an erratum and keep original numerical receipt intact. A new synthetic Phase B PASS can be issued **only after real-format tests execute with success**.

## 5. Safe research-stage decision

**Keep PR #38 open temporarily, despite the runner's proposal to close**; its goal as a synthetic format parser has a definite standard mismatch. The branch also contains **44 commits / 70 changed files, still DRAFT and UNMERGED into main** as observed at review. Closing without merge leaves these important artifacts off main. This is reversible but not a good scientific handoff until the current acceptance standard is honestly satisfied and the preservation/merge strategy is explicit.

Once the small correction passes, choose consciously between:
- merging the validated synthetic code/results to main (if this broad PR is suitable and reviewed), OR
- closing this draft PR as an archived experiment while preserving its source branch/ref, its precise head SHA and a README pointer from a tracked integration branch.

The next real-material gate still requires a **separately authorized VASP/Wannier producer** with actual \`hr/r/spn/U\` and true common-source provenance. Nothing here authorizes that work.

**Reviewer scope:** verified actual GitHub source lines and official Wannier90 file order; did not execute this project's recent synthetic test suite in its private WSL environment. Source-reported Phase A and response numbers are accepted at that level of evidence.
