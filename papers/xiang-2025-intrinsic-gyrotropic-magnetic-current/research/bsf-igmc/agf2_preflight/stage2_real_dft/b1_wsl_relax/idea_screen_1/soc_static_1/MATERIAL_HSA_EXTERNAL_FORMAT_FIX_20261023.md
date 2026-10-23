# MATERIAL_HSA_EXTERNAL_FORMAT_FIX — 2026-10-23

Verdict: **PASS_SYNTHETIC_FORMAT_PARSER** (fix round; NOT `MATERIAL_AUTHENTICATED`)

Task: `HARNESS_W90_EXTERNAL_FORMAT_FINAL_FIX_20261023.md`. Machine-readable
receipt: `MATERIAL_HSA_EXTERNAL_FORMAT_FIX_QA.json` (written by
`scripts/test_import_hsa_operator_bundle.py`, exit code **0**). The prior
1a70a23 receipts are untouched; this dated note does not overwrite them.
Only existing Python/NumPy on the four-band synthetic model; **no
VASP/PBS/Wannier executable, no new AgF₂ calculation**
(`no_wannier_executed: true`). Local CPU ≈ 200 s (import audit) + gate
re-run. The PR stays open; no merge/close until independent review.

---

## FIX 1 — external U-family files are column-major (rows fastest)

`seedname_u.mat`, `seedname_u_dis.mat`, `seedname_v.mat` store matrix
elements column-major (rows change fastest, then columns; official docs,
https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/).
Both the production writer and the parser previously iterated row-major,
which silently transposes/misassigns a genuine complex nonsymmetric U.

- `scripts/wannier90_synthetic_fixtures.py` `mat_block` now loops
  `for c: for r:`; `scripts/import_hsa_operator_bundle.py`
  `parse_kpoint_matrix_file` now fills `m[r, c]` with `c` outer, `r` inner.
- **Independently authored golden fixture** (hand-written literals recorded
  in `test_import_hsa_operator_bundle.py` BEFORE any project writer runs):
  a complex nonsymmetric 4×4 U with nontrivial imaginary phases at two
  k-points, and a rectangular 6×4 U_dis (header `1 4 6`), each with the
  physical (row,column) indices documented:
  - `u_golden_exact: true` — column-major text parses EXACTLY to the
    recorded matrix,
  - `u_rowmajor_rejected: true` — a row-major reordering of the same values
    does NOT match the golden (strict comparison),
  - `u_dis_golden_exact: true` — the rectangular block parses exactly,
  - `writer_column_major: true` — the production writer's on-disk block
    equals the column-major flat of the in-memory U (independent
    reconstruction from `eigh(H)`).
- The old "wrong U orientation" (per-k transposed u.mat) negative still
  fails closed: `H(k) mismatch vs eig at k=0: 6.629e-01 (tol 1e-08)`.
  The transposed file is now written and re-read in the corrected
  column-major convention.

## FIX 2 — real-space spin Hermiticity is pairwise only

`read_hsa_carrier` no longer requires every off-site S_j(R) to be
individually Hermitian. The physical condition is S_j(−R)=S_j(R)† for the
R/−R pair; only R=0 is individually Hermitian; the reconstructed S_j(k)
must be Hermitian with a physical Pauli norm bound (checked on a k sample,
with a documented 0.05 slack for legitimate off-site terms; the tight
‖S(0)‖ ≤ 1+1e-9 check stays in the pinned `validate()`). Numerical norm
bounds are necessary-not-sufficient evidence and are **not** provenance.

- Positive fixture `offsite_spin_positive`: complex off-site S(+R)
  non-Hermitian (imaginary-symmetric perturbation δ, δ†=−δ, preserving
  S(0)) with S(−R)=S(+R)† — **imports** (was rejected before the fix);
  off-site non-Hermiticity and the pairing are asserted on the re-read
  bundle.
- Negative `mismatched_spin_pair`: S(−R) ≠ S(+R)† — raises
  `R=(-1, 0) missing Hermitian partner -R`.

## FIX 3 — k-list validation and honest spn limitation

- The parser now keeps the k-coordinate arrays of `u.mat` / `v.mat` /
  `u_dis.mat` (no `[1]` discard): each must be a finite, unique
  fractional k list (`_check_k_list_consistent`) and match the `seed.win`
  explicit kpoint block and one another at 1e-8 (orientation, count,
  values). Violations raise.
- When both `v.mat` and `u_dis.mat` exist, the reconstruction verifies
  `V = U_dis U` per k (≤1e-8) instead of silently choosing v; the
  consistency value is recorded in the receipt (`V_UdisU_consistency`).
- Formatted `spn` has k blocks but **no independent per-block k labels**:
  a permuted spn block file is structurally legal and no parser can
  fail-closed on its k mapping from the file alone. The importer now
  labels the bundle provenance
  `spin_k_provenance: "SPN_K_PROVENANCE_UNVERIFIED"`, and the permuted-spn
  control in the audit is **relabelled honestly** as
  `RESPONSE_LEVEL_CONTROL_NOT_INTRINSIC_REJECTION` (response drift
  8.1e-6 ≫ the 1e-11 numerical floor, but detectable only against a known
  synthetic reference).
- `MATERIAL_SOURCE_AUTHENTICATED` remains off and will additionally require
  a trusted shared source-producer manifest with file hashes, an ordered-k
  checksum and an independently validated spin/k-texture oracle before it
  can ever be used; for unsupported external input the importer says
  `SPN_K_PROVENANCE_UNVERIFIED`, never successful authentication.
- Gauge-law comments corrected to the plus-i convention:
  A' = W†AW + i W†∂_k W under D_k H = ∂_k H − i[A,H] (the numeric writer
  already added the positive analytic term; the response invariance
  survives, see Phase B).

## Results (actual exit codes)

`python3 scripts/test_import_hsa_operator_bundle.py` → **exit 0** (200.7 s):

| item | value |
|---|---|
| external golden | u exact, row-major rejected, u_dis exact, writer column-major — all true |
| Phase A | SYNTHETIC_CARRIER_PASS; carriers bitwise (0.0 recovery) incl. off-site-spin positive; 9/9 negatives raise |
| v2_positive 41² | direct +0.0342937475, Wilson +0.0341517445 (reproduced) |
| Phase B | SYNTHETIC_FORMAT_PARSER_PASS; B1 41² op-recovery ≤1e-10 (spn 16-dec floor), H-vs-eig 2.3e-12; B2 smooth 21² direct 0.0342937475 (connection law), H-vs-eig 7.7e-10; B3 disentangle 21² direct 0.0342937475, V†V−I 3.3e-15; negatives: wrong-U raises 6.629e-01, permuted-spn relabelled response-level |
| Phase C | WILSON_WINDING_FIXED; ramp winding 1 closure 0.0; Chern 2-band ±1; v2_positive winding 0 |
| k validation | U-family k lists matched win kpoints and each other (atol 1e-8) |

`python3 scripts/test_gate_regressions.py` → **exit 0** (7/7) and
`python3 scripts/test_multiband_igmc_gate.py` → **exit 0** (PASS; direct vs
Wilson 0.41% at 41², `chern_winding_flagged: false`) — the column-major and
provenance changes leave the hardened gate green.

## Honest limitations (unchanged scope)

- `MATERIAL_AUTHENTICATED` is never claimed: synthetic authoring only;
  upstream has LWAVE=.FALSE., so no real WAVECAR/W90 files exist.
- Smooth/disentangled final gauges carry the documented 361-R window
  Fourier-truncation error (H-vs-eig ~5–8e-10; direct response d ~1e-10).
- The spn k-mapping cannot be proven from the file alone
  (`SPN_K_PROVENANCE_UNVERIFIED`).

## Deliverables in PR #38

- `scripts/import_hsa_operator_bundle.py` — column-major parse, pairwise
  spin rule, k-list validation, V=U_dis U check, spn provenance label,
  k-sample spin checks.
- `scripts/wannier90_synthetic_fixtures.py` — column-major write.
- `scripts/test_import_hsa_operator_bundle.py` — golden literals + FIX 1–3
  tests; QA JSON.
- `MATERIAL_HSA_EXTERNAL_FORMAT_FIX_20261023.md` (this file) and
  `MATERIAL_HSA_EXTERNAL_FORMAT_FIX_QA.json` (new SHA below).

## Next step (single)

Independent final review of this fix round; only after its explicit sign-off
should PR #38 be archived — as an **archived, unmerged experiment** with the
branch/ref preserved (no automatic merge, delete-branch or compute action) —
or merged after a separate broad-PR review.
