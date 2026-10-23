# Stage8k — DeepSeek harness one-file READ-ONLY assignment

**Only this repository document is required.** Continue from newest GitHub head of research/bsf-nmee-stage7-sw-gate-20261020. This is a **no-new-calculation** validation of old Stage8i spinor output on Japan IMR. Do NOT submit PBS jobs, start VASP, run Wannier, or modify the frozen PRL reproduction. The user should not be asked to copy terminal commands, send WAVECAR files or manually collect output.

## Inputs already available

- Independent physical error correction: STAGE8J_INDEPENDENT_VALLEY_SPIN_ERRATUM_20261022.md
- Stage8k independent physics: STAGE8K_NONDEGENERATE_SPIN_SCIENCE_REVIEW.md
- New machine-readable target selection: STAGE8K_NONDEGENERATE_SPIN_CANDIDATES_20261022.json
- Deterministic zero-job filter: stage8k_existing_spin_prescreen.py
- Original historical evidence and hashes: STAGE8I_AUDIT_RAW_20261022.json, STAGE8I_SPIN_TEXTURE_20261022.json, STAGE8I_ISING_SPLIT_20261022.json
- Real old PBS IDs: AB 491186.cmspbs1, half-slide 491187.cmspbs1. Discover actual existing output paths from your authorized Japan PBS/SSH environment and persisted resume state; do not guess or overwrite them.

## Autonomous work required

1. Fetch repository branch to recent HEAD. Read source documents and run the offline standard-library Stage8k filter yourself. Confirm it selects nondegenerate AB pseudospin at actual k=(-1/3,+1/3) (NOT Kprime), (0.476190,0) near M, and (0.142857,0.142857) on Γ-K; emphasize AB band35 at these three k and the true stage8h midpoint comparison. Treat glide midpoint old falsely named Kp band35/36 as exactly degenerate, not a valid single-eigenvector spin indicator.
2. Use existing authorized SSH session to READ ONLY access original Stage8i EIGENVAL/WAVECAR, verify files' SHA and record lengths. Do not request POTCAR bytes, SSH keys, licenses or HPC site metadata in chat or GitHub. Confirm EIGENVAL energies and wavefunction band indexing are from the SAME calculation. The old reader incorrectly parsed WAVECAR raw band energies. Write a second, independent reader or perform a documented independent library/raw-format check, and show that its reconstructed eigenvalues agree with EIGENVAL, or explicitly mark BAD_WAVECAR_PARSER and stop.
3. For each selected nondegenerate state, verify k coordinates, BOTH spinor channels and G-vector mapping, phase-invariant spin expectations, local band isolation, time reversal spin reversal at k and -k, and consistent normalization under the correct pseudo/PAW overlap. If full physical PAW spin cannot be reconstructed reliably (PROCAR truncated), say so precisely and only report pseudized values with their uncertainty.
4. For AB and glide midpoint same k, compare actual nondegenerate state spin and energy splitting; DO NOT identify same band number across stacks with a single identical orbital character without wavefunction/orbital overlap tracking. Check the extra glide/mirror-imposed band degeneracies; single spin expectation of a degenerate pair is gauge arbitrary.
5. Correct old time-reversal valley mistake: real K=(1/3,1/3), real Kprime=(-1/3,-1/3), old mislabeled Kp=(-1/3,+1/3) is at 2/3 Γ-M. Compare full EIGENVAL sorted eigenvalue vectors at true K and -K, and spinor spin vectors where trustworthy. Do not re-label same-spin 87.54/35.31 meV band35–36 gaps as an Ising splitting.
6. Return a strict evidence table per key band/k: actual gap to nearby bands, validated spinor pseudo Sx/Sy/Sz, whether PAW augmentation was included, WAVECAR read validation, exact physical evidence vs unverified proxy, and supported conclusion on [H,Sz] nonconservation. This conclusion **must** separate physical spin mixing from full material static Gamma: the latter needs physical layer spin, velocity vertices, complete occupied-empty transitions, full BZ and Stage7j gauge.
7. Commit ONLY sanitized re-analysis scripts, checked derived spin matrices/energies, source checksums and concise Chinese physical report on the INDEPENDENT BSF research branch. Keep original external report with linked erratum and PR #3 Draft; do not rewrite prior raw external audit files or the Jia original reproduction. No POTCAR, WAVECAR, PROCAR, SSH secrets or large raw outputs in GitHub.

## Hard stop / reporting

- If valid physical PAW spin is unavailable, report BLOCKED_PAW_SPIN_OPERATOR. If independent WAVECAR reader fails, report BAD_WAVECAR_PARSER. Retain already known EIGENVAL sampled positive gaps but reject all unreliable spin claims.
- If a correctly reconstructed isolated band has demonstrably nonzero **physical** Sx or Sy (beyond numerical tolerance), you may conclude [H,Sz]≠0 at that k (hence not globally conserved); **do NOT conclude the interband Kubo response or NMEE is nonzero**.
- If this can only be seen in pseudo-wavefunction spin with incomplete PAW, report UNVERIFIED_SPIN_MIXING_PROXY even when numerical values look compelling.
- Do not ask the user to transfer files. If automatic authorized SSH access fails, report the single actual access blocker and preserved resume instructions.
- **Authorized new VASP/PBS jobs: ZERO.** No new NEB, Wannier, static NMEE or 12-point sweep. No new physics results may be fabricated.
