# Independent science/code review of adcd75e: multiband IGMC gate

**Date: 2026-10-23. Scope: synthetic Python evidence only. Verdict: MAIN_POSITIVE_CONTROL_CREDIBLE / NEGATIVE_CONTROLS_NEED_FIX / MATERIAL-INTERFACE_GO_CONDITIONAL. No VASP, PBS or Wannier executed.**

## 1. Accepted evidence (from source, not independently re-running all 73 seconds)

The committed JSON, report and executable consistently record the four-band v2_positive direct Dxx/omega 0.0342937475, rank-2 Wilson slope 0.0341517445 at 41² and a **0.4141%** relative difference. Its mesh sequence decreases 1.572% → 0.724% → 0.414%; varying the Zeeman energy b from 1e-4 to 1e-3 eV barely changes the 41² result. Joint constant H/S/A gauge rotation leaves direct and Wilson responses invariant to ~2.8e-17 and ~1.1e-12, respectively. These are useful synthetic positive controls and **not** a real AgF₂ material response.

The source pinned Cui Round26 synthetic H(R)/S(R)/A(R) contract is reusable, but still only tests *model* joint gauge consistency. The reporter has not proven an authentic Wannier90 hr/spn/mmn/r conversion path or the real-material spin operator.

## 2. A definite NC4 programming/physics error: the 'no SOC, conserved spin' claim is invalid as coded

In scripts/test_multiband_igmc_gate.py, function build_no_soc:

- Actual `h_r[i] = np.kron(np.eye(2), h2[i])`.
- Actual `s_r[2,j] = np.kron(np.eye(2), [sigma_x,sigma_y,sigma_z][j])`.
- Both h2 and the declared spin Pauli matrices act on the **same second factor**, so they **do not generally commute**. This is not the stated `H=h2⊗I_spin`, `S=I_orb⊗sigma` spin-conserving construction.
- At k=0, h2 = beta sigma_z with beta=0.6 eV; the x spin operator is I⊗sigma_x. Consequently the 4×4 commutator has spectral norm **||[H,Sx]||₂ = 2 beta = 1.2 eV**. I independently computed 1.2 with NumPy, versus 0 for the correctly ordered tensor factors.

**Conclusion:** NC4's `D≈0` is a real result for its supplied fixture, but it **does not validate the stated no-SOC/conserved-spin zero theorem**. It might be an accidental symmetry cancellation. Fix by making H = kron(h2, I_spin) while S = kron(I_orb, sigma); explicitly assert max_k ||[H(k), Sx(k)]|| < 1e-10 eV before describing NC4 as a zero-spin-conserving control. Ensure that H/spin basis ordering is also consistent for ALL components.

This is a mandatory correction before treating the gate as an accepted test suite.

## 3. A reporting inversion in the NC2 omit-position-connection control

The source JSON has `NC2_omit_A.response_cov = -0.02535718098` and `response_bare = +0.03429374753` for **v2_spinZ_pos**. The Markdown table reverses these labels, stating 'cov +0.03429 → omit-A −0.02536'. The **absolute difference 0.05965093** is correct but the assignment of signs to physically correct and incorrect operators is reversed. Use the source code and JSON as authoritative, issue a dated erratum and do not silently rewrite historical recorded numbers. The positive v2_positive fixture has +0.03429, but that is a *different position connection*.

## 4. Two other checks are weaker than their description

1. **Direct gap vs true global insulating gap:** `gap_scan` checks min over k of epsilon_3(k)−epsilon_2(k). This is the *minimum direct gap*, not automatically the **indirect global gap** min_k epsilon_3(k)−max_k epsilon_2(k) required for a single global Fermi level at half filling. Calculate both; 'global insulator' only if indirect gap >0 for each tested b. If negative, label response a rank-two occupied-band subspace formal exercise, not a bulk band insulator, and revise the fixture explicitly before promoting claims.
2. **Wilson branch continuity:** the 'fail_closed' code presently checks gap floor and |det W|≈0.982 (actually it builds an **unembedded**, b=0 Wilson at N=21); it sets `no_crossing_asserted=true` **without explicitly testing Wilson-phase branch winding or finite-b branch jumps**. A determinant modulus near 1 is *not* a proof of phase continuity. Add wrapped phase differences across neighboring ky and ±b, and check the embedded Wilson actually used by the response; fail closed on unresolved phase branch jumps.

## 5. Two program-quality improvements

- Add assertions that the 2-band anchor reproduces the frozen numbers to a declared tolerance, that Wilson b-step convergence satisfies its own bound, and that the program exits nonzero for failed core physical checks. Presently the anchor values are saved, but the verdict logic does not include them; b-step results are likewise printed without a hard pass/fail criterion.
- The reported ad hoc orbital-Liouville referee (`+0.034304` vs direct `+0.034294`) has **no saved reproducible implementation in the committed script**. Preserve this as development-only evidence; do not promote to an independent third oracle until a standalone routine/test and traceable data are committed. 'Orbital Liouville' here refers to a model density-matrix approach, **not a calculation of physical orbital magnetic coupling**.

## 6. Scientific go/no-go after this review

**Keep the 4-band main positive result; downgrade full-suite gate GO to _CONDITIONAL_GO_PENDING_TEST_HARDENING_.** There is no reason to run DFT yet. The fastest useful next action is to repair NC4 and NC2, add a genuine global-gap/phase-branch assertion, and freeze the development Liouville referee as a reproducible test, using small-mesh NumPy only.

Once the corrected full suite actually passes, move to a **fail-closed material-operator interface specification** (not a material compute): file schemas, common gauge/provenance for H(R), physical projected S(R), A(R), Wannier centers, spinor SOC convention and counterexamples (mismatched spin gauge, position center double count, missing operator matrices, wrong units). No new VASP/Wannier without separate authorization.

**This report does not alter the original adcd75e files or the user's prior receipts.**
