# Multiband same-gauge IGMC gate (four-band synthetic, zero VASP)

**2026-10-23 ｜ Read-only synthetic Python/NumPy validation. Verdict: GO (for developing an authenticated material operator interface — NOT for costly DFT/Wannier; no material IGMC value implied).** Executable: `scripts/test_multiband_igmc_gate.py` (+ pinned reference `scripts/wannier_operator_contract_pinned26.py`); machine-readable results: `MULTIBAND_SAME_GAUGE_IGMC_GATE.json`. Wall time 73 s ≪ 15-min cap.

## 1. Reused source (pinned, not on this branch)

The existing four-band same-Wannier-gauge operator contract of the same repository, pinned at commit `cdac6418166466b5243e73d8eb291d6dea8d4772`:

- `papers/cui-2026-curvature-flux-zeeman-qgt/ROUND26_SAME_GAUGE_WANNIER_OPERATORS.md`
- `.../src/flux_lab/wannier_operator_contract.py` — **extracted verbatim** into `scripts/wannier_operator_contract_pinned26.py` (provenance header added; no edits to the module body)
- `.../scripts/audit_wannier_operator_contract.py`, `.../results/frozen/round26_same_gauge_operator_audit.json` (read as reference)

The contract provides an in-memory common-basis H(R), physical Pauli S(R) and full position connection A(R) with analytic/covariant derivatives and positive/negative gauge controls. No Wannier file importer, no material operators. The pinned fixture `synthetic_soc_bundle()` (4-band SOC toy, spin-1/2 × two orbitals) is insulating at half filling but its Dxx response is an **identically zero (degenerate symmetry) baseline** — it is used as an additional zero control (NC5), exactly as the task anticipates ("if source fixture is gapless, make a separately named, explicit synthetic modification" — here the fixture is not gapless but response-degenerate).

## 2. Synthetic fixtures (all explicitly named; NOT material data)

| fixture | construction | role |
| --- | --- | --- |
| `pinned_round26` | `synthetic_soc_bundle()` verbatim | contract validation + zero-D baseline (NC5) |
| `v2_positive` | pinned + on-site `0.08 kron(σ_y,σ_y)` + y-hopping `0.05 kron(σ_x,σ_x)` + position y-hoppings `0.02 I` | positive multiband Dxx case |
| `v2_spinZ_pos` | v2 + on-site position off-diagonal `0.3 kron(σ_x,σ_z)` | omit-A response-level control (NC2) |
| `no_soc_control` | spin-diagonal H = kron(h₂(q₁), I_spin), [H, S^x] = 0 | spin-conserving SOC-free zero-displacement control (NC4) |

## 3. Conventions (locked to the frozen two-band benchmark)

- Reduced units: ħ = e = 1; energy in eV, length in Å; Zeeman coupling H_b = H − b·S^x, b in eV.
- Direct (Eq. (12) shape, per ω): `Dxx/ω = mean_k Σ_{n occ, m empty} 2 Re(i ⟨u_n|D_k H|u_m⟩ ⟨u_m|S^x|u_n⟩) / (ε_m − ε_n)²`, with D_k H the **covariant** derivative (includes A via −i[A, H]); mean over the uniform fractional-q grid = BZ measure.
- Wilson polarization: row phase `γ(ky) = −Im ln det W`, `W = ∏_j M_j`, `M_j = ⟨u_occ(k_j)|T(k_mid)|u_occ(k_{j+1})⟩`, `T = exp(−i dk A^x)` the **A-embedded** parallel transporter; `P_x = (a₁/2π) mean_ky γ`; slope `= [P(+b) − P(−b)]/(2b)`. Link orientation ⟨u(k)|u(k′)⟩ matches the frozen benchmark; the transposed orientation was diagnosed as a sign error during development (see §7).
- **No fitted scalar factor.** All factors (a₁/2π, 2b) are analytic. Mesh and b-step are separately converged.

## 4. Results

**Two-band anchor (self-contained reimplementation of the Eq. (17) Dirac, a = 1 Å):**

| quantity | this gate (41²) | frozen round3 / review |
| --- | --- | --- |
| direct Dxx/ω | 0.015236647006 | 0.0152366470 |
| Wilson slope | 0.015182770723 | 0.0151827707 |
| relative | 0.354 % | −0.354 % |

The anchor reproduces the frozen benchmark at printed precision — the normalization and sign conventions are locked before any multiband claim.

**Main gate (`v2_positive`, b-step 1e-4 eV):**

| N | direct Dxx/ω (Å/eV) | Wilson slope (Å/eV) | relative discrepancy |
| --- | --- | --- | --- |
| 21² | +0.03429375 | +0.03375450 | 1.57 % |
| 31² | +0.03429375 | +0.03404560 | 0.72 % |
| 41² | +0.03429375 | +0.03415174 | **0.41 %** |

Converges with mesh like the frozen benchmark (4.4 % → 0.35 %); b-step 1e-4/5e-4/1e-3 eV at 41² give Wilson slopes 0.0341517/0.0341518/0.0341519 (stable, discrepancy 0.41 %). Gap: min 0.985 eV, max 1.855 eV (all sampled k, b ∈ {0, ±1e-4, ±1e-3}). An ad hoc orbital-Liouville referee (development check, N = 21/31) gives Im K(ω)/ω = +0.034304 vs direct +0.034294 (0.03 %) — the current-response route agrees with the direct geometric sum.

**Joint gauge invariance (constant unitary on H/S/A, seed 26):** direct residual 2.8e-17, Wilson residual 1.1e-12 — invariant to numerical precision.

**Negative controls (each must produce a detectable incorrect result):**

| control | correct | wrong | detection |
| --- | --- | --- | --- |
| NC1 rotate H/A but not S | +0.03429 | +0.17040 | |Δ| = 0.136 ≫ 0.02 ✓ |
| NC2 omit A (bare velocity) | operator max|∂H−D_kH| = 0.799 eV·Å (>0.005) ✓ | response on `v2_spinZ_pos`: cov +0.03429 → omit-A −0.02536, |Δ| = 0.0597 (>0.01) ✓ |
| NC3 orbital pseudospin | +0.03429 | +0.31591 | |Δ| = 0.282 (≫ 5 %) ✓ |
| NC4 no-SOC spin-conserving | D = −5.1e-18, slope = 0.0 | — | exactly zero ✓ |
| NC5 pinned round26 | D = −2.6e-17 | — | exactly zero ✓ |

Note: on `v2_positive` the covariant-vs-bare velocity response is model-insensitive (Â response channel pointwise zero for that position structure); the spin-Z position connection fixture `v2_spinZ_pos` provides the response-level omission control.

**Fail-closed:** gap floor 0.02 eV enforced on all grids/fields (raises on violation); |det W| ∈ (0.9817, 0.9823) — no ambiguous/non-unitary Wilson branch; no band crossing.

## 5. Verdict

**GO** — for developing an **authenticated material operator interface**:

1. direct multiband Dxx/ω ↔ finite-field rank-2 Wilson slope match at 0.41 % (41²) with proper mesh and b-step convergence, at the same normalization as the frozen two-band benchmark;
2. joint gauge invariance at machine precision;
3. all negative controls produce detectable incorrect results;
4. fail-closed criteria hold.

Per the task's go criterion this is **not** GO for costly DFT/Wannier and **not** evidence of a nonzero real AgF₂ IGMC. The material operator chain for AgF₂ (SOC spinor Wannier H(k), same-gauge S_b with g·μB, orbital −B·L, multiband Wilson cross-check on a matched model) remains missing per the previous operator feasibility matrix.

## 6. One recommended next step (highest value)

**Build the matched-model falsification pair in the same framework**: finite-field non-Abelian Wilson polarization derivative vs static low-ω displacement kernel (orbital Liouville / Kubo) on the SAME synthetic fixture after fixing SI units — the Liouville referee (+0.034304 vs direct +0.034294 on `v2_positive`, 0.03 %) is already implemented ad hoc and should be made a permanent gate check; then re-verify the `-42'm'` SM row (1/116 residual) against a fresh SM page when the pinned SM PDF becomes available. No new DFT authorization is requested or implied.

## 7. Development notes (fail-closed honesty)

- The first det-Wilson implementation used the transposed link orientation ⟨u(k′)|T|u(k)⟩ and produced −0.03415 vs direct +0.03429 — a sign error diagnosed against the frozen benchmark's ⟨u(k)|u(k′)⟩ orientation (single-link phase anatomy documented in the working notes). The corrected orientation matches.
- The two-band anchor initially used a [0, 2π) grid and a mis-indexed einsum; both fixed (grid [−π, π)/a per `dirac_grid`; einsum `…ai,…ab,…bj`). Final anchor reproduces the frozen CSV at printed precision.
- `-42'm'` supplemental row status unchanged: not resolved without the primary SM PDF.
