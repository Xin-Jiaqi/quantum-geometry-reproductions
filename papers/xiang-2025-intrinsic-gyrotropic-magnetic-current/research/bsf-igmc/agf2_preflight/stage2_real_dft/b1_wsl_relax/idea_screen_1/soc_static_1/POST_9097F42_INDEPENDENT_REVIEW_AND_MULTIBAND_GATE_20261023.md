# After commit 9097f42: independent audit and highest-value next physics gate

**2026-10-23 | PR #38 | Read-only independent review; no VASP/PBS/Wannier or paid compute.**

**Decision: ACCEPT PARTIAL as the correct global project status. SOURCE CONVENTION is strongly supported but not publisher-SM-fully-verified. Do not rebuild existing two-band Wilson tests. Next useful step is a *multiband same-gauge mixed Zeeman-response falsification gate*.**

## 1. What was actually verified

The four new PR files are `SOC_MXPRIME_INDEPENDENT_RECHECK.json`, `IGMC_EQ14_SOURCE_CONVENTION_AUDIT.md`, `IGMC_RESPONSE_OPERATOR_FEASIBILITY_MATRIX.md`, and `BSF_IGMC_GO_NO_GO_AFTER_SOC.md`, reported in commit 9097f428d6ab61917a723e5d28e1970e22baafa6.

- The 12-site \(m'_x\) magnetic pairing with spin residual ≤0.001 μB and orbital residual 0.000 μB agrees with the independent prior audit. The common weak \(M_z\) and local Ag \(L_z\) are **allowed**, not forbidden, by \((M_x,M_y,M_z)\to(-M_x,+M_y,+M_z)\). Full spin-density magnetic space group, ground state, easy axis and response still not certified.
- The main paper's **printed Eq. (14)** uses one antiunitary \(\eta_T\) for joint \(\sigma^{C/D}\), while the existing source inventory/SM-table transcription and model parity strongly favor **\(\eta_C=+1,\eta_D=-1\)** for the separate real quadratures. Under the \(m'_x\) B1 hypothesis, \(D_{xx}\) is allowed under **both** compared conventions. Do **not** promote “115/116” to source-perfect confirmation: it compares against a **transcribed** SM and has one unresolved \(-42'm'\) D-row; the originally pinned SM PDF was not retrieved. A public arXiv v2 main text is independently available (arXiv:2312.04093), but not identical in file hash to the formerly pinned publisher article.
- \(C\) is the clean in-gap-vanishing Fermi-surface contribution. A nonzero \(\sigma^D_{xx}/\omega\) is **not** implied by symmetry or by the SOC 1.3616 eV occupation-gap screen. The same-gauge physical-spin/Berry-connection matrices remain missing *for actual AgF₂*.
- No real-material IGMC, source-validated full B1/B2 stability, or material-specific response amplitude exists.

## 2. Independent model calculation: avoid repeating already closed tests

I independently implemented the **same analytical two-band Dirac Hamiltonian** and direct velocity-spin geometric kernel with a **separate finite-field occupied-band Wilson loop** in a local NumPy calculation. This is a manual independent implementation/check of the frozen equations, **not a full rerun of the actual repository Python modules or CI**. Parameters: a=1 Å, β=0.25 eV, v=0.825 eV Å, M=0.1 eV Å², μ=0, \(H_b=H-b\sigma_x\), central field step \(b=10^{-4}\,\mathrm{eV}\), real/frequency-free \(D_{xx}/(\hbar\omega)\) convention as implemented in the source (reduced units). Results:

| N×N grid | direct geometric \(\sigma^D_{xx}/\omega\) | independent finite-\(b\) Wilson derivative | relative error |
| --- | ---: | ---: | ---: |
| 21² | 0.0161234760 | 0.0154084315 | −4.435% |
| 31² | 0.0152809661 | 0.0151633699 | −0.770% |
| 41² | 0.0152366470 | 0.0151827707 | −0.354% |
| 61² | 0.0152340859 | 0.0152103618 | −0.156% |

The 41² Wilson number agrees with the already frozen repository [Round-3 convergence CSV](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/research/agf2-b1-wsl-relax-20261021/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/results/raw/round3_wilson_convergence.csv) at the printed precision; that repository itself documents further 71²→301² convergence to a relative Wilson/direct difference below 0.01%. This is **existing tested functionality**, not a newly discovered numerical result for AgF₂ or BSF.

Moreover the *BSF toy eight-band model* already has a documented four-occupied-band **non-Abelian** Wilson finite-field check for its own toy \(D_{xx}/\omega\), with finite-mesh error scaling. Thus simply asking for “finite-field Wilson vs \(D_{xx}\) on the same model” duplicates **two already delivered baselines** and is **not the next highest-information experiment**.

## 3. Found useful existing infrastructure outside this PR's branch

There is a directly relevant, **different line of work in the same repository** at pinned git commit

`cdac6418166466b5243e73d8eb291d6dea8d4772`:

- `papers/cui-2026-curvature-flux-zeeman-qgt/ROUND26_SAME_GAUGE_WANNIER_OPERATORS.md`
- `papers/cui-2026-curvature-flux-zeeman-qgt/src/flux_lab/wannier_operator_contract.py`
- `papers/cui-2026-curvature-flux-zeeman-qgt/scripts/audit_wannier_operator_contract.py`
- tests in the corresponding Cui package.

**Important scope:** that is an **in-memory, synthetic four-band H(R), physical projected Pauli spin S(R), position A(R) joint operator contract**, with covariant \(D_kH=\partial_kH-i[A,H]\), controlled unitary-gauge invariance, Wilson finite-difference checks and explicitly failing *wrong-spin-gauge / missing-A* negative controls. It reports local QGT residual 2.09×10⁻¹⁵ under a joint constant rotation, and a large error under mismatched spin-operator gauge. **It does not import authentic Wannier90 .hr/.spn/.mmn/.r files or authenticate material-spin provenance; it does not deliver an AgF₂ IGMC coefficient.**

Because this component was visible at a different pinned commit rather than this PR's current base, **do not assume these paths are present under PR #38 or merge/rebase whole unrelated work without permission**. Treat it as an upstream reusable component with explicit pinned provenance.

**Updated operator readiness language:**
- **Synthetic four-band same-gauge H/S/A input + derivative/gauge negative controls: EXISTS at pinned Cui snapshot**.
- **Actual AgF₂ spinor Wannier operators, certified source gauge, physical orbital driving and quantitative material C/D integration: STILL MISSING**.
- **Full multiband, frequency-consistent IGMC response validation on the Cui synthetic bundle: NOT YET DEMONSTRATED**.

## 4. Concrete next test (zero VASP)

**Single best next physics experiment:** using the **pinned, synthetic Cui four-band operator bundle**, add an isolated read-only/few-second numerical test of *insulating, low-frequency spin-Zeeman displacement response*:

1. First show the candidate synthetic four-band model is gapped at fixed filling over the selected 2D k mesh and at \(\pm b\), and that its projected physical-spin operator is the one declared in the bundle. If not, adjust **only an explicitly documented synthetic toy parameter** or create a fresh minimal 4-band bundle with a known gap, never silently repair real materials.
2. From the **same bundle** and same gauge construct \(H(k)\), \(S_x(k)\), covariant \(D_{k_x}H\) and Berry/position connections. Compute multiband \(D_{xx}/\omega\) using occupied–empty geometric terms, with clear denominators, k units, and Zeeman-field energy \(b\). Independently compute a **rank-2 occupied non-Abelian Wilson-loop charge-polarization derivative** at \(\pm b\). Check sign, normalization, mesh and field-step convergence.
3. Apply one **constant joint** random unitary to all H/S/A: physical \(D_{xx}\) and Wilson slope must be invariant. Then deliberately rotate H/A *without* S and deliberately omit the connection A; both negative controls must produce a significant detectable error rather than accidentally pass.
4. Recover the existing two-band Dirac benchmark at a suitable reduction, and add one clean **no SOC / conserved spin** control where the clamped-texture spin-only displacement is zero. Different metallic/degenerate cases must fail closed.
5. **Output:** one small JSON/CSV and a concise Markdown go/no-go, with pinned operator source and reproducible code/test command in PR #38. No external file import, no new DFT, no material IGMC claim, no `qsub` or Wannier. Do not silently replace prior paper-source Eq. (14) ambiguity.

If this multiband joint-gauge **positive** control and all **negative** controls pass, the remaining engineering task is authentic material H/S/A export and its verification. If they fail, stop and repair the operator/gauge theory **before** any material-Wannier calculation.

## 5. Research narrative after audit

**Most defensible hypothesis:** Type-I sliding connects a pair in which a finite insulating displacement \(D_{xx}\) changes sign while an out-of-plane axial magnetic signature (weak \(M_z\), orbital \(L_z\), possible AHE) remains even. Strained Type-II would invert the partner parity pattern. This is a symmetry fingerprint **conditional on actual partner states**, and remains unmeasured.

The **minimum publication-quality nontrivial evidence** is not another symmetry matrix or SOC single-point, but a **material finite, converged \(D_{xx}\) in two physically grounded stacking partners**, with an honest spin/orbital/induced-electric-field distinction and feasible detection estimates. The single B1 SOC screen does not establish this.

**Decision:** `PARTIAL → LOW_COST_MULTIBAND_OPERATOR_GATE` (no materials compute). Never assert a new go-for-expensive-DFT solely from the 115/116 comparison or from an allowed \(D_{xx}\).
