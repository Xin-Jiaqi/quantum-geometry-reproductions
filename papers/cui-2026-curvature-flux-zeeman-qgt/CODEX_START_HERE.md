# CODEX START HERE | Cui v3 next independent research shift
**ONE-ENTRY LOCAL CODEX HANDOFF — updated 2026-10-23.** Do not ask the user to paste prior rounds or raw datasets again.

**Read first:** [FULL_REVIEW_2026-10-23.md](FULL_REVIEW_2026-10-23.md) (scientific status, false positives, research ideas). Then original source arXiv:2604.09725v3, [Issue #105](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/issues/105), open [PR #107](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/107), and [independent review](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/107#issuecomment-6364242247).

## Mission and physical target

**Complete a correct, independently tested v3-specific theory replication, prioritizing finite-frequency antiunitary Table I + Dirac topology + exact Kubo/GHE/KME signs.** Do NOT work on new EPW, DFPT or a large material calculation to reproduce this specific *theoretical* article. Previous validated operator/Wannier/thermal models are useful optional extensions, not source claims.

Target model \(H=\mathbf d(\mathbf k)\cdot\boldsymbol\sigma\). Core observables \(r^a_{nm}\sigma^b_{mn}\), four symmetric/antisymmetric real/imag Zeeman sectors, Dirac winding and Hodge dual radial flux, clean interband Zeeman-current GHE / reciprocal electric-field magnetization KME, *actual* 2D magnetic-point-group selection of the finite-frequency tensor components.

Do **not** declare a source erratum unless the author's **v3 TeX** equation, our operator definitions, the antiunitary Hamiltonian, and an independent microscopic calculation all support the comparison.

## Current baselines and branch hygiene

- **Accepted main:** SHA **ac9b0cb9d5dce19baf55e6fe95d00428bab5ee53**, Python 3.11/3.12 **371/371 PASS**, postmerge GitHub Actions [#39916655754](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39916655754). All Round17–34 historic audits green.
- **UNMERGED PR #107:** SHA **c97dee0b0a39747b17451b53f70652dfbdf2dc7e** at review, Python **3.12 383 PASS / 3.11 FAIL** in Actions [#39942431072](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39942431072), due to Python 3.11 syntax error at zeeman_v3.py:741. The 3.11 run showed **372 tests + 1 import error**, **NOT 383 tests passing**. *Inspect CURRENT HEAD before acting; another agent may have updated it.*
- **Local Round35 Fe:** user-reported official W90 spinor postprocess; git-ignored receipt and raw data are not visible here. DO NOT delete, reset or upload them.
- **Other collaborators** are active in the same repository on unrelated papers and PRs. **Never execute git reset --hard, git clean, checkout that would overwrite the user's work, force push or raw upload**.

Recommended independent worktree (only if git and the directory are available; otherwise adapt paths and stop before destructive changes):
~~~bash
git fetch origin main research/cui-round36-zeeman-v3-20261023
git status --short
git worktree add ../cui-v3-codex-review -b research/cui-round36-fix-verified origin/research/cui-round36-zeeman-v3-20261023
cd ../cui-v3-codex-review/papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
~~~
If a similarly named branch or worktree already exists, choose a fresh name; do not remove/overwrite it. Do not assume this handoff PR is already merged until main shows the files.

## P0 | Correct the rejected evidence BEFORE adding new features

**Gate 0 — exact version/CI compatibility.** Fix the Python-3.11-unparsable interpolated f-string:
~~~python
# Incorrect for Python 3.11:
# rule[f'{name}{"'" if theta else ""}'] = {...}
label = name + ("'" if theta else "")
rule[label] = {...}
~~~
Compile/import and run the full test suite under **BOTH** 3.11 and 3.12. Add a reproducible `scripts/audit_zeeman_v3.py` that re-computes and compares the *corrected* Round36 frozen JSON, and call it in the exact Cui GitHub Actions workflow. A numeric frozen file that no CI step reads is not standalone scientific acceptance.

**Gate 1 — physical antiunitary group action, not a static sign guess.** For every antiunitary source/symmetry model, independently validate the **Hamiltonian** before trusting local tensors:
\[
Q=(-1)^\theta R,\quad \Lambda=(-1)^\theta\det(R)\operatorname{diag}(R,1),\quad
 U_A H(\mathbf k)^* U_A^\dagger=H(Q\mathbf k).
\]
\(\theta=1\) for antiunitary; use existing **src/flux_lab/antiunitary_response.py** and **src/flux_lab/group_microscopic.py** from merged R15/16 as cross-checks. Physical spin transforms as an **axial TR-odd** vector and the complex frequency needs the correct \(-z\) *advanced* branch under conjugate antiunitary covariance.

For **m′=m_x T**, \(Q=\mathrm{diag}(1,-1)\) and \(\Lambda=\mathrm{diag}(-1,1,1)\). The unmerged source's claimed counterexample \(d=(k_x,k_x k_y,m+k^2)\) **DOES NOT satisfy symmetry**, so its nonzero \(\langle g^N_{xy}/\epsilon\rangle\) is **not** a counterexample to manuscript Table I. Add the actual matrix-level negative test; add a symmetry-valid control \(d=(k_y,k_x,m+k^2)\). Repeat with several symmetry-valid warping terms if searching for nonzero transverse \(g^N\). Any conclusion about m′ remains **UNRESOLVED**, not 'SOURCE CONTRADICTION', unless a legitimate invariant response establishes it. Reclassify test/docs/frozen output accordingly.

**Gate 2 — remove universal C3 radial theorem.** C3 covariance does not require \(d_+=zf(|z|^2)\). Add invariant tangential test \(d=(-k_y,k_x,m)\) (whose \(\Omega_z^N=n(r^2+2m^2)/(2(r^2+m^2)^{3/2})\) is nonzero) and trigonal warping \(d_+=z+\lambda\bar z^2\). Ensure the code checks actual \(H(Rk)=U_RH(k)U_R^\dagger\) at random k and integral/quadrature outputs independently. Withdraw any claim that all Cn≥3 2-band models have zero normal curvature from one chosen radial model. Distinguish 'this representative is zero', 'symmetry-forbidden tensor', 'all 2-band models zero', and 'general multiband forbidden'.

**Gate 3 — physically reconcile Eq. exp/exp1 sign.** Make a concise source-convention ledger: \(H_B\), \(H_E\), Fourier \(e^{-izt}\) vs original omega, \(\epsilon_{nm}\), \(r^a_{nm}\), spin Pauli versus \(\hbar\sigma/2\), static \(\omega^0\) subtraction, \(\eta\to0^+\), Onsager reciprocal definition. Derive O(ω)/O(ω²) GHE and KME from both full complex spectral Kubo and independent orbital Liouville. Match actual **v3**, not just old v2 equation numbers. Report sign disagreement only if it survives all convention conversions.

**Gate 4 — actual Table I output.** Group by unitary versus antiunitary. Apply exact complex-frequency transformation to observables **current under magnetic B** and **magnetization under electric E**, not hand-assigned ±1 local sector signs; *track odd/even frequency separately*. Run complete 31-group **group closure and local Kubo covariance** from accepted R16; use a small number of genuinely group-invariant 2-/4-band periodic models to test full BZ, including C3/hexagon from R17. A symmetry-allowed *rank-2* response does not imply an allowed *symmetric transverse xy* component. Print source Table-I row versus actual realized sector with `PASS/CONDITIONAL/BLOCKED`, never 'source error' from a non-invariant test.

## P1 | Physics ideas only after P0 passes

**A. Geometry stability:** anisotropic velocities, mass, trigonal-warping and multiple-node annihilation; compare local flux/winding to Berry curvature and Wilson loop. IR regularization and cutoff dependence must be explicit; distinguish topological contour winding from smooth nonuniversal curvature.

**B. Independent Kubo–Liouville reciprocity:** source GHE four channels \(g^N,g^A,\Omega^N,\Omega^A\), KME four complementary channels; verify actual small-ω slopes with properly subtracted static terms, frequency-sign negative controls and occupied-band partitions. This is the **primary source paper**. Preserve old equilibrium static Ward (R23) as a distinct theorem: no bulk equilibrium charge current, not a ban on driven finite-ω KME.

**C. Multiband group test (optional, separate PR):** 4-band physical spin × orbital fully invariant magnetic models to distinguish 2-band accidental vanishing from symmetry-forbidden; many-body thermal/collision and EPW extensions are lower priority.

## Deliverables required, without further user handholding

1. One **corrected research PR** based on a fresh worktree/branch, or coordinated commits to existing #107 only if ownership and HEAD are known. Avoid force push and shared-directory writes.
2. Short **scientific acceptance report** (1–3 pages) with a table: v3 exact claim/equation/table, invariant model(s), independent math/algorithm, positive/negative controls, numerical residual, genuine status, source-vs-our-convention caveat.
3. Executable **Round36 audit script** in the dual-version workflow, numeric output frozen *from a passing independent run*. Preserve all previous frozen results. Include explicit invalid m′ and false C3 radial assertion as **rejection tests**.
4. GitHub check: exact PR SHA -> 3.11+3.12 fully green -> review negative controls -> safe merge ONLY after explicit scientific acceptance -> postmerge *same merged SHA* green. Never claim green merely because local 383 passes.
5. Final single concise report with **what is proven, what changed, remaining mathematical blocker, PR/CI links and a one-line best next idea**. Whole original Cui status stays **PARTIAL** until source Eq./Table validations truly complete.

### Existing code and source pointers

- **Merged:** [main README](README.md), **REPRODUCTION_PROGRESS.md**, **SCIENTIFIC_ACCEPTANCE.md**, **ROUND14_VORTICITY_AND_MAGNETIC_SYMMETRY.md**, **ROUND15_ANTIUNITARY_MICROSCOPIC.md**, **ROUND16_ALL31_MICROSCOPIC.md**, **ROUND23_EQUILIBRIUM_STATIC_WARD.md**, **ROUND27_FINITE_Q_GAUGE_WARD.md**, plus **ROUND34_NATIVE_SPN_SOURCE_CHAIN.md**.
- **Unmerged #107:** **ROUND36_SOURCE_SCOPE_AUDIT.md**, **ROUND36_CUI_V3_SCIENTIFIC_ACCEPTANCE.md**, **src/flux_lab/zeeman_v3.py**, **tests/test_table1_zeeman_v3.py**, **results/frozen/round36_zeeman_v3_table1_audit.json**.
- **Source:** https://arxiv.org/abs/2604.09725v3 ; confirm actual revised title/version (last revised 2026-09-17). The source v3 TeX used by the local harness is git-ignored, not a checked-in public file.
- **Review thread:** https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/107#issuecomment-6364242247 .

**User preference:** work independently, do not request many copy/pastes, keep the report concise but do not hide wrong physics or failed CI.
