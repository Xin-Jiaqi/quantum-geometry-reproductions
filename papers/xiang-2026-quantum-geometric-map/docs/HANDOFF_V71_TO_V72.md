# v71 → v72 | Source-Observable Boundary and Exact Contact Ward Handoff

Date 2026-10-22. Repository Xin-Jiaqi/quantum-geometry-reproductions, Xiang PRB 113 L201406 (2026). **Scientific status PARTIAL**.

## What changed and why this matters

v60–v70 focused on source-independent finite orbital/charge models, scalar gamma Liouville, microscopic semi-infinite lead Keldysh, independent full-Hermitian finite-lead Lehmann and L/eta protocol. v71 demonstrates why those successful finite **device polarization-current** responses CANNOT be identified with the author's bulk spin PHE/orbital MNHE without a separate physical current/field-limit mapping.

Given Ptilde=P_device direct-sum zero leads, the **FULL** finite-system Heisenberg derivative is i[H_device+H_leads+H_contact-E Pb,Ptilde] = J_device + i[H_contact,Ptilde]. The second contribution was absent from the *definition* of the prior device-only observable. It is not zero in general.

At L=22 leads, T=.37, mu=.13, B=.23, isothermal static dE<J_device>=+0.002996376694233, dE<J_contact>=-0.002996376694233, and dE<J_full> approximately zero. Exact finite-size thermodynamic identity Tr(f(H(E))[H(E),Ptilde])=0; full independent E rediagonalization confirms the current and its derivatives. This is NOT zero Landauer conductance nor zero infinite-device transport. At omega=.29 eta=.025, complete finite-frequency J_full response is .030134799972-.017968308910i, versus J_device -.005310644134-.007389851857i, and a nonzero missing contact .035445444107-.010578457054i. Full finite spectral Ward residual ~1e-16 in three magnetic/orbital frames.

SM source distinctions, taken from primary-source SHA-pinned v49/v53 published repository audits (not a new PDF reading): spin PHE couples to spin Zeeman; orbital MNHE to orbital magnetic positional shift and interband geometry; neither is the finite-device internal Peierls polarization-current coefficient. Eq20 direct velocity contraction = +tau*sum Z*df whereas printed Eq31 = -tau*sum Z*df for the specified source epsilon/v conventions. Hermitian random 4-band independent check gives direct -0.2039672645 and printed +0.2039672645. No author correction inferred. Main Eq10 vs SM Eq62 helicity switch must also rotate physical spin. v53's independently found orbital Eq49/52 relative metric-factor issue and authors' omitted three-band interference are separate and remain unsolved.

## v71 acceptance

- **VERIFIED:** exact finite bounded polarization contact current identity; isothermal static cancellation, independent finite-E/B re-diagonalization, optical Lehmann Ward, three-frame gauge covariance; literal SM Eq20/31 indexed contraction and spin-helicity basis rotation.
- **CONDITIONAL:** 4-site supplied finite graph with two finite leads, device-only bounded P, finite T, omega/eta and arbitrary Hermitian spin-index test. No bulk response, physical lead-terminal DC, or source solution inferred.
- **UNRESOLVED:** printed source sign/factor/three-band issues, author relaxation pole and full source tensor, material Wannier position and magnetic embedding, bulk finite-q/omega/L limit, Fig2 SI source normalization.
- 42 physics gates + one independent frozen numerical regeneration, 43/43 local PASS. Exact full remote GitHub CI and merge SHA should be recorded only after real Actions success.

## Physical capability integration status

A presence-only *fail-closed* input inventory now distinguishes finite_device, spin_phe, and orbital_mnhe. It demands actual Zeeman spin for spin_phe; the orbital pathway requires full Wannier/interband positions, orbital moment and magnetic positional shift plus BZ/units. No input inventory grants solver authorization or proof of source reproduction. Do not fill missing orbital data with diagonal Wannier-center approximations.

## Progress, management estimates not completeness claims

GitHub reproducibility engineering ~90–95% (v71 PR, main CI and portable archive still release gates); finite model physical API ~75–80%; entire original PRB source paper ~30–40%. These are subjective **planning assessments** not physics-accuracy percentages. v71 materially sharpens a no-go condition but does NOT raise the whole-paper completion figure.

## v72 recommended project

1. Start from the actual original source mixed tensor's **Bloch-band current definition**, Zeeman spin versus orbital B, finite-q/omega and selected eta_a pole. Build source-observable interface without inferring it from finite device J.
2. Derive a **full periodic bulk finite-q source response**, keeping position/magnetization versus transport and contact subtraction explicit. Ask: when does the source spin PHE reduce to a steady state, and can any orbital MNHE device limit be equated with the published expression?
3. Supply real Wannier90 Hamiltonian, intercell position and spin/orbital magnetic matrices. Require actual field derivatives and unit conversion; fail closed for missing projection.
4. Return to printed SM Eq20/31, Eq49/52 (factor 2 vs 1), omitted m != l three-band contribution. Preserve literal original branches as simultaneous acceptance tests; no fitted signs or factors.
5. Independently validate macroscopic BZ-normalized response and published Fig2 under documented SI conditions. Preserve source original **PARTIAL** unless all specified gates close.

Paths: src/qgmap/source_observable_bridge_v71.py, tests/test_source_observable_bridge_v71.py, tests/test_source_observable_bridge_v71_frozen.py, scripts/reproduce_v71.py, results/frozen/v71_source_observable_bridge.json, docs/SOURCE_OBSERVABLE_BOUNDARY_AND_WANNIER_GATE_V71.md. Frozen source PDF original hashes are pinned in v49; copyrighted PDFs are not redistributed.

**GitHub release proof:** append exact PR number, final HEAD, PR CI, merge SHA and post-merge main CI after verified; do not invent results.


## v71 final GitHub release certification (2026-10-22)

- [PR #79](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/79) **squash merged** as `e946d841cdb32c849a66f401a853767895f0e9b7`; final PR HEAD `7d8540be79fec5f3995d47f81c07475b4e6265e4`.
- Final PR-head [Actions 39757723955](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39757723955): **969/969 PASS**, `Ran 969 tests in 154.511s`, completed/success, original source scientific status explicitly `PARTIAL`.
- Independent merged **main** [Actions 39758266016](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39758266016): **969/969 PASS**, `Ran 969 tests in 198.909s`, completed/success on exact merge SHA `e946d841cdb32c849a66f401a853767895f0e9b7`, original source `PARTIAL`.
- 42 v71 scientific checks plus 1 independently reconstructed frozen-value gate; independent local and fresh-extracted **43/43 PASS**. Earlier v60–v68 standalone baseline **306/306 PASS** in a bounded batch and v70 independent standalone 10/10 PASS. A single local all-versions combined run exceeded its execution cap, so **do not claim local 100% suite PASS**; GitHub complete 969/969 core suite is the release authority.
- Standalone archive (after certification regeneration): 96 files, 95/95 SHA256 verified; v71 physical source/test/reproduce/frozen-test Git blob SHA matched independently tested local versions; GitHub frozen JSON has same numeric values but different ordering/serialization.
- **VERIFIED** bounded finite closed-system J_full=J_dev+J_contact, isothermal cancellation and optical polarization Ward, plus source-authentic SM Eq20/31 signed index contraction and main/SM spin helicity. **CONDITIONAL** finite supplied Peierls/position system. **UNRESOLVED** physical bulk PHE/orbital MNHE response map, full magnetic Wannier embedding, original source orbital factor/3band and SI concerns; original Xiang PRB remains **PARTIAL**.
- Final release-documentation-only commits after merged CI modify no source code or tests.
