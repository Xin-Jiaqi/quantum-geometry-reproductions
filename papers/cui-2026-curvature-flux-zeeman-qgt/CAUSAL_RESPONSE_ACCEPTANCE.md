# Cui 2026: causal interband response, v2 scientific acceptance

**Status: named MODEL/METHOD VALIDATED; full paper remains PARTIAL.**

**Source:** Cui et al., *Curvature-Flux Representation of Local Dirac-Node Topology from Non-Hermitian Zeeman Quantum Geometry*, arXiv:2604.09725v3, main Eqs.(16)–(19), Appendix B Eqs.(90)–(126). Independently checked against the local 15-page paper; no publisher PDF/figures included.

## Physical input and conventions

- H(N,N), two D_k H(2,N,N), and three *physical dimensionless Pauli spin* matrices S(3,N,N) are supplied in the SAME fixed, orthonormal orbital basis. A missing/mismatched Wannier spin or position matrix is **not inferred** from hr.dat. All operators must be Hermitian; unresolved degeneracy fails closed.
- H_B=-κ B·S; H_E=+E·r; J_a=Tr[(D_a H)ρ]; M_b=κ Tr(S_bρ), ħ=1, *natural units*. This M normalization is declared model convention, **not automatically the paper's physical SI magnetization**.
- Only n≠m interband position r_nm=V_nm/[i(ε_n−ε_m)]. Intraband electric acceleration, Fermi-surface terms, collisions, full periodic position, contacts, degenerate bands, SI conversions and real-material operators are excluded.
- Source Appendix B uses F[e^{-i ωβ t}]=δ(ω+ωβ); hence **ω_source = −ω_drive**. source_four_sector_* takes ω_source, Kubo/Liouville takes ω_drive.

## Independent physical gates

| Gate | Independent comparisons | Result |
|---|---|---|
| Finite-frequency J/B, M/E | orbital-basis N² Liouville vs directed eigenband Kubo sums | 2.85e−14 and 7.91e−15 max on selected 4-band benchmark |
| Original Eq. (16), O(ω²) | explicit four gN,gA,ΩN,ΩA contributions vs direct complex T_nm band pairs; both vs source-frequency Kubo | four-sector match at floating precision, remaining error scales O(ω³) |
| Local interband reciprocity | independent Kubo χ_JB(z) and α_ME(−z); static subtraction at same η | χ_JB(z)−χ_JB(z0)=i[zα_ME(−z)−z0α_ME(−z0)]ᵀ, z=ω+iη, z0=iη; max residual below 1e−12 |
| Wrong protocol controls | delete static KME term; substitute retarded alpha(−ω+iη); mismatch spin orbital gauge | demonstrably nonzero residual |
| Fixed number comparison | frozen Round8 values and independent whole-BZ synthetic midpoint sum | local Im Δχ_yz=−0.0310000213902597, Re=0.0002776323787547; 48² BZ Im Δχ_yz=4.0616172620068635e−5 |

The **full interband α_ME** on the reciprocal right side must retain its static term; removing the zero-frequency part from *both* channels and asserting the same simple equality is not justified for this particular local observable construction. At η>0, alpha(−z) uses the *advanced* analytic branch. Neither statement is an assertion that the source paper's Eq.(19) is wrong, or proof of a full experimental Onsager relation.

## Runnable acceptance

~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/demo_causal.py --kmesh 48
~~~

**Total 39 tests**: original 21 geometry and 18 new causal/sector/negative-controls. CI runs Python 3.11 and 3.12. The 11-round (185-test) private/local research archive covers *different, broader* model protocols and is not uploaded or counted as this public verification.

## Remaining gates

Original Eq.(19) exact tensor-by-tensor author-normalized KME, full metallic diagonal electric response, many-band degeneracies, a physically conserving collision operator, electromagnetic contact and magnetization currents, experimental sheet/SI KME, and true Wannier spin/position matrices are still open. This release is therefore **a reproducible mathematical interband capability, not a finished bulk KME solver**.
