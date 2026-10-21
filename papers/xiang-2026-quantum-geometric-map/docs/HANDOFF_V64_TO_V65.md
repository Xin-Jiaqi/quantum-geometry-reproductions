# v64 -> v65 scientific handoff

Scientific audit of Xiang et al. PRB remains **PARTIAL**. v64 addresses a **source-independent** QWZ edge-device method: the positive-weak-coupling virtual sensor limit, not the complete source orbital mixed E/B response.

**Starting point:** v63 verified that top/bottom floating probe couplings scaled by lambda give lambda² selfenergies and finite-lambda local current rearrangement. In v64, do NOT simply set lambda=0 in four-terminal Landauer–Büttiker, because then sensor voltages are undefined.

**New code and method:** use the two-drive-only G0, with left/right physical QWZ strip selfenergies, and compute unit-amplitude sensor weights w_pL/R=Tr[Gamma_p G0 Gamma_L/R G0†]. The lambda→0+ voltage is (w_pL VL+w_pR VR)/(w_pL+w_pR), and the Hall resistance is voltage difference divided by the unperturbed two-lead Caroli current. No sensor selfenergy enters G0. Test positive spectral weights/current before reporting a voltage.

**Bounded result:** m=-1, L=7, eta=1e-4, fixed width 3 gives Rxy(0+)=1.0021574603585706 h/q² and I_L=0.9956990025136869 q²/h under unit bias. Finite lambda=1e-4 differs by about -3.568654e-10 h/q², with O(lambda²) approach. L=5,7,9,11 gives 1.000984,1.002157,1.002729,1.002957, NOT a controlled thermodynamic extrapolation. Probe offsets still change Rxy even though I_L stays fixed. For width 2 versus 7, spectral sensor strength changes substantially but symmetric chiral voltage ratio barely changes, a *model-specific* result. Trivial m=3 Rxy is deliberately undefined; metallic mu=1.5 is a negative control.

**Acceptance:** 26 new independent gates, 137 standalone total v60–v64. Frozen JSON and reproducer are in qgmap. Report: QWZ_VIRTUAL_PROBE_ORDERED_LIMITS_V64.md. Release is complete only after final PR-head and post-merge main full GitHub Actions logs are verified; a passing local suite alone is insufficient.

**Next research:** if pursuing size scaling, keep a fixed physical contact fraction and explicitly order lambda→0+, eta→0+, L→∞. For real bulk vs NEGF, choose the same transport observable, include sigma_xx, scattering and conserving current vertices. Prioritize the original missing Wannier embedded orbital magnetic mixed E/B, with operator/contact/bath complete independent Kubo–Streda protocol. Do not convert charge Hall topology into evidence for original orbital response.
