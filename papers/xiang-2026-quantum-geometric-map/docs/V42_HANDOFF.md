# v42 handoff | Xiang et al. PRB 113 L201406 (2026)

Scientific acceptance remains **PARTIAL**. v42 adds the full **exact finite-disorder, finite-volume linear Hall Smrčka–Středa split** using exactly the v41 torus Hamiltonian/currents/disorder, explicitly calculating I(a) and II, and verifying their sum against the independent full Lehmann Hall as external \(\eta\to0\). The Hall-odd I(b) cancels. These I/II components have an overlap and are not uniquely separable physical mechanisms. Independent dense matrix inversions/energy quadrature check the eigenstate primitive; clean W=0 also matches v38.

Only in the fixed finite model and zero external regulator limit is this a **correct linear Hall reference**. It is not a conserving disorder-averaged SCBA Středa response; no static nonlinear EB PHE, Fig.2(b), multiband orbital completion, or real-material claims.

See `docs/EXACT_STREDA_V42.md`, ten tests, script CLI, frozen gap/metal/ensemble JSON. `L=4,W=.12,seed=1947,theta=(.11,.33),mu=1.1`: Středa total at eta=.003 `-0.306993647`, exact v41 `-0.307005888` e²/h. Gap mu=0 at same eta has error `8.33e-8`.

Next v43: derive a **single conserving disorder average protocol** including energy-dependent Sigma^R,A, RA/RR current and density vertex structures, explicit Středa I/II, contacts; compare against ensemble v42 with size/seed/broadening scaling. Do not splice v39 sea and v40 RA without deriving response prefactors/vertex/limit identities. Original PRB second-order EB remains OPEN. Do not contact authors or modify other paper modules. Report science/code/remaining rounds each turn.
