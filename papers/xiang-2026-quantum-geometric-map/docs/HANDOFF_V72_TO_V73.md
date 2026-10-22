# v72 → v73 Source-grounded physical research handoff | 2026-10-22

Original: Xiang et al., PRB 113 L201406 (2026); source paper status **PARTIAL**. Source-specific citations refer to earlier original SHA-pinned main and SM audits v49, v53 and v71. These were not newly re-read in v72.

## New result

v71 found a missing lead-contact term for the complete bounded finite-device Heisenberg polarization current and proved it cannot be directly identified with bulk DC transport. v72 **changes the geometry and driving field**: a periodic 2D two-band Bloch lattice, finite commensurate q_y scalar potential and a physical-model spin Zeeman source B sigma_x. It probes χ_{jx,rho}(q,omega;B), normalized by E_y=-i q_y φ_q, and its B derivative. It does NOT substitute orbital Peierls field or assert source-identical steady-state spin PHE.

Exact longitudinal density Ward uses J_y=[H(k+q)-H(k)]/q and the inter-momentum density vertex W=U(k+q)†U(k). For periodic q=2πs/N the integrated f(H(k+q))-f(H(k)) trace cancels, yielding (omega+i eta)χ_rhorho-qχ_Jy,rho=0. The B-derivative also vanishes to machine precision. Spin Zeeman Fermi preparation, resolvent and current sources are independent in the analytic derivative; ∂_B J=0 here solely because Zeeman sigma_x is k independent in the fixed physical orbital frame. NOT true for orbital-B transport.

## Numerical and negative-control evidence

Model N=24, q=2π/24, omega=.27, eta=.045, B=.11, mu=.32,T=.12, beta=.24, alpha=.13, m=.32, w=.35:

- finite-q transverse electric-normalized response: +.006895888683488+.003128753730087i.
- mixed spin Zeeman response d_B χ_{jx,rho}/(-iq): **−.005924308951739−.007974278348438i**.
- Dense per-k 4×4 Sylvester (separate implementation) vs eigenstate spectral response machine-precision agreement; analytic Fermi-Fréchet+Liouvillian B derivative vs full 5-point separately rediagonalized response ~2e-13 absolute difference.
- Longitudinal charge Ward and its Zeeman mixed derivative residuals <1e-17. A smooth momentum-dependent orbital basis changes all k vertices but none of the complete observables.
- Fake intermomentum density vertex I instead of U_+†U_- breaks Ward by ~8e-4; replacing exact finite-q longitudinal current with naive averaged v_y violates Ward by ~1.66e-5. These are REJECTION tests.
- 37 new physical tests plus 1 frozen-data regeneration, local **38/38 PASS**. Separate remote CI and merged-main evidence required before claiming release.

## Scope and open core

VERIFIED finite-q lattice scalar-potential/spin-Zeeman method, Ward, derivative, and gauge covariance. CONDITIONAL physical transverse vertex, finite q and eta, noninteracting toy two-band lattice, charge=1 model units. UNRESOLVED original SM Eq15 electric tau pole, separate eta_a/eta_b and BOTH E/B insertion orders, Eq20/31 sign, main Eq10/SM Eq62 helicity, spin PHE symmetric tensor, full orbital MNHE + Eq49/52 factor-2 + omitted third-band sectors, magnetization vs transport, complete intercell Wannier magnetic shift and experimental SI/Fig2. Source Xiang remains **PARTIAL**.

## v73 research priority

(1) Derive source-specific Bloch two-color EB Liouville current from SM Eqs4-20, not a single scalar φ and frozen B: preserve poles eta_e/eta_b and all orderings, with negative controls for incorrect regulator/order. (2) Show what combination/limit of finite-q χ_jx,rho corresponds to E uniform in the source length gauge, if any, with gauge and magnetization-contact subtractions. (3) Keep literal Eq20 vs Eq31 signs, and main/SM helicities in parallel rather than fit a sign. (4) For orbital MNHE use full magnetic Wannier P, velocity, moment and position shift plus all three-band terms; do not infer them from Pauli spin. (5) Confirm SI normalization and published figures ONLY after source operators and filling matched.

## Files and release

src/qgmap/bulk_finite_q_spin_v72.py; tests/test_bulk_finite_q_spin_v72.py (37); tests/test_bulk_finite_q_spin_v72_frozen.py (1); scripts/reproduce_v72.py; results/frozen/v72_bulk_finite_q_spin.json; docs/BULK_FINITE_Q_SPIN_ZEEMAN_V72.md. Github PR, final SHA, actual Actions CI, merge SHA/postmerge run, and archive SHA256: append **only after independent verification**.

Project planning prior estimate unchanged: GitHub reproducibility ~90–95%, finite-model method integration ~75–80%, **complete original-paper science ~30–40%**. v72 strengthens *bulk periodic method* but not source-identical full-paper completion; percentages are NOT numerical measures of physical correctness.
