# Stage 3 reproducibility log | 2026-10-20

- Prior stage-2 GitHub commit: `c429336f1a3e4827559c89a75425ced7af69c696`, preserved unchanged.
- First discovered a gapped-sign reversal in the 20-meV-coupled 8-band sliding model; 28-meV coupling is **not** a valid insulating scan through `tau=0` with the frozen mu.
- Sign reversal robust to 41→61 valley-integrator refinement and 0.25→0.125 T field step; root near `u=0.1974`, not symmetry enforced.
- K and K′ each switch sign and contribute equally; no intervalley cancellation.
- Comparison `same-sublattice-only` versus all orbital channels survives BZ-averaged tunneling-Frobenius-norm matching at two sliding points.
- Independent unit tests added. Exact connection-operator and previous stage 2 source files are untouched.
- Physical protocol remains zero-temperature gapped *static* interband finite-B diagnostic, not the disputed finite-T driven DC NMEE.
- OPEN: mechanism attribution beyond a model correlation; real ferroic energetics and material relevance; author's reply on Eq. (9) vs final Eqs. (13)-(14).
