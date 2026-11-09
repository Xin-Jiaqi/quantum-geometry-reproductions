# InSe homobilayer NMEE — running research checkpoint

The five baseline records (mono at 0 and +/-0.10 V/nm; A and B at zero gate)
are complete. The two bilayer gate SCFs, electrostatic/channel controls and
an independent denser-mesh check are not yet complete at this checkpoint.
Do not infer an electrical gain or a hybridization mechanism from these data alone.

BiDB 862 is the parent monolayer; 864 is the parallel shifted bilayer A.
B is its m001-related sliding partner with the same lower-layer orientation.
The isolated layer uses the exact constituent geometry and common cell.
All science runs in Japan PBS compute allocations; no frontend computation.

Pre-result hypotheses and protocol: [PLAN.md](PLAN.md). Raw point records:
`mono0`, `monom`, `monop`, `A0`, `B0`. The `summary.json` completion flag remains
false until all seven preplanned points have completed. C is Pauli-spin density
response to Ex Bz, in 1e-9 (nm V T)^-1. Gate F is a separate perpendicular field.
Only the total entry of each mono record denotes the physical monolayer; its
raw lower/upper entries bisect that single layer and are not two constituent layers.

Initial results: isolated-layer C is numerically zero without a gate, and
changes sign under gate reversal. A has approximately (-1053,+1052) lower/upper
responses, with a much smaller net value -1.37. B swaps and reverses the two
layers and has net +1.38. These are 6x6 exploratory values, not converged claims
of giant response. The unexpectedly large cancellation motivated the 8x8 check.

The original bilayer allocation completed A/B; its newly started Am was stopped
to avoid the 2h batch limit. That incomplete attempt is preserved remotely and
its metadata appears in `audit`. Am and Ap are now independent allocations;
interrupted data are never used as a result. See `jobs.json` for provenance.

Prior InSe spin-charge conversion from mirror-breaking Rashba SOC and a Mexican
hat dispersion is already established theoretically:
[Zhou et al., PRB 99,155402 (2019)](https://doi.org/10.1103/PhysRevB.99.155402).
The present question is the neutral insulating EB layer response and whether
interface mixing additionally changes its screened gate sensitivity.
