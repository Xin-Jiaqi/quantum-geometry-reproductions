# MoSSe: screened common-channel test

2026-11-09. Reuses existing zero/±0.10 V/nm gate-SCF data; no new DFT or material.
PBS 559934, csuper001, 1 min 52 s, Exit_status=0. All numerical and wavefunction
processing ran on a compute node. Local work only prepared scripts and summarized
the returned small result files.

## Result and interpretation

The selected cross-layer valence–empty channel supports BOTH the negative
lower-layer zero-gate response and the strong upper-layer gate sensitivity.
Removing it at the unchanged original screened potentials restores positive
lower-layer C and greatly reduces upper selectivity. The two responses are
effects of a common electronic channel; local sign reversal does not cause, or
form a necessary threshold for, gate selectivity.

Representative 112-state, frozen-8-empty definition:

| Retained VC amplitude alpha | Zero-gate lower C | Lower G | Upper G | Net G |
|---|---:|---:|---:|---:|
| 0 | +12.383 | −0.410 | −0.399 | −0.809 |
| 0.5 | +6.110 | −0.202 | −0.785 | −0.987 |
| 0.75 | +0.702 | −0.082 | −1.078 | −1.160 |
| 1, physical bilayer | −6.193 | +0.057 | −1.427 | −1.370 |

C is Pauli-spin density in delta s_x=C E_x B_z, units
1e-9 (nm V T)^-1. G=[C(+0.10)-C(-0.10)]/0.20, units 1e-9 (V^2 T)^-1.
The perpendicular gate is distinct from the Ex/Bz probe.

Across all six 96/112-state and 4/8/12-empty projection settings:

- Physical lower C = −6.92…−6.19, no-VC lower C = +9.08…+13.92.
- Physical net |G| = 1.363…1.370; no-VC |G| = 0.734…0.834.
  Existing matched isolated-pair |G| is 0.921. Thus the selected removal removes
  the relative gain in this conditional comparison; this is not an additive
  percentage of the physical response.
- Upper share of absolute layer changes falls from 95–96% to 43–64%.
- At alpha=.75 lower C remains positive in ALL settings, yet upper share is
  already 84–95%. Strong selectivity does not require same-layer reversal.
- Removing lower-valence/upper-empty mixing alone restores positive lower C
  (+2.86…+5.00) and changes gate dominance to the lower layer: G_lower
  −1.32…−1.20, G_upper +0.089…+0.207. This identifies a particularly useful
  channel to investigate, not a unique additive microscopic contribution.
- Removing VV+CC leaves lower C negative and upper gate dominance intact
  (87–97%). Removing the reverse VC direction has projection-dependent local
  signs; do not present a unique sign rule for that control.

## Protocol and checks

Channels are defined in a 44-state d+p subspace at zero gate, embedded in the
full KS response. The valence projectors of its separated layer blocks are
fixed once at zero gate. At each field, the same physical selected orbitals
are represented by <u_F|u_0> PAW overlaps and polar orthonormalization. Scaling
the projected Hamiltonian's selected VC channel at EVERY field consistently
scales its field dependence too. All unselected states and physical spin,
Zeeman and Bloch-link vertices remain present. No band is fitted or discarded
from the full response to manufacture a favorable outcome.

Same-k overlaps include pseudo-wave and PAW augmentation terms. Equal atoms,
cell, potentials and reciprocal wave sets are asserted. Maximum metric error
2.69e-6; transported selected-space singular-value deviation <1.50e-4.
Synthetic basis-covariance/orientation and Hermiticity checks passed locally.
The alpha=1 results recover the independent gate-SCF outputs exactly at both
cutoffs. All 126 points remain insulating, smallest sampled gap 0.795 eV.
Halving the magnetic probe step changes C by at most 3.90e-5 in reported units.

The inherited 6x6 mesh and finite KS cutoffs are not final convergence. Raw
B-even spin-electric terms are retained, including a larger residual in the
112-state representation (up to 0.073 in the internal normalized chi units,
already present in the full baseline; 96-state values much smaller). The
reported coefficient uses the symmetric B-odd derivative, not chi(B)/B.
Cross-cutoff agreement of the mechanism does not prove all time-reversal or
empty-band artifacts are resolved. Final quantitative claims require that audit.

This is a frozen-electrostatics Hamiltonian intervention on independently
screened reference states. The hypothetical modified Hamiltonians are NOT
rescreened self-consistently. No-VC is not an isolated layer or a unique
electrostatic-only material. It supports a common-channel mechanism conditional
on the specified decomposition; it does not prove a universal causal law,
complete experimental feasibility, or PRL novelty by itself.

## Material exploration consequence

Do not demand a perfect MoSSe reduction before exploring. Use a small set chosen
to change the physics: stronger-SOC non-Janus homobilayer; a heterobilayer with
different layer-energy alignment; later a different lattice/symmetry. Selection
must use orbital compatibility, momentum and energy alignment, not small gap
alone. Do not presume that type-II alignment implies strong valence–empty
hybridization or a particular NMEE sign.

MoSe2/WS2 is an illustrative candidate because gate-controlled type-I/type-II
alignment has been observed [Nat. Commun. 15, 4075 (2024)](https://doi.org/10.1038/s41467-024-48321-1).
That is optical/electronic feasibility evidence, NOT an NMEE result. Lattice
mismatch/moire versus a strained commensurate approximation needs evaluation
before selecting a concrete calculation. MoSe2/WSe2 is another alignment
comparison candidate. No new materials were submitted in this round.

## Reproduction

`PLAN.md` precedes results. `export_overlaps.py` and `screened_channels.py` run
through `run.pbs` using the existing remote source caches and licensed wavefiles.
`overlaps.json` identifies source-cache hashes. Large wave data, POTCAR and
cross-field matrices remain private on the compute filesystem; no credentials
are included. `result.json` retains every point/control/window, `summary.csv`
the zero-field values and gate slopes, and `audit.json` the checks.
`summarize.py` uses the sibling archived gate-SCF results for exact baseline
comparison. Existing `analyze.py` and `mixing_channels.py` are unchanged helpers.
