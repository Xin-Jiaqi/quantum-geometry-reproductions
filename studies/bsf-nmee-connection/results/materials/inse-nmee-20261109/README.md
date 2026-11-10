# InSe homobilayer: interface response and screened electrical control

The seven preplanned gate-SCF states, electrostatic/channel controls, independent
8x8 zero-gate SCF, and Gamma character diagnostic are complete. All allocations
have ended. These are exploratory material results, not final amplitude
convergence or experimental validation.

## Physics result

InSe extends the test from Mo d / chalcogen p materials to an In s,p / Se p
system. The isolated D3h layer has zero zero-gate Cxx. Its parallel shifted
bilayer has opposite layer responses and a much smaller nonzero net response;
the m001-related sliding partner exchanges their magnitudes and reverses the
net sign. The predeclared symmetry expectations are met.

At the initial 6x6 resolution, the net gate sensitivity is about **1.4-1.6
times** that of two geometrically matched isolated layers over tested band
cutoffs, including self-consistent gate screening. However, the two layers
change in opposite directions and partially cancel. The upper supplies only
61-62% of the sum of absolute layer changes, compared with 95-96% in the prior
MoSSe S/S calculation. Similar net gains do not imply similar layer selection.

Positive gate field makes the negative lower coefficient less negative, and
the positive upper coefficient less positive. Both magnitudes decrease, the
upper more strongly, so the net coefficient becomes more negative. Neither
layer nor the net spin response reverses within the tested +/-0.10 V/nm range.
The conditional hypothesis of same-sign layer increments was not observed;
its pre-result statement and outcome remain in [PLAN.md](PLAN.md).

## Structure and definition

BiDB parent **862**, parallel shifted bilayer **864** (A): identity stacking
plus (1/3,2/3), C3v. B is its m001 partner plus a common in-plane translation,
equivalent to (2/3,1/3), with the lower chemical orientation unchanged. Both
use the same database layer geometry. The database calls A sliding-stable;
inversion structure 863 is slightly lower among sampled structures. This is
not a unique-ground-state, phonon, or finite-polarization claim.

Observable: delta s_x,l = C_l Ex Bz, using Pauli-spin density per area.
C units: **1e-9 (nm V T)^-1**. The perpendicular gate F, in V/nm, is a
separate control, with VASP EFIELD=-F/10 eV/Angstrom. G is the central gate
sensitivity [C(+0.1)-C(-0.1)]/0.2, in **1e-9 (V^2 T)^-1**.
Only a mono record's **total** denotes the physical isolated-layer response;
its raw lower/upper entries bisect that single layer, not a bilayer.

Fixed database ions, neutral PBE+SOC slab; gate-SCF followed by a full-KS,
frozen-potential spin-only EB derivative. This is not a self-consistent probe
magnetic response, orbital magnetization, or a relaxed device calculation.

## Initial gate results

Representative bilayer 160 / monolayer 80 spinor bands; all occupied states
retained. Alternative cutoffs are in `points.csv` and `slopes.csv`.

| State | F (V/nm) | C lower | C upper | C net |
|---|---:|---:|---:|---:|
| Matched isolated layer | 0 | — | — | 0.000020 |
| A | -0.10 | -1053.944 | +1052.993 | -0.951 |
| A | 0 | -1053.124 | +1051.755 | -1.369 |
| A | +0.10 | -1052.527 | +1050.731 | -1.796 |
| B | 0 | -1051.752 | +1053.130 | +1.379 |

| Gate sensitivity | G lower | G upper | G net |
|---|---:|---:|---:|
| Matched isolated pair | -1.314 | -1.314 | -2.628 |
| Specified electrostatic reference | -1.191 | -1.190 | -2.382 |
| Full A bilayer | +7.087 | -11.311 | -4.224 |

All initial states remain insulating; the sampled zero-gate SCF gaps are
1.414 eV (mono) and 0.937 eV (bilayer). The large individual bilayer C values
strongly compensate. They motivate the independent mesh check; they must not
be advertised as a giant net signal. Gate-slope convergence is separate from
zero-gate C convergence.

The independent 8x8 SCF retains the layer signs and gives
(-1061.649, +1060.758) for lower/upper C at 160 bands, changes of only
0.81%/0.86% relative to 6x6. The small net C changes from -1.369 to -0.891,
about **35%**. The sampled gap is 0.900 eV. Thus opposite local responses are
stable in this check, whereas their small remainder is quantitatively more
sensitive. The 1.4-1.6 gate-gain range is an initial **6x6 estimate across band
cutoffs**, not a mesh-converged enhancement interval. Finite-gate 8x8 states
have not been calculated. See `mesh_summary.json` for all cutoffs and maps.

The symmetric +/-B derivative is used throughout. Halving its step changes
the large layer C by about 0.014 or less at the representative full-space
cutoff, much less than its magnitude. Nevertheless, finite-space B-even
susceptibility residuals remain (e.g. about -0.017 for the A0 net value in the
internal chi normalization), and the A/B net mirror mismatch reaches 0.010
at 160 bands. These are retained numerical limitations, not claims of final
time-reversal/empty-state precision. See the raw `even_chi` entries.

## What the model controls establish

The electrostatic reference transfers the gate-dependent smooth bilayer-minus-
monolayer planar potential onto a matched isolated layer at the same field.
Its two taper choices and band cutoffs both give nearly equal, negative layer
gate sensitivities and positive-lower/negative-upper static C. This specified
reference cannot explain the full bilayer's signs or unequal electrical
changes. It is not a unique additive partition or a newly self-consistent
material without hybridization.

In the full KS response space, use a selected In s,p / Se p subspace to change
interlayer Hamiltonian blocks. Zero-gate channel projectors are transported
to both screened gate states by PAW overlaps. The original screened potential
and the physical response operators remain; these artificial controls are
not new self-consistent structures and not separately measurable contributions.

| Retained/removed mixing (160 bands, 8 frozen empty states) | C lower at F=0 | G lower | G upper | G net |
|---|---:|---:|---:|---:|
| Full material | -1053.12 | +7.09 | -11.31 | -4.22 |
| Remove selected occupied-empty mixing | -226.01 | -0.99 | -0.99 | -1.98 |
| Retain 50% selected occupied-empty mixing | -807.43 | +0.72 | -3.91 | -3.19 |
| Remove occupied-occupied and empty-empty mixing | -215.94 | +0.79 | -2.42 | -1.63 |
| Remove lower-occupied to upper-empty mixing | -581.86 | -38.96 | +35.88 | -3.08 |
| Remove upper-occupied to lower-empty mixing | -594.75 | +42.26 | -45.39 | -3.13 |

Across all nine band/subspace choices, removing the lower-occupied to
upper-empty block changes the signs of both **gate sensitivities**, while
the static layer-spin signs remain unchanged. Thus changing how a gate tunes
a layer is distinct from reversing that layer's spin response.

Important counterevidence: after removing all selected occupied-empty mixing,
the lower static C remains negative. The residual lower gate-slope sign is
window-sensitive (-1.09 to +0.69), and some choices retain net gain. Removing
occupied-occupied/empty-empty mixing also weakens the gate response substantially.
Therefore the MoSSe statement that a selected occupied-empty channel alone
controls same-layer reversal does **not** transfer unchanged to InSe.
The useful lesson is competing channel control, not a universal one-channel
mechanism or an additive percentage decomposition.

The initial momentum map concentrates near Gamma. At radius 0.4 Gamma-K,
seven of 36 points carry about 93% of the absolute lower-C integrand and 86%
of the absolute net-G integrand. Signed regional sums can exceed the net
because other regions cancel them. These percentages are mesh-dependent.
The 8x8 check retains Gamma concentration: the same radius contains about
92% of the absolute layer-C integrand (13 of 64 points). The Gamma character
diagnostic finds the top valence pair (bands 75/76) to have 74-77% Se pz weight
in the selected orthonormal trial representation; the lowest empty pair
(77/78) has about 53% In s and 45% Se pz. Nearby occupied pairs 69/70 and 73/74,
about 0.63-0.69 eV below the sampled VBM, are predominantly Se px/py.
These support the choice of an s,p orbital contrast to the TMDs, but are
descriptive projected weights, not a unique pz or band-edge NMEE mechanism.
The directional channel intervention includes more states than these pairs.

## Relevance to the PRL candidate

The transfer to InSe supports examining how an interface jointly determines
the signed response of a layer and its electrical tunability. It also rejects
the shortcut "net gain means selecting one layer." It does not yet supply a
predictive rule for which layer a real gate strengthens or weakens, or an
experimentally accessible operation that changes the relevant channel balance.
Those remain physics questions; a larger fitted minimal model is not the goal.

A concrete **next prediction, not an additional calculated result**, follows
from the verified A/B mirror relation. At zero gate,
G_lower(B)=G_upper(A) and G_upper(B)=G_lower(A), while G_net(B)=G_net(A).
Thus, in a small positive gate, A should weaken both opposite layer-spin
responses, whereas B should strengthen both. The net gate-induced change has
the same sign in both domains. B at finite gate has not been independently
calculated in this seven-state set. This is a useful sliding-domain prediction
and a warning that net measurements alone miss the local tuning pattern;
the mirror identity itself is already standard layer-tensor symmetry, not a
new universal law. Substrate/gate asymmetry can invalidate the exact relation.

Prior InSe mirror-breaking Rashba/Mexican-hat spin-charge conversion is already
established theoretically ([PRB 99, 155402](https://doi.org/10.1103/PhysRevB.99.155402)).
It is not evidence that the present neutral insulating EB gate response was
measured. Likewise, generic gate selection and layer-tensor symmetry are
background, not by themselves the new result claimed here.

## Data, figure, and audit

- `mono0/monom/monop/A0/B0/Am/Ap/result.json`: completed independent points.
- `summary.json`, `points.csv`, `slopes.csv`: initial seven-state comparison.
- `reference.json`, `planar.json`, `reference_summary.json`: electrostatic reference.
- `channels.json`, `overlaps.json`, `channel_summary.json`: common-channel controls.
- `A0_mesh8/result.json`, `mesh_summary.json`: separate zero-gate mesh audit.
- `orbital_snapshot.json`, `orbital_summary.json`: Gamma pair characters and representation ranges.
- `inse_gate_comparison.pdf/.svg/.png`, `figure_source.csv`, `figure_qa.md`:
  discussion figure, source data, caption, and numerical-scope caveats.
- `jobs.json`: scheduler provenance, including the intentionally interrupted
  original Am attempt, whose incomplete data are not used.

All VASP, wavefunction processing and remote numerical postprocessing execute
inside Japan PBS compute allocations. The local small-data summaries and
figure are reproducible with the adjacent Python scripts. Raw licensed PAW
data, WAVECARs and large wavefunction caches remain on the compute storage;
they are not committed to Git. Python export imports the existing PAW/wave
helpers from the prior MoSSe workflow, so this archive is a research provenance
bundle, not an independently packaged end-user installation.

Follow-up: [independent mirror-domain gate validation](../inse-mirror-gate-20261109/README.md) confirms local tuning directions; net gain is not quantitatively converged.
