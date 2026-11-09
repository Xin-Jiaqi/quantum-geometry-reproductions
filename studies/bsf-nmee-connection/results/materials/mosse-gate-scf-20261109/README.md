# MoSSe S/S: independent self-consistent gate test

Run date: 2026-11-09. Same relaxed zero-strain BiDB 7269 S/S bilayer as the
preceding model study. No ionic relaxation, strain, added carrier doping,
new material, or same-layer reversal target.

## Result: relative gate gain and upper-layer dominance survive screening

All six gate-SCF points completed and remain insulating. At the same **applied**
field interval, the bilayer net sensitivity is **48–49% larger** than the
sum of the two matched isolated layers (64-state references; 96/112-state
bilayer). Including all 48/56/64 reference settings gives 48–52%.

| Gate-self-consistent response | Lower G | Upper G | Net G |
|---|---:|---:|---:|
| Bilayer, 96 states | +0.07849 | −1.44161 | −1.36312 |
| Bilayer, 112 states | +0.05701 | −1.42653 | −1.36952 |
| Exact isolated layers, 64 states each | −0.45790 | −0.46278 | −0.92068 |

G units are `1e-9 (V^2 T)^-1`. The upper layer accounts for 95–96% of the sum
of the **absolute layer changes**, not 95% of the net signed change: the small
lower change opposes the upper one. The bilayer upper layer is about three
times as gate-sensitive as the actual isolated upper layer. Both bilayer
layer responses retain their original signs (lower negative, upper positive).
For positive physical gate, the positive upper response weakens substantially;
the negative lower response changes very little and becomes slightly less
negative. The net therefore moves negative.

Absolute gain must be distinguished from relative advantage. Screening
reduces the bilayer sensitivity from about 15.7 to 1.37, and the isolated pair
from about 11.9 to 0.92. Thus absolute electrical control is approximately an
order of magnitude weaker than the unscreened model, even though the bilayer's
relative advantage survives. The Mo-plane averaged potential contrast is
0.09178 eV per V/nm, versus the bare linear 0.60389: about 15%. This supports
strong screening but is not a unique dielectric constant.

The useful materials-physics observation is: **nearly compensated layer
responses need not have equal gate sensitivities**. At this interface, gate
control is concentrated in the upper layer and the net signal changes more
rapidly than the sum of isolated-layer controls. This gives independent
self-consistent materials support to the candidate interface/gate-control
story. It does not make a 50% enhancement alone a completed PRL novelty claim,
nor uniquely attribute the SCF gain to hybridization instead of differential
screening. The previous channel interventions supply mechanism evidence;
this round verifies survival under gate screening.

No additional fields or zero crossings were computed. Compensation-bias
estimates must be recalibrated with the screened response; do not reuse the
bare-model field scale or extrapolate a confirmed compensation point from
this two-field test.

## Figure and checks

`gate_scf.png/pdf/svg`: (a) Absolute net weak-field sensitivity with and without
gate screening; (b) signed layer and net sensitivities after gate SCF. Open
points and spans retain all tested response-state cutoffs; filled points use
the highest cutoff. They are not statistical errors. Both panels use the same
±0.10 V/nm difference interval, with the model's imposed slope distinguished
from the DFT applied field. The earlier model uses the near-mirrored upper
reference; this SCF round uses its exact geometry. The exact upper-reference
slope differs from the mirrored-lower estimate by only 0.53% of their sum at
the highest cutoff, much smaller than the enhancement.

All actual input hashes and isolated-layer geometries were checked. Minimum
sampled probe gap is 0.79536 eV; the maximum physical response change upon
halving the magnetic probe step is 2.02e-5 in C units. Vacuum slopes confirm
the imposed sign and magnitude to better than 1e-4 eV/Å per V/nm. No carrier
number change occurs. These checks are not a substitute for final mesh and
experimental readout validation. Raw isolated half-cell responses are kept
only in the original result JSON; physical CSV comparisons use their totals.

Both PBS jobs ended normally: 559806 on csuper001 (32m39s) and 559807 on
csuper003 (22m16s). All six SCFs, wavefunction exports and response evaluations
ran on compute nodes. No project job remains active. `scheduler.txt`,
`jobs.json`, `audit.json`, `inputs/`, `points.csv`, `slopes.csv` and
`summary.json` preserve the evidence; `figure_source.csv` preserves plotted
numbers. Figure QA passed; the only advisory is no TIFF submission raster,
with PNG preview and editable PDF/SVG available.

## Protocol and observable

At physical perpendicular fields F = −0.10 and +0.10 V/nm, independently
converge the bilayer and each of its two isolated constituent layers. Both
single layers use the actual bilayer intralayer geometry; the upper reference
is **not** replaced by an approximate mirror of the lower. Initial zero-field
wavefunctions/density are reused for bilayer/lower, but the gate calculation
updates the density self-consistently (ICHARG=1, not a fixed-charge run).
Upper references begin from atomic densities.

The gate sign is conventional physical electric field, corresponding to
electron potential energy +eFz. VASP EFIELD = −F/10 in eV/Å because its sign is
opposite to conventional field; see the [official EFIELD documentation](https://vasp.at/wiki/EFIELD).
[LDIPOL](https://vasp.at/wiki/LDIPOL), IDIPOL=3 and the same vacuum cell are used.
Measured vacuum potential slopes provide a separate sign/magnitude check.

The observable is the existing static spin-only coefficient C in
`delta s_x = C E_x B_z`, using Pauli-spin density per area. The gate-SCF
Hamiltonian and its full KS links/physical spin vertices are re-exported at
each bias. The **perpendicular gate** is self-consistent; the subsequent
Ex/Bz probe response remains the established frozen-potential interband
definition, without new ionic or probe-field exchange-feedback contributions.
This distinction is essential; do not call this a full self-consistent
three-field magnetoelectric calculation.

Gate sensitivity is the finite weak-field secant
`G = [C(+0.10)-C(−0.10)]/0.20`, not a percentage relative to a nearly
compensated C(0). C units: `1e-9 (nm V T)^-1`; G units: `1e-9 (V^2 T)^-1`.
For isolated controls compare the **sum of total responses** of the actual
lower and actual upper layers. Their artificial half-cell partitions in raw
response output are not separate layers and are not used as physical data.

SOC, PBE-D3(0), ENCUT450 eV, full 6×6 mesh, EDIFF1e-9, 112 bilayer/64
single-layer KS bands and fixed atomic positions match the prior quick
materials test. Response cutoffs 96/112 and 48/56/64, plus halved B-probe
increments, test numerical sensitivity. This is not a final k-mesh/device
convergence study. Existing zero-gate results are reused for the bilayer and
lower layer; central slopes do not require an upper zero-gate calculation.

## Sources and reproduction

Inputs derive from `/work/sunq/xjq/mosse-strain-20261108/e000`. New calculations
are separate at `/work/sunq/xjq/mosse-gate-scf-20261109`; source data are
unchanged. `point.py` constructs fixed-geometry inputs, runs VASP, verifies
SCF completion and exports full-KS response data inside PBS. `point.pbs`
accepts `CASE=fm010,FIELD=-0.1` or `CASE=fp010,FIELD=0.1` and refuses login-node
execution. `export.py` follows the prior validated PAW/link implementation,
with the single-layer occupation check extended to the exact upper layer.
`analyze.py` is the unchanged response implementation.

`prior.json` and `PLAN.md` record the pre-result prediction and decision rule.
Each final result JSON includes geometry, input hashes, SCF convergence,
export checks, planar potential/charge profiles, and responses at all state
cutoffs. `zero_reference.json` and `model_reference.json` retain prior context;
`model_slopes.csv` allows the same ±0.10 interval to be compared. No licensed
POTCAR, WAVECAR or credentials belong in the public archive.

Once both final JSONs are collected, run `summarize.py` and `plot_result.py`
locally. Field screening is diagnosed from central potential differences,
including the potential contrast averaged within ±0.25 Å of each Mo plane.
That contrast is an explicit electrostatic diagnostic, not a unique
macroscopic dielectric constant or a complete microscopic mechanism split.

All six gate-SCF calculations and their scientific exports run on PBS compute
nodes, not the Japan frontend. Local work only prepares small inputs and
analyzes returned small JSON/CSV. This round tests whether the **phenomenon**
survives screening; by itself it does not assign the entire self-consistent
gain uniquely to hybridization rather than environment-dependent screening.
