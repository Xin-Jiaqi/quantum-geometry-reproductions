# Jia et al. (2026): Nonlinear Magnetoelectric Edelstein Effect

**正文主图：** [图片、数据与重算说明](figures/README.md)

**Reference:** Jinxiong Jia, Longjun Xiang, Zhenhua Qiao, and Jian Wang, *Nonlinear Magnetoelectric Edelstein Effect*, Physical Review Letters **136**, 226302 (2026). DOI: [10.1103/r3x2-kb93](https://doi.org/10.1103/r3x2-kb93).

This folder is the independently tested **paper reproduction module**, not an implementation endorsed by the original authors. It reproduces the Kane–Mele and CuMnAs tight-binding mechanisms, separates *published* and *static* response conventions, and offers a checked input interface for an existing Wannier Stage-1 pipeline. It does **not** claim perfect agreement with published curves or universal metallic driven-DC physics.

## Core access (2026-11-01)

The module now has three explicit, callable paths:

| Input / question | Entry | Supported scope |
|---|---|---|
| Nondegenerate model, published two/three-band formula | `evaluate_clean(..., CleanSpec(branch="printed"))` | Existing printed and full-static branches retained |
| PT-degenerate / general multiband model, mixed AC E/B response | `evaluate_dynamic(model, k, DynamicSpec(...))` | Full finite-frequency SM (6)/(7), both time orderings, finite-T Fermi derivative, physical position/spin operators |
| Gapped T=0 static interband response | `evaluate_insulating_static(...)` | Occupied/unoccupied groups, including degeneracy inside each group |

CuMnAs Eq. (2) Fermi-surface response remains available via `scripts/reproduce.py --case fig1d-printed`; the separate BGK path retains its stated collision/preparation protocol. These capabilities can be used without closing every paper discrepancy.

```python
from nmee_lab import DynamicSpec, evaluate_dynamic, cumnas_dynamic_model
from nmee_lab.ptbackend.cumnas import Model, points

p = Model(a_nm=1.)  # Explicit model length; not measured CuMnAs data.
model = cumnas_dynamic_model(theta_rad=.15, params=p)
spec = DynamicSpec(
    electric_energy_eV=.043+.017j, magnetic_energy_eV=.071+.023j,
    mu_eV=.03, temperature_K=100., observable="A_z", magnetic="x", beta=0,
)
result = evaluate_dynamic(model, points(481)/p.a_nm, spec)
print(result["response"], result["parts"])
```

`Re(z)=hbar*omega` and `Im(z)>0` are supplied in eV; output is at `z_E+z_B`.
The response is per **complex Fourier amplitude** E*B, in `(nm V T)^-1` for
`s/mu_B`. For two cosine drives, each signed pair carries **1/4** of this
coefficient; retain its output frequency when summing the pairs. Imaginary
energies specify adiabatic switching, not a fitted scattering time. `T>0` is
required to resolve Fermi derivatives; the T=0 static path is separate.

For another model, supply `CovariantBloch2D`: H(k), partial H, position
connection A, observable/Zeeman spin and partial Zeeman spin, all in one
orbital gauge. The covariant derivatives include both `[A,H]` and `[A,S]`.
No physical spin or position connection is inferred from HR alone.
Equal-weight BZ integration is chunked to allow dense meshes.

`python examples/use_dynamic.py` reproduces the **specified AC benchmark**:
481² -> 721² changes the complex result by **5.84e-7 relative**. The reported
result is approximately `0.000797145299 + 0.000757237073j`; this is not a
replacement for the published intrinsic DC Fig. 1(c). See
[verification](docs/FINAL_VERIFICATION.md) and the
[source-sign audit](docs/AUTHOR_DISCUSSION_BRIEF.md).

## Source article, Supplemental Material, and study notes

- [Published PRL article (official DOI)](https://doi.org/10.1103/r3x2-kb93)
- [Official Supplemental Material](https://link.aps.org/supplemental/10.1103/r3x2-kb93)
- [PDF provenance, SHA256 and redistribution policy](references/README.md)
- [Researcher's PRL-specific Chinese Markdown notes, preserved verbatim](docs/RESEARCH_NOTES.md)
- [Final verification matrix and the open scientific acceptance gate](docs/FINAL_VERIFICATION.md)
- [Frozen offline reproducibility gates](docs/RELEASE_GATES.md)
- [One technical question for future author discussion (not sent)](docs/AUTHOR_DISCUSSION_BRIEF.md)

Publisher-formatted article/SM PDFs were supplied for **offline research**, but are **not mirrored publicly** without evidence of a license or publisher's permission. The uploaded full research notebook also includes unpublished cross-paper research plans; only its PRL-specific Section 3.3 is released here.

## Reproduce the verified model paths

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/verify_release.py
python scripts/reproduce.py --case fig2-printed-vs-static --mesh 451 --output results/fig2_km_451.json
python scripts/reproduce.py --case fig1b --mesh 101 --temperature 100
python scripts/reproduce.py --case fig1c-static-diagnostic --mesh 81 --temperature 100 --mu 0.03 --theta 0
python scripts/reproduce.py --case fig1d-printed --mesh 151 --temperature 100 --mu 0.03 --theta 0.15 --tau 50 --a-nm 1
python scripts/reproduce.py --case fig1-bgk --mesh 81 --temperature 100 --mu 0.03 --theta 0.15
```

### What the cases mean

| Command | Physics | Validation level |
|---|---|---|
| `fig2-printed-vs-static` | PRL Eq. (1) branch and independent SM Eq. (9) full-static band response, including their difference | **Verified model kernels**, different protocols clearly named; published full curves not strictly matched |
| `fig1b` | CuMnAs Eq. (5) two exactly double-degenerate dispersions | Analytic band structure, no full high-symmetry-path rendering in CLI |
| `fig1c-static-diagnostic` | Sublattice-resolved **frozen-occupation static** finite-B derivative | Staggered-sign check, **not** a completed paper-defined degenerate driven-DC Eq. (1) calculation |
| `fig1d-printed` | PRL Eq. (2), PT-doublet projector spin metric/Fermi-surface relaxation-time expression | Formula/symmetry/scaling verified; real CuMnAs length/absolute amplitude unfixed |
| `fig1-bgk` | **Model extension for verification only:** charge-conserving BGK with specified magnetic-field preparation | Independent numerical checks; **not a claim that PRL uses BGK** |

All results print machine-readable JSON containing a protocol label, units, and known limitations. `results/` contains small examples computed from this repository; they are **not** original paper data.

## Wannier workflow compatibility

The *read-only* preflight accepts the same signed Stage-1 top-K package contract already used by the existing NLO Stage-2 workflow:

```bash
python integration/stage1_handoff.py /path/to/rank-01_candidate --structure L1
```

It checks `package-manifest.json`, four file hashes, quality PASS, SOC, 2D dimensionality, Wannier dimensions and occupation. It deliberately returns `BLOCKED_PHYSICAL_OPERATORS` until a same-gauge Berry/position connection and projected spin operators are supplied. **This is an actual package-interface test, not a demonstration that a real material NMEE value can be calculated from HR alone.** See [method integration](docs/WORKFLOW_INTEGRATION.md).

## Organization

- `nmee_lab/`: calculation code preserved from the previous 32 rounds plus a corrected package-relative import.
- `scripts/reproduce.py`: reproducible Fig. 1/Fig. 2 model entry point with explicit protocol tags.
- `integration/stage1_handoff.py`: non-mutating signed top-K package preflight.
- `tests/`: old 20 tests plus regression tests for paper cases and integration.
- `results/`: example outputs (synthetic/model inputs only).
- `docs/`: physics response definitions, acceptance matrix, provenance, open questions.
- `validation/`: gauge-covariance synthetic benchmark and validation source.

## Historical records

The preceding 1–32 research rounds, including 33 archived ZIPs, their SHA256 index, original comparison charts, and physics reports, are kept in a **separate offline research archive**, not embedded in this public-facing source tree. Retaining full internal research history is different from licensing original publisher figure images or depositing unpublished private data.
