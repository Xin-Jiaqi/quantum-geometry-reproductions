# DBC follow-up: finite-width optical-response design before DFT

Independent theory/model study based on Zou et al., PRL 135, 226901 (2025), DOI **10.1103/sqwy-xhdf**. This is not the authors' implementation. The parent paper reproduction/SHG leading-cancellation module remains separate.

**Scope: MODEL_READY for the explicitly declared 2D insulating three-orbital model and conditional local design protocol. NOT MATERIAL_READY.** No microscopic lifetime, DFT, exciton, generic degenerate-band optical engine, full SHG spectrum or device efficiency claim.

## Run

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/reproduce_core.py --output outputs/check.json
python scripts/reproduce_core.py --extended --output outputs/extended.json
```

The runner recalculates both line-shape optima, three perturbed TB families, a full-BZ competing-feature benchmark, and an independent projector-derivative oracle. `--extended` adds skew scans and denser full-BZ comparisons. It does not just load old results. Use the package via `PYTHONPATH=src`.

## APIs

- `optical_design.critical`, `predict`: local Hamiltonian extraction and line-shape-aware candidate selection.
- `optical_design.integrate_point`: both optical pairs and the complete S10 kernel, energy-aligned finite-BZ quadrature; specified model families only.
- `optical_design.predict_two_feature`: the validated quartic feature plus competing quadratic threshold in the strong-mixing model (Gaussian).
- `optical_design.path.path_stationarity`: coupled-parameter chain rule after optimizing photon energy.
- `optical_design.closure`: mixed normal-form correction, symmetry-breaking skew integral, signed quadratic edge for both shapes, conditional error bounds and a material-input contract.
- `optical_design.projectors.projector_kernel`: independent phase-free finite-difference reference; isolated exactly-degenerate scalar-energy subspaces supported, unequal-energy clusters rejected.

## Read

[Theory and its assumptions](docs/THEORY.md) · [Fresh verification](RESULTS.md) · [DFT handoff](docs/DFT_HANDOFF.md) · [Capability contract](capability.json).

K-space and response units are model-normalized. Gaussian width is sigma; Lorentzian width is HWHM. The small conduction-band denominators are never clipped to generate the reported peaks. Accurate peak quadrature does not certify a complete spectrum or a global optimum.

All new source is independently implemented. Publisher PDFs and private NLO source are excluded. Local audit records distinguish previous frozen results from this run's recomputation. Repository license policies apply; no publisher or other submodule license is silently extended here.
