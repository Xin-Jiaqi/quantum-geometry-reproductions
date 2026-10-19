# Existing Wannier workflow handoff (read-only)

The `integration/stage1_handoff.py` adapter deliberately matches the existing `W90Flow` Stage-1 → NLO Stage-2 export contract, without modifying either private production workflow.

```text
W90Flow quality-ranked top-K candidate
  ├─ package-manifest.json (num_wann, nocc_model, physics.soc, dimensionality)
  └─ L1/ (or another named structure)
      ├─ POSCAR
      ├─ wannier90_hr.dat
      ├─ wannier90-wcenter.dat
      └─ quality.json
           ↓
  NMEE signed-package preflight (SHA256 + NOCC + SOC + dimensionality + quality + HR/centres)
           ↓
  BLOCKED_PHYSICAL_OPERATORS (until operator-gauge inputs supplied)
           ↓
  after external data added: covariant T=0 insulating response interface
```

## Commands

```sh
python integration/stage1_handoff.py /path/to/rank-01_candidate --structure L1 \
  --output /path/to/results/nmee-input-preflight.json
python -m unittest discover -s tests -v
```

The adapter does not copy or alter private input. The model code can be imported as `nmee_lab` after installing the package. It offers:

- `read_wannier90_hr(path, lattice_nm)` for 2D, with lattice in nm and dimensionality explicitly validated.
- `to_clean_model` only with the user's **explicit** zero-connection exploratory opt-in (not for final DFT quantities).
- `covariant_from_wannier` requiring same-gauge `A_beta(k), S_alpha(k), dS_gamma/dk` and `evaluate_insulating_static` (T=0, gapped only).
- CuMnAs PT collision backend only for the validated four-band model, **not** arbitrary degenerate Wannier Hamiltonians.

Production-quality real-material NMEE beyond the scope of this PRL model reproduction requires additional operators and independent physical benchmarks; those are not fabricated here.

This is Stage-1 **integration contract acceptance** and model calculation end-to-end acceptance, not a claim that Stage 2 NMEE has already computed a real material.
