# BSF × NMEE — Stage 5: Weak-coupling channel competition

Core physics: `docs/STAGE5_PHYSICS.md`.

**Main result:** the leading two-channel mixed coefficient `C_sc` stays negative at the three checked slides (`u=0.1,0.2,0.4`). The total static NMEE reverses because positive even-in-cross-hopping contributions overtake the negative odd-in-cross-hopping contribution. The prior stage-4 inclusion-exclusion quantity contains *higher-order* cross-channel terms, so its sign change should not be called a sign change of the *leading interference coefficient*.

Additionally at a **fixed stacking `u=0.2`**, increasing tunneling from 15 to 20 meV reverses the static response without a sampled insulating-gap closure; the refined crossing is `t_perp ≈ 19.63385 meV` (toy-model parameter, not measured).

Important restrictions: this is a **synthetic, gapped 0 K static** response, *not* the disputed PRL driven DC NMEE, not a material prediction, and not evidence of actual ferroelectric switching.

New independent code: `bsf_nmee/channel_perturbation.py`; existing stage-1–4 kernels and frozen tests are left unchanged.

### Reproduce (from this study directory)

```bash
OPENBLAS_NUM_THREADS=1 python scripts/run_stage5_weak_coupling.py --output results/sliding_stage5_weak_coupling.json
python scripts/plot_stage5.py
python scripts/verify_stage5_global_gap.py  # optional, grid + local optimization; SciPy
python -m unittest discover -s tests -v
```

Requires installing the parent `papers/jia-2026-nonlin-magnetoelectric-edelstein/` package, and the pinned `group-theory-operations-toolkit` for the separate group-theory tests.

New results: `sliding_stage5_weak_coupling.json`, `sliding_stage5_fixed_shift_transition.json`, `sliding_stage5_global_gap_screen.json`, and two PNG figures. The fixed-slide root JSON includes separate coarse- and fine-grid bisection checks; the primary weak-coupling script does not recompute that bisection.
