# BSF + NMEE: orbital interference physics (stage 4)

**Key result:** at small positive sliding, a negative total static response appears only when the same-sublattice and cross-sublattice hopping channels are combined. Removing either channel produces a nonnegative signal in the studied points. This nonadditivity is a model-specific interference clue, *not yet a universal quantum-geometric rule*.

- Theory, counterfactual controls and limitations: [`docs/STAGE4_PHYSICS.md`](docs/STAGE4_PHYSICS.md)
- Source: [`bsf_nmee/response_mechanism.py`](bsf_nmee/response_mechanism.py)
- Frozen independently rerunnable calculations: [`results/sliding_stage4_mechanism.json`](results/sliding_stage4_mechanism.json)
- Stage 3 sign reversal and original PRL reproduction: **unchanged**.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/run_stage4_mechanism.py --quick --output /tmp/sliding_stage4.json
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
```

This is **zero-temperature, insulating static** NMEE-like susceptibility, *not* the unresolved driven DC coefficient in the original PRL. No BSF structure stability is claimed.
