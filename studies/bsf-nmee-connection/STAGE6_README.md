# Stage 6: symmetry-enforced sliding response orders

**Scientific result:** for a C3-symmetric bilayer reference stacking with a layer-exchanging horizontal mirror, the static in-plane `Gamma_(i,z;j)` response tensor is odd under reversal of sliding. Its rotational **trace** has no linear sliding term and begins at cubic order if analytic; the **traceless nematic** part can begin linearly. This distinguishes two measurable response channels even before specifying microscopic hopping parameters.

See [the complete proof and validation](docs/STAGE6_PHYSICS.md). Calculations for our eight-band *synthetic* bilayer reproduce near-origin cubic-versus-linear scaling and the C3/mirror covariance of independently evaluated Hamiltonians.

## Reproduce

```bash
# From studies/bsf-nmee-connection with existing NMEE installed:
python scripts/run_stage6_symmetry.py --quick --output results/stage6_quick_local.json
python -m unittest discover -s tests -p 'test_stage6_symmetry.py' -v
python scripts/plot_stage6_symmetry.py
```

The tests require the existing NMEE code and use the **actual** `group-theory-operations-toolkit` in GitHub Actions. A toolkit-only test is skipped locally if it is not installed, and must not be counted as a successful external integration test. The original paper reproduction is not modified.

**Important:** this is a **0 K, gapped static model response**, not a resolved driven-DC NMEE or a physical ferroic-material calculation. Orbital-channel even/odd pieces are comparisons of different explicitly chosen Hamiltonians, not gauge-independent experimental observables.
