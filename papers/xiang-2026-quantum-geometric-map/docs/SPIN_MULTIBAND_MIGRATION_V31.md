# v31: Four-band spin-MNHE method migration

Migrated from the frozen offline `v30` research checkpoint of Xiang et al., PRB **113**, L201406 (2026). The implementation retains the original independent four-band density-matrix algorithm (SM Eqs. 16–19), with only its basic grid/Fermi/Pauli helpers factored into `_spin_grid.py`.

## Executable protocol
```python
from qgmap.spin_multiband_v5 import four_band_audit
d=four_band_audit(N=23)
print(d["j_rho_terms"], d["total_rho"], d["direct_SM29"])
```

- The four separately computed density-matrix currents add to the independently evaluated SM Eq. (29), with representative (N=23) values (0.000789812871384695) and (0.000789812871385176); discrepancy less than (10^{-12}).
- The four-band model is periodic and deliberately nondegenerate. Bloch U(1) phase changes and Hall index reversal are tested.
- Model units: (e=\hbar=E=1), BZ-mean normalized response to (ar B=-g\mu_B B).
- **Not implemented:** exact PT degenerate subspaces, conserving impurity vertex, driven-DC RTA/thermodynamic protocol matching, or quantitative Wannier material calculations.
- v30 offline has additional tests; GitHub CI certifies this migrated subset alone.

Next migration gates: `orbital_phe_v4`, `orbital_mnhe_v5`, then `microscopic_ward_v11` and SCBA/Kubo paths, each with dependencies and explicit limitations.
