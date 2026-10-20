# v32 orbital-response migration: printed versus independently derived formulas

The source implementations are from the offline v30 reproduction of Xiang et al., PRB 113 L201406 (2026), with the small periodic Bloch model support extracted into `_orbital_grid.py`.

Independent callable functions:
- `orbital_phe_v4.orbital_phe_eq56`: separate literal first line, literal second line, corrected-sign integration by parts, and finite-B Boltzmann derivative for 2D Bz.
- `orbital_mnhe_v5.orbital_mnhe_2d`: strict 2D squared Berry curvature and finite-field derivative.
- `orbital_mnhe_v5.orbital_mnhe_3d`: two-band 3D positional-shift curl versus literal SM Eq.(52) and Eq.(49)-derived factor-two branch.

All calculations are in normalized e=ħ=τ=1 units; no material-scale SI prefactor is inferred. **SM Eq.(56) internal sign and Eq.(52) factor-two disagreement remain OPEN.** No omitted three-band processes or generic disorder self-energy is implemented here. CI tests preserve the disagreement and check independent finite-field paths, rather than disguising it with fitted factors.

For testing run `python -m unittest discover -s tests -v` in this PRB module after `pip install -e .`.
