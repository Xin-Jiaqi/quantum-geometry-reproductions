# Core capability verification (2026-11-01)

**Core access is usable within the named protocols. Whole-paper status remains
PARTIAL.** This is sufficient to use the verified calculations in later work;
removing PARTIAL is not the development target.

| Capability | Evidence | Boundary |
|---|---|---|
| Nondegenerate printed/static band response, including two/three-band terms | Existing Kane–Mele 451² frozen data; independent static/gauge tests | Printed and full-static branches remain distinct |
| CuMnAs bands and Eq. (2) Fermi-surface response | Existing five-case release verification; tau/length/Neel scaling | Model length and absolute experimental scale require inputs |
| Full finite-frequency mixed E/B response | New `evaluate_dynamic`; orbital-superoperator + five-point k derivative, independent CuMnAs doublet projectors | Adiabatic response, not isolated intrinsic DC |
| General multiband / exact degeneracy | Six-band, three-doublet test with thermal occupation derivatives | T>0; finite positive switching rates |
| Same-gauge position and spin operators | k-dependent unitary test with nonzero A and dS; negative control omitting dS | Physical matrices must be supplied |
| Dense BZ quadrature | Chunking invariant; 241²/481²/721² CuMnAs AC example, last change 5.8311e-7 | Convergence established for stated parameters only |
| Source diagonal-term contraction | Direct SM p3 Eq. (9) -> Eq. (11) algebra and independent 2/3/5-band matrix contractions | Confirms a local sign inconsistency, not a corrected full paper |

The dynamic kernel solves, with `R_z=(z-[H,.])^-1` and `D=partial-i[A,.]`:

```
rhoE = i R_E D f(H)
rhoB = R_B [S,f(H)]
D rhoB = R_B ([D H,rhoB] + [D S,f(H)] + [S,D f(H)])
rhoEB = R_(E+B) (i D rhoB + [S,rhoE])
```

This is the full mixed equation in SM (6)/(7) in a common orbital gauge. It
retains both time orderings and exact degenerate blocks. The independent
oracle solves the N²-dimensional orbital linear equation and differentiates
the actual first-order response at neighboring k points; it does not reuse
the production spectral resolvent or the analytic derivative above.

## Verification commands

From the paper directory after `pip install -e .`:

```
python -m unittest discover -s tests -v
python scripts/verify_release.py
python examples/use_dynamic.py
python validation/check_diagonal_source.py
```

38 tests (30 prior + 8 new). The release manifest tracks changed/new source
files explicitly; all five prior frozen response files remain unchanged.
New outputs are `results/dynamic_core.json` and
`results/diagonal_source_audit.json`.

## Remaining source questions, separate from usable capabilities

CuMnAs Fig. 1(c)'s finite-T intrinsic DC term is still not equated with a
static derivative, a BGK steady state, or a Laurent constant of arbitrary
adiabatic rates. The printed diagonal sign conflict is now a definite local
algebraic finding; deciding a corrected total response still requires all
terms and the selected limit. See `AUTHOR_DISCUSSION_BRIEF.md`.

No new real material, complete microscopic scattering theory or generic
Kubo–NEGF theorem is required to use the capabilities in the table. True
material predictions still need physical same-gauge operator input.
