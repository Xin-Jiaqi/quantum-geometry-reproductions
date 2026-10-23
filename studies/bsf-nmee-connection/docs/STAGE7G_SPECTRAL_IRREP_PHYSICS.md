# Stage 7g: Hidden virtual-band cancellation behind valley-localized NMEE

**Protocol:** synthetic 8-band sliding Kane–Mele homobilayer; T=0, insulating, equilibrium **static** finite-Zeeman NMEE. Not the paper's unresolved driven-DC coefficient, a stable sliding ferroelectric, or a DFT/material result. All source and response operators follow the existing Bloch embedding. The earlier Stage 1–7f and paper-reproduction source/results remain unchanged.

## Question

Stage 7f found a positive, approximately cubic **scalar** coefficient \(\mathsf B\) from AB/BA-only hopping, competing with the negative, approximately linear **nematic** \(\mathsf C\) from AA/BB–AB/BA interference. Stage 7b also showed that remote occupied/empty bands matter. Is the net *valley-locality* sufficient to justify discarding remote bands? And are both tensor irreps carried by the same virtual transitions?

## Exact spectral-group audit

At each k, and at both \(\pm B_z\), partition the complete eight-band interband static curvature into four **whole band-group** blocks (occupied→empty): LL (middle four only), RL (remote occupied→low empty), LR (low occupied→remote empty), and RR (remote→remote). Each block sums over its 2x2 transitions; as long as these band groups are spectrally isolated, the sum is invariant under mixing inside each degenerate group. There is **no** arbitrary single-band attribution.

For a spectral block \(\mathcal S\),
\[
M_{ij}^{\mathcal S}=\int_{D_K\cup D_{K'}}\frac{d^2k}{(2\pi)^2}
\,\partial_{B_z}\left[-2\mathrm{Im}\sum_{v,c\in\mathcal S}
\frac{\langle v|S_i|c\rangle\langle c|\partial_{k_j}H|v\rangle}
     {(\epsilon_v-\epsilon_c)^2}\right].
\]
From independently recomputed cases \((s,c)=(0,t),(0,0),(t,t),(t,-t)\), take
\(\mathsf B=[M(0,t)-M(0,0)]/t^2\), \(\mathsf C=[M(t,t)-M(t,-t)]/(2t^2)\).
At finite \(t\), these contain higher-order corrections and approach the exact weak-hopping projector derivatives as \(t\to0\).

## Central physical result: a valley-local total hides enormous annulus redistribution

For u=0.4, t=0.15625 meV, \(\pm0.25\) T and 40×96 radial/angle quadrature **per valley** (values of \(B_{\rm scalar}=\mathrm{Tr}\mathsf B/2\)):

| Valley radius R (nm^-1) | LL | RL+LR+RR | Full scalar B |
|--:|--:|--:|--:|
| 0.10 | +0.25489591 | -0.08550753 | +0.16938838 |
| 0.20 | -0.00983792 | +0.18674740 | +0.17690947 |
| 0.35 | -0.34066345 | +0.50998053 | +0.16931708 |
| 0.50 | -0.40306206 | +0.57129689 | +0.16823483 |
| 0.75 | **-0.42038993** | **+0.58843725** | **+0.16804732** |
| 1.00 | -0.42343622 | +0.59146875 | +0.16803252 |

In the annulus 0.10 < |q| < 0.75 nm^-1, \(\Delta B_{\rm LL}=-0.67528584\) and \(\Delta B_{\rm remote}=+0.67394478\), whereas the complete integral changes only \(\Delta B_{\rm full}=-0.00134106\). An informative **cancellation condition number** is
\[
\kappa=\frac{|\Delta B_{\rm LL}|+|\Delta B_{\rm remote}|}
{|\Delta B_{\rm full}|}\approx 1006.
\]
Each internal redistribution is ~503× the net change. Thus the total response can appear converged to <1% although *individual virtual-band sectors are vastly nonconverged and change sign*.

This is a genuine eight-band static Kubo spectral identity, not a claim that a bare four-band model equals LL or that virtual bands can be assigned independent measurable currents.

## Tensor-sector distinction

At u=.4, R=.75 nm^-1, same weak-t setting:

| Complete-band spectral group | B scalar | C nematic \(N_{xx}\) |
|--|--:|--:|
| LL | -0.42038993 | -0.06649126 |
| RL | +0.19857064 | -0.04133402 |
| LR | +0.60936337 | -0.01424803 |
| RR | -0.21949676 | +0.03616709 |
| **Total** | **+0.16804732** | **-0.08590622** |

The **positive cubic scalar** requires remote-band virtual transitions to overcome the negative LL contribution. The **negative linear nematic** has a negative LL contribution that accounts for ~77.4% of its net value, but remote pathways remain substantial and themselves strongly cancel. Distinct irreducible tensor channels have different spectral origins.

This complements, rather than replaces, Stage7e's physical **intralayer velocity/state dressing versus explicit AB/BA interlayer velocity** partition. The spectral-sector partition is a different cut through the same source-response expression.

## Validation and limits

- Decreasing t from .625→.3125→.15625 meV converges the R=.75 scalar B as .16800002→.16803786→.16804732, compared with zero-hopping projector B_scalar≈.16804370. Nematic C converges to -.08590622 vs projector -.08590333; differences are a few parts in 10^5, including finite-Zeeman effects.
- Zeeman ±0.125,.25,.5 T maintains sector and full **signs** at R=.1 and R=.75. Mirror u↔-u flips the grouped scalar B and nematic C consistently.
- Sector ratios are **not symmetry constants**: at fixed R=.75, varying Rashba 16→20→24 meV changes B_LL / B_full from about -7.998→-2.502→+0.117. The 24 meV run has a sampled µ-to-band gap of ~4.67 meV, which warns of increased sensitivity.
- All results are **two finite valley disks, not full BZ bounds**. A global claim about sector capture outside R=1.0 nm^-1 is not proven. The cross/same A/B orbital-channel basis is declared; under arbitrary k-dependent orbital reembeddings one needs connection-covariant observables.

## Consequence for building a minimal valley theory

A proposed reduced model must reproduce **more than the valley-integrated net NMEE**. At a minimum it must:
1. preserve C3/mirror scalar-cubic and nematic-linear selection rules;
2. include the self-energy and source-dressed **spin, velocity, Zeeman and geometric/contact** contributions of remote states;
3. match tensor components and their **momentum-radius dependence**, not only one numerical \(\Gamma_{xz;x}\) or its sign.

Until these checks succeed, the validated eight-band projector response remains the benchmark. The next analytical task is to derive which remote-band spin/velocity excitation loops contribute the compensating scalar response, ideally in a symmetry-adapted resolvent without near-degenerate bare denominators.

**Study-only files:** \`bsf_nmee/spectral_irrep_radius.py\`, \`tests/test_stage7g_spectral_irrep.py\`, \`scripts/run_stage7g_spectral_irrep.py\`, \`results/sliding_stage7g_radial_spectral.json\`. An independent offline archive contains the 40×96 data, field-step and coupling tests, plots, source and checksum manifest.
