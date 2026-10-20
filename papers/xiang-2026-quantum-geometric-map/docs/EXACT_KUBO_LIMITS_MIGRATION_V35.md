# v35: Exact Kubo continuity, contact term and SCBA order-of-limits migration

Source: the original v20 scientific-method subset from the offline Xiang PRB 2026 reproduction, **not** an external source of materials data. The two files `real_frequency_continuity_v20.py` and `limit_order_v20.py` and the original 18 method tests are now present alongside v17–v19 in the public PRB module.

## Two independent numerical protocols

**(1) Exact finite-system real-frequency Kubo** (`real_frequency_continuity_v20`):
A finite (L\times L) periodic two-orbital lattice is exactly Fourier-transformed and combined with onsite scalar disorder. Full Lehmann eigenstate sums yield density–current response and the equal-time commutator. By defining the longitudinal lattice current (J_q=-[H,\rho_q]/Q(q)), with (Q(q)=2\sin(q/2)), the matrix relation

[
(\omega+i\delta)\chi_{\rho B}-Q(q)\chi_{JB}
=\langle[\rho_q,B]\rangle
]

holds for density and current choices of (B). The contact term for the current source is large and essential; missing it fails dramatically.

**(2) SCBA regulator-order scan** (`limit_order_v20`):
Keep (q=0.06\) and (omega=0.04) fixed and decrease external (eta) in the v19 finite-((q,\omega)) resolvent identity. The regulator leakage decreases, but at (q=omega=0) the density vertex grows so that the leakage need not vanish. This is a **noncrossing approximate resolvent** study, independent of the exact finite-system Kubo calculation.

## Reproducible usage

Run from the PRB module after `pip install -e .`:

```bash
python scripts/reproduce_advanced.py --case exact-kubo-v20
python scripts/reproduce_advanced.py --case limit-order-v20
python scripts/reproduce_advanced.py --case exact-kubo-v20 --mesh 4
python -m unittest discover -s tests -v
```

`--mesh` denotes finite system size (L) for `exact-kubo-v20` (3–12); for `limit-order-v20` it is the momentum-grid size (n). Output includes a protocol label and PARTIAL science status. No routine infers real-material SI PHE coefficients.

## Frozen independent sanity numbers

For the default (L=5), (omega=0.27), (delta=0.035):
- Max longitudinal density/current Ward errors approximately (10^{-15}).
- Missing-contact current-residual magnitude approximately **19.6950466983** in the (y) direction.
- Uniform dynamical density response approximately zero.
- Grand-canonical thermodynamic compressibility **3.7449214344**.

For (n=41) at (eta=.02,.005,.001), the fixed-nonzero-((q,omega)) regulator-leak magnitudes are approximately (0.05805908,0.01719904,0.00378760). The (q=omega=0) leakage at (eta=.001) remains approximately (0.03562275).

These are model numbers and are **not** measurements or direct Fig. 2 PHE reproductions.

## Scientific limitations

- The exact benchmark enforces longitudinal lattice charge continuity by operator algebra; it is not a calculation of transverse mixed (EB) magnetotransport.
- The parameter (delta) in the retarded Lehmann resolvent and the parameter (eta) in SCBA should not be treated as a proven relaxation time (	au).
- A dynamic uniform density source and a quasi-equilibrium chemical-potential shift are different protocols; their zero-frequency limits cannot be interchanged without proof.
- Full electromagnetic gauge response, physical Kubo–Středa/Bastin current, contact vertices at finite field, and published nonlinear (EB) steady-state PHE remain **OPEN**.
