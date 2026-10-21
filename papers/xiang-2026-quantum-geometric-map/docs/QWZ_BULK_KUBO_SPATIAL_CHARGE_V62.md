# v62 | Independent periodic charge Kubo and spatial four-terminal QWZ current

**Scientific acceptance:** VERIFIED named periodic interband Kubo and analytic Pauli-band charge Hall for the *same* QWZ Hamiltonian, and finite-QWZ device Caroli + local spectral/bond-current identities. Finite Hall voltage, metal and contact-dependent spatial transport are **CONDITIONAL**. A generic, correctly matched bulk conductivity-tensor → mesoscopic resistance limit and original intrinsic multiband orbital mixed E/B response remain **UNRESOLVED**. Xiang source paper remains **PARTIAL**.

## 1. Do not confuse observables or protocols

The microscopic Hamiltonian is shared: H(k)=sin(kx)σx+sin(ky)σy+(m+cos(kx)+cos(ky))σz, with real-space Tx=(σz−iσx)/2, Ty=(σz−iσy)/2.

* **Periodic intrinsic charge σxy** in q²/h: independently (a) analytic Pauli d·(∂xd×∂yd)/(2|d|³) weighted by f(ε−)−f(ε+) and (b) nondegenerate interband Kubo velocity matrix elements with Fermi occupation, integrated over identical Brillouin-zone grids. Explicit μ and T, where T is energy (k_B=1). This is the *intrinsic band contribution* only, not σxx or an experimental rho_xy at finite disorder.
* **Finite four-terminal Rxy** in h/q²: zero-temperature Fermi-surface Caroli–Landauer observable with VL−VR=1, IT=IB=0, and retarded semi-infinite QWZ lead η=10^-4. The leads have their own *real transverse strip width* and physical QWZ bond matrices. Only full-side-width electrodes continue all boundary bonds exactly; narrower contacts introduce explicit physical edges and scattering (not identical infinite 2D bulk).
* **Local current**: rho_bias = Σ_a V_a Gr Γ_a Ga (without common 1/2π); J(i→j) = −2 Im Tr[H_ij (rho_bias)_ji] in q²/h per unit voltage. All bonds are oriented +x/+y and evaluated from a commutator, not by fitting an edge-current amplitude. The local spectral A=i(Gr−Ga) must equal Σ_a Gr Γ_a Ga. Current continuity applies to interior sites without a lead; open-boundary contacts can exchange current locally even when the *net* floating-probe current is zero.

Do **not** equate finite-T bulk Kubo and T=0 device; do not use Rxy≈1/σxy outside the clean chiral plateau (including metallic filling).

## 2. Frozen evidence and meaningful negative controls

All energies/hoppings and η in dimensionless tight-binding units. Periodic σxy uses 31×31 k sampling.

| m | μ | T | periodic σxy [q²/h] | Notes |
| --- | --- | --- | ---: | --- |
| -1 | 0 | 0 | +0.9999999981 | clean topological charge Kubo |
| +1 | 0 | 0 | -1.0000000019 | reversal |
| +3 | 0 | 0 | ~0 | trivial insulator |
| -1 | 0 | 0.2 | +0.9892200797 | explicitly thermal, NOT device T |
| -1 | 0 | 1.0 | +0.4619335085 | thermal suppression |
| -1 | 1.5 | 0 | **-0.0774126575** | metallic filling: NOT quantized |

Periodic analytic Pauli and independent interband velocity/Kubo oracles agree in each case (numerical precision), including in the metal.

**Spatial probe/contact control:** at L=7, m=-1, μ=0, width=3, Rxy=+0.9999978305, I_L=0.9975552616 (unit voltage) while four central x-cut bond-current sums are **not all equal**. Cuts from left to right: 0.9975553, 0.9975553, **0.0228959**, **0.0228959**, 0.9975553, 0.9975553. This is *not* continuity failure; floating top/bottom probes redistribute local charge injection/extraction while keeping net IT=IB=0. Exact spectral identity residual <2×10^-14; bulk-site current-divergence residual <2×10^-16. Edge local density of states (x=y=0) 0.360964 vs center 0.000679, at that specific chosen contact geometry. Moving top/bottom strip contacts changes interior current-profile while Rxy remains near the plateau.

**Width effect:** L=7, width=5: Rxy=+0.9999999993, I_L=0.9999879091; middle cut ~0.0018273, further demonstrating voltage-probe invasiveness and finite geometry. Full-width L=5 reproduces v61 transmission and Rxy to numerical precision.

**Opposite magnetic chirality:** L=7, width=3, m=+1 yields Rxy=-0.9999978305, I_L=0.9975552616.

**Fail-closed trivial state:** L=7, width=3, m=+3, I_L≈3.54×10^-13, return Rxy=None rather than divide by tiny current.

**Metallic negative control:** L=5, width=3, m=-1, μ=1.5: Rxy=+0.01001052, I_L=1.5993993, while intrinsic bulk σxy=-0.07741266. Naively computing 1/σxy would predict ≈−12.92 (h/q²), plainly wrong. **No general equality is asserted; metallic σxx and lead/contact/bath physics have not been matched.**

## 3. Scientific tests and validation boundaries

30 **new** v62 tests, plus previously standalone v60 21 and v61 28 = 79 standalone independent tests. Full GitHub qgmap regression must additionally pass on the PR head (CI is not presumed from local results).

Independent physical controls include real finite-width QWZ principal layers; exact L=5 full-contact v61 limit; independent analytic Pauli vs velocity Kubo across μ/T/m; local bond-current commutator and interior continuity; spectral A vs Σ injected spectral densities; local orbital gauge transformation **including all contact couplings**; offset/probe-width change; sign and trivial-insulator controls; deliberately invalid inverse σxy in metallic regime.

## 4. What remains open for v63

1. Coherent spatial bond flux is not yet an independently implemented full Kubo two-terminal/four-terminal **same-observable** conductivity in the thermodynamic limit; no disorder scattering, inelastic/dephasing, physical reservoir model, or controlled contact order-of-limits proof. Avoid claiming this correspondence.
2. Repeat with explicitly **noninvasive** (weak and narrow) probes; compare probe-strength→0, width→0, L→∞, η→0 in controlled ordering. Distinguish probe-induced current rearrangement from bulk response.
3. Higher scientific priority: full position-operator-embedded Wannier orbital B perturbation and conserving Kubo–Středa mixed E/B with contact / Ward identity. Original Xiang source audit remains **PARTIAL**, not rehabilitated by charge-topology controls.

## Reproduce

From the paper module: `python -m pip install -e .`, `python -m unittest discover -s tests -v`, then `python scripts/reproduce_v62.py`. Frozen results: `results/frozen/v62_qwz_bulk_spatial_charge.json`. Use the independent standalone ZIP for the 79-test v60–v62 subset and hash verification.
