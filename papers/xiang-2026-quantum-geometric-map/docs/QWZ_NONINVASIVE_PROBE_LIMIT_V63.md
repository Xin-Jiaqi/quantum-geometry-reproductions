# v63 | Non-invasive QWZ voltage probes: coupling, width, position and finite-size limits

**Scientific outcome:** VERIFIED for the stated finite coherent QWZ scattering problem and its local spectral/current identities. The extrapolation toward weak voltage-probe coupling is **CONDITIONAL** on this Hamiltonian, geometry, Fermi energy, and fixed retarded regulator. The generic thermodynamic Kubo–NEGF limit and the original Xiang orbital mixed E/B response remain **UNRESOLVED**; the source paper remains **PARTIAL**.

## 1. Protocol: change only the voltage-probe coupling

The finite device and all lead principal-layer blocks continue to use the same QWZ Hamiltonian
[
H(mathbf{k})=sin k_xsigma_x+sin k_ysigma_y+
[m+cos k_x+cos k_y]sigma_z.
]
Drive leads are left/right; top/bottom leads are floating voltage probes. The drive-contact width (w_d) and voltage-probe width (w_p) are independently controlled. All leads have real QWZ transverse blocks and matching hopping matrices, not wide-band self-energies.

Only the device-to-probe coupling matrices are rescaled:
[
C_{T,B}(lambda)=lambda C_{T,B}(1),qquad
Sigma_{T,B}(lambda)=C g^r C^dagger,quad
Gamma_{T,B}(lambda)=i(Sigma-Sigma^dagger).
]
Thus (Sigma,Gammaproptolambda^2); left/right drive leads remain fixed. (lambda=0) is **not** reported as a voltage measurement because the floating-probe equations become singular. Instead, finite positive (lambda) values are reduced and the observable limit is checked numerically.

Unless otherwise specified: (T=0), (m=-1), (mu=0), (L=7), (w_d=w_p=3), (eta=10^{-4}), (V_L-V_R=1), (I_T=I_B=0). The retarded (eta) is an explicit numerical regulator; it is not claimed to be a physical inelastic lifetime.

## 2. Weak-coupling sweep: local invasiveness disappears as (lambda^2)

The local x-cut current is the sum of all oriented +x bond currents crossing that cut. Define
[
Delta J_{m cut}=max_xleft|I_x^{m cut}-I_Light|.
]
The lead current (I_L) is the unit-voltage drive response in (q^2/h); the cut current uses the same normalization.

| (lambda) | (|Gamma_T|_F) | (R_{xy};(h/q^2)) | (I_L;(q^2/h)) | (Delta J_{m cut}) |
|---:|---:|---:|---:|---:|
| 1 | 5.21823437 | 0.999997831 | 0.997555262 | 0.974659367 |
| 0.3 | 0.469641094 | 1.000614457 | 0.997164159 | 0.432670464 |
| 0.1 | 0.052182344 | 1.001835491 | 0.996008323 | 0.071283895 |
| 0.03 | 0.004696411 | 1.002125648 | 0.995729604 | 0.006785592 |
| 0.01 | 0.000521823 | 1.002153896 | 0.995702432 | 0.000757794 |
| 0.003 | 0.000046964 | 1.002157139 | 0.995699311 | 0.000068241 |
| 0.001 | 0.000005218 | 1.002157425 | 0.995699037 | 0.000007583 |
| 0.0001 | 0.000000052 | 1.002157460 | 0.995699003 | 0.000000076 |

The measured top/bottom floating voltages stay finite, while the probe contribution to the local spectral injection/extraction vanishes. At (lambda=1), the two central cuts carry only about (0.0228959) while (I_L=0.9975553). At (lambda=10^{-4}), every x-cut differs from (I_L) by less than (7.6	imes10^{-8}). The floating leads have **zero net current at every scale**, but only the weak-coupling limit suppresses their *local* current exchange.

At fixed L, the resistance converges from the strong-probe value (0.999997831) to the finite-device weak-probe value (1.002157460). Do **not** round the latter to an exact quantum Hall plateau or claim that strong/weak probes give identical voltages. Here the Hall-resistance change is modest while the spatial current rearrangement is dramatic.

## 3. Width and position: two independent geometric controls

Keep drive-contact width (w_d=3) and vary probe width (w_p) independently.

| (w_p) | (Delta J_{m cut}), (lambda=1) | (Delta J_{m cut}), (lambda=10^{-4}) | (R_{xy}), (lambda=10^{-4}) |
|---:|---:|---:|---:|
| 2 | 0.957229156 | (4.33	imes10^{-8}) | 1.002157458 |
| 3 | 0.974659367 | (7.58	imes10^{-8}) | 1.002157460 |
| 5 | 0.996209415 | (1.93	imes10^{-7}) | 1.002157456 |
| 7 | 0.996835097 | (4.32	imes10^{-7}) | 1.002157457 |

At unit coupling, wider floating contacts produce a stronger local redistribution; at (lambda=10^{-4}), that redistribution vanishes. The weak-coupling Hall voltage is nearly width-independent in this particular clean model/geometry, but is not position-independent: shifting the top and bottom contacts in opposite directions changes (R_{xy}) from (1.002157460) to (1.001921570). This is a finite-device spatial-potential effect, not a universal value.

## 4. Explicit regulator and size controls

### Retarded regulator (eta)

At (L=7), (w_d=w_p=3), (lambda=10^{-3}):
- (eta=10^{-3}): (R_{xy}=1.002165695)
- (eta=10^{-4}): (R_{xy}=1.002157425)
- (eta=10^{-5}): (R_{xy}=1.002156599)
- (eta=10^{-6}): (R_{xy}=1.002156517)

The range is about (9.2	imes10^{-6}) in this sweep. This is a bounded regulator sensitivity check for the chosen device, not proof that every observable has a universal (eta	o0) limit.

### Device size (L)

At (lambda=10^{-4}), (w_d=w_p=3), (eta=10^{-4}):
- (L=5): (R_{xy}=1.000983988)
- (L=7): (R_{xy}=1.002157460)
- (L=9): (R_{xy}=1.002728517)
- (L=11): (R_{xy}=1.002957110)

The successive increments decrease over these sizes, but the values remain size-dependent and are **not** converged to a thermodynamic limit. This finite-size sequence is deliberately retained as a caveat rather than extrapolated into an exact result.

## 5. Independent gates and fail-closed controls

The 32 new v63 tests verify:
- (Sigma,Gammaproptolambda^2) while drive-lead couplings/broadenings stay fixed;
- the (lambda=1) and full-width calculations reproduce v62/v61 values;
- cut-current nonuniformity shrinks monotonically with probe coupling;
- finite-(lambda) voltage measurement is distinguished from invalid (lambda=0);
- current conservation, floating-probe zero-net-current, interior continuity, and (A=i(G^r-G^a)=sum_aG^rGamma_aG^a);
- independent drive/probe widths, probe position, (eta), and size controls;
- opposite Chern/chirality, trivial-insulator current threshold, and metallic (R_{xy}
e1/sigma_{xy}) negative control.

Local standalone v60–v63: **111/111 tests passed** (79 prior + 32 v63). Full GitHub regression is separately required; its test count must be read from actual remote logs.

## 6. What remains unresolved

1. The weak-probe extrapolation is a finite coherent QWZ result, not a proof of a general Kubo–NEGF equality or a fully ordered ((lambda,eta,L)) thermodynamic limit.
2. The finite-size Hall resistance is still measurably size- and position-dependent. Larger devices, multiple geometries, and explicit order-of-limits studies remain necessary.
3. This work computes **charge** transport only. It does not supply the original paper's full multiband orbital (B)-perturbation, Wannier position/embedding completion, conserving mixed (E/B) Kubo–Středa response, or Ward/contact/bath terms. Original source audit remains **PARTIAL**.

## Reproduction

```bash
unzip Quantum_Geometry_v63_to_v64_Standalone.zip -d v63-check
cd v63-check
python -m pip install numpy scipy
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python -m unittest discover -s tests -p 'test*.py' -v
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python scripts/reproduce_v63.py --output results/reproduced_v63.json
```

Check `SHA256SUMS.json` and compare `results/reproduced_v63.json` byte-for-byte with `results/frozen/v63_noninvasive_probe_limit.json`.
