# Stage 3: electronic-only sliding-energy gate and symmetry-protected nonpolar response

**Scope:** additional *research* on the separate BSF × IGMC eight-band model. This is **not** part of Xiang et al.'s PRL reproduction, and **not** a DFT total-energy / experimentally realizable ferroelectric claim.

## 1. What was checked

The earlier stage established a nonzero, insulating reduced displacement IGMC kernel at two *chosen* sliding configurations `u=±0.19`, without imposing a layer-odd onsite potential. However, those states were **not established as energy minima**. Here the missing eligibility gate is examined rather than asserting a sliding-ferroelectric material.

Recompute the sum of the four occupied tight-binding band energies per model unit cell,

\[
E_{\mathrm{band}}(u)= \frac1{N_k^2}\sum_{\mathbf k}\sum_{n=1}^{4}\varepsilon_{n\mathbf k}(u),
\]

for a complete translation period `u/a ∈ [0,1]`, 41 uniformly spaced points (`Δu=0.025`), `13×13` midpoint k mesh, real-space overlap radius `shell=3` and the identical eight-band parameters from Stage 2 (`mass=0.85`). The FM-like branch uses the same frozen model and magnetic order, not the AF parameters. The band sum is **not** the ionic or DFT total energy; it is used only as a diagnostic *within the declared frozen toy model*.

## 2. Result: no electronic-only BSF double well in this path

| Quantity | AF-like (`m_z` partners) | FM-like (`m'_z` partners) |
| --- | ---: | ---: |
| `E_band(u=0)` (eV/model cell) | −3.957643858254 | −3.939304821194 |
| `E_band(u=0.5)` (eV/model cell) | **−4.122542078219** | **−4.087254354277** |
| `E_band(u=1)` (eV/model cell) | −3.957643858261 | −3.939304821200 |
| `E_band(u=0.19) − E_band(u=0.5)` (eV/model cell) | +0.106008199234 | +0.098753437273 |
| Minimum on the **sampled 1D grid** | `u=0.5` | `u=0.5` |
| Monotonic 0→0.5→1, **on sampled grid** | decrease→increase | decrease→increase |

The difference between `u=0` and `u=1` is <`7e-12 eV` with `shell=3`; the smaller `shell=2` approximation yielded ~`2e−5 eV` residual because of finite overlap cutoff. A `21×21` mesh independently gives the same ordering of `u=0`, `u=0.19`, and `u=0.5` for both branches.

**Interpretation:** `u=±0.19` do **not** form local minima of the scanned electronic-only model. The apparent minimum at `u=0.5` and its `u=−0.5` partner are related by **one full lattice translation**, so they are the same stacking class, *not* two nonvolatile opposite-polarization states. This one-dimensional sampled band sum cannot exclude double wells after realistic ion-ion, electrostatic, magnetic and lattice-relaxation contributions. Equally, it cannot demonstrate them.

## 3. Separate, symmetry-relevant observation at translation-identified `u=0.5`

On the `13×13` mesh and `shell=3`:

| Response at `u=0.5` | AF-like | FM-like |
| --- | ---: | ---: |
| Layer population difference `〈τ_z〉` | `≈0` | `≈0` |
| Reduced `D_xx / ω` | `≈0` | **+0.006771784140** |
| Indirect gap (eV) | +0.444149 | +0.300144 |

There is thus a **nonpolar magnetic model configuration with a nonzero displacement IGMC geometric kernel** in the FM-like branch, even though its layer-charge difference vanishes. This is allowed by the `m'_z` antiunitary partner constraint; the corresponding unitary `m_z` branch forces the AF-like `D_xx` to vanish at the translation-identified configuration. This observation distinguishes *magnetic* symmetry control of Zeeman geometry from ferroelectric polarization, and specifically **does not** demonstrate a polarization-controlled switch at `u=0.5`.

The numbers refer to the reduced, dimensionless-BZ averaged *spin-Zeeman* kernel (model energy conventions), **not** an SI-unit material magnetocurrent. FM-like and AF-like Hamiltonians have different frozen exchange structures, so their numeric magnitudes must not be interpreted as the measured effect of changing an experimental magnetic domain.

## 4. Next physics gate: real BSF, not additional arbitrary toy scans

1. Identify genuine **relaxed structural partners** in one candidate bilayer, beginning with AgF₂ only as a literature-guided *hypothesis*. Fix exchange-correlation / vdW / U / SOC / magnetic-order protocols before comparing.
2. Confirm **two local minima plus a finite energy barrier** on the full ion-relaxed sliding landscape. Measure charge and Berry-phase polarization differences in consistent axes; do not infer ferroelectricity from `〈τ_z〉` alone.
3. Resolve the **actual magnetic layer group** of each state and their mapping with the existing `group-theory-operations-toolkit`. Check internal allowed `D_xx` and partner transformation separately. Include SOC magnetic anisotropy / domain consistency.
4. Only for a passing material: export same-gauge spin and position operators along with the Wannier Hamiltonian, freeze signed manifests, and compute gauge-covariant multiband IGMC. `hr.dat` and centers alone are insufficient.

**Go/no-go:** the eight-band *model mechanism* is validated, but this electronic-only energy landscape **does not pass** the BSF-bistability gate. No DFT and no real-material response has been computed.

## 5. Reproduction

```bash
python -m pip install -e .
OPENBLAS_NUM_THREADS=1 python scripts/audit_sliding_energy.py --mesh 13 --shell 3 --intervals 40 --outdir results
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 python scripts/verify_stage2.py
```

Outputs are [`results/stage3_electronic_energy_gate.json`](../results/stage3_electronic_energy_gate.json) and [`results/stage3_energy_scan_one_period.csv`](../results/stage3_energy_scan_one_period.csv). This Stage-3 audit does not mutate the original Stage-2 frozen checksums or any of the Xiang reproduction files.
