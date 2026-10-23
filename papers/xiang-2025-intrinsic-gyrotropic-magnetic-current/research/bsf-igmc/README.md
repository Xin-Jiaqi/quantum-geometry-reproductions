# BSF × IGMC: eight-band microscopic-hopping study

**Independent follow-up study; NOT the Xiang et al. (2025) paper reproduction.**

Scientific status: **MODEL-LEVEL RESULT VERIFIED / AgF₂ B1 NON-SOC AND SOC SCREENING PASSED / MATERIAL IGMC AND FERROELECTRIC SWITCHING NOT DEMONSTRATED**.

The existing reproduction under `../../` remains immutable. The present study asks whether sliding-modulated real-space interlayer hopping, **without inserting any phenomenological odd-in-sliding onsite polarization term**, can generate a nonzero insulating Zeeman-geometry displacement response, with different transformation laws under a unitary and an antiunitary stacking connector.

The answer is **yes for this declared eight-band toy model**. This is a proof of model-level *possibility*, **not** proof of switchable ferroelectricity, of a new material phase, or a DFT prediction. A switchable *pair* of Hamiltonians is not a stable ferroelectric pair.

## What is new versus the previous four-band pilot?

- Hilbert space: **2 layers × 2 spatially displaced orbitals × physical spin** (eight bands; half filling = 4 occupied).
- Two orbitals at in-plane sites `r_A=(0,0)` and `r_B=(0,0.37)` (lattice units).
- Layer displacement `u` enters explicitly through **real-space two-center interlayer hopping** `t_ab(R,u)`, including orbital-diagonal and orbital-offdiagonal hopping.
- The code **rejects** layer-odd onsite potentials, removing the previous pilot's manually prescribed `delta(u) * tau_z`.
- Analytical derivative `∂kx H` is checked independently by numerical difference quotients. Direct interband Kubo geometry and independently diagonalized finite-field non-Abelian Wilson polarization agree as the Wilson mesh is refined.
- Both positive- and negative-u partners exhibit a **nonzero layer charge imbalance** which reverses sign without any imposed layer-odd onsite energy. This is a layer-dipole proxy, **not a modern-polarization or ferroelectric hysteresis calculation**.

## Results (frozen model parameters, not physical units)

At `u = ±0.19`, `31×31` midpoint mesh, `T=0`, `|spin|=1`, `ℏ=e=a=1`:

| Connection branch | Chemical potential (eV) | +u: Dxx/ω | −u: Dxx/ω | Indirect gap at +u (eV) | 〈τz〉 at +u |
| --- | ---: | ---: | ---: | ---: | ---: |
| AF-like, unitary `m_z` | 0.00 | +0.003085332739 | −0.003085332739 | 0.23390069 | −0.11361271 |
| FM-like, antiunitary `m'_z` | −0.03 | +0.019375728336 | +0.019375728336 | 0.09104118 | −0.10350517 |

The **geometric kernel** `Dxx/ω` uses a dimensionless Brillouin-zone average and eV-band parameters; it is **not** `A/(m·T)`. The actual AC current is proportional to drive frequency (and depends on unit conversion). No SI material claim is made.

The 21² → 31² geometric-kernel change is ≈ `2.7e-9` (AF) and `5.5e-9` (FM). Non-Abelian Wilson polarization derivatives at 21² are `0.00304023` (AF) and `0.01940637` (FM), approaching the geometric coefficients with mesh refinement. The Wilson deviation is larger than the direct BZ quadrature deviation and is **not** silently described as exact agreement.

AF-like `u=0`: `Dxx/ω≈0`, `〈τz〉≈0`. In FM-like branch, the antiunitary connector does not require `Dxx/ω` to vanish at zero shift even though `〈τz〉=0`.

**Mechanism control:** For the gapped AF-like model, removing an allowed spin-dependent interlayer hop reduces the response from `0.00308533` to `0.00005396` (≈57×), without closing the global gap. This is a controlled *model* sensitivity, not a quantified material effect.

The occupied-band sum **decreases** over `u=0 → 0.48`; `±0.19` are not demonstrated metastable minima. The study supplies **no total-energy double well, barrier, Born effective charges, magnetic ground state, or DFT structure**. This is the decisive outstanding gate for claiming BSF.

## Reproduce

From this study directory, in a fresh environment (Python ≥3.10):

```bash
python -m pip install -e .
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 python scripts/verify_stage2.py
# Optional full benchmarks; regenerates research results (not paper results):
OPENBLAS_NUM_THREADS=1 python scripts/run_stage2.py --branch AF --output results/new_af.json
OPENBLAS_NUM_THREADS=1 python scripts/run_stage2.py --branch FM --output results/new_fm.json
# Optional figures require matplotlib; frozen CSV/PNG are already included:
python scripts/plot_shift.py
```

Review [physical definitions](docs/MODEL_AND_VALIDATION.md), [scientific judgment](docs/SCIENTIFIC_ASSESSMENT.md), and [DFT gates](docs/NEXT_DFT_GATES.md). Frozen result bundles: `results/stage2_af.json`, `results/stage2_fm.json`, plus `results/shift_scan_17.csv` and two generated plots.

## Group-theory provenance and scientific separation

Magnetic operation matrices are attributed to [`Xin-Jiaqi/group-theory-operations-toolkit`](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit), pinned at v0.18.0 (`551a0daf8d4ba5a5f24b2811694bfe438616de48`). A 2.8kB subset of its generator data and license notice is included under `data/`; **no claim is made that the upstream live API was installed or called**. Internal single-state `m'_x` and partner `m_z` / `m'_z` identities are also validated directly on the computed 8×8 Hamiltonian.

This directory does not modify `src/igmc_lab`, the original PRL-related datasets, baseline CI gates or research history. The first four-band pilot remains an offline antecedent (with a phenomenological `delta(u)`), not rewritten as an eight-band result.

## Stage 3: electronic energy and nonpolar IGMC gate

The Stage-2 `±0.19` responses were **not** established as stable BSF states. [Stage 3](docs/STAGE3_ELECTRONIC_ENERGY_GATE.md) audits the electronic occupied-band sum across a complete translation period: with the frozen eight-band parameters the **sampled 1D energy minimum is at `u=0.5`**, translation-equivalent to `u=-0.5`, whereas `±0.19` are **not local minima** in this band-energy model. A realistic structural total-energy double well remains unverified.

At this translation-identified nonpolar configuration, the AF-like `m_z` branch requires `Dxx/ω≈0`, while the FM-like `m'_z` branch retains **nonzero** `Dxx/ω≈0.00677` despite `〈τ_z〉≈0` (mesh 13², reduced model units). This is a **magnetic-symmetry model observation**, not ferroelectric switching, and not a physical SI response.

```bash
OPENBLAS_NUM_THREADS=1 python scripts/audit_sliding_energy.py --mesh 13 --shell 3 --intervals 40 --outdir results
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
```

## Stage 4: first real-material candidate, AgF₂ (source-grounded preflight)

See [`agf2_preflight/`](agf2_preflight/) for the independent AgF₂ Type-I/strained Type-II magnetic-connector tensor predictions, the published BSF structural evidence and the **fail-closed DFT input gate**. This follow-up does **not** modify or promote the eight-band model's artificial energy surface to a ferroelectric double well. A C2DB-derived AgF₂ B1 candidate has now passed a **non-SOC, low-cost structural/magnetic screen** in PR #38; no SOC-on material response or AgF₂ IGMC amplitude has been computed.


## Stage 5: real-material B1 screen and next gate (2026-10-22)

The single authorized WSL screening run is documented in [PR #38](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38): 11/11 electronic steps converged, VASP reached the screening force criterion at 0.0782 eV/Å, the Type-I Ag moment signs and a preliminary 1.3875 eV gap were retained, and the nearest interlayer F–F distance remained 3.002 Å. This is a **screening candidate**, not a verified local minimum, magnetic ground state or FE double well.

The previously proposed one-point SOC-on static screen in [B1_SOC_STATIC_SCREENING_TASK_20261022.md](agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/B1_SOC_STATIC_SCREENING_TASK_20261022.md) **has now run** (commit fee4e62): 64 SCF iterations converged, preliminary occupation gap 1.3616 eV and Type-I spin pattern preserved. Actual runtime 4.658h **exceeded** the authorized 4h limit and is disclosed in the original receipt. No material IGMC or second FE minimum was computed. The corrected channel-resolved response and literature gates are documented in [the tensor erratum](agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/B1_IGMC_TENSOR_ERRATUM_20261022.md) and [the novelty/story assessment](agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/B1_BSF_IGMC_NOVELTY_STORY_GATE_20261022.md).


## Stage 6: SOC source evidence and response-theory gate (2026-10-23)

The SOC B1 site-moment vectors are consistent with the conditional magnetic operation m'_x at the printed site precision: the antiunitary mirror sends axial time-odd (Mx,My,Mz) to (-Mx,+My,+Mz), so weak common Mz and PAW-projected Ag Lz are **allowed**, not symmetry violations. The complete magnetic density and the actual easy axis remain unverified. The [independent site/magnetic-response review](agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/soc_static_1/B1_SOC_MAGNETIC_AND_RESPONSE_FEASIBILITY_AUDIT_20261023.md) corrects the earlier canting interpretation and records the resource overrun. The [operator feasibility matrix](agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/soc_static_1/B1_RESPONSE_OPERATOR_FEASIBILITY_MATRIX_20261023.md) distinguishes validated toy kernels from missing AgF₂ H/spin/position matrices. The [next task](agf2_preflight/stage2_real_dft/b1_wsl_relax/idea_screen_1/soc_static_1/B1_SOC_RESPONSE_THEORY_READONLY_HARNESS_TASK_20261023.md) is **read-only**, focused on primary-source Eq. (14) convention and operator readiness, with no further VASP authorization.
