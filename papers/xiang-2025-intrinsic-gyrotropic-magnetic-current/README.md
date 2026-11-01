# Xiang et al. (2025): intrinsic gyrotropic magnetic current

**Reference:** Xiang, Jia, Xu, Qiao, Wang, *Intrinsic Gyrotropic Magnetic Current from Zeeman Quantum Geometry*, Physical Review Letters **134**, 116301 (2025), [DOI](https://doi.org/10.1103/PhysRevLett.134.116301). This module is an **independent model calculation**, not endorsed by the authors.

**2026-10-20 status: PARTIAL.** The two Hamiltonian model mechanisms, Rashba Fig.1 AC SI, the bulk of Dirac Fig.2, and the independently checked Eq.(12) displacement factor are reproducible. The **SM Table I–III group-symmetry coverage now includes all 122 magnetic point groups**, using real data from the user's group-theory toolkit. The published Fig.2(d) absolute ordinate and source Fig.2(d) plotting provenance **remain unresolved**; we do not rescale the validated physical response to fit them.

**Latest Round 8:** `REPORT_ROUND8.md` and `docs/ROUND8_SEMICLASSICAL_GREEN_CLOSURE.md` provide direct independent SM semiclassical wavepacket and lesser Green methods, 71 unit tests (70 pass, 1 optional external-toolkit skip), and two additional source-expression findings. The original Fig.2(d) amplitude and Eq.(14) convention remain OPEN. No GitHub upload.

## Run

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/reproduce_all.py
python scripts/verify_release.py
```

Outputs: `results/figures/fig1_rashba_model.png`, `fig2_dirac_model.png`; raw JSON, convergence CSV, frozen reference. `verify_release.py` verifies the baseline frozen model cases, the independent Round-3 Wilson / Liouville checks and the Round-4 122-MPG symmetry snapshot. Previous Round-4 checkpoint: **40 tests discovered, 39 passed, 1 expected skipped** (the optional external-toolkit import). Its PASS is not a full-paper physics claim.

Read the audit chain: [`docs/SOURCE_MAP.md`](docs/SOURCE_MAP.md) → [`THEORY_PROTOCOLS.md`](docs/THEORY_PROTOCOLS.md) → [`FIGURE_ACCEPTANCE.md`](docs/FIGURE_ACCEPTANCE.md) → [`INDEPENDENT_THEORY_AUDIT.md`](docs/INDEPENDENT_THEORY_AUDIT.md) → [`NUMERICAL_CONVERGENCE.md`](docs/NUMERICAL_CONVERGENCE.md). Continue from [`docs/HANDOFF.md`](docs/HANDOFF.md).

## Physics protocols and uncertainty

- `paper-defined-low-frequency` computes the conduction `ΣfZ` and displacement `Ω Σ(2fQ/ε)` channels; unlike static NMEE, the latter multiplies a sinusoidal AC current.
- `independent-full-frequency` computes the causal clean interband Kubo susceptibility for Ω below the direct gap, with no arbitrary η. At small Ω, `K(Ω)≈C+iD`.
- Main Eq. (14) seems to contain a **time-reversal/tensor-character ambiguity** for the conduction channel: treat the printed form and separately implied C/D symmetry behavior as different pending protocols, not a silently corrected formula. See audit.
- Round-1 overview figures plot natural-unit response kernels. Round-2 SI audit separately converts to A/(m·T) and uses the 1 THz cyclic-frequency convention supported by SM §VI. Main Fig. 2(d) remains quantitatively unmatched; do not treat the overview figure as a digitized agreement.
- Pauli spin matrices and projected same-gauge velocity/spin are essential. `wannier90_hr.dat` alone cannot reproduce this response. Existing NMEE code is kept separate and untouched.

## Licensing and provenance

Publisher article/SM PDFs and image scans are not included. Original source hashes and remote NMEE baseline commit are documented in `docs/SOURCE_MAP.md`. A public source license has not been chosen; none is assumed.

## Round 2 update (2026-10-19)

`REPORT_ROUND2.md` and `docs/ROUND2_TIME_REVERSAL_AND_UNITS.md` document independent antiunitary C/D parity checks, SI conversion with `g=2`, explicitly disambiguated 1 THz source frequency, and a **quantitatively unresolved Fig. 2(d) amplitude mismatch**. Model has **23 tests** and adds `scripts/audit_round2.py` plus a symmetry/SI summary plot. `python scripts/verify_release.py` checks both the original six frozen cases and new round-2 checks. This is still a **partial physical reproduction**, not a claimed publisher-figure match.

## Round 3: independent displacement factor verification

Run `python scripts/audit_round3.py` and `python scripts/verify_release.py` after installing the package. Two *independent* checks have been added in `src/igmc_lab/independent_checks.py`: Berry-phase finite-field Wilson loops and orbital-basis Liouville inversion. These support **main Eq.(12)**, while **SM p.6 Eq.(34)** contains an intermediate factor-of-two normalization tension detailed in `docs/ROUND3_INDEPENDENT_DISPLACEMENT_AUDIT.md`. Never use its literal ×2 as a validated physical rescaling to match Fig.2(d). Source graph remains partially unmatched; no GitHub publication claimed.

Public third-party prior-art is documented in `docs/GITHUB_PRIOR_ART_AUDIT.md`; it corroborates the paper's main-text Zeeman geometry definitions but does not provide the original Fig.2(d) code/data. This module contains independently written Python checks, not imported GPL source.

## Round 4: full magnetic point-group classification from the user's toolkit

- `src/igmc_lab/toolkit_bridge.py` uses a pinned minimal operation snapshot from `Xin-Jiaqi/group-theory-operations-toolkit` v0.18.0 at Git SHA `551a0daf...`. It **independently solves** the tensor Neumann nullspaces, preserving original physical C/D parity; it does not copy the toolkit's implementation.
- `python scripts/audit_round4.py` computes all **122** magnetic point groups and **244** C/D tensor sectors, and *independently* checks each Neumann-SVD dimension against a full-group character average. It freezes counts/dimensions in `results/raw/round4_mpg_122.csv` and canonical tensor structures in `results/raw/round4_mpg_tensor_shapes.json`. See `docs/ROUND4_TOOLKIT_SYMMETRY_AUDIT.md`.
- Optional direct API comparison (not run in this offline environment): install the pinned toolkit and rerun `python scripts/audit_round4.py --require-live-toolkit`. Original toolkit's source repo is unchanged.
- Main Eq. (14)'s shared antiunitary sign changes the **number of C degrees of freedom for 58/122 MPG** relative to channel-resolved C/D rules; this remains labelled an author-definition dispute, not a verified full-frequency error.

**Toolkit provenance / license:** the frozen minimal generator data are derived from the user's BSD-3-Clause repository. Its full relevant upstream license is retained at `data/TOOLKIT_BSD3_LICENSE.txt`. This notice is specific to the borrowed generator data and is not an assertion that the entire IGMC reproduction repository has the same license.

### Fourth-round independent character verification

For a 3×3 polar-current response to an axial Zeeman field, the 9-dimensional tensor representation has character `det(R) [tr(R)]²`. The magnetic group average with `eta_C=+1` and `eta_D=-1` produces the allowed tensor-space dimension. Its **244 results agree with all 244 independently computed Neumann/SVD nullspaces**, without reusing the tensor-transformation implementation. All original rounds 1–3 frozen values are preserved. See the v2 checkpoint `results/frozen/round4_mpg_reference_v2.json`.

## Round 5: strict printed-SM Tables I/II audit

`src/igmc_lab/paper_table_source.py` contains an independently and visually transcribed source table from SM pp.8–10 (23 rows; 116 MPG-channel allowed entries). `python scripts/audit_round5.py` compares all 244 source-versus-derived **nine-component tensor spaces**, not only their dimensions. Result: **243/244 exact subspace matches**; the one exception is SM Table II's `-42'm'`, for which the printed shared row omits the additional relation `D_xy = D_yx`. The allowed/forbidden category does not change. See [`docs/ROUND5_SM_TABLE_SOURCE_AUDIT.md`](docs/ROUND5_SM_TABLE_SOURCE_AUDIT.md) and [`REPORT_ROUND5.md`](REPORT_ROUND5.md). Earlier frozen physics and source references are preserved; native full toolkit API remains optional/not locally installed. The paper remains PARTIAL because Fig.2(d)'s absolute magnitude is unresolved.


## Round 6: independent two-band audit (2026-10-20)

- **Eigenvector-free Bloch-sphere proof:** `Z=-∂_k n_beta` and `2Q/Δ=d·(∂_k d × e_beta)/(2|d|³)` for the paper's Zeeman convention. Source: `docs/ROUND6_DVECTOR_FREQUENCY_AUDIT.md`.
- **Full complex frequency and TR:** at 1 and 3 THz, the six complex response components satisfy `K[TR(H)](ω)=K[H](ω)*` in the declared real-phase convention. Real `C` is T-even; imaginary `D` is T-odd. Printed Eq. (14) convention remains open.
- **Independent Gauss–Legendre BZ integral:** order 128² gives `D/ω=0.01523407801563 Å⁻¹ eV⁻¹`, agreeing with old midpoint 301² within `2.01e−11`. Absolute in-gap result remains `0.01775389 A/(m T)` at 1 THz and `g=2`; published Fig.2(d) amplitude is not quantitatively matched.
- **Scientific correspondence:** `docs/AUTHOR_DISCUSSION_BRIEF.md` contains *unsent* questions about Eq. (14), SM Eq. (34) and plotting provenance. No source factor was fitted or altered.

```bash
PYTHONPATH=src python scripts/audit_round6.py --no-write
PYTHONPATH=src python scripts/verify_release.py
```

**Overall status: PARTIAL.** This additional model validation is not a substitute for authors' original plotting data.



## Round 7: metallic sharp Fermi surfaces

```bash
PYTHONPATH=src python scripts/audit_round7.py --no-write
PYTHONPATH=src python scripts/verify_release.py
```

The exact `cos(a*kx)` contour partition, independent Eq.(13) Fermi-line identity, separate 1 THz SI conversion and final acceptance matrix are documented in `docs/ROUND7_METALLIC_FERMI_CONTOUR.md` and `docs/FINAL_ACCEPTANCE_MATRIX.md`. The published Fig.2(d) absolute amplitude is **still unresolved**; no phenomenological ×2 normalization is accepted.

## Round 8: finite-frequency wavepacket and lesser-Green checks

From this directory, `PYTHONPATH=src python scripts/audit_round8.py --no-write` computes the direct orbital Green-resolvent and wavepacket matrix-element responses (no reused quantum-geometry pair kernel). `PYTHONPATH=src python scripts/verify_round8.py` checks SHA256 and frozen values; the global `scripts/verify_release.py` calls this gate. See the standalone Round-8 report for SM Eqs. (34),(46) source inconsistencies and the subleading main Eq. (10) sign question. **Never interpret source-local factors as authorization to rescale Fig. 2(d).**
