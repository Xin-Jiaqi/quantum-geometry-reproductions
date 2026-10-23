# Round32 | Real-material EPW scattering input contract: strict provenance without counterfeit physics

**Whole Cui source PAPER PARTIAL. Real EPW/Wannier spin+position source data ABSENT. This round validates an *import and eligibility contract*, not actual ab-initio KME.**

## Why the new stage is not another toy scattering theory

Round26 already defined a rigorous *in-memory* same-gauge real-space H/S/A contract, and Round31 showed coherent finite-frequency Ward responses for an **assumed Ohmic** bath. Neither gives material-specific first-principles \(g_{mn\nu}(\mathbf k,\mathbf q)\). Filling that gap starts by preventing misinterpreted electronic gauges, fake spin operators, dimensionless couplings, unconverged k/q meshes and substituted material labels from passing a "real material" gate.

### Verified external references and boundaries

EPW's official theory and interpolation tutorials and the first-principles EPW review describe \(g_{mn\nu}(\mathbf k,\mathbf q)\) as a *complex*, gauge-dependent scattering amplitude, and EPW's interpolation uses electron and phonon Wannier basis transformations. EPW examples contain intermediate \`*.epmatwp\` and \`*.ukk\` files; **neither is a universal exported \`npz\` scatter tensor**, and a printed degenerate-averaged \`|g_sym|\` cannot restore the complex phase/coherent vertex.

- [EPW theory](https://docs.epw-code.org/Theory.html)
- [Official EPW electron–phonon interpolation tutorial](https://docs.epw-code.org/tutorials/tutorial_01/index.html)
- [EPW first-principles Wannier/interpolation review](https://www.nature.com/articles/s41524-023-01107-3)
- [2026 EPW interpolation exposition](https://www.nature.com/articles/s41524-026-02042-9)

**Implementation is explicitly a *normalized external export* reader, NOT a parser for EPW's binary \`epmatwp\`, \`epmatwp1\`, \`epf\` records or any particular QE/EPW version.** A separate verified extraction workflow must specify both original source-version-specific mappings and physical gauge corrections.

## Required manifest and arrays

The loader reads a directory containing **both** \`normalized.npz\` and \`manifest.json\`. It verifies the NPZ SHA256 before loading with \`allow_pickle=False\`.

Required arrays, all same electronic Wannier *declared* gauge:

| Normalized field | Shape / physical meaning |
| --- | --- |
| \`k_frac\`, \`q_frac\` | \`(Nk,2)\`, \`(Nq,2)\`, fractional reciprocal coordinates in [0,1) |
| \`lattice_angstrom\` | \`(2,3)\`, independent real-space primitive vectors, Å |
| \`h_eV\` | \`(Nk,Nb,Nb)\`, Hermitian electronic spinor Hamiltonian, eV |
| \`g_eV\` | \`(Nk,Nq,Nmode,Nb,Nb)\`, complex g, output **m at k+q**, input **n at k**, eV |
| \`phonon_eV\` | \`(Nq,Nmode)\`, phonon energies **ħω**, eV |
| \`phonon_active\` | \`(Nq,Nmode)\`, Boolean mask for validated nonzero channels; inactive modes have zero energy and zero g |
| \`spin_pauli\` | \`(Nk,3,Nb,Nb)\`, **projected physical spin Pauli**, not orbital pseudospin, Hermitian, norm ≤1 |
| \`connection_angstrom\` | \`(Nk,2,Nb,Nb)\`, projected position connection including centers **exactly once**, Å |

Manifest schema \`cui-epw-normalized-scattering-v1\` requires explicit eV/Å/Pauli units, all **H/g/physical spin/position gauge IDs identical**, declared physical spin and position convention, matched \(\pm q\) **phonon branch phase convention**, \`g\` final/input index order and polar-long-range treatment.

\(\mathbf k+\mathbf q\) **must appear exactly** on supplied finite electronic mesh, and \(-\mathbf q\) on the q mesh. For **the explicitly selected phase-paired convention** the importer checks

\[
g_{mn\nu}(\mathbf k,\mathbf q)
=\big[g_{nm\nu}(\mathbf{k+q},-\mathbf q)\big]^*,
\]
with matching partner phonon energies and active-mode masks. **General EPW raw Bloch/phonon gauges need not satisfy this literal identity until phases/branches are re-matched.** Missing or degenerate/mixed phonon branches require a separate defensible mode-space unitary/degenerate-subspace adapter; silently assuming mode index equality is unacceptable.

The validator checks tensor shapes, Hermiticity, Pauli norm, finite values, mesh closure and inverse pairing. It does *not* authenticate whether a declared spin matrix is physically \(\sigma\), whether a tensor came from DFPT or whether operator IDs honestly refer to the same underlying Wannier gauge; these require **independent raw source inspection**. The k-grid position connection also needs the **inhomogeneous derivative term** under nonconstant electronic Bloch gauge transformations; the implemented independent-k gauge-invariant **rate-only** oracle rotates H and g but *never falsely treats an untransformed connection as valid*.

## Source evidence and anti-spoofing semantics

\`stage="synthetic_fixture"\` is always **SYNTHETIC_FIXTURE_ONLY**, even if someone changes \`material_id\` to a real material. Synthetic files may not attach fake evidence. For \`stage="ab_initio_candidate"\`, seven separate **within-directory** sources, each with exact SHA256, are mandatory:

1. Electronic DFT/SCF ground state (\`electronic_scf\`)
2. DFPT phonons (\`phonon_dfpt\`)
3. Original complex electron–phonon data (\`electron_phonon_complex\`)
4. Same-gauge Wannier Hamiltonian (\`wannier_hamiltonian\`)
5. **Physical** Wannier spin (\`wannier_physical_spin\`)
6. Wannier position connection (\`wannier_position\`)
7. Exact export mapping and tool version (\`exporter_mapping\`)

Sources must be readable relative files and must not escape their base folder by a path or symlink; bytes and hashes must match. Even if all seven files exist and hashes match, the only possible result is **CANDIDATE_BYTES_VALIDATED_PHYSICAL_AUTHENTICITY_UNVERIFIED** and \`real_material_verified=False\`. To enforce this scientific distinction, the test suite explicitly hashes **FAKE_SYNTHETIC_EVIDENCE** bytes for all roles: they pass the byte gate but **do not** pass real-physics authenticity. Tampering even one source byte fails.

This prevents source readiness or **paper completion** from being inferred merely from matching filenames, IDs, calculated hashes or a pasted \`material_id\`.

## Limited finite-width scattering diagnostic (not conductivity)

For a validated normalized bundle, diagonalize \(H(\mathbf k)\) and transform the matrix amplitude using independent \(\mathbf k\)-dependent electronic eigenvector unitaries:
\[
g^{\mathrm{band}}_{mn\nu}(\mathbf k,\mathbf q)
=U^\dagger_{\mathbf{k+q}}\,g^{\mathrm{Wannier}}_\nu(\mathbf k,\mathbf q)\,
 U_{\mathbf k}.
\]
For initial \(n\mathbf k\), final \(m,\mathbf{k+q}\), the displayed *per-event* Gaussian-broadened golden-rule rates (s⁻¹) are
\[
\Gamma^{\rm abs/em}_{mn\nu}
=\frac{2\pi}{\hbar}|g_{mn\nu}|^2
\begin{cases}
n_B(\omega)\,\delta_\sigma(\epsilon_{m,k+q}-\epsilon_{n,k}-\omega),\\
(n_B(\omega)+1)\,\delta_\sigma(\epsilon_{m,k+q}-\epsilon_{n,k}+\omega),
\end{cases}
\]
and Pauli-blocked outgoing event rates multiply \(1-f(\epsilon_{m,k+q})\). Units: \(g,\epsilon,\omega,\sigma\) in eV, \(\beta\) in eV⁻¹, \(\hbar=6.582119569\times10^{-16}\) eV·s. These are **not** a normalized BZ integration, conductivity, momentum-relaxation rate, self-energy, screened electron–phonon band basis, or a conserving Ward–Takahashi vertex. A finite Gaussian broadening is **off-shell** and cannot be asserted to satisfy exact energy-resolved detailed balance. Zero-energy acoustic singularities are **masked**, not silently divided by zero; a proper acoustic long-wavelength treatment remains an external dependency.

Positive test: under an independent \(U(\mathbf k)\), transform H and g jointly and re-diagonalize H: the bare and Pauli-blocked rates agree. **Deliberate negative**: transform H alone but not the complex g and the rate changes appreciably. This rate oracle does **not** validate the inhomogeneous position connection and is strictly separate from Round26's full H/S/A physical gauge gate.

## Outcomes and exact next real-material dependency

**New capability accepted only after remote dual-version CI:** an auditable, fail-closed **normalized input readiness stage**, finite-k/q/g and spin/position semantic checks, positive/negative electronic gauge-oracle and finite-width Pauli-blocked rate diagnostic. Existing model work Round17–31 and frozen results are untouched.

**Blocked for real Cui KME:** a real DFPT/EPW calculation's complete original archive and its verified complex g extraction; raw physical-spin and position matrix provenance; lattice and phonon gauge alignment; acoustic and polar long-range treatment; k/q and BZ convergence; actual momentum-conserving collision integral/self-energy/vertex; correct metallic intraband and SI normalization. Passing the automated *candidate* evidence gate can NEVER alone promote the original paper above **PARTIAL**.

Run:
~~~bash
cd papers/cui-2026-curvature-flux-zeeman-qgt
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/audit_epw_material_scattering_contract.py
~~~

The standalone audit produces *synthetic-only* numerical results and a **fake-hashed-candidate non-authentication** negative control. Any honest material stage needs true source files separately supplied and independently reviewed.

## First independent GitHub numerical acceptance and frozen test

The initial [GitHub Actions run 39908180114](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/actions/runs/39908180114) independently passed **326/326 Cui tests** in Python **3.11 and 3.12**, including all prior Round17–31 audits and the new Round32 audit. A separate push-run Python 3.11 log (39908137196) supplied the original machine-measured [Round32 frozen synthetic JSON](results/frozen/round32_epw_material_scattering_contract.json). The final candidate adds **one mandatory frozen numerical test**, raising the target to **327**, to be revalidated in both interpreters on the exact final HEAD and again after merge.

Frozen **synthetic-only**, `Nk=3,Nq=3,Nband=2,Nmode=1`:
- Inverse-q paired complex e-ph g Hermiticity residual **0**, exact k+q closure.
- Broadening `sigma=0.08 eV` produces max *per-event model diagnostic* **4.182276705766e12 s^-1**. Blocking effect max **4.171935511536e12 s^-1**. These are NOT physically converged lifetimes, rates for a material, or KME.
- Independently k-dependent **matched H/g** gauge rotation alters the maximum rate by **0.00244140625 s^-1**, a relative error about **5.84e-16** against the reference max.
- Deliberately **unmatched H-only** rotation produces error **2.715703434731e13 s^-1**. This is the negative control proving that phase-resolved complex g matters.
- The deliberately fake seven-source file bundle is successfully SHA256 checked but remains **CANDIDATE_BYTES_VALIDATED_PHYSICAL_AUTHENTICITY_UNVERIFIED**, `real_material_verified=false`; one tampered source file is rejected. Named real material without evidence likewise remains SYNTHETIC_FIXTURE_ONLY.
- **No authentic real EPW matrix was discovered or ingested.** The frozen file contains explicit `real_material_verified=false`, `physical_spin_provenance_authenticated=false`, `raw_binary_epw_reader_claim=false`, `real_material_kme_si_claim=false`.

The final **327** and merger state require actual exact-commit CI, and are not inferred from the initial 326-result.
