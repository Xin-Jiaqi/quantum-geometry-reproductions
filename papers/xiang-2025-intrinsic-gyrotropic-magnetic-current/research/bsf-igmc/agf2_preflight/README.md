# AgF₂: symmetry-first DFT preflight for BSF × IGMC

**Research extension, not a reproduction of Xiang et al. (2025).** This directory is independent of the frozen IGMC paper model code. Current status: **SYMMETRY PASS / LITERATURE BSF SUPPORT / DFT INPUTS BLOCKED / MATERIAL IGMC UNCOMPUTED**.

This first real-material candidate follows the **AgF₂ bilayer** reported by Zhu et al., *PRL* **135**, 056801 (2025), DOI [10.1103/dmzg-ck2t](https://doi.org/10.1103/dmzg-ck2t), rather than interpreting the preceding eight-band toy Hamiltonian as AgF₂.

## Main result: crossed switching parities

Use experimental/lab axes `x=[100], y=[010], z=[001]`; let the polar AC current obey `J_i = chi_ij B_j` with an axial AC magnetic field. Within the *explicitly declared real-phase* spin-Zeeman IGMC C/D protocol, the unitary / antiunitary connector acts as

\[
\chi_X^{(B_2)}=\eta_X^\theta\det(R)R\chi_X^{(B_1)}R^T,\quad\eta_C=+1,\quad\eta_D=-1.
\]

For the AgF₂ **Type-I** SOC ground-state magnetic point group `m'_x`, single-state allowed components are

\[
C=\begin{pmatrix}0&c_1&c_2\\c_3&0&0\\c_4&0&0\end{pmatrix},\qquad
D=\begin{pmatrix}d_1&0&0\\0&d_2&d_3\\0&d_4&d_5\end{pmatrix}.
\]

The two structural FE partners are connected by `m_z`. The magnetic connector is `m_z` for Type-I (zero strain), and `m'_z` for **Type-II** (literature: 4% biaxial tensile strain). Source-defined connector relations imply:

| Observable across B₁→B₂ | Type-I: `m_z` | Type-II: `m'_z` |
|---|---|---|
| Structural `P_z` | − | − |
| `C_xy` (if metallic) | − | − |
| Insulating-eligible `D_xx` | **−** | **+** |
| AHE `sigma_xy` (if nonzero) | **+** | **−** |

**Candidate physics:** an orthogonal *AHE/Displacement-IGMC parity fingerprint* of the two interlayer magnetic orders. This is a **selection-rule prediction**, not a calculated current magnitude or a demonstrated measurement; a symmetry-allowed component can vanish accidentally. Type-II must also be checked against its actual magnetic symmetry after relaxation, with correct local Cartesian setting.

## Source-grounded material qualification

Zhu et al. report distinct stable B₁/B₂ AgF₂ bilayer stacking translations `(0.5,−0.16)` and `(0.5,+0.16)` (fractional in-plane), a polarization pair of `±1.66 pC/m`, and a barrier of `17.1 meV/f.u.`. Their DFT establishes an AgF₂ bilayer BSF energy barrier **in that paper's computational protocol**. This does **not** prove that our previous synthetic eight-band Hamiltonian has a double well, nor is it independent verification here.

They report Type-I below Type-II by `0.26 meV/f.u.` at zero strain and a strained Type-II branch around 4% biaxial strain. A sub-meV magnetic-energy difference demands an independently converged and fully specified first-principles comparison before promoting a given magnetic ordering to a robust device proposal.

## Run the source-constrained preflight

From this directory (Python ≥3.10; standard library only):

```bash
python -m unittest discover -s tests -v
python agf2_preflight.py --out results/agf2_preflight.json
```

The **DFT input gate intentionally returns `BLOCKED_INPUTS`** when using the example contract. The original paper provides stacking translations and polarization but these are **not complete POSCAR atomic coordinates**, and we have not located/validated the author's released pair of relaxed structures or their frozen functional/U/POTCAR protocol. No guessed POSCAR, Hubbard U or pseudopotential was created; no DFT was run.

Read [materials/DFT readiness](docs/AGF2_DFT_READINESS.md), [analytic symmetry check](docs/AGF2_SYMMETRY_DERIVATION.md), [sample DFT contract](inputs/agf2_dft_contract.example.json), and [machine-readable gate output](results/agf2_preflight.json).

## Provenance and separation

- Zhu et al., *PRL* 135, 056801 (2025), DOI `10.1103/dmzg-ck2t`: structural and magnetic-source facts; the numbered data above are **literature reports**, not new DFT results.
- Xiang et al., *PRL* 134, 116301 (2025), DOI `10.1103/PhysRevLett.134.116301`: IGMC mechanisms under independent reproduction; its Fig. 2(d) normalization ambiguity stays open.
- [`Xin-Jiaqi/group-theory-operations-toolkit`](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit) v0.18.0 is the background source for layer/magnetic-group capabilities; **this small gate derives the elementary diagonal-mirror tensors independently**, rather than claiming a live invocation of the toolkit's API. Revalidate against the native toolkit in the DFT structure's *actual* Hall/layer setting.
- No published figure, supplemental PDF, POTCAR, or unpublished atomic data is redistributed.
