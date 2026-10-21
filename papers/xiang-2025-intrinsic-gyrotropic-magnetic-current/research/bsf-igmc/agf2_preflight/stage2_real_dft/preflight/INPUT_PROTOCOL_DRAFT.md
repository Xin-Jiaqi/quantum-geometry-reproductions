# AgF₂ Stage-2 · Input Protocol Draft (G0: prepared, NOT submitted)

Status: **BLOCKED_PROTOCOL / NOT_RUN** — inputs are ready for human approval;
no job was submitted (`actual_jobs_submitted=0`, `submitted_job_ids=[]`).

## 1. Campaign

First campaign (only): three B1 initial candidates from
`stage2_geometry/candidates/B1_g2.6|g3.0|g3.4` (unrelaxed, C2DB-derived,
12 atoms Ag₄F₈, species order `Ag F`, ~20 Å vacuum; **not** author-relaxed
structures). Physics: **0% strain, Type-I interlayer Néel OPPOSITE, collinear
spin-polarized, non-SOC**, fixed in-plane cell, ions free (sliding not
artificially constrained).

| case | POSCAR sha256 (full) | atoms | cell (a,b,c Å) | vacuum (Å) |
| --- | --- | --- | --- | --- |
| B1_g2.6 | 966174568594104e… (see manifest below) | 12 (4 Ag, 8 F) | 3.8270, 5.8797, 26.5385 | ~20.00 |
| B1_g3.0 | 2d323c3fe2353ef9… | 12 (4 Ag, 8 F) | 3.8270, 5.8797, 26.9385 | ~20.00 |
| B1_g3.4 | 04282d8a1756333a… | 12 (4 Ag, 8 F) | 3.8270, 5.8797, 27.3385 | ~20.00 |

Full SHA-256 and provenance: see `stage2_geometry/candidates/GEOMETRY_MANIFEST.json`
and the per-file hashes below. B1_g3.0 fractional coordinates are **identical to
the Stage-1 frozen B1 candidate** (max fractional residual 0.0) — cross-checked.

### Full SHA-256 of the three input POSCARs
- B1_g2.6: `966174568594104ebf8faadfd02ce3025498bc7b5f890d4fab677d5d0a6439b2`
- B1_g3.0: `2d323c3fe2353ef984330b3b99fc24a75a85337883bbaf36d08b20ed45e6d60a`
- B1_g3.4: `04282d8a1756333aaa6d431ca049b801e53a70697a9afee53f5539fe56a0aa48`

(All reproducible via `sha256sum stage2_geometry/candidates/B1_g<g>/POSCAR`.)

## 2. Per-case private input files (prepared in a private workdir, not committed)

Each case directory contains: `POSCAR` (copy of candidate), `INCAR`, `KPOINTS`,
`run.pbs` (unexecuted PBS template), `POTCAR_checksum_record.md`.

| file | B1_g2.6 sha256 | B1_g3.0 sha256 | B1_g3.4 sha256 |
| --- | --- | --- | --- |
| INCAR | 32118d04… | 31c677cb… | ea9fb7ff… |
| KPOINTS | ff823701… | ff823701… | ff823701… |
| run.pbs | 43e733ff… | 43e733ff… | 43e733ff… |

## 3. Parameter classification (per CALCULATION_GATE / Stage-1 review)

### 3.1 Frozen from literature (Zhu et al. PRL 135, 056801 (2025), SM)
- XC: PBE (PAW, LEXCH=PE)
- Hubbard: Dudarev `U_eff(Ag-4d)=4 eV` → `LDAU=.TRUE. LDAUTYPE=2 LDAUL=2 -1 LDAUU=4 0 LDAUJ=0 0`
- Dispersion: DFT-D3 family (Grimme); **`IVDW` damping variant must be documented**。根据文献所引 Grimme 2010 D3 方法及 VASP 手册，首选核对 `IVDW=11`（D3 zero damping）；`IVDW=12` 是不同的 D3(BJ) 方法，不能静默混用。需结合原作者 SM/计算输入和所选 VASP 实际输出确认。
- `ENCUT=500 eV` (≥ both PAW ENMAX: Ag 249.8, F 400.0 — factual from TITEL/ENMAX headers)
- Vacuum ≈ 20 Å (verified in all three POSCARs)
- `EDIFF=1E-6 eV`, max residual force `< 0.01 eV/Å` (`EDIFFG=-0.01`)
- Magnetic input: **Type-I = interlayer Néel OPPOSITE** (Stage-1 review
  conclusion; the connector `{m_001|t₁}` flips the in-plane axial spin at the
  top images)

### 3.2 Operator choices (explicit, not author parameters — each to be recorded in the convergence record)
- k mesh: **12×8×1 Gamma-centered** (screening candidate); 14×9×1 and 16×10×1
  reserved for convergence checks. **Not claimed converged.**
- **科学审查不予批准原候选：** 私有输入草稿曾建议 `ISMEAR=1, SIGMA=0.1`，并错误称为“metallic Ag”。AgF₂ 不是元素 Ag 金属；对于待检验的磁性绝缘体/未知带隙，应先采用 `ISMEAR=0, SIGMA=0.05 eV` 的低成本单点与弛豫候选，并检查占据、展宽和带隙。**此处仅更新公开审查结论，私有 INCAR 仍须在 G1 前实质改写、重新计算 SHA256 并核验。**
- `MAGMOM` seed magnitude m=1.0 μB per Ag — TO CONFIRM (see §4)
- `ISIF=2` (ions relax, cell fixed; keeps the in-plane lattice fixed, also
  fixes c) — chosen for the pilot; ISIF=4 alternative to be discussed
- `ISYM=0` (avoid symmetry-op spin handling in the collinear pilot) — TO CONFIRM
- `LDIPOL=.FALSE.` (pilot) — TO CONFIRM
- `PREC=Accurate`, `ALGO=Normal`, `LREAL=Auto`, `NCORE=8` (32-core run) — TO CONFIRM
- `LWAVE=.FALSE.` screening; `LCHARG=.TRUE.`
- PAW dataset revision: user private PBE.64 set (Ag 02Apr2005, F 08Apr2002) vs
  site shared older set (Ag 06Sep2000) — TO FINALIZE with user
- VASP cpu build: candidate `vasp.6.6.1_cpu/bin/vasp_std` (versions 6.4.3 /
  6.5.1 / 6.6.0 / 6.6.1 cpu builds observed; none executed in G0) — TO CONFIRM
- PBS shape: `-q CP_001 -l select=1:ncpus=32:mpiprocs=32:mem=120gb
  -l walltime=24:00:00` (site-verified syntax; within queue maxima
  56 ncpus / 230gb / 72h) — TO FINALIZE with user + approved budget

### 3.3 Missing / BLOCKED (kept null or BLOCKED_*, not fabricated)
- Compute-budget approval from the user: **false** (default per task)
- VASP executable version **not execution-verified** (running VASP is not
  permitted in G0)
- IVDW D3 variant on site: unverified
- POTCAR content: **not uploaded/committed** (private; element order Ag→F;
  reference combined sha256 `8c0da262849433cad092946af66c1b584d457e46a70d4c1eb768ea68305f0fff`
  recorded from the private PBE.64 set)

## 4. Type-I initial magnetic-moment ordering (MAGMOM) and per-atom mapping

B1 POSCAR Ag order (e.g. B1_g3.0):
| idx | frac (x,y) | layer | sublattice | moment (Type-I) |
| --- | --- | --- | --- | --- |
| Ag0 | (0.5,0.5) | **upper-z reference layer** | B | −m |
| Ag1 | (0,0) | **upper-z reference layer** | A | +m |
| Ag2 | (0,0.34) | **lower-z mirror-image layer** | image of B | +m (flipped) |
| Ag3 | (0.5,−0.16) | **lower-z mirror-image layer** | image of A | −m (flipped) |

`MAGMOM = -1 1 1 -1 0 0 0 0 0 0 0 0` (F atoms zero). Basis: Stage-1 review —
monolayer altermagnet Ag(0,0)=+x̂, Ag(½,½)=−x̂; the stacking mirror
{m_001|t₁=(0.5,−0.16)} maps bottom→top; physical spin is an axial vector so
the unitary connector flips the in-plane S_x; therefore Type-I (Néel OPPOSITE)
mirror images carry the flipped sign. **Important coordinate convention:** in the frozen B1 POSCAR, Ag0/Ag1 have the larger physical z coordinate and Ag2/Ag3 the smaller z; earlier narrative labels “bottom/top” incorrectly referred to construction roles, not geometric height. The sign sequence itself is unchanged. Non-SOC collinear input: direction is a seed, not a magnetocrystalline-anisotropy result.

## 5. PBS template (site-verified syntax, unexecuted)

```bash
#PBS -N agf2_B1_typeI
#PBS -q CP_001
#PBS -l select=1:ncpus=32:mpiprocs=32:mem=120gb
#PBS -l walltime=24:00:00
#PBS -j oe
source <site oneAPI/MPI loader>      # TODO: confirm exact loader with user
export OMP_NUM_THREADS=1; ulimit -s unlimited
cd "$PBS_O_WORKDIR"
VASP_BIN=<candidate vasp.6.6.1_cpu/bin/vasp_std>   # TO CONFIRM
# POTCAR assembled privately on cluster at run time: cat <Ag> <F> > POTCAR
exec mpirun -np 32 -ppn 32 -hostfile "$PBS_NODEFILE" "$VASP_BIN"
```
The template in the private input dir carries a `BLOCKED_PROTOCOL` header and
TODO placeholders for the loader path and VASP binary; it must not be run
until the user approves budget/queue/binary/PAW set.

## 6. Independent physics audit (2026-10-21)

This report is a **prepared-input snapshot**, not a ready-to-run scientific protocol. Changes to comments here do **not** automatically modify private INCAR files or their existing SHA-256 values. Before any PBS submission: (i) regenerate the actual private INCARs with verified Gaussian smearing, (ii) freeze the selected D3 variant and complete the VASP/PBS interface, (iii) compare **effective VASP tags** across all three cases (their draft INCAR hashes differ), and (iv) record the resulting updated SHA-256 hashes. Preserve `BLOCKED_PROTOCOL`/`NOT_RUN` until user budget approval and an actual successful calculation.
