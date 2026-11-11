> **Independent scientific ERRATUM (2026-10-22):** This original external-agent Stage8i report is retained for provenance, but its ~1.3-eV “K/K′ valley offset”, 87.5/35.3-meV “Ising spin splitting”, and inference of Sz nonconservation from the Γ-point Kramers spin texture **are not scientifically supported as labeled**. K′ must be −K; the old Kp (−1/3,+1/3) is Γ→M interior, and the K band35/36 spin polarizations have the SAME Sz sign. Read the detailed [independent Stage8j correction](STAGE8J_INDEPENDENT_VALLEY_SPIN_ERRATUM_20261022.md) before citing any spin/valley result. Original external observations are preserved below; sampled positive gap and SCF electronic convergence remain provisionally supported, with MPI exit-255 and incomplete PROCAR retained as limitations.

# Stage8i Review Conclusions — WSe2 SOC Electronic-Structure Preflight (AB s=0.00, Glide s=0.50)

Date: 2026-10-22 · Branch: research/bsf-nmee-stage7-sw-gate-20261020
Scope: two 21×21×1 noncollinear (SOC) VASP SCF runs, per TASK_CONTRACT (electronic
feasibility gate only; static NMEE and Kubo response NOT computed).

## 1. Execution facts

| job | PBS id | binary | queue | TOTEN (eV) | SCF | OUTCAR sha256 |
|---|---|---|---|---|---|---|
| soc_AB_s0.00 | 491186.cmspbs1 | vasp_ncl 6.4.3 | CP_001 16c/32G/06:00 | −46.07647873 | EDIFF reached (DAV, last dE −3.2e−17) | 8f88851d…c534b |
| soc_GLIDE_s0.50 | 491187.cmspbs1 | vasp_ncl 6.4.3 | CP_001 16c/32G/06:00 | −46.05232372 | EDIFF reached | 0b13aa2a…e3acb |

- Protocol: PBE+D3(BJ), ENCUT=520 (≥1.3×ENMAX both species), 21×21×1 Gamma, NBANDS=96,
  NSW=0, ISYM=0, EDIFF=1e−7, LORBIT=11, LDIPOL/IDIPOL=3, MAGMOM=18×0, nonmagnetic spinor SOC.
- Inputs immutable: both POSCARs byte-identical to Stage8h relaxed CONTCAR (SHA-verified by
  prepare_soc.py fail-closed gate); POTCAR (W 08Apr2002, Se 06Sep2000) sha 7ef8b236…a2a5c.
- ΔE(glide − AB) = +24.16 meV/cell (SOC) — consistent with +24.12 meV/cell non-SOC (Stage8g), cross-check passed.

## 2. Sampled gap & band-edge results (441 k, 96 spinor bands, NELECT=36)

| quantity | AB s=0.00 | Glide s=0.50 |
|---|---|---|
| sampled indirect gap | 1.0881 eV | 1.1490 eV |
| sampled min direct gap | 1.1527 eV | 1.2124 eV |
| VBM (band36) | −0.1448 eV @ Γ | −0.2025 eV @ K |
| CBM (band37) | 0.9433 eV @ λ≈(1/6,1/6) | 0.9465 eV @ λ′≈(−0.19,0.38) |
| K-valley VBM | −0.1765 eV (Γ higher by 32 meV) | −0.2025 eV (highest) |
| K-valley valence spin split (Ising) | 87.5 meV | 35.3 meV |
| K′ valley VBM | −1.4756 eV | −1.4469 eV |
| K vs K′ valley offset | 1.30 eV | 1.24 eV |

Sampled-gap caveat: k-grid sampling only; continuous-BZ gap NOT claimed.

## 3. Symmetry checks (EIGENVAL)

- Time reversal: band-for-band E(K)=E(−K) verified ((1/3,1/3) vs (−1/3,−1/3) rows identical).
- C3 rotation (hexagonal fractional map (u,v)→(−(u+v),u)): K is C3-invariant; the K′ family
  (three equivalent vertices) is internally C3-consistent.
- K vs K′ splitting (~1.3 eV) is a real 3R-bilayer valley effect: under C3v no operation links
  K and K′ (no C2′ / σh, and −K belongs to the K family), so unequal K/K′ energies are allowed
  and are consistent with symmetry. Earlier apparent “C3 violation” traced to an incorrect
  reciprocal-space rotation map, not to the calculation.

## 4. Spin texture (plane-wave spinor coefficients from WAVECAR, pymatgen ncl parse)

| point | AB s=0.00 | Glide s=0.50 |
|---|---|---|
| K valley (split 87.5/35.3 meV) | |Sz|=0.498–0.500, Sx=Sy=0 → pure Ising | |Sz|≈0.499, Sx=Sy=0 → pure Ising |
| Γ VBM Kramers pair | S=(−0.181,0.441,0.080); partner reversed; strong in-plane (non-Ising) mixing | S=(−0.036,−0.460,0.135); in-plane mixing |
| Γ CBM pair | Sz=±0.436, Sx≈Sy≈0 (Ising-like) | mixed Sz=±0.30, Sx,Sy small |
| K′ deep valley | in-plane (Sz≈0) | in-plane (Sz≈0) |

Interpretation: the band-edge K valley is Ising (spin quantized along z), while the Γ
valence-band-edge Kramers pair shows genuine in-plane spin components — a signature that
[H,Sz]≠0 locally. Because the Γ pair is degenerate, its texture carries gauge freedom and is a
proxy (plane-wave/pseudized part only, PAW and gauge dependent); it is NOT proof of nonzero
interband spin-mixing matrix elements.

## 5. Degraded outputs & warnings

- Both jobs exited with Exit_status=255 at MPI finalize; SCF had converged (EDIFF reached),
  EIGENVAL/OUTCAR/CHGCAR/WAVECAR coefficient body complete and mutually consistent
  (spin textures physical, TR pairs opposite). OUTCAR ends without “Total CPU time used”
  (write stopped inside LORBIT=11 projection output).
- PROCAR: 21 bytes (header line only) — writing interrupted; PROCAR-based spin projections unusable.
- WAVECAR header record is nonstandard (recl=0-style first words; only 3 floats at offsets 4/12/20);
  pymatgen’s Wavecar parsed it (vasp_type=ncl, nk=441, nb=96, encut=520) and the coefficient body
  reads correctly; pymatgen band_energy read is unreliable for this file → energies taken from EIGENVAL.
- Magnetization (x/y/z) all ~0 → nonmagnetic as intended.

## 6. Gate conclusion

Level: ELECTRONIC_PREFLIGHT_PASS — STATIC_NMEE_NOT_CALCULATED
(equiv. SCF_COMPLETE_BUT_FEASIBILITY_AMBIGUOUS for NMEE).

Supported: converged spinor SOC SCF, positive sampled gap (indirect ≈1.09/1.15 eV,
direct ≥1.15/1.21 eV), K-valley Ising texture (|Sz|≈1/2, Sx=Sy=0), coexisting in-plane
texture at Γ (gauge-dependent), ΔE(glide−AB) consistent with non-SOC.

NOT supported / not claimed: continuous-BZ gap; nonzero interband off-diagonal spin matrix
elements (needs Stage7j Wannier/operator gate); any static NMEE response; Kramers-pair texture
as quantitative physics without gauge fixing.

## 7. Recommended next physics gate

1. Stage7j: Wannierization of the two SOC band windows (Γ+K edge) and explicit evaluation of
   off-diagonal spin/velocity matrix elements ⟨c|ŝi|v⟩ — the decisive Kubo go/no-go.
2. If matrix elements are nonzero and gaps survive: spinor Wannier static-response calculation
   (NMEE) for s=0.00 and s=0.50 — NOT authorized yet.
3. Revisit with denser k-grid (e.g., 33×33×1) only if a new task is authorized.

## 8. Compliance

- Exactly 2 SOC SCF jobs submitted (one per structure); no further jobs started.
- No POTCAR/credentials/site.json/raw outputs committed; evidence files are derived JSONs + this MD.
- Static NMEE not run; 12-job/NEB/Wannier/NMEE items remain unauthorized.
