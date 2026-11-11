# Stage 8h FINAL: Japan PBS raw-output symmetry gate passed (external harness audited)

Date: 2026-10-21. **Independent original BSF × static-NMEE study; NOT Jia et al. (2026) driven-DC reproduction.**

## Provenance / what was actually checked

The user's local DeepSeek harness **read the original Japan IMR PBS OUTCAR/CONTCAR bytes onsite**, compared source/input/script hashes, inspected final VASP ionic/electronic convergence, and committed **sanitized SHA256-plus-derived-value evidence**:
- [STAGE8H_AUDIT_RAW_20261021.json](../hpc-stage8g/STAGE8H_AUDIT_RAW_20261021.json)
- [STAGE8H_REVIEW_CONCLUSIONS_20261021.md](../hpc-stage8g/STAGE8H_REVIEW_CONCLUSIONS_20261021.md)

ChatGPT inspected these committed reports and sanitized evidence, **but did not independently fetch the original remote OUTCAR/CONTCAR bytes**. This is a verified external-agent audit trail, distinct from byte-level verification by ChatGPT. Do not convert it to first-person claim of having directly inspected the raw files.

Three tasks (all PBS Exit_status=0, user-harness reported):
| s | PBS job ID | Final TOTEN eV/2 f.u. cell | ΔE meV/2 f.u. cell | max |Fz| eV/Å | final group / op order |
|:--:|---|---:|---:|---:|---|
|0 | 490621.cmspbs1 | -45.66480588 | 0 | 0.006172 | C3v, 6 |
| .25 | 490622.cmspbs1 | -45.65117241 | **+13.63347** | 0.006209 | Cs, 2 |
| .50 | 490623.cmspbs1 | -45.64068876 | **+24.11712** | 0.009974 | glide-containing order 4 |

The 15×15×1 PBE+D3(BJ), ENCUT 520 eV, non-SOC, ISYM=0, fixed-cell **F F T / z-only** protocol yielded final six-atom geometries with no lateral drift or cell relaxation. This is NOT fully relaxed sliding ferroelectricity, a CI-NEB, SOC electronic structure or a material static spin response.

## New material-specific geometric result

At **actual z-relaxed** s=.5, glide {m_z|t_half} survived, symmetry matching stable down to 10^-5 Å in the local audit. It maps the entire six-atom layer motif onto itself, including layer exchange. As a *selection rule*, for the gapped equilibrium/static source response M^±_ij=∂Bz χ(S^±_i,E_j), i,j=x,y:

- s=0: **C3v** forces M^+ = γ I_2 (one independent scalar).
- s=.25: **Cs** leaves two independent tensor entries.
- s=.5: total-spin tensor **rank of allowed invariant subspace ZERO**, while physically layer-staggered spin tensor has **allowed invariant subspace dimension TWO**.

The total-spin null constraint applies to the symmetry-preserving constrained geometry, regardless of whether microscopic SOC causes a finite layer-staggered response. There is STILL NO actual WSe2 NMEE coefficient, no evidence that M^- ≠ 0, and no claim the s=.5 elevated-energy reference is a stable third FE state.

## Stationarity but not a rigorously established switching barrier

Summed lateral forces on upper layer (W2,Se5,Se6) in original force tables:
- s=0: (0,0) eV/Å;
- s=.25: (-0.031622,+0.018257) eV/Å;
- s=.5: (0,0) eV/Å.

For arm A Δd=(a1-a2)/3≈(1.659503,-.958115) Å, the actual **rigid-slide energy derivative** at s=.25 obeys dE/ds = -F_upper·Δd ≈ +0.0700 eV per unit s (≈ +70 meV per cell per full slider coordinate). This is compatible with rising E(s). s=.5 is stationary by glide symmetry. The second finite energy difference E(.5)−2E(.25)+E(0)=−3.14982 meV/cell is **not** a saddle Hessian. A stationary high-energy point on the prescribed line may not lie on the relaxed minimum-energy switching path, and the transverse curvature remains unknown.

## Numerical and theory gates still open

1. Midpoint final max |Fz| 0.009974 eV/Å barely satisfies current 0.01 criterion; do not advertise convergence of ~meV response energies without k-mesh/ENCUT/vacuum/dispersion tests. Those are separate new job proposals.
2. Exactly protecting layer-exchange symmetry can fail when full ionic/lattice relaxation, external perpendicular electric field, substrate, twist or domain wall strain is included. The present geometric no-go is **not proof of universal robustness to all those conditions**.
3. Stage9's exact Ising [H,Sz]=0 no-go for transverse static spin response remains untested for real WSe2. Full-BZ spinor-SOC matrix elements and physical layer spin operators, not strong atomic SOC alone, decide whether M^- can be finite.
4. Current equilibrium/static mixed-Berry quantity is NOT automatically a current-driven layer Edelstein response and cannot inherit the source PRL's unresolved driven-DC claim.

**Decision:** real structural/symmetry gate PASS. Next priority should be **minimal SOC electronic feasibility / spin-mixing check** and a separately budgeted constrained-PES convergence gate, rather than automatic expansion of 12 path jobs or immediate full NEB/Wannier. **No new jobs authorized here.**
