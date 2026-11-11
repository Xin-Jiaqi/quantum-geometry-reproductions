# Stage 8i proposal: minimal SOC feasibility before spending on complete slide/NEB

**Status: DESIGN ONLY, NOT AUTHORIZED, NOT SUBMITTED.** Applies to user-harness-confirmed Stage8h Japan IMR VASP WSe2 z-only relaxed source geometries. Exact scientific results of Stage8h are stored in [STAGE8H_FINAL_SCIENTIFIC_ACCEPTANCE_20261021.md](STAGE8H_FINAL_SCIENTIFIC_ACCEPTANCE_20261021.md). Maintains separate original research; frozen Jia et al. PRL reproduction untouched.

## Why the next scientific calculation is NOT 12 more non-SOC slide jobs

We established a nontrivial symmetry contrast for real z-only constrained geometry:
- AB (C3v) allows an in-plane SCALAR total static spin-magnetoelectric coefficient;
- half-slide keeps a layer-exchanging glide so the TOTAL tensor is strictly zero, while the LAYER-ODD tensor has a two-dimensional allowed space.

This is a **symmetry theorem, not a response amplitude**. The main potential falsification is Stage9:
\[
[H_{\rm SOC}(\mathbf k),S_z]=0
\Longrightarrow \Gamma^{\rm stat}_{xz;x}=\Gamma^{\rm stat,\;odd}_{xz;x}=0
\]
for transverse spin sources and physical layer operators commuting with Sz, even if W atomic spin-orbit splitting is large. We must learn whether the real WSe2 bands have *sufficient non-Ising spin mixing*, whether they remain insulating over the FULL BZ, and whether a physical layer-resolved spin operator can be formed.

## Efficient staged resource proposal (distinct authorizations)

### Phase A: ZERO additional PBS jobs — already admissible now

1. Freeze stage8h evidence, validate the original six-atom final CONTCAR symmetries and hashes, record the tested tolerance and the upper-layer force as a rigid-slide gradient; done by local DeepSeek.
2. Locate existing SOC/DFT WSe2 results in the current authorized research archive or literature without claiming they substitute for the exact s=.5 stack. Diagnose W, Se PAW spinor channel and code compatibility in read-only mode.
3. Draft one or two identical-protocol **nonmagnetic spinor-SOC electronic feasibility** inputs using existing relaxed AB and s=.5 geometries, with all inputs versioned/hashes frozen. The necessary VASP tags include LSORBIT=.TRUE. (and noncollinear mode), a consistent SAXIS axis, ISYM=0 initially, high-quality charge/spin density and recoverable k-resolved spin projections; VASP may require a site-specific executable and k-point restrictions. Preserve the same licensed PBE-PAW species W→Se and D3/dipole/cell protocol. Do not silently convert this to an energy-barrier calculation. Do NOT submit until separately approved.

### Phase B: suggested first NEW DFT jobs **only after approval**

**Primary science option: 2 spinor-SOC single-point SCF calculations** on the actual stage8h CONTCARs:
- AB s=0 (C3v), candidate for finite longitudinal \(B_z E_x\) static response;
- glide midpoint s=.50, where total response must be zero for retained geometric symmetry but physical staggered-spin response may be finite.

Analyze full-BZ sampled gaps, K/K', Gamma, and Lambda band edges; use sufficiently dense mesh for a genuine first-pass 2D gap estimate but do not claim an exhaustive/full-BZ gap from a few high-symmetry path plots. Inspect physical spin textures and *off-diagonal spin-mixing proxies/operators* in the valence→conduction window, especially whether \(S_z\) is approximately conserved. K-resolved \(\langle S_{x,y,z}\rangle\) from PROCAR is only a diagnostic: vanishing diagonal \(\langle S_{x,y}\rangle\) does NOT prove \([H,S_z]=0\), and nonzero spin textures do NOT alone give the required cross-gap Kubo products. Do not claim \(\Gamma\) from band plots or orbital projections.

**Secondary convergence option: 2 non-SOC z-only comparison jobs**: 21×21×1 for s=0 and s=.50, otherwise identical parameters, with same source geometry/updated residual-force threshold and actual post-relaxation symmetry check. Separate k-mesh numerical changes from tighter EDIFFG; do not mix protocols when estimating energy changes. This tests the current +24.11712 meV per 2-f.u. cell, NOT NMEE. Could be performed after SOC feasibility rather than as the primary scientific deliverable.

**Optional** one midpoint tighten-only z re-relaxation to EDIFFG=-.005 eV/Å, then repeat symmetry/relative energy audit. The current 0.009974 residual makes this useful but does not supersede the SOC material-physics go/no-go.

**Do not yet submit** s=.75, s=1, the other two C3 arms, full CI-NEB, a 12-job batch or spinor Wannier/NMEE response workflow merely because they have generated input files. Such work needs a verified finite observable and physically consistent spinor operator inputs.

## Stop/go acceptance

**STOP material static-NMEE push** if: the two representative structures turn metallic over a reliably sampled full 2D BZ; required transverse spin mixing is vanishing/negligible in a consistently converged spinor analysis; or layer-spin projectors are not physically defensible and stable against reasonable spatial partition changes. A merely nonzero Ising spin splitting does NOT pass.

**GO to Wannier feasibility** if: a stable occupied gap and significant spin-nonconserving interband pathways exist across actual spinors, preserved exact stack symmetry is confirmed, and physical-layer spin operators can be implemented in the same orbital Hilbert space. Wannier requires the proper position/Berry connection and remote-band completeness; HR alone is not an acceptable numerical response basis.

## Publication decision

The physically distinct falsifiable signature under C3v is ordinary nonequilibrium current-driven spin Edelstein \(\delta\mathbf S_\parallel\propto\hat z\times\mathbf E\) versus the allowed **longitudinal, Bz-odd** equilibrium/static susceptibility \(\delta\mathbf S_\parallel\propto B_z\mathbf E_\parallel\). Its angular rule alone is group theory, NOT novel quantitative materials physics. The potential original mechanism requires a real finite static spin response, a finite (not merely symmetry allowed) hidden layer-staggered response at the glide midpoint, and a microscopic SOC/remote-band/quantum-geometry analysis.

**Execution boundary:** this document requests an explicit computational plan review, authorizing NO new HPC resources. Do not mark existing Stage8h outputs as Wannier results or compute \(\Gamma\) from unavailable data.
