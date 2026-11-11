# Stage 8h | Scientific review of REAL Japan PBS Stage8g z-only WSe2 pilot

**Evidence level:** \`USER_REPORTED_HARNESS_CONVERGED_REVIEWED / RAW_OUTPUT_NOT_YET_INDEPENDENTLY_OPENED\`. The numbers below were furnished in the user's full 2026-10-21 local DeepSeek-harness report. We have NOT independently read the original remote OUTCAR/CONTCAR or attached resume-state.json in this ChatGPT session. The next onsite read-only audit is distinct. This is independent BSF×NMEE research, NOT a claim from Jia et al. PRL, and the frozen paper reproduction remains untouched.

## Confirmed by the user's reported HPC audit (not rerun here)

Execution: Japan Tohoku IMR PBS Pro \`cmspbs1\`, directory \`/work/sunq/xjq/runs-stage8g-20261020\`, queue CP_001, 16 cores/32GB/12h; Git commit \`52f126345ee3d7a89d3b6c4ea199a359f35bd094\`. VASP 6.4.3, oneAPI MPI 2023.2, locally licensed PBE-PAW W/Se (\`potpaw_PBE.54\`, W ENMAX 223.057 eV, Se ENMAX 211.555 eV, W then Se), POTCAR SHA256 \`7ef8b236d47cb52384021b0d5255b9dd2cd6fb636f1b7e0a093eff353b3a2a5c\`; no POTCAR or original binary is committed. ENCUT 520 eV, 15×15×1 Gamma, PBE+D3(BJ) (\`IVDW=12\`), no SOC; fixed cell, fixed x/y atomic fractional positions, ion-only z relaxation; EDIFF=1e−7 eV, EDIFFG=−.01 eV/Å, dipole correction along slab normal. WSL checked as a local Ubuntu agent host, not a separately available VASP scheduler.

| Real PBS Job ID | Slide fraction s | Reported final TOTEN (eV/cell) | ΔE relative AB (meV/cell) | ΔE relative AB (meV/f.u., 2 f.u. per cell) | final max |Fz| eV/Å | ionic steps |
|---|---:|---:|---:|---:|---:|---:|
| 490621.cmspbs1 | 0 | −45.664806 | 0 | 0 | .00617 | 16 |
| 490622.cmspbs1 | .25 | −45.651172 | +13.634 | +6.817 | .00621 | 13 |
| 490623.cmspbs1 | .50 | −45.640689 | **+24.117** | **+12.0585** | **.00997** | 6 |

The three reported outputs possess OUTCAR/CONTCAR, ionic \`reached required accuracy\`, final SCF EDIFF convergence and frozen-xy/cell integrity. **This is real external DFT evidence as reported by the harness**, not a calculated NMEE tensor. Remaining 12 generated inputs were not submitted.

### Physically justified inference from the actual report

1. Along this particular **constrained z-only PBE+D3(BJ) AB→BA shortest-arm** sample, energy increased monotonically at the sampled half-path nodes (0, .25, .50); s=.5 exceeds s=.25 by **10.483 meV/cell**. The sampled finite second difference E(.5)−2E(.25)+E(0)=**−3.151 meV/cell**, but this is NOT the analytical saddle curvature or a continuous-path barrier.
2. Stage8f predicts mirror-related \`E(s)=E(1−s)\` and \`M_total(s)=−M_total(1−s)\`, \`M_staggered(s)=+M_staggered(1−s)\` **if the actual relaxed geometry preserves the full atomic-motif exchange symmetry**. At s=.5 the ideal nonsymmorphic \`{m_z|t}\` forbids uniform total-spin **equilibrium/static** NMEE, but permits nonzero layer-staggered static NMEE. This conclusion is a symmetry selection rule only: it does NOT establish real WSe2 spin response or its magnitude.
3. A structurally stationary s=.5 under the symmetry can be a plausible switching-barrier candidate. However, three sampled constrained energies are not a CI-NEB, do not locate the true minimum-energy path, and do not prove that the midpoint is the **continuous** path's highest-energy saddle. Off-path relaxation or other stacking routes can reduce the true barrier.
4. \`|Fz|=.00997\` eV/Å at s=.5 is very close to the .01 threshold. The 15×15×1 k-mesh, XC/dispersion choice, vacuum and dipole treatment and ENMAX criterion pass only **initial setup checks**, not physical meV-scale convergence.
5. The actual relaxed **s=.5 CONTCAR layer-glide symmetry was NOT explicitly demonstrated by the supplied report**. With selective F F T relaxations performed at ISYM=0, tiny or spontaneously symmetry-breaking z displacements can remove the glide no-go even though initial atomic input was glide-symmetric. This is the immediate highest-priority test. At s=.25, the remaining vertical mirror allows a restricted nematic response; at AB C3v, static M_ij must be scalar if symmetry survives.
6. x/y forces are **not expected to vanish atom-by-atom** under an F F T constraint. The physically useful sliding-stationarity diagnostic is the **sum of x/y Hellmann–Feynman forces over all upper-layer atoms**. Stage8h read-only checker reports these sums instead of incorrectly treating individual forbidden-direction forces as ionic convergence failures. Use the projected sum along the actual Cartesian A-arm as derivative −∂E/∂s.

## Stage 8h: ZERO new PBS jobs, an independent raw-evidence check

Run the new \`hpc-stage8g/review_phase1.py\` from the Git research branch against existing authorized VASP run directories on Japan/WSL. It reads all three actual OUTCAR/CONTCAR **locally at the HPC**, enforces original manifest input hashes, re-extracts independent final forces/cell/xy integrity, then enumerates whole six-atom point operations at a tolerant 0.001 Å z threshold and reports final permitted **TOTAL and LAYER-STAGGERED** static tensor dimensions. It separately reports upper-layer net lateral force and recomputes ΔE from actual final TOTEN.

Expected if the relaxed coordinates preserve Stage8f:
- AB s=0: C3v order 6, total tensor dim 1;
- s=.25: vertical mirror group order 2, total tensor dim 2;
- **s=.50: glide-containing order 4, total tensor dim 0, layer-staggered dim 2.**

The checker can produce \`RELAXED_SYMMETRY_NOT_CONFIRMED_STOP_BEFORE_NMEE\` even with electronically/ionically converged VASP output. That would be a scientifically significant **failed symmetry gate**, not a bug to paper over. It only certifies actual geometry/symmetry plus local force checks, NOT convergence of energy differences or actual NMEE.

**Security:** run on the authorized HPC and retain VASP original outputs locally; only a sanitized evidence JSON with input/OUTCAR hashes and derived quantitative checks may be handed back or archived after user-controlled disclosure. Do not leak POTCAR, site credentials, raw large logs or private records.

## Next compute decision AFTER read-only recheck (not authorized now)

Prioritize explicit **comparison-quality gates**, not an unreviewed 12-job batch:
- tighten s=.5 Fz tolerance if needed (ideally ≤.005 eV/Å); preserve and compare symmetry after relaxation;
- repeat **only** s=0 and s=.5 at a denser in-plane k mesh (e.g. 21×21×1) with *otherwise identical* PBE-D3/POTCAR/ENCUT/dipole protocol and careful source-version tracking, to assess the ~24 meV/cell energy rise;
- for mirror relation independently compute/compare s=.75 and BA endpoint s=1 on this arm only **if new jobs are explicitly approved**, rather than assuming mirror pair energies from theory;
- only with convergence and symmetry accepted, select a complete path, refined 1D PES, full-ionic endpoints and CI-NEB. SOC full-BZ insulating-gap checks and spinor Wannier static/layer-resolved NMEE remain later independent gates.

**No additional PBS jobs are authorized by this document.** Its report is a user-reported real-computation record with a pending independent raw-output symmetry audit, not a published prediction.
