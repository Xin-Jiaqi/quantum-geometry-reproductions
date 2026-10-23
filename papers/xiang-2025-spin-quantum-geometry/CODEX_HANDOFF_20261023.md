# spinQG / layer-Edelstein Codex handoff (2026-10-23)

**Authoritative GitHub repository:** `Xin-Jiaqi/quantum-geometry-reproductions`. **This handoff branch:** `handoff/spinqg-codex-20261023`. **Goal:** local Codex can orient and begin *read-only scientific verification* without asking the user to relay assistant messages.

## What is on GitHub, and what is not

| Area | Actual source of truth | State |
|---|---|---|
| **Xiang 2025 PRL spin QG** reproducible paper subset | `papers/xiang-2025-spin-quantum-geometry/`, `docs/ACCEPTANCE.md`, `methods/`, `tests/`, `scripts/reproduce.py`, `results/` | **PRESENT in main**, source science **PARTIAL**. V32: 70 scoped tests reported in repo, *not rerun by this ChatGPT environment*. |
| Offline full V1–V27 Xiang research | User's offline archives | **NOT fully uploaded**. No claim of full original-paper recreation. |
| A/B/C/D research & MnPS3 jobs | `research/brainstorm-20261021/` `tasks/CURRENT.md` | **PRESENT** but its `CURRENT` is mainly MnPS3, R1 closed. Do not interpret as latest layer Edelstein task. |
| BD31 layer-registry exact symmetry | `research/brainstorm-20261021/docs/BD31_LOCAL_EQUIVARIANT_RESPONSE_FIELD_RULE_20261022.md`, group theory [Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32) | **PRESENT** idealized horizontal-mirror group and stabilizer test; no first-principles result. |
| BD32–BD35 PtSe2 model/theory and numerical audits | Independent snapshots recovered from conversation artifacts, **portable source ZIP** (listed below) | **NOT previously uploaded to GitHub**; full snapshots are provided as local downloadable archive, awaiting direct Git import. Summary and exact scientific caveats are written below so Codex does not rely on conversation. |
| BD36 BiDB PtSe2 and para preflight | [PR #99](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/99), `research/bd36-ptse2-bidb/` in that **separate** branch | **PRESENT in PR #99 branch**, not main. True para preflight FAIL: no .spn-capable VASP, no account Pt/Se POTCAR; zero jobs, zero core-hours. Latest `FAST_IDEA_GO_NOGO_20261023.md` is the proper next physical-gate plan. |

### Portable BD32–BD35 research source package (one-time import)

The ChatGPT conversation has a verified portable archive `spinqg_handoff_sources_bd32_bd35_20261023.zip` containing:
- `scripts/BD32_ptse2_inversion_registry_gate.py`, `BD32_PT_trace_model_gate.py`, `BD32_torus_index_gate.py`;
- `scripts/BD33_periodic_PT_spin_response_gate.py`, `BD33_plot_response.py`;
- `scripts/BD34_interlayer_Se_sublayer_PT_gate.py`, `BD34_make_figure_and_tables.py`;
- `scripts/BD35_independent_audit.py`;
- four reports `reports/BD32...md` through `BD35...md`;
- `reference/BD34_interlayer_channel_controls.json`, `BD34_interlayer_scaling.csv`, `BD35_audit_results.json`.

**ZIP SHA256:** `605843ccf73d64d29f05da7181b915429167337a852f2b025f6adfc1c66d1602` (43, - see size from actual archive: 42998 bytes). The file is **not on GitHub yet**; do not claim it was pushed. After obtaining it, stage files in `papers/xiang-2025-spin-quantum-geometry/research/lee-ptse2-bd32-bd35/` as *separate exploratory research*, preserving their source hashes. Scripts BD34/BD35 currently write JSON to hard-coded `/mnt/data/`; make output path configurable when importing and preserve an original snapshot.

**ChatGPT's direct spot checks** on locally extracted source on 2026-10-23:
```bash
python scripts/BD32_ptse2_inversion_registry_gate.py  # PASS: 6 high-symmetry site stabilizers; 144-grid even/odd
python scripts/BD32_torus_index_gate.py               # PASS: 3(+1)+3(-1) mathematical Fourier witness
```
These two scripts were independently executed in ChatGPT's container and exited with code 0. BD33/34/35 main loops have **not been freshly reexecuted** in this handoff pass; frozen BD35 results are from the earlier audit, not a fresh independent acceptance here.

## Essential physics for continuation: separate actual theorem from toy output

Observable: `S_+^z=S_top^z+S_bottom^z`, `S_-^z=S_top^z-S_bottom^z`, response `δS_±^z(u)=F_±(u)·j_xy`. `j_xy` is **charge current**, not spin current, and S is **nonequilibrium layer-resolved spin density**, not spin current.

**BD31 ideal horizontal mirror family** (e.g. idealized Pd2Te2 geometry with pure m_z at aligned registry): *conditional* leading local invariants `f_- ~ i γ z` (local index +1 if γ≠0) and `f_+ ~ β (z̄)^2` (index −2 if β≠0). Do **not** transplant these to PtSe2.

**BD32 actual identical parallel centrosymmetric 1T-PtSe2 geometry**: a moving inversion center for every in-plane displacement `u` ensures `F_+(u)=0` for the entire rigid-sliding family. The C3 fixed registries G=(0,0), M=(1/3,2/3), N=(2/3,1/3) force `F_-=0`; A=(1/2,0), B=(0,1/2), C=(1/2,1/2) have at most a 1D allowed response, **not forced zero**. This is an ideal rigid model and requires correctly tracking atom motifs and Seitz operations. For isolated nondegenerate zeros, if the three C3 zeros each have index +1, torus Euler index zero forces other zeros with total index −3, with a mirror-orbit parity constraint.

**BD32 Fourier witness**: periodic potential `V=cos(2πx)+cos(2πy)+cos(2π(y−x))+.16[sin(2πx)−sin(2πy)+sin(2π(y−x))]` intentionally places zeros at G/M/N and three satellites (~[.101003,.550502],[.449498,.550502],[.449498,.898997]); indices +1 and −1 respectively. **These are not PtSe2 predictions.** Hex Cartesian `g=A^{-T}∇_fractional V`, not direct fractional gradient.

**BD33 4-band toy**: PT doublet trace Boltzmann response can be nonzero; however isolated monolayer limit is spin-active, unlike an inversion-symmetric PtSe2 monolayer. Existence witness only.

**BD34 8-band toy**: sublayer `η`, whole monolayer `τ`, spin `σ`; hidden-spin + two distinct assumed interlayer channels; PT degeneracy, `S_+` zero and layer-odd nonzero in a selected parameter set. The registry field still comes from *prescribed Fourier potential*, not PtSe2 hopping fits. Parameters/energies are arbitrary dimensionless; not material calibrated.

**BD35 independent correction (critical, do not omit):**
`F_-(u)=C(|g(u)|²; μ,other params) J g(u)` with `J=[[0,-1],[1,0]]` from the model's artificially enlarged O(2) covariance. Four probed residuals ~1e-17; BD34 zero coordinates inherited from `g`, **not independent microphysical prediction**.
At `μ=.49835`, `C` switches sign along a path with `g≠0`: extra zero responses (often a continuous level set) can appear, so six-isolated-zeros is **not parameter-universal**. A tested root shifts ~(.2953,.4456)→(.3004,.4517) across k grids. Also `t=0, λ≠0` retains interlayer spinful hopping, so `t=0` is **not** the disconnected-bilayer limit; only `t=λ=0` is.

**BD36 actual BiDB**: 2252 monolayer, 2253/2255 parallel centrosymmetric bilayers, 2254 distinct-orientation bilayer; derived ±(.08,.19) rigid slides are not new BiDB entries. BiDB PBE gaps 2253 .219 eV, 2255 .614 eV are **database metadata**, no new SOC/Wannier results. para platform preflight: no current usable Pt/Se POTCAR and no `.spn` writer; do not reuse obsolete task. PR #99 now recommends *physics-first* G0 platform/read-only, then one generic-slide SOC spin-texture check with site-PROCAR doublet traces, then current-response proxy, then only if GO Wannier. G0/G1/G2/G3 are **not run**.

## Direct local Codex first run (without network compute)

```bash
git status --short --branch
python -m pip install -r papers/xiang-2025-spin-quantum-geometry/requirements.txt
cd papers/xiang-2025-spin-quantum-geometry
python -m unittest discover -s tests -v
python scripts/reproduce.py
```

Open `research/brainstorm-20261021/tasks/CODEX_LEE_CURRENT_20261023.md` for the independent exploratory track. Do not overwrite the historical `tasks/CURRENT.md` (MnPS3 read-only state). When importing the portable ZIP, verify SHA256, extract **only into the named exploratory directory**, run BD32 scripts, then inspect/repair BD33–BD35 portability separately. Log actual status and commit; no remote jobs.

## Code-vs-science definition of done

- CODEX_READS_CONTEXT: file paths and exact current tasks identified.
- SOURCE_PRESENT: original source archive imported to Git branch and checksum-verified.
- TESTS_RUN: independent clean-branch test results captured, no inferred PASS.
- BD32_SYMMETRY_VALIDATED: exact Seitz gate passes on local Python.
- BD35_CRITIQUE_PRESERVED: artificial response-direction locking and extra zeros explicitly retained.
- MATERIAL_RESPONSE_VERIFIED: **NOT RUN** until actual structure/hopping and layer-spin protocol are validated.
- NOVELTY_CLEARANCE: **OPEN**, direct prior-art includes 2026 LEE PRL and same-team sliding-REE companion not fully reviewed.

**User support needed, only if necessary:** one-time download/placement of the BD32–BD35 portable ZIP into local checkout; no repeated AI-to-AI copy-pasting. Afterwards Codex should update its own GitHub PR with script migration and test evidence, not ask the user to relay reports.