# B+D Real-material Gate R0 | Bilayer MnPS₃ — evidence-first preflight

**Date:** 2026-10-21. **Status: LITERATURE PRE-FLIGHT / NOT SUBMITTED / NO MATERIAL RESULTS FROM THIS PROJECT.**  
**Purpose:** transition from synthetic B+D toy models to a published sliding-polar magnetic bilayer while using the existing GPT → GitHub → local DeepSeek harness → WSL/Japan HPC workflow. This is **independent research**, **not** reproduction/validation of Xiang–Jin–Wang PRL 135, 146303 (2025), and not a frozen DFT campaign.

## 1. What is established in the papers — with exact source boundaries

**Primary structural and symmetry source**

Sun, Wang, Yang, Huang & Cheng, *A unified symmetry framework for spin–ferroelectric coupling in altermagnetic multiferroics*, **Nature Communications 17, 3101 (2026)**, [DOI](https://doi.org/10.1038/s41467-026-69635-2) ([publisher full text](https://www.nature.com/articles/s41467-026-69635-2)). Its DFT-constructed bilayer MnPS₃ uses an oppositely oriented monolayer and lateral translation, then finds six degenerate global extrema in its **reported stacking energy map**:
`τ_o = ±(a,0)/3, ±(0,b)/3, ±(a,b)/3` in its notation. **Do not assume the raw coordinates have been retrieved from the SI or that a different cell uses the same fractional translations.**

It defines `P_n^±(L_m)`: `n=1,2,3` are three in-plane polarization/mirror orientations, `±` is out-of-plane polarization, `m=1,2` labels *two distinct magnetic interlayer configurations*. The authors report **L₂** as their lower-energy magnetic configuration, with an inter-magnetic-state energy difference ~0.22 meV (check unit/context in original text). For the **L₂ path-1 pair**:
- `P_1^+(L_2)` and `P_1^-(L_2)`;
- their **published connection** is `[E || M_z M_y]` (Nature Commun. Eq. 13), not a generic assumed `m_z T`;
- this connection is reported to give a **Type-II pseudo-time-reversal relation for spin bands**, whereas the L₁ path-1 connection `[E || M_z]` (Eq. 11) is Type I. It does **not** automatically imply that *every component of a SOC-included complex response function* simply flips sign: derive tensor and antiunitary frequency transformation for the actual observable.
- path 1 retains the in-plane polar/mirror orientation, switches the out-of-plane component and has the reported barrier **11.8 meV/f.u.** Paths 2/3 modify the in-plane orientation and report **9.2 meV/f.u.** barriers each.

Published computational *reference settings* (not yet a run recipe): VASP PAW, PBE, Dudarev `U_eff(Mn d)=4 eV`, DFT-D3, cutoff **500 eV**, centered **9×9×1** k mesh, vacuum **>15 Å**, and relaxed energy/force targets **10⁻⁶ eV / 10⁻² eV Å⁻¹**. Polarization evaluated by Berry-phase method; path barrier by climbing-image NEB. Published electrical transport used Boltzmann/Wannier constant relaxation-time spin-resolved charge conductivities, **not the proposed layer-resolved spin-density/local-torque AC Kubo**. **Before using these settings, validate the exact magnetization orientation, SOC setup, D3 variant, slab electrostatics, actual POSCAR and magnetic orders from the SI.**

**Other, overlapping REAL-MATERIAL work — DO NOT TREAT ALL AS INDEPENDENT VALIDATIONS**

1. Sun et al., *Symmetry-locked six-state control of altermagnetism via sliding ferroelectricity*, **Science Advances 12, eaec5229 (2026)**, [DOI](https://doi.org/10.1126/sciadv.aec5229). **Same core author group** as above, very similar MnPS₃ six-state setting. It treats spin-splitting texture and BCD-driven nonlinear Hall fingerprints. It is *separate publication, not automatically independent experimental confirmation*.
2. González, Brumme, Morell & León, *Engineering altermagnetism via layer shifts and spin order in bilayer MnPS₃*, **npj 2D Materials and Applications 10, 11 (2026)**, [DOI](https://doi.org/10.1038/s41699-025-00648-z). **Different authors, distinct stacking/ordering search**; AA-series lowest minimum is reported **7.67 meV/cell below AA′** series. Hence a reported local polar minimum in one layer-orientation family cannot be claimed *globally preferred among all possible MnPS₃ stackings*. Source also distinguishes parallel `↑↓/↑↓` and antiparallel `↑↓/↓↑` order.
3. Yang, Xu & Jin, *Symmetry-Programmed Pure Spin Photocurrents in Sliding Ferroelectric Altermagnets*, **ACS Nano (2026)**, [DOI](https://doi.org/10.1021/acsnano.6c06799). Already predicts symmetry-controlled **pure nonlinear spin photocurrent in bilayer MnPS₃** with SOC and other systems. Avoid calling generic switching of spin photocurrent new; exact peer equation-level comparison still required.
4. Tang, Zhang & Cheng, *Néel spin-orbit torque in antiferromagnetic quantum spin and anomalous Hall insulators*, **Nature Communications 16, 7790 (2025)**, [DOI](https://doi.org/10.1038/s41467-025-63171-1). Prior art for insulating E-induced staggered spin torque and AC magnetic resonance.
5. Zhu et al., *Sliding Ferroelectric Control of Unconventional Magnetism in Stacked Bilayers*, **PRL 135, 056801 (2025)**, [DOI](https://doi.org/10.1103/dmzg-ck2t). Prior art for structural/magnetic connecting-operator methodology.

**Important status of structural data:** publisher states relevant data are in the article and supplemental information, but **I have not verified/reconstructed a complete author-provided relaxed pair of POSCAR/CIF files or its magnetic moment ordering and numerical cell settings**. No structure is an accepted DFT input until the source coordinates and layer/site mapping are independently checked. *Do not copy a schematic into guessed fractional coordinates.*

## 2. Actual candidate, not another toy Hamiltonian

**One initial candidate only:**
`BD-MNPS3-L2-P1-plus-minus`:
published `P_1^+(L_2)` and `P_1^-(L_2)`, **path 1**, held in the **same documented magnetic-order convention**. Do not mix the `L_1` and `L_2` magnetic groups. Optional **next**, not part of the initial campaign: same `L_2` and path 2 (`P_1^+(L_2)→P_2^-(L_2)`), to compare geometry-only vs spin-label changes; do not launch it first.

Physics to falsify:

> The full complex **subgap or near-gap** electric-field-induced **layer-staggered spin-density** and **interlayer-electron-spin-transfer** spectra may be very different across *published FE switching paths*, even when a static total Edelstein/magnetoelectric measure or a Boltzmann spin-current texture cannot resolve the microscopic channel. A potentially large interlayer-transfer torque with local exchange/SOC torque cancellation must be demonstrated as a **material-specific dynamical observable**; standard operator Ward identities alone are not novelty.

Choose layer-resolved spin-density vector (S_{ell}^alpha) and (mathcal T_{-,mathrm{inter}}^alpha=i[H_perp,S_-^alpha]) with **unambiguous operator normalization**; split local torque into exchange and atomic SOC before inferring what drives magnetic order. The initial response tensor component (chi_{S_-^alpha,E_eta}(omega)) must be selected *after* the **SOC magnetic point group** and the exact `L_2` moment axis are fixed. **Do not assume our honeycomb toy's `S_-^x,E_x` component is symmetry allowed in the real cell.**

Primary scientific eliminators:

- Evaluate full causal (chi_{S_-,E}(omega)), (chi_{mathcal T_{mathrm{inter}},E}(omega)), (chi_{mathcal T_{mathrm{exchange}},E}(omega)), (chi_{mathcal T_{mathrm{SOC}},E}(omega)) in **one same-gauge** orbital/position/spin/layer representation; check electronic continuity at each omega.
- Test physically paired FE states with exact connection transformation, not an arbitrary chosen (pmdelta) or (pm lambda) sign. Compare the **whole complex response** with correct antiunitary causal transformation; simple plus/minus parity is not assumed without proof.
- Controls: `t_perp→0` (same monolayer state/orientation), `SOC→0`, `Néel→-Néel`, spectral cutoff, k-mesh, broadening, and full occupied/unoccupied completeness.
- Compare with **known mixed (k,Zeeman) Berry-curvature static response, Néel torque and ordinary layer Edelstein**. Do not rename these familiar mechanisms as original. If no independent torque-dynamical prediction survives, stop.

**Success is not a precondition:** the first scientific result may be a clean, source-grounded **NO-GO**.

## 3. READY vs BLOCKED, before anyone touches WSL/Japan jobs

| Stage | Current | Acceptance evidence required |
| --- | --- | --- |
| Primary peer identification, L₂ path-1 state names, structural barriers, reference VASP protocol | **SUPPORTED FROM PAPERS** | The DOI citations above and their SI |
| Two exact published **atomic structures** `P_1^±(L_2)` and matching magnetic order, relaxation provenance, license review | **BLOCKED / NOT EXTRACTED** | Author SI/archival coordinates OR independent source-driven construction from an authoritative monolayer structure, with explicit `M_z`, translation, magnetic layer mapping; hashes |
| True switching pair and actual SOC magnetic symmetry / allowed response tensor | **NOT VERIFIED FOR THIS PROJECT** | Two POSCAR/manifest hashes, magnetic moments and symmetry operation closures; symmetry check over tolerances |
| Same-gauge **spin, layer, position and electric vertices** for Wannier response | **NOT YET A VERIFIED CAPABILITY** | Validated operator export/embedding beyond HR alone; benchmark versus explicit DFT or minimal exact model; `physics/README.md` limitations apply |
| DFT setup on WSL or Japanese HPC | **NOT SUBMITTED** | Dedicated local DeepSeek task, platform validated, explicit user approval, PIN code/input commit and convergence standards |
| Real-material spin-density, torque, frequency spectra, energy/polarization, novelty | **NO RESULTS** | Native converged output, accepted Wannier/response metrics and reviewed peer equation matches |

### Why the present repos can be reused without pretending that everything is ready

- `materials-structure-core`, `extension-to-BSF`, `batch-symmetry-checker`, `group-theory-operations-toolkit`: **candidate structure/geometry/symmetry preparation**; do not invent magnetic order from spglib.
- Private `wannier-workflow`: hashed input, source reuse, dry-run, real fitting, stage-1 top-K. Its successful band fit does **not** certify spin/layer/position-operator export.
- Private `nlo-response-workflow`: actual implemented Stage 2 is shift current/SHG, **NOT the requested retarded layer spin/torque solver**. The new response kernel needs its own scientific oracle/contract; do not route through an unrelated validated NLO tensor.
- Private `wannauto-research-archive` tracks deployed versions/actual evidence and must not turn historical `RUNNING` into a live job status.
- `local-codex-bridge` includes a separate `hpc` instance template, but has **not** been connected/validated for our WSL/Japan scheduler in this chat.

## 4. Next concrete DeepSeek task — preflight only, no submission

**Task** `BD-MNPS3-L2-P1-STRUCTURE-R0` to be sent by user through their already configured local DeepSeek harness (or copied into its own private worktree).

```text
Goal: On real source material, prepare only an auditable, zero-DFT-cost
structure and magnetic-symmetry preflight for bilayer MnPS3.

Primary source: Sun et al., Nat Commun 17, 3101 (2026),
DOI 10.1038/s41467-026-69635-2; distinguish L1 vs L2.
Target states: P_1^+(L_2) and P_1^-(L_2), path 1, author Eq.13.

1. Obtain the original monolayer/reference structure (cite file/source),
   verify A'/antiparallel layer orientation, atomic species, 2D lattice,
   fractional coordinates, all Mn moments and site-to-layer mapping.
   Do not construct POSCAR from a screenshot.
2. If the published relaxed FE structures are downloadable, use them
   with exact provenance; otherwise build two unrelaxed candidate
   bilayers using extension-to-BSF from a validated monolayer,
   applying the author's mirror/orientation operation BEFORE shift.
   Mark generated structures GENERATED_NOT_DFT_RELAXED.
3. Verify that source-to-daughter mapping and connection operator agree
   with author Eq.13, with magnetic SOC point symmetry measured on the
   actual structures. Record every tolerance.
4. Write one private machine-readable manifest: POSCAR SHA256 for both
   states, magnetism, source DOI/SI references, code Git SHAs, symmetry
   verdict, missing metadata and human approval gates.
5. Produce a SANITIZED Markdown handoff without POTCAR, WAVECAR,
   CHGCAR, private raw data or login/queue credentials.
6. Do NOT submit VASP, do NOT queue Japanese HPC, do NOT fabricate
   coordinates or declare DFT convergence. If input files are missing,
   return BLOCKED with exact provenance gap and recommended next source.
```

**Second-stage task, only after preflight passes and the user explicitly approves costs:** a cheap paired magnetic/structural calculation on WSL first; confirm actual deployed W90Flow commit and strict convergence; then consider expensive SOC Wannier and HPC. No broad high-throughput scan. The hard question is whether this is a **distinct response** beyond published MnPS3 magnetic FE spin-band/Boltzmann/NLH/spin-photocurrent literature.

## 5. What GPT may subsequently update

This file, a small independently verified symmetry/response component, peer-equation comparisons, and sanitized acceptance reports may be updated in the **public** `research/brainstorm-20261021/` directory. Unpublished structures, complete Wannier HR/spin operators, remote paths and licensed software reside in private storage. If a material-specific result cannot safely be disclosed, the GitHub note only records its evidence level, sanitized summary and SHA/provenance pointer.

**Explicit provenance separation:** TOY_MODEL evidence lives in Gate 1–5 notes; MnPS₃ `R0` above is **PUBLISHED_CANDIDATE / REAL_COMPUTATION_NOT_STARTED**. Do not promote either into full Xiang original-paper reproduction.
