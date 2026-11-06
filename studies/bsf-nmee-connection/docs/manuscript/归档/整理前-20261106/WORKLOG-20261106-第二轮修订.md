# WORKLOG — 2026-11-06 second manuscript revision

## Request and outcome

Replace the generic quarter displacement by verified AB/BA where supported; simplify geometry; separate captions from physical argument; repair Markdown table; identify valence band groups; inspect the user's pentagonal-chain Fig.4; unify sign conventions without changing data.

AB/BA are horizontal-mirror partners. Main Fig.1 now uses independent AB/BA calculations, not relabelled quarter-displacement data. The quarter displacement and its pure-compensation gate path are retained in Fig.S4.

## New physics / correction to an earlier extrapolation

At m=45 meV, AB lower C changes from negative to positive between t=10 and12 meV at U=0. At t=10 meV, changing U from0 to2 meV also reverses the lower response. The AB net C stays positive over the computed U=0..12 meV range. Thus the old quarter-path gate result cannot be directly transferred to AB. Independent BA respects layer exchange and sign reversal. Do not claim that current AB reproduces the opposite compensation fields found for MoSSe; its shared results are local reversal and mirror relation.

Decoupled bias control: U=−12,0,+12 meV does not change a layer response while occupations remain insulating. This follows because each bias block is a scalar. Finite-t gate reversal therefore requires the interface, not a rigid isolated-layer energy shift.

Unified physical statement: forming an interface can reverse a specified layer's response relative to its same-orientation isolated reference. R=C/C_ref is invariant to a common overall spin sign convention. Keep laboratory-coordinate raw C values. A projected alpha0 Hamiltonian is not a fresh self-consistent monolayer.

## Calculations and verification

- 111 local KM evaluations archived in seven CSV files (pilot, bounded mass44/46/48 screen, main AB, independent BA, decoupled controls). Existing run_fixed_stack_layer_cxx.py and AuditBilayer were used without physics-code changes.
- Main AB: fraction0, mass45meV, grid81x80x120, step0.03125T; BA fraction1, same setup. Pilots grid41x45x72, step0.0625T; screen grid23x24x36. All evaluated outcomes, including nonmatching gate trends, retained.
- Hamiltonian mirror at three generic k points: max~5.3e−18eV. Matching integrated layer-mirror residual~1.22e−15. Decoupled bias variation~5.8e−16. Minimum common main-series gap~3.93meV. See AB_BA_symmetry.json and revision_checks.json.
- No new DFT, no Japanese login-node or compute-node use, no remote job submitted.

## Material source / figure clarifications

Main2b adds same-mesh independent monolayer reference11.3148 in 1e−9 units (oddE,oddB from four local result.json files). This replaces using the more precise-stage13..14 range as an exact plot baseline; the sign is unchanged.
Main2c groups occupied bands1–36 and37–48. Main2d displays actual Gamma/K levels21–56; bands1–20 remain in the response sum. Orbital S5 uses all archived bands39–56 at four k points, two separations, both-layer PAW sphere fractions. No fabricated continuous fatband or interpolated DFT path. Raw levels come from the positive prepared Zeeman field; captions state this.
Correction to previous theoretical note: vc_by_valence uses Q=1−Pocc, including the unoccupied complement; only vc_pairs is finite explicit empty-band pairing. The band-group sum differs from direct density response by less than0.3%. Energy-ranked groups are not asserted to be unique orbital channels or globally isolated subspaces.

## Writing and literature

Research question precedes the model; each result paragraph states its physical implication. Eight figures have distinct captions. Five-column table repaired by using LaTeX mid rather than raw pipe inside stacking operators. Previously accepted Hamiltonian definitions retained. Original-paper missing diagonal term acknowledged once in the theoretical note, not used as a main-text controversy.
Re-read model/material transition passages in band-edge mixing, layer disorder and LEE. Read the user's PRMaterials7,074001(2023) Fig.4 and explanation from Zotero. Source paper PDF/text render stays private, excluded from repository snapshot.

## QA and next step

Three main / five supplementary figures, PNG320dpi and editable PDF/SVG; source and data archived. Main changes and new supplemental figures visually inspected; unchanged environmental panels retain prior QA. Static preflight12pass,2draft-exportwarnings(noTIFF,320ratherthan600dpi),0fail. Markdown image links, math delimiters, table columns and source bounds checked. Original annotated file unchanged: SHA256 a14a7cd5397acde40a5b62a87be06a8487c5cc518c361937f39eb3dc5f3b7ff5.

Next proposed work: disentangle energy denominators and field-dependent spin/electric matrix elements in AB, then test the same decomposition in MoSSe. Search for stable band-subspace/energy alignment conditions and predict an unused point. Do not equate PDOS weight with response sign or claim PRL novelty established. No next-round scan started.

Saved in Obsidian with first-round backup; local Git snapshot only, no push in this round.
