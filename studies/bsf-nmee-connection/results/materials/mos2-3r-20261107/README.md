# Non-Janus 3R-MoS2 NMEE control

Status: unbiased monolayer, gated monolayer controls, and unbiased stacking A are complete. Mirror B and gated bilayers are pending; the smaller-probe check is submitted. See PLAN.md for the physics question, independent symmetry predictions, database-based selection and decision paths. Twelve-state model work is paused.

BiDB 1929 is the isolated layer; BiDB 1930 is the sliding-stable parallel 3R bilayer. B is the exact horizontal mirror partner of A, not a second unrelated chemical orientation. Slabs are centered and share the same supercell height. A source and B partner have equal zero-field energy by symmetry, not a new independent stability calculation.

Fields: conventional gate Fz in V/nm, with VASP EFIELD=-Fz/10. Probe Ex is EFIELD_PEAD and probe Bz is spin-only BEXT. Total/layer coefficients use Pauli spin density, actual area, and fixed ions/neutral insulating filling. Monolayer half-spaces are not treated as layers.

response_initial_four_signs.py is the exact source used by the three initial zero-gate jobs. response.py additionally permits the established nonmagnetic time-reversal reduction for gate points. A missing negative-B check is recorded as null, never zero. POTCAR, WAVECAR, CHGCAR and transport credentials are not distributed. Reproduction requires authorized VASP/PAW access and a compatible compute-node PBS environment.

No new PRL novelty, universal mechanism, experimental switch or quantitative convergence claim at this stage.

## Prediction sealed before gated bilayer results

At 2026-11-07 19:50 Asia/Shanghai, the isolated-layer-addition hypothesis predicts net A responses -0.8221 and -0.2642 at +0.25 and -0.25 V/nm, respectively. Each prediction adds twice the isolated-monolayer field increment to the independently calculated zero-gate A baseline. Both remain negative. See gate_prediction.json for layer predictions and exact inputs. A disagreement would reject this simple addition model, without by itself isolating hopping from self-consistent interface electrostatics.

Metadata correction: finite-PEAD empty-state EIGENVAL entries are placeholders, not physical gaps. Raw records are retained under results/raw. Only seed_gaps_eV in canonical results is a reported gap; profile.py performs this correction and validates the spatial integrals on a PBS compute node. Initial input_hashes.json refers to the setup snapshot in commit fc43d5e, not subsequent analysis edits.
