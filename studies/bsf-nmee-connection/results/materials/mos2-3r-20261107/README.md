# Non-Janus 3R-MoS2 NMEE control

Status: calculations submitted; no response values available at this checkpoint. See PLAN.md for the physics question, independent symmetry predictions, database-based selection and decision paths. Twelve-state model work is paused.

BiDB 1929 is the isolated layer; BiDB 1930 is the sliding-stable parallel 3R bilayer. B is the exact horizontal mirror partner of A, not a second unrelated chemical orientation. Slabs are centered and share the same supercell height. A source and B partner have equal zero-field energy by symmetry, not a new independent stability calculation.

Fields: conventional gate Fz in V/nm, with VASP EFIELD=-Fz/10. Probe Ex is EFIELD_PEAD and probe Bz is spin-only BEXT. Total/layer coefficients use Pauli spin density, actual area, and fixed ions/neutral insulating filling. Monolayer half-spaces are not treated as layers.

response_initial_four_signs.py is the exact source used by the three initial zero-gate jobs. response.py additionally permits the established nonmagnetic time-reversal reduction for gate points. A missing negative-B check is recorded as null, never zero. POTCAR, WAVECAR, CHGCAR and transport credentials are not distributed. Reproduction requires authorized VASP/PAW access and a compatible compute-node PBS environment.

No new PRL novelty, universal mechanism, experimental switch or quantitative convergence claim at this stage.
