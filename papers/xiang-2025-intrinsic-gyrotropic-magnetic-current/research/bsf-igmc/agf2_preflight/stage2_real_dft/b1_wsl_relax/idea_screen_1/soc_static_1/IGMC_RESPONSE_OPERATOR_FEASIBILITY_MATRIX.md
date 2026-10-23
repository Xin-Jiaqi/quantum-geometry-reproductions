# IGMC response-operator feasibility matrix — AgF₂ B1 (SOC screen to Wannier response)

**2026-10-23 ｜ Read-only inventory of `src/igmc_lab` and the SOC screening outputs. Statuses: EXISTS_AND_TESTED / EXISTS_UNVERIFIED / MISSING / NOT_ACCESSIBLE.**

| # | operator / validation needed for a material IGMC calculation | repo artifact | status | evidence / note |
| --- | --- | --- | --- | --- |
| 1 | 2-band Zeeman QGT kernels Z (curvature) and Q (metric), z = i·v·σ/ε | `core.py: pair_kernels` | EXISTS_AND_TESTED | toy 2-band only; non-Abelian multiband version absent |
| 2 | σ^C (Eq. 11), σ^D (Eq. 12), Kubo cross-check (Eq. 8) | `core.py: kernel_response` | EXISTS_AND_TESTED | ran this session (Rashba/Dirac) |
| 3 | Eq. (16) Rashba analytic + Eq. (17) Dirac lattice + Eq. (13) Fermi-surface route | `core.py: rashba_analytic/rashba_integrate, dirac_response, fermi_surface_response` | EXISTS_AND_TESTED | |
| 4 | finite-field Wilson-loop polarization derivative (insulating) | `independent_checks.py: finite_field_wilson_derivative` | EXISTS_AND_TESTED | ran: 0.01335 (toy units, 2π-normalized) |
| 5 | finite-ω orbital Liouville Kubo kernel | `independent_checks.py: orbital_liouville_kernel` | EXISTS_AND_TESTED | ran: Re ≈ 1e-16 (in-gap), Im = 2.75e-5 |
| 6 | MPG tensor symmetry protocols (channel-resolved / paper-printed) + SM Tables I/II transcription + 122-group toolkit snapshot | `tensor_symmetry.py`, `paper_table_source.py`, `toolkit_bridge.py` | EXISTS_AND_TESTED | this session: channel_resolved 115/116 vs SM transcription; paper_printed 75/116 (41 mismatches) |
| 7 | physical SI conversion of σ (A/(m·T)) | `core.py: physical_response_si` | EXISTS_UNVERIFIED | docstring: "Physical SI current conversion unverified"; g-factor documented as input |
| 8 | SOC spinor Wannier H(k) for AgF₂ B1 | — | MISSING | SOC run saved no WAVECAR/CHGCAR (LWAVE=.FALSE., LCHARG=.FALSE.); no Wannier run |
| 9 | same-gauge spin/Zeeman matrix S_b with g·μB in the Wannier basis | toy Pauli SX/SY/SZ only | MISSING | `wannier90_hr.dat` alone insufficient (SOURCE_MAP); gauge-consistent spin operator required |
| 10 | covariant velocity / Berry position matrix for the material (8-band, Wannier-embedded, degeneracy-safe) | 2-band analytic v only | MISSING | |
| 11 | multiband non-Abelian Wilson polarization cross-check on a matched model | 2-band toy only | MISSING | the review's proposed independent cross-check |
| 12 | orbital magnetic-field coupling (gauge-covariant orbital Zeeman −B·L) | LORBMOM diagnostic only | MISSING | LORBMOM = PAW-sphere projection, NOT the orbital-Zeeman operator |
| 13 | k-convergence + degeneracy handling for non-Abelian bands | `gap_tol` raise-only | MISSING | no symmetry-adaptive degeneracy treatment |
| 14 | AgF₂ SOC electronic structure for operator construction | SOC single point (this PR) | EXISTS (screen only) | no WAVECAR; moments/gap usable as screening evidence only |

**Shortest implementation route (no new DFT):**
1. Build the missing operator chain on a **matched tight-binding/Wannier-style spinor model** of the B1 band structure (same-gauge spin matrix S_b = ∂H/∂b with g·μB, covariant velocity v = (1/ħ)∇_kH, interband position) — all in the existing 2-band kernel framework first.
2. Validate units and gauge: reproduce σ^C/σ^D on the toy Rashba/Dirac in the new same-gauge code, matching `core.py` numbers; fix SI conversion explicitly.
3. Only then, with fresh authorization, connect a Wannier spinor H(k) (requires a new DFT/Wannier run with WAVECAR retained — currently not authorized).

**One falsification test for the missing chain:** the finite-field **non-Abelian Wilson polarization derivative** (already implemented for the toy) must agree with the static low-ω displacement kernel σ^D on the **same** matched model, after frequency factors and field units are fixed. If they disagree, stop and re-audit the operator gauge before any material work.
