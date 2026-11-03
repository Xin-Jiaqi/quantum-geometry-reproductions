# BSF IGMC GO / PARTIAL / NO-GO after SOC screen + Eq. (14) source audit

**2026-10-23 ｜ Decision: PARTIAL (GO to the next low-compute validation stage; NO new VASP/PBS/Wannier compute authorized by this audit).**

## Three-line summary

1. The SOC-on electronic screen **passes** its magnetic-symmetry gate: Ag Type-I moments (±0.568 μB along x, net 0.038 μB) and the m'_x moment mapping hold at printed precision (max spin pair residual 0.001 μB, orbital 0.000), on a structure whose {m_x|(1/2,1/2,0)} glide is exact to 1.28e-5 Å; gap 1.3616 eV (screening grid).
2. The Eq. (14) convention conflict is **resolved to the channel-resolved protocol (η_C = +1, η_D = −1)** using the paper's own SM Tables I/II (transcription) and its own model results (T-invariant Rashba σ^C_xy ≠ 0, σ^D ≡ 0; TR-breaking Dirac σ^D flips under TR partner); the literal main-text sentence is documented as shorthand, not silently fixed. Under m'_x, **D_xx (insulating displacement) is allowed under BOTH protocols**; C_xy allowedness is protocol-sensitive (a falsifiable discriminator).
3. But **no material response operator chain exists**: no Wannier spinor H(k), no same-gauge spin/Zeeman or covariant velocity matrices for AgF₂ (SOC run saved no WAVECAR), no quantitative σ value, no lower bound. So the project is PARTIAL, not GO-for-costly-work.

## Prioritized science story (rank order)

- **Primary:** Type-I B1 → B2 sliding switches a finite insulating **D_xx** (protocol-robust allowed under m'_x) while the weak out-of-plane axial response is even — connector-selective displacement IGMC switching.
- **Alternative (strained Type-II):** D_xx even, M_z/AHE odd — the parity fingerprint flips with the connector class.
- **Fail path:** the quantitative D slope vanishes once a real response calculation is done, or the source operator chain proves unbuildable (unit/gauge mismatch) — then demote AgF₂ and spend nothing further.

## Real blocker (single)

A source-validated, unit-tested **Wannier-operator response chain** (spinor H(k), same-gauge S_b with g·μB, covariant v, orbital −B·L, multiband Wilson cross-check) does not exist and cannot be built from the current outputs (no WAVECAR saved; `wannier90_hr.dat` alone insufficient). Additionally the pinned SM PDF is unreachable, capping SM evidence at INDEPENDENT transcription (1 residual on `-42'm'`).

## Low-compute validation recommendation (no new DFT)

1. Matched-model falsification: finite-field Wilson polarization derivative vs static low-ω σ^D on the same model after fixing units — the single decisive free check.
2. Re-verify the `-42'm'` SM row (1/116 residual) against a fresh SM page image when the pinned SM is available.
3. Re-run the protocol-vs-SM projector comparison in CI on every push.

## Boundary statement

No VASP, qsub, PBS, B2 optimization, MAE, Wannier execution, or material IGMC integration was performed or is authorized by this audit; the single SOC-on WSL run remains the only DFT compute in this thread. The resource-cap breach of that run (4.66 h vs 4 h) stays documented in RUN_RECEIPT.json and is not normalized.
