# Open physics questions · do not silently resolve

1. **High priority: Eq. (14) antiunitary sign.** Printed Eq. (14) uses `η_T=-1` for *both* C/D. However Table I assigns Z T-even/Q T-odd; Rashba is T invariant and has nonzero C, and SM §V uses distinct eV² vs aeV². Is Eq. (14)'s antiunitary sign intended only for D, or is an AC phase/time-convention implicit? Check directly before any symmetry-based exhaustive classification. Treat `paper_printed_eq14` and `tensor_channel_parity` as **separate named branches**.
2. **Frequency nomenclature and SI normalization.** Fig. 2 caption states `ω=10^12 Hz` although theory uses `cos(ωt)` and ℏ=1. Is `10^12` a cyclic Hz or angular s^-1 numerical value? Expect possible 2π in absolute D, must not fudge.
3. **Very low ω in a clean metal.** Finite ω, η→0, and ω→0 order can affect intraband preparation; main split keeps only intrinsic interband contribution. Distinguish the published clean prescription from disorder-induced gyrotropic response.
4. **Bare momentum quantum distance vs projected QGT.** SM §I explicitly discards `m=n` as gauge-dependent instead of deriving the phase-invariant Fubini-Study distance. Independent geometrical derivation should specify the projection that makes this legitimate.
5. **Induction fields.** SM §VI ideal 2D geometry assumes induced E perpendicular to the sheet and `r^z_nm=0`; actual sample geometry/finite thickness and far-field B amplitude require their own experimental regime. Not a blocker for bare model calculations.
6. **Rashba k=0 degeneracy.** Zeeman geometry diverges locally; response integral finite only with its exact occupation cut and appropriate integration. Do not insert ad hoc `η` to conceal underresolved integrals.

No author contacted. No resolution assumed.

## Round-2 status note (2026-10-19)

- Antiunitary **channel-resolved** sign is now independently VERIFIED for nondegenerate clean models via exact TR partners and two SM point-group examples. **Original printed Eq.(14) remains unreconciled as an author-intent issue.**
- SM §VI's `ℏω≈4×10^-3 eV` supports `f=1THz -> omega=2πf`, but does not prove the author's Fig2 numeric settings.
- **New highest-priority question: Figure 2(d) absolute amplitude.** Unit-restored Eq.(12) gives in-gap `0.017754 A/(mT)` vs visually ~`0.04`, despite C curves and Rashba SI checking correctly. See `ROUND2_TIME_REVERSAL_AND_UNITS.md`; no fabricated factor-of-two cure.

## Round-3 focused update (2026-10-19)

- **CONFIRMED source-local algebra discrepancy:** SM p.6 Eq.(34) uses `Q_SM=2Re(rσ)` in an intermediate cos term, despite main Eq.(5) `Q_main=Re(rσ)` and despite the preceding SM Eq.(34) first line giving an exact coefficient proportional to `Re(rσ)` with the full `ω±` bracket. A literal factor two results; the main Eq.(12) is independently supported by finite-field Wilson phase and orbital Liouville.
- **Fig2(d) raw numerical provenance still OPEN:** the source figure is close to the *literal* ×2 branch, but neither PDF nor SM provides author code or raw plot values. This is not a basis for fitting or changing the validated physical formula. See `ROUND3_INDEPENDENT_DISPLACEMENT_AUDIT.md`.

## Round-4 MPG symmetry gate

- The user's group-theory toolkit 0.18.0 magnetic generators allow a full 122-point-group audit; 244 allowed tensor spaces recomputed with independent Neumann SVD solver. Conduction C supports 58 groups (24 only diagonal, 18 only Hall, 16 both) while displacement D supports 58 (23, 19, 16). 64 groups each prohibited.
- `paper_printed_eq14` **C tensor dimension** differs from channel-resolved convention for 58 groups, not just one Rashba case. Still limited to real in-phase/quadrature coefficients and does not resolve general finite-complex-frequency antiunitary definitions; see `ROUND4_TOOLKIT_SYMMETRY_AUDIT.md`.
- Toolkit vs source monoclinic setting difference corrected with y/z coordinate permutation, not treated as an error.
- Complete 122×2 live API invocation is not yet executed (toolkit unavailable in offline container); pinned source-generator independent implementation and 122 group order cross-check do pass. Fig.2(d) amplitude problem remains the key scientific block.

## Round-5 source table relation (2026-10-20)

- SM Table II (p.9) lists `-42'm'` in a shared row with `2'2'2` and `mm2`, showing `D_xy` and `D_yx` with separate symbols. For `-42'm'` alone, the unprimed improper fourfold forces `D_xy=D_yx`, leaving **one** independent coefficient, as confirmed independently by its `C2xT` generator; this is a *missing equality*, not an incorrect allowed-component position. Full 244-space audit: **243 source subspaces exact**, this one printed 2D super-space differs from the actual 1D space. `docs/ROUND5_SM_TABLE_SOURCE_AUDIT.md`.
- Remains separate from the global Eq. (14) antiunitary sign ambiguity and SM Eq. (34) normalization issue.

## Round-6 (2026-10-20) scientific status

- New independent d-vector/Pauli trace formulation `Z=-∂n`, `2Q/Δ=d·(∂d×e_β)/(2|d|³)`, and a completely different BZ Gauss-Legendre integration support main Eq. (12) prefactor at the 2×10⁻¹¹ integrand level. Historical Wilson-phase/Liouville tests also support main Eq. (12).
- Full complex nonresonant `K(ω)` contains even-real/odd-imag TR-domain quadratures. Printed Eq. (14)'s single real antiunitary sign does not transparently describe these distinct quadratures; this remains a source-convention question, not a proved author error.
- Fig. 2(d) original plotted data/code have not been supplied. It is scientifically unsupported to apply a 2× factor to the independently validated signal solely because the literal SM intermediate matches the approximate figure height.
- `docs/AUTHOR_DISCUSSION_BRIEF.md` holds two focused clarification requests (unsent). The standalone SM Table II `-42'm'` issue is lower priority.


## Round 7 update · metal integration no longer explanatory

The ±1 eV Fig2(d) metallic ordinate has been independently converged by exact Fermi-boundary integration, removing uniform-mesh resolution as a plausible explanation for the factor-of-two source-plot mismatch. `µ=±1 eV` gives `D=−0.06417844743 A/(m T)` on main Eq12. The missing author graph normalization/data remain an unresolved *source-provenance* problem; do not promote the SM Eq34 ×2 diagnostic. See `docs/ROUND7_METALLIC_FERMI_CONTOUR.md`.
