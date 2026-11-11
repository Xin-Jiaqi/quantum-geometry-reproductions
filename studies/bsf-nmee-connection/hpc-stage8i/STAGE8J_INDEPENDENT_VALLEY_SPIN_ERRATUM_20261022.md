# Stage 8j | Independent scientific erratum for Stage8i SOC evidence (2026-10-22)

**Original external DeepSeek scientific report from commit 6ffef840ca5ca5045e9073dbf781c20164b1a071 is retained verbatim as historical evidence. Its valley and spin interpretations are partly incorrect, as independently shown below.** This work examines committed EIGENVAL and WAVECAR-derived JSON, not the original remote VASP binary outputs. No new VASP/PBS jobs and no modifications to frozen published PRL reproduction.

## 1. Retain the sound, but conditional, computational evidence

Two noncollinear WSe2 spinor SOC SCFs ran on Japanese IMR PBS, with incomplete MPI shutdown (both Exit_status=255; PROCAR contained a 21-byte header only). External-agent evidence reports completed electronic iteration and usable EIGENVAL/WAVECAR data, although neither run has a normal complete OUTCAR termination. The two 441-point sampled indirect gaps are **1.088127 eV (AB)** and **1.148975 eV (half-slide)**; min sampled direct gaps **1.152668, 1.212372 eV**, with no claim of continuously positive gap at all k.

Recorded SOC energies are E_AB=−46.07647873 and E_mid=−46.05232372 eV / two-formula-unit cell, ΔE=**+24.15501 meV/cell**, near the non-SOC ΔE=+24.11712 meV/cell. HOWEVER both SOC and sampling mesh changed (15×15 to 21×21), so their 0.03789 meV difference is **not solely the SOC correction**.

## 2. High-impact reciprocal-space label error: the reported Kp is NOT K′

Original source hexagonal direct vectors: a1=(a,0), a2=(−a/2,√3 a/2), with a=3.319006624 Å. The authors' C3 reciprocal fraction operation (u,v)→(−u−v,u) fixes the corner K=(1/3,1/3) modulo integers. Time reversal sends that true K to its opposing **K′=(−1/3,−1/3)≡(2/3,2/3)**.

The existing DeepSeek evidence mislabeled (−1/3,+1/3) as Kp. This point lies on the **Γ→M line at 2/3 of M=(−1/2,+1/2)**: |k_false|=0.72865150 Å⁻¹, |k_trueK|=1.26206142 Å⁻¹, |k_M|=1.09297725 Å⁻¹. It is NOT a corner or time-reversal partner of K. Its different valence band energy (~1.3 eV below K) is simply a dispersion comparison between distinct k points, **NOT valley splitting**. The archive separately asserts it already tested bandwise **E(K)=E(−K)**, so labeling a 1.3-eV K-K′ asymmetry a real nonmagnetic valley effect contradicts its own time-reversal check.

Required correction: **delete from scientific conclusions** any interpretation that C3v permits energy-inequivalent time-reversal K and K′ or that −K belongs to the K family. Preserve the raw energy sample and rename its location Γ–M (2/3 M). Test actual −K against K from the existing EIGENVAL rows on HPC. No SCF rerun is needed.

## 3. 87.5 and 35.3 meV are SAME-spin neighboring band separations, not identified Ising splittings

From the archived K-point WAVECAR pseudo-spin expectations:
| Stack | E_band36−E_band35 | Sz band35 | Sz band36 |
|---|---:|---:|---:|
| AB | 87.54 meV | −0.49961 | −0.49788 |
| Half-slide | 35.31 meV | −0.49872 | −0.49886 |

**Both band 35 and 36 have nearly the same spin direction.** Their energy difference is not an established splitting between opposite-spin partner branches. Rename as **energy separation between adjacent K-point occupied valence bands**, pending orbital-matched and opposite-spin state pairing. Other K-point bands 33/34 have opposite spin orientation, but no PAW/orbital matched pair was proven (PROCAR incomplete). Pure Ising *orientation at K* remains a qualitative proxy; quantitative SOC splitting is **undetermined**.

## 4. Kramers-pair spin texture does not demonstrate broken Sz conservation

At Γ, valence bands 35/36 are exactly degenerate. Within a degenerate Kramers manifold, unitary basis rotations alter the spin expectation on each nominal band. Exact counterexample: H=ε I_2 obeys [H,Sz]=0, while |+x>=(|↑>+|↓>)/√2 has nonzero Sx. Thus finite Sx/Sy of Γ band35 cannot by itself prove **[H,Sz]≠0**.

More relevant existing-data check: evaluate *nondegenerate* generic-k spinors away from high-symmetry (including the former mislabeled Kp, now Γ–M interior) after verifying pseudized WAVECAR spinor normalization, record shape, and PAW limitations. Finite off-axis spin of a reliably nondegenerate state would indicate local Sz nonconservation; it is still NOT proof of finite occupied–empty spin×velocity Kubo matrix elements.

Caution: The original spin-texture JSON shows unreliable WAVECAR raw energy fields, deviations of computed spinor norms from unity, only pseudized plane-wave coefficients, and MPI shutdown abnormality. No quantitative band spin or matrix-element claim should proceed without independently validating the WAVECAR reader and physical spin operator.

## 5. Stop/go without any additional HPC jobs

**Corrected status: SOC_SAMPLED_GAP_PASS__SPIN_MIXING_NOT_YET_ESTABLISHED__OUTPUT_COMPLETION_WARNING.**

Keep original four Stage8i evidence files and commit unmodified as provenance. A read-only follow-up by the local DeepSeek harness should: (i) correct reciprocal K/K′ labels and compare all sorted band energies K versus true −K from EIGENVAL; (ii) track opposite-spin AND orbital-matched pairs before naming any SOC spin splitting; (iii) validate generic nondegenerate-k spinors and PAW/gauge limitations; (iv) inspect WAVECAR record lengths/energies and show the effect of Kramers gauge at Γ; (v) produce corrected sanitized evidence and scientific report without new qsub.

**We do not yet recommend an expensive full spinor Wannier/NMEE run.** If existing wavefunctions cannot establish trustworthy spin-mixing indicators, a narrowly scoped operator-complete Stage7j feasibility prototype may then require separate scientific and compute approval. The original spin-only equilibrium/static Gamma is NOT the published driven-DC NMEE; actual WSe2 material gamma remains uncomputed.
