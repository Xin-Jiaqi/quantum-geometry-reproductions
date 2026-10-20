# Scientific decision: eight-band hopping-only BSF × IGMC toy

**Evidence level:** symmetry plus two independent electronic calculation routes, **model only**.

| Check | AF-like (`m_z` partner) | FM-like (`m'_z` partner) |
| --- | --- | --- |
| Input: layer-odd onsite field | exactly zero, rejected by API | exactly zero, rejected by API |
| Dxx/ω at `u=+0.19`, 31² | +0.003085332739 | +0.019375728336 |
| Dxx/ω at `u=-0.19`, 31² | −0.003085332739 | +0.019375728336 |
| Indirect gap, 31² (eV) | 0.2339006858 | 0.0910411766 |
| 〈τz〉 at `u=+0.19`, 31² | −0.1136127099 | −0.1035051697 |
| 〈τz〉 at `u=-0.19`, 31² | +0.1136127099 | +0.1035051697 |
| Wilson finite B, 21² | +0.0030402337 | +0.0194063714 |
| Single-state symmetry | `m'_x` | `m'_x` |

The Wilson finite-mesh discrepancy relative to converged geometric response is about 1.46% (AF) and 0.16% (FM) at 21². Wilson 31² for AF was separately calculated as ≈0.00306465 (0.67% residual), and Wilson 41² as ≈0.00307352 (0.38% residual). Treat this as controlled finite-link convergence, not exact equality at finite N.

**Mechanism beyond a fitted polarization potential:** A nonzero **layer charge imbalance** and nonzero insulating spin-Zeeman geometrical response arise despite not applying a hand-picked odd-u onsite term. Their ±u parity is dictated by the verified partner symmetry, but their amplitudes reflect explicitly declared microscopic two-center hopping and spin couplings. AF-like spin-dependent interlayer tunneling suppression tests show >50× lower response without gap closure.

**What cannot be claimed:** stable BSF, new material, large experimental current, full spin+orbital+Faraday response, or any quantitative AgF₂ / Fe₂MoSe₄ result. Moreover, the electronic occupied-band sum falls monotonically in the sampled positive-u window; the selected `±0.19` are **not proven minima**. This is a major negative scientific gate, not a detail to hide.

**Go/no-go decision:** sufficient to store as a **separate research-stage model** and justify a model-driven DFT structure/magnetism feasibility study; insufficient for a publication claim about switchable ferroelectric control. A convincing follow-on result needs a physically relaxed bistable bilayer or other demonstrable nonvolatile structure pair, consistent magnetic symmetry, and a validated spin/position Wannier operator basis.
