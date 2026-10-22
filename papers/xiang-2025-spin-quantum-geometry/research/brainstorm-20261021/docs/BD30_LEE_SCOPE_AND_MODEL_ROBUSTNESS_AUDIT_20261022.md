# BD30 | LEE scope audit + independent rerun + registry-winding robustness gate

**Date:** 2026-10-22  
**Decision:** Do not claim generic LEE, stacking-based LEE selection rules, or a universal finite-radius winding as new. Keep the operator-resolved registry-field idea OPEN, with direct-competitor review still pending.

## 1. What the published LEE work already owns

Primary sources:
- Zhou et al., *Layer Edelstein Effect*, PRL **137**, 076401 (2026), DOI [10.1103/qhvs-pfkq](https://doi.org/10.1103/qhvs-pfkq), [arXiv:2607.17209](https://arxiv.org/abs/2607.17209).
- The user's local full-paper/SM reading note, Layer Edelstein Effect.md, especially the tensor rules and SM Table S1 (80 layer-group stacking routes).

The article explicitly says its model-independent stacking symmetry criterion applies to **all 80 layer groups and their stacking configurations**. The main text derives a second-rank axial spin/current tensor, divides bilayer symmetries into layer-preserving and layer-exchanging operations, and classifies Type-I/Type-II LEE. Its SM Table S1 enumerates initial layer group + stacking operator + resulting bilayer point group and whether the LEE is allowed.

Therefore the following are **already covered / direct prior art**, not novelty claims:
1. Generic current-induced layer-resolved spin magnetization in nonmagnetic bilayers.
2. A general bilayer stacking framework for LEE.
3. Neumann tensor-allowed/forbidden components by bilayer point group.
4. Classification over the 80 layer groups and their stacking routes.
5. Field-induced activation and reversal of layer-opposite spin response.

## 2. Narrower possible gap, with a serious competitor caveat

The public main text and the local SM reading note are organized around **point-group tensor feasibility at a stacking**, plus Type-I/Type-II activation. A targeted text search in the public arXiv manuscript finds no discussion of winding, and the available text does not formulate a global vector field \(\boldsymbol\chi_\pm(\mathbf u)\) over the complete continuous registry torus with operator-resolved zero sets and local indices.

However, the LEE article cites a same-group companion work, explicitly labelled **submitted**, titled *Control of the Rashba-Edelstein effect via sliding ferroelectricity in bilayer systems* (reference 105 in the arXiv version). A public search did not locate a freely accessible full text as of 2026-10-22. This is a **potential direct competitor**, not merely Introduction context. Until its manuscript is obtained or its exact scope is clarified, the continuous-sliding part of the novelty assessment remains **UNRESOLVED**.

Working hypothesis worth testing, not yet a novelty claim:
- For a fixed-current response vector \(\boldsymbol\chi_\pm(\mathbf u)\) of layer-even/odd \(S_z\), the output-operator parity under layer exchange changes the equivariance law over registry space.
- The isotropy group at each special registry can force different zeros in the two fields.
- Near the high-symmetry G registry, a generic nonzero linear term for the layer-odd field gives local index +1; a generic nonzero quadratic term for the layer-even field gives local index −2.
- Standard equivariant-map/index theory already exists. Any paper contribution must be the **new response-specific classification and physically validated consequences**, not “discovering winding” or renaming Neumann's principle.

## 3. Independent reproduction of the current toy model

I reran the exact research-branch script against the user-supplied local BiDB SQLite (PRAGMA quick_check=ok, 11,184 systems) in a read-only connection. No DFT/HPC. The raw database was not uploaded.

Script: [pd2te2_te_p_soc_boltzmann_gate.py](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/pd2te2_te_p_soc_boltzmann_gate.py).

At \(u_x=+0.10\) Å, \(\mu=1.5\) eV, \(T=0.13\) eV and illustrative \(\lambda_{\rm SOC}=0.48\) eV:

| \(N_k\) | layer-odd \(\chi_{-,j_y}\) | layer-even \(\chi_{+,j_x}\) |
|---:|---:|---:|
| 8 | −0.0065826711 | −0.0000623071 |
| 12 | −0.0065247867 | −0.0000857695 |
| 16 | −0.0065113140 | −0.0001019749 |
| 20 | −0.0065073858 | −0.0001133953 |
| 24 | −0.0065060434 | −0.0001201083 |

All are **arbitrary unfitted toy-model units**. C3/vertical mirror/layer-swap Hamiltonian covariance and time reversal reproduce at machine precision. At the six registry representatives, I independently reproduced the same zero/nonzero pattern: layer-odd response vanishes at G/A/B/C/M/N; layer-even vanishes at G/M/N and is nonzero at A/B/C, with three equal-length C3-related vectors. SOC-off and interlayer-hopping-off controls give zero.

## 4. Crucial robustness result: the finite-radius winding is NOT universal

The earlier baseline scan around G used radius \(r=0.125\) Å and \(\mu=1.5\) eV, with 12 angular points. For \(\lambda=0.48\) eV it gives odd winding +1 and even winding −2 at \(N_k=12\) and 16. But a parameter stress test reveals that the **finite-radius even winding can change**, because additional response zeros can cross that chosen contour:

### \(r=0.125\) Å, \(N_k=8\), 12 angular points

| \(\mu\) (eV) | \(\lambda_{\rm SOC}\) (eV) | odd index | even index | minimum even-vector norm |
|---:|---:|---:|---:|---:|
| 1.0 | 0.24 | +1 | −2 | \(1.20\times10^{-4}\) |
| 1.0 | 0.48 | +1 | −2 | \(1.89\times10^{-4}\) |
| 1.0 | 0.72 | +1 | −2 | \(3.19\times10^{-4}\) |
| 1.5 | 0.24 | +1 | −2 | \(2.17\times10^{-5}\) |
| 1.5 | 0.48 | +1 | −2 | \(1.53\times10^{-5}\) |
| 1.5 | 0.72 | +1 | **+4** | \(8.87\times10^{-5}\) |
| 2.0 | 0.24 | +1 | −2 | \(3.98\times10^{-5}\) |
| 2.0 | 0.48 | +1 | −2 | \(7.15\times10^{-5}\) |
| 2.0 | 0.72 | +1 | −2 | \(5.79\times10^{-5}\) |

For the \(\mu=1.5\) eV, \(\lambda=0.72\) eV case, the \(r=0.125\) Å even-field index +4 is reproduced at \(N_k=12\), with minimum sampled norm \(3.33\times10^{-5}\). In contrast, at **smaller circles** the same parameter sets recover the expected local even index −2:

| \(\lambda_{\rm SOC}\) (eV) | radius (Å) | \(N_k\) | odd index | even index | minimum even norm |
|---:|---:|---:|---:|---:|---:|
| 0.48 | 0.02 | 8 | +1 | −2 | \(2.33\times10^{-6}\) |
| 0.48 | 0.02 | 12 | +1 | −2 | \(3.38\times10^{-6}\) |
| 0.48 | 0.05 | 12 | +1 | −2 | \(2.47\times10^{-5}\) |
| 0.72 | 0.02 | 8 | +1 | −2 | \(3.26\times10^{-6}\) |
| 0.72 | 0.05 | 8 | +1 | −2 | \(1.44\times10^{-5}\) |
| 0.72 | 0.125 | 8 | +1 | +4 | \(8.87\times10^{-5}\) |
| 0.72 | 0.125 | 12 | +1 | +4 | \(3.33\times10^{-5}\) |

**Interpretation:** the symmetry-allowed leading local index is conditionally +1/−2 if its leading coefficient is nonzero and the contour is sufficiently local. A finite-radius winding is a property of the *whole enclosed field*, not just the Taylor coefficient at G. The even response is tiny and sensitive to \(\mu\), SOC and mesh, so a universal finite-radius −2 claim is a **NO-GO**. The layer-odd +1 is much more stable in these tests. A non-Abelian/gauge-safe response formulation remains a separate numerical-method requirement near degeneracies.

## 5. Updated scientific decision

**Keep as OPEN, but downgrade the claim:**
- Robust toy result: symmetry-distinct global zero sets at six special registries, with the layer-even channel nonzero at A/B/C; nonzero SOC-dependent layer-odd response requiring interlayer hybridization.
- Plausible local theorem: operator-parity changes the leading equivariant Taylor representation, giving conditional local index +1 vs −2 near G.
- Not established: a universal global winding; an actual Pd2Te2 material response; experimental observability; independent originality against the unpublished/submitted sliding-REE companion.

**Next highest-value actions:**
1. Obtain/check the companion submitted manuscript *Control of the Rashba-Edelstein effect via sliding ferroelectricity in bilayer systems*. This is the highest-priority prior-art gate.
2. Formalize a response-specific local equivariant-index theorem, explicitly separate symmetry-forced zeros from coefficient-dependent indices and accidental zeros.
3. Replace the band-diagonal Boltzmann toy with a gauge-covariant/degeneracy-safe response calculation and material-grounded Pd-d + Te-p Hamiltonian before treating finite-radius winding as physical.
4. Do not start DFT/HPC until the first two gates survive.

The public repository contains only this analysis and reproducible source code, not the user-provided BiDB SQLite. R1 MnPS3 remains CLOSED; no new external compute was started.
