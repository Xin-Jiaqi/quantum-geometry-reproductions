# v64 | Independent virtual voltage probe and ordered limits

**Scope:** VERIFIED model-level zero-positive-coupling probe formula, CONDITIONAL finite QWZ coherent Hall reading, UNRESOLVED bulk thermodynamic Kubo/NEGF correspondence and Xiang original orbital mixed E/B. Xiang source-paper audit remains **PARTIAL**.

## Mathematical observable and independent method

The existing v63 physical voltage-probe calculation scales only top and bottom lead-device hopping amplitudes by lambda>0, with Sigma_p,Gamma_p proportional to lambda², leaving drive leads unchanged.

A fully detached probe at lambda=0 has NO measured voltage (floating equations become singular). We instead compute the **lambda→0+ limiting reading** from the unperturbed *two-drive-terminal* Green function, without either top or bottom probe selfenergy:

    G0^r = [E-H_QWZ-Sigma_L-Sigma_R]^{-1}
    w_pL = Tr[Gamma_p G0^r Gamma_L G0^a]
    w_pR = Tr[Gamma_p G0^r Gamma_R G0^a]
    V_p(0+) = (w_pL V_L + w_pR V_R)/(w_pL+w_pR)
    V_L=+1/2, V_R=-1/2
    I_L(0+) = Tr[Gamma_L G0^r Gamma_R G0^a]*(V_L-V_R)
    Rxy(0+) = (V_top(0+)-V_bottom(0+))/I_L(0+).

Here Gamma_p is the *unscaled* imaginary selfenergy of a physical semi-infinite QWZ strip; finite probe-to-drive transmissions scale as lambda²*w_pL/R, and probe-to-probe terms as lambda^4. The limit is reported only if both sensor spectral weights and the two-terminal current clear explicit thresholds. Units of Rxy are h/q²; drive current is q²/h for unit bias. This is not bulk conductivity or the source's orbital E×B response.

## Main bounded numerical observations

QWZ m=-1, mu=0, eta=1e-4, left/right and proposed T/B probe strip widths=3:

| L | virtual Rxy(0+) [h/q²] | I_L(0+) [q²/h] |
| --- | ---: | ---: |
| 5 | 1.000983988533 | 0.998034885976 |
| 7 | 1.002157460359 | 0.995699002514 |
| 9 | 1.002728517058 | 0.994565219087 |
| 11 | 1.002957109982 | 0.994111910474 |

The finite L=5..11 values drift. Contact widths are kept **fixed in lattice units**, not at constant fraction of L. This is NOT a thermodynamic convergence demonstration and these data must NOT be fit to claim one.

At L=7, finite-lambda four-probe Rxy differs from the new two-terminal virtual limit by approximately -3.56e-6 at lambda=.01, -3.57e-8 at .001, and -3.57e-10 at .0001: the predicted O(lambda²) convergence. The maximum probe-induced x-cut redistribution correspondingly decreases from ~7.58e-4 through ~7.58e-6 to ~7.58e-8 (q²/h under unit bias). The virtual limit was **not fitted to the finite-lambda points**.

## Necessary physical negative controls

- Sensor offsets (0,+1,0,-1) versus (0,-1,0,+1) give Rxy(0+) **1.00192157** versus **1.00238115**, even though the two-drive-terminal current is unchanged. A noninvasive sensor is not guaranteed position independent.
- Changing sensor width from 2 to 7 lattice sites changes virtual total spectral weight approximately **8.20 to 86.71**, but in this *particular symmetric chiral sample* changes the Rxy ratio by less than 2e-9. This is NOT a general width universality theorem.
- eta=1e-3, 1e-4, 1e-5 yields Rxy about 1.0021657308, 1.0021574604, 1.0021566349 respectively; these are finite-regulator results, not eta→0 proof.
- Trivial m=3 has current below reporting resolution: Rxy=**None**. Metallic mu=1.5 virtual Rxy ~0.143616; do not apply Rxy=1/sigma_xy in metal where sigma_xx, scattering and geometry matter.
- A virtual reading is not reported if the proposed probe's unscaled spectral transmission is below the declared threshold, even if algebraic ratios could formally be computed.

## Independent scientific acceptance and next steps

**26 new tests** independently check two-drive-only resolvent, four weights and floating potential formula, spectral sum rule, left/right Caroli reciprocity, lambda² finite-probe resistance and internal-current errors, chirality, width and position, eta/size, metallic/trivial controls and fail-closed inputs. Frozen output in results/frozen/v64_virtual_probe_ordered_limits.json is produced byte-for-byte with scripts/reproduce_v64.py.

Combined v60–v64 standalone local suite: 137 tests. Remote GitHub CI is the separate release gate.

Next: design a contact-fraction-preserving size family and explicitly fix the (lambda,eta,L) order of limits. For a **same-observable** Kubo/NEGF comparison, supply longitudinal conductivity, scattering and conserving current vertices. The main unresolved science remains Wannier embedding and complete orbital magnetic E/B response, not more charge Hall topology. Original source-paper state **PARTIAL**.
