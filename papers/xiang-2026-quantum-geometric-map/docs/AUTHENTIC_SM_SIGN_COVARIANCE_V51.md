# v51 | Original SM Eq20→31 sign and spin-basis covariance

**Source:** Original user-provided published 8-page PRB 113 L201406 and original 16-page SM, read directly. Main PDF SHA256 b87e0970f25dd45c5048f8cbbb8691b53ea215040935e9b21a4cde80c69b2730; SM SHA256 9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4. Source PDFs remain private; original-paper science **PARTIAL**.

## SM p5 Eq20 vs p8 Eq31

SM Eq20 states rho_mn^(EB;5) = -tau sigma_mn^c * d_b(f_n-f_m) * E^b Bbar^c/(eps_m-eps_n). SM below Eq31 explicitly prints v_nm^a = i(eps_n-eps_m) r_nm^a and Z_nm^ac = i(r_nm^a sigma_mn^c - r_mn^a sigma_nm^c), with Bbar=-g muB B.

**Literal substitution in j_a = sum rho_mn v_nm:**

\[
j_a^{(EB;5)}/(E^b\bar B^c)=+\tau\sum_{nm}Z_{nm}^{ac}\partial_b f_n.
\]

**Printed Eq31:** MINUS tau times that very sum. Renaming the band indices and converting Bbar to physical B flip neither *relative* sign. The independently solved original two-frequency Liouville SM Eq4/9 in v50 also approaches the Eq20 branch, not printed Eq31. 120 **random Hermitian two-band** velocity and spin matrices independently confirm this sign identity at machine precision, with arbitrary band U(1) gauge phases.

## Main Eq10 versus SM Eq62 helicity

\[
H_{\rm main}=\beta k_x I+k_x\sigma_y-k_y\sigma_x,\quad
H_{\rm SM}=\beta k_x I+k_y\sigma_x-k_x\sigma_y.
\]

They are related by \(U=\sigma_z\): H_SM=U H_main U^dagger. But also U sigma_x U^dagger = -sigma_x, and U transforms both current operators. **A physically consistent spin-basis rotation must transform H, current, and Zeeman source together; it leaves both responses unchanged and cannot fix their Eq20/31 relative sign.** Holding sigma_x fixed after switching helicity changes the physical spin-field coupling and changes both responses' signs.

Numerical results beta=1/8, mu=.3, T=.03, 80 radial x112 angular Gauss mesh, UV=2.5, per E*Bbar (in e=hbar=vF=tau=gmuB=1 units):

| Hamiltonian and spin source | Direct Eq20 symmetric | Printed Eq31 symmetric |
| --- | ---: | ---: |
| Main Eq10, +sigma_x | +0.002506163610 | -0.002506163610 |
| Main Eq10, -sigma_x | -0.002506163610 | +0.002506163610 |
| SM Eq62, +sigma_x | -0.002506163610 | +0.002506163610 |
| SM Eq62, -sigma_x | +0.002506163610 | -0.002506163610 |

All physical-B coefficients would carry the **same extra minus sign**. The properly basis-paired cases (main,+sigma_x) and (SM,-sigma_x) are identical. Prior apparently equal signed main Eq11 outputs are achieved through *different* formula and magnetic-source branches; numerical agreement alone does not prove the entire source sign chain.

## SM Eq33 independent verification

SM p8 Eq33 uses sum_{m!=n} Z_nm^ac = d_a <n|sigma_c|n>. Independently numerical-differentiating eigenstate spin expectation confirms the local Zeeman curvature identity to about **8e-11**, across both k directions, both helicities and both spin signs. Separately checked Eq33 **integrated xy/yx symmetry**, up to bounded finite quadrature error. These checks support Eq33 itself and rule it out as an explanation of the Eq20/31 sign mismatch.

## Acceptance boundary

**12 new source-level tests PASS locally**: generic random Hermitian index algebra; independent U(1) phases; full spin-basis SU2 covariance; all four source/hamiltonian combinations; doping and tilt sweeps; Eq33 pointwise and integrated checks; v50 multi-frequency source sign consistency; invalid inputs. Overall paper remains **PARTIAL**: original Eq20/31 contradiction, published absolute Fig2(b) SI ordinate (~2pi audit pending), actual bath/contacts, and orbital higher-band/three-band effects remain OPEN. This finding is a **local original-formula algebra discrepancy**, not a claim that the experimental response is zero or that the full paper is invalid.

Reproduce from isolated PRB paper module:

    python -m unittest tests.test_sm_sign_covariance_v51 -v
    python scripts/reproduce_sm_sign_covariance_v51.py --beta .125 --mu .3

The CLI optionally SHA256-verifies both original private PDFs. Data: results/frozen/v51_sm_sign_covariance_summary.json. Next v52 should test published Fig2(b) absolute SI conversion using original figure units and all physical constants, without fitting 2pi or flipping signed source branches. Other paper modules unchanged.
