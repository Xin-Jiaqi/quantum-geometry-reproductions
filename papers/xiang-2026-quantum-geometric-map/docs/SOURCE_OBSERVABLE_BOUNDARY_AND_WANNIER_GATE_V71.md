# v71 | Source-observable bridge: bounded device current, contact Ward and Xiang source spin signs

**Scope and provenance.** Xiang et al., *Quantum Geometric Map of Magnetotransport*, PRB 113 L201406 (2026). This report relies on the repository's previously original-PDF/Supplemental-Material-verified [v49](ORIGINAL_SM_SOURCE_AUDIT_V49.md) and [v53](ORIGINAL_ORBITAL_SM_EQ49_52_56_V53.md) source audits; it does NOT claim fresh independent access to the original copyrighted PDFs. Original source remains **PARTIAL**. Published SM SHA256 is 9b0b8603d13962e100cc059bf4d96b91aa37062360d13decd1cb2eda495332c4.

## 1. Exact finite-device observable separation

For the same explicitly Hermitian four-orbital device plus finite tight-binding leads of v69, define the bounded, **device-only supported** polarization Ptilde_a = P_a direct-sum 0_leads. The current measured by v68–v70 is the **internal device polarization current**:

\[
J_\mathrm{device}(E)=i[H_\mathrm{device}-E P_b,P_a]\oplus0 .
\]

It includes the electric current contact -i[P_b,P_a]. The Heisenberg derivative of exactly that same bounded polarization in the **complete** Hamiltonian satisfies

\[
J_\mathrm{full}(E)=i[H_\mathrm{device+leads}-E\tilde P_b,\tilde P_a]
=J_\mathrm{device}(E)+i[H_\mathrm{contacts},\tilde P_a].
\]

The second term is the **missing lead-contact current**. The bounded projector polarization is not the physical all-space position of lead carriers, and Jfull is not the terminal Landauer flux.

At isothermal equilibrium f(H(E)-mu) commutes with the complete finite Hamiltonian. Hence the finite equilibrium expectation Tr(f(H(E)) Jfull(E)) is exactly zero, as are its E and B derivatives. This does NOT imply zero DC transport current in an infinite biased electrode system.

Numerical independent full-E re-diagonalization versus finite-T Fermi Frechet for four device orbitals and two 22-site leads (B=.23, mu=.13, T=.37) gives:

| Static isothermal E derivative | Dimensionless value |
| --- | ---: |
| Device-only J | +0.002996376694233 |
| Lead-contact contribution | -0.002996376694233 |
| Complete finite-system Jfull | <1e-15 |

Symmetric, Landau and noncommuting orbital frames agree. Different static E and B controls preserve the cancellation. Contact-operator norm ~0.847293.

At finite omega=.29 and eta=.025, three frames give sigma_device=-0.005310644134-0.007389851857i, sigma_contact=+0.035445444107-0.010578457054i, and sigma_full=+0.030134799972-0.017968308910i. This satisfies the exact retarded Lehmann kinematic Ward relation

\[
\sigma_\mathrm{full}=-i(\omega+i\eta)\,\chi_{\tilde P_a,\tilde P_b}
\]

to about 1e-16. The dynamic omega=0, eta>0 regulator response is not the isothermal thermodynamic limit and need not vanish for fixed nonzero eta. The finite bounded-polarization result is therefore an **exact negative bridge**: reproducing the device-only internal current does NOT establish the original bulk spin PHE/orbital MNHE transport tensor.

## 2. Source-specific spin PHE sign branch (previous v49 discrepancy reaffirmed)

The earlier SHA256-pinned SM source audit records literal Eq20:

\[
\rho_{mn}^{(EB;5)}
=-\tau\sigma^c_{mn}(\partial_b f_n-\partial_b f_m)/\varepsilon_{mn},
\quad \varepsilon_{nm}=\varepsilon_n-\varepsilon_m,
\quad v_{nm}=i\varepsilon_{nm}r_{nm}.
\]

With Z_nm=i(r_nm sigma_mn-r_mn sigma_nm), the **literal index contraction** gives PLUS tau sum_nm Z_nm partial_b f_n. The source's printed Eq31 gives MINUS tau times the same quantity. In an unrelated seeded Hermitian 4-band control at tau=1:

| Indexed source branch | Value |
| --- | ---: |
| Direct Eq20 velocity contraction | -0.203967264526674 |
| Independently inferred Z contraction | -0.203967264526674 |
| Printed Eq31 sign | +0.203967264526674 |

Numerical direct-minus-inferred residual ~1.1e-16, and the distinction survives U(1) phases and band relabeling. This is a **literal source-level convention/algebra discrepancy**, not an author-certified erratum or a judgement of measured physics. The actual spin PHE couples to a physical **spin Zeeman operator**, unlike the orbital/Peierls finite device.

Likewise, main Eq10 and SM Eq62 helicities differ at fixed sigma_x. A sigma_z rotation takes main H to SM H but also transforms sigma_x to MINUS sigma_x; it is not a physical equivalence when the spin operator is held fixed. The exact 2x2 control verifies this.

The **orbital** SM Eq49→Eq52 coefficient and explicitly omitted m != ell three-band interference were already independently audited in v53, and remain unresolved in the full published tensor. SM Eq56's two integrated printed lines passed the named integration-by-parts check. This v71 work does NOT falsely equate the separate spin and orbital source branches.

## 3. Fail-closed physical capability inputs

Source-specific inventory gates distinguish finite_device, spin_phe and orbital_mnhe. Spin PHE must supply Bloch Hamiltonian and velocity, spin Zeeman, filling/preparation, the two-color pole convention, magnetic units and tensor definition. Orbital MNHE requires the **full Wannier position and interband connections, orbital magnetic moment and magnetic positional shift**, BZ measure, filling, units and tensor definition. Even a complete input inventory is *not* a solver, matrix validator or proof of bulk response. Incomplete/missing fields are explicitly rejected; no invisible conversion from orbital centers to full magnetic Wannier embedding is allowed.

Scientific classification: **VERIFIED** bounded finite closed-system current/contact Ward and source-index algebra, **CONDITIONAL** finite supplied graph and seed spin test, **UNRESOLVED** source-prescribed bulk spin/orbital response, real magnetic Wannier embedding, original SI panels and Eq49/52/three-band issues. Overall Xiang PRB **PARTIAL**.

## Reproducibility and next priorities

- Code: src/qgmap/source_observable_bridge_v71.py
- 42 named physics tests: tests/test_source_observable_bridge_v71.py
- 1 frozen-data regeneration gate: tests/test_source_observable_bridge_v71_frozen.py
- scripts/reproduce_v71.py and results/frozen/v71_source_observable_bridge.json
- Initial independent local v71 acceptance: **43/43 PASS**.

v72 should prioritize true bulk finite-q source-specific PHE/MNHE with a defined current, full magnetic Wannier embedding and independently checked source signs/factors. Do not rebrand the finite bounded-polarization Heisenberg derivative as DC Hall, or infer author intent from the printed mismatch.
