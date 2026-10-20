# v34: SCBA, ladder, spectral Kubo and finite-(q,omega) Ward migration

This migration copies the audited v11/v17/v18/v19 model kernels and their 48 original regression tests from the frozen v30 research tree. These methods are independent periodic two-band regulator audits, not a quantitative fit to the continuum surface or the published paper's PHE.

## Newly callable kernels

- \`qgmap.microscopic_ward_v11.scba\`: self-consistent retarded Born self-energy \(\Sigma^R=W\langle G^R\rangle\), fixed-point residual and causality.
- \`qgmap.microscopic_ward_v11.ladder_vertex\`: independently solves the four-dimensional Pauli-basis noncrossing Bethe-Salpeter vertex for a specified pair of Green functions.
- \`qgmap.warped_disorder_v17.ra_observable\`: finite-eta RA Fermi-surface bubble, with bare and ladder-dressed branches and an exact momentum-shift control.
- \`qgmap.spectral_kubo_v18.spectral_bubble\`: symmetric spectral Kubo-Greenwood combination \(RA+AR-RR-AA\), with a distinct ladder solve in each channel.
- \`qgmap.finite_qomega_ward_v19.finite_ward\`: finite-q, finite-omega matrix resolvent identity using separately converged retarded/advanced self-energies and the exact finite-difference longitudinal insertion \(\Delta H/q\).

## Reproduction commands

Run from this directory after \`python -m pip install -e .\`:

\`\`\`bash
python -m unittest discover -s tests -v
python scripts/reproduce_advanced.py --case scba-ra
python scripts/reproduce_advanced.py --case kubo-spectral
python scripts/reproduce_advanced.py --case finite-qomega-ward
\`\`\`

The original v17–v19 suites contribute 48 tests and are part of the ordinary GitHub CI test discovery.

## Physical interpretation and limits

1. SCBA and the ladder use the same momentum-independent scalar disorder variance \(W\); the vertex is not fitted separately.
2. RA alone is not a full dissipative spectral conductivity. The spectral implementation uses all four \(RA, AR, RR, AA\) channels but is still not complete Kubo–Středa: Hall-sea terms with self-energy/vertex frequency dependence and crossed diagrams are absent.
3. The finite-(q,omega) Ward identity includes \(\omega+2i\eta\). The extra \(2i\eta\) is an external resolvent regulator, not a physical scattering rate or proof of charge conservation at finite \(\eta\).
4. The finite-q longitudinal bare vertex is \([H(k+q/2)-H(k-q/2)]/q\), not the midpoint velocity. Replacing it by \(v(k)\) fails the exact matrix identity at finite q.
5. The periodic lattice completion has extra Dirac points and regulator-specific symmetries. Results are dimensionless BZ averages; \(b\) is a model energy or momentum shift, not Tesla. They do not reproduce the paper's two-color Liouville/RTA dc PHE, and do not settle the order of the physical \(q,\omega,\eta\) limits.

## Next migration gate

v20's exact finite-system real-frequency Kubo continuity/contact-term benchmark and the v20 limit-order scan remain offline. They should be migrated next as a separate protocol: their exact finite-system Kubo identity is independent of the SCBA/ladder resolvent Ward test above.
