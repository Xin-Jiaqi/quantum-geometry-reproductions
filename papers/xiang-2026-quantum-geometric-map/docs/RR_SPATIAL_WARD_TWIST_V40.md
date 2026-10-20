# Independent v40b RR spatial Ward validation by a uniform vector-potential twist

**Scope:** Same synthetic QWZ + scalar white-noise SCBA model used in \`current_vertex_v40\`. This file adds an independent **finite-difference** verification of the RR current ladder without copying the ladder derivation.

For \(H_{\mathbf A}(\mathbf k)=H_0(\mathbf k+\mathbf A)\),

\[
\Sigma^R(\mathbf A)=W\langle G^R_{\mathbf A}\rangle,\qquad
\Gamma^{RR}_a=v_a+\partial_{A_a}\Sigma^R,\qquad
\partial_{A_a}\Sigma^R=W\langle G^R\Gamma_a^{RR}G^R\rangle.
\]

The code separately self-consistently solves \(\Sigma^R(A_a=+h)\) and \(\Sigma^R(A_a=-h)\) on the **same finite momentum mesh** and checks that the resulting numerical derivative agrees with the 4-Pauli RR current Bethe–Salpeter correction. It also independently verifies \(\partial_{A_a}G^R=G^R\Gamma_a^{RR}G^R\) and \(\partial_{A_a}H=v_a\) in both x and y directions.

The *RA* current correction is computed separately with exactly the same \(W\), and can be much larger than the RR spatial correction; the RR **density** vertex from v39 has yet another bare source. The distinction is essential: a large RA antisymmetric bubble does not by itself give a physical Hall conductivity, and \((RA+AR-RR-AA)\) symmetric spectral response does not replace the full Kubo–Středa I/II response.

Run:

\`\`\`bash
python scripts/reproduce_rr_ward_v40.py --mesh 21 --W .12 --energy 1.1 --eta .12
python -m unittest discover -s tests -v
\`\`\`

This method uses e=ħ=1, an independent periodic lattice regulator and external finite η. **It is not a complete charge-conserving dc magnetotransport response, and it does not evaluate the original paper's nonlinear EB conductivity.**
