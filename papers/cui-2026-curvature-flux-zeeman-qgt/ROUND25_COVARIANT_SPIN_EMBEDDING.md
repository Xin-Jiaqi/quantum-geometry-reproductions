# Round 25 | Covariant static Ward in a k-dependent Bloch basis

**Acceptance before CI: PENDING.** This adds an independent equilibrium-static mathematical and numerical method, not complete Cui et al. KME reproduction.

For smooth Bloch-gauge \(W(k)\), define the anti-Hermitian connection \(A_a=W^\dagger \partial_a W\) and covariant derivative \(D_a X=\partial_a X+[A_a,X]\). The explicit k-dependent Zeeman operator is \(H_B(k)=H(k)-\kappa B_b S_b(k)\). Its **physical covariant current** is \(D_a H_B\), which is B dependent whenever \(D_a S_b\ne0\).

At T>0 the complete equilibrium response is

\[
\chi_{JB}^{ab}(k;0) =
-\kappa\operatorname{Tr}\{(D_a H) Df_H[S_b]\}
-\kappa\operatorname{Tr}\{f(H)D_aS_b\}
=-\kappa\partial_{k_a}\operatorname{Tr}[f(H)S_b].
\]

The second contribution is the **explicit Zeeman current-vertex contact**, absent from a population-only response. \(Df_H\) is the full Fréchet divided difference of the Fermi matrix, including off-diagonal entries at exact Kramers degeneracies from verified Round24; no band flattening.

## Independent method controls

- Construct a smooth and exactly reciprocal-periodic unitary Bloch rotation of the Round24 Kramers-compatible synthetic four-band H and physical Pauli spin.
- Verify \(D_a \tilde H=W^\dagger (\partial_a H)W\) and, for purely gauge-transformed constant spin, \(D_a \tilde S=0\). A naive raw \(\partial_a \tilde H\) need not equal the physical velocity vertex.
- Independently insert finite uniform B in **both the Hamiltonian and the current operator**; compare to analytic response.
- Independently finite-difference the equilibrium spin density in k, including explicit S(k).
- Add a controlled k-dependent **effective projected spin** test vertex to make the contact nonzero. This is not an actual Wannier/PAW spin operator; it deliberately checks the general algebra.
- Test exact Kramers doublets, arbitrary fixed global gauges, reciprocal periodicity, T sensitivity and invalid operator/connection inputs.

## Scope

Only a static equilibrium Ward and a synthetic projected-spin/Bloch-gauge test. The complete periodic-BZ integral of the full derivative vanishes, so this does not generate DC equilibrium bulk transport. No finite-frequency, conserving collisions, Peierls phase at finite q, diamagnetic/transverse Maxwell contacts, realistic orbital centers or SI KME.

The full original paper remains **PARTIAL**. No source published PDF or private Wannier files are included. Prior 203 accepted tests remain unchanged. Target: **213 tests total** (10 new), to be counted as verified only after Python 3.11/3.12 remote CI and runnable audit.


## Independent initial CI acceptance and frozen numerical reference

GitHub Actions run **39877965441** (Python 3.11/3.12) passed **213/213** original + 10 new tests and executed the numerical audit. These actual CI values are frozen in `results/frozen/round25_covariant_spin_ward_audit.json` with one additional frozen-reference regression test (final target **214**).

At (k=(0.36,-0.29)), (T=0.22,mu=0.15,kappa=0.7), synthetic (k)-dependent projected spin amplitude 0.16:

- Maximum basis-gauge invariance discrepancy: **2.37e-16**;
- Complete Ward vs independent finite uniform B: **2.65e-11**;
- Complete Ward vs independent finite k spin-density gradient: **4.16e-11**;
- Maximum nonzero explicit spin-vertex contact: **0.0421665**;
- Omitting that contact while retaining correct physical covariant velocity gives an erroneous response by **0.0421665**;
- Raw basis derivative vs physical covariant velocity difference: **0.425737** (operator maximum norm).

When spin variation is exactly zero, the effective spin covariant derivative/contact vanishes despite the same nontrivial k-dependent Bloch basis. At exact Kramers point k=0, complete Fréchet theory also remains regular. All figures are natural-unit synthetic-model results, NOT a material/Wannier contact coefficient.

The frozen test adds **1** more test and must run successfully in both remote Python versions before final merge; full paper remains **PARTIAL**.
