# Stage 7j: Wannier gauge covariance across momentum and spin source

**Original follow-up research, NOT a Jia et al. PRL 136 226302 (2026) reproduction result.** Protocol: synthetic gapped eight-band sliding Kane–Mele model, \(T=0\) **equilibrium/static** spin–electric Zeeman derivative. Not an actual DFT result, ferroic-stability calculation, or driven-DC NMEE claim.

## 1. General physical theorem: basis covariance must extend across all response parameters

For physical \(H(\mathbf k,\lambda,B)\), let an arbitrary smooth **representation-only** unitary \(U(\mathbf k,\lambda)\) give \(H'=U^\dagger H U\) and \(Y'=U^\dagger Y\). The physical occupied-subspace Wilson overlap between parameter corners \(a,b\) is

\[
\boxed{\mathcal L_{ab}=Y_a'{}^\dagger\,\mathcal B_{ab}\,Y'_b,\qquad
\mathcal B_{ab}=U_a^\dagger U_b.}\tag{1}
\]

The loop \(\prod_\square\det\mathcal L_{ab}\), normalized to unit phase, is independent both of the orbital basis and of an arbitrary independent occupied \(U(N_{\rm occ})\) rotation at each vertex. **Using just \(Y_a'{}^\dagger Y'_b\) is wrong unless the two coordinate sets refer to the same fixed physical basis.** For distinct Wannier subspaces, real cross-parameter basis overlaps must be computed in the common Hilbert space; independent hr.dat models cannot supply them.

With \(\Omega_\alpha=U^\dagger\partial_\alpha U\),

\[
\boxed{D_\alpha H'=\partial_\alpha H'+[\Omega_\alpha,H']
=U^\dagger(\partial_\alpha H)U,\quad
\alpha=k_x,k_y,\lambda.}\tag{2}
\]

For \(\mathcal A_\alpha=i\Omega_\alpha\), this is equivalently \(D_\alpha H'=\partial_\alpha H'-i[\mathcal A_\alpha,H']\). The same connection rule applies to the physical spin-source derivative, not only to velocities. These are mathematical **covariance identities**, not newly discovered transport effects.

## 2. Exact gauge counterfactual, with unchanged physical spectrum

Introduce a *synthetic*, diagonally phase-shifted spinor orbital frame:

\[
U(\mathbf k,\lambda)=
\operatorname{diag}\exp[i(k_xr_{a,x}+k_yr_{a,y}+\lambda\xi_a)].
\]

The arbitrary positions and source slopes are specified in \`bsf_nmee/wannier_gauge_links.py\`; they are deliberately NOT fitted Wannier centers. The same physical \(H\), \(S_{x,y,z}\), magnetic Zeeman coupling and occupied projector are used. At first K plus \(q=(.05,.02)\,\mathrm{nm}^{-1}\), \(\delta k=.006\ \mathrm{nm}^{-1}\), \(\delta\lambda=.00015\) eV, \(B=\pm.25\) T:

| u | t_perp | Original Wilson local static kernel | Correct \(\mathcal B_{ab}\) links | Naive eigenvector coordinate links |
|--:|--:|--:|--:|--:|
| .1 | 5 meV | -0.032808558314 | -0.032808558175 | **+0.20141988734** |
| .4 | 5 meV | -0.086926169511 | -0.086926169141 | **+0.14175483318** |
| .4 | 20 meV | -1.701673546616 | -1.701673546280 | -1.66408557228 |

The first two **spurious sign flips** result exclusively from the missing basis-overlap factors; there is no physical symmetry breaking or band-gap closure. Corrected and original phases agree at \(\sim10^{-10}\) absolute local-kernel units. These are **single-k tests, not integrated BZ or real-material predictions**.

A separate, independently evaluated **k-only** gauge test at \(u=.1,t=5\) meV yields original Wilson \(-.032808558314\), covariant Wilson \(-.032808557990\), and naive Wilson \(-.054691128530\), an artificial response change of ~67%. The independent Kubo derivative is original \(-.032844952934\), correctly transformed \(-.032844952933\), and naive bare-\(\partial_kH'\) \(-.054703355637\).

At this k, naive \(k_x\) velocity matrix error \(\|\partial_{k_x}H'-U^\dagger v_xU\|_{\max}\approx0.07799886\) eV nm; adding the basis commutator reduces the error to \(\sim1.3\cdot10^{-11}\) eV nm. The naive source derivative \(\partial_\lambda H'\) incurs a maximum spin-source matrix-element error \(\approx2.7000\) (Pauli units), corrected to \(\sim9.1\cdot10^{-10}\).

## 3. Why this matters for the scalar–nematic research mechanism

Stages 7e–g validated a symmetry-constrained scalar-cubic AB/BA response, a nematic-linear same/cross interference and exceptionally large cancellation of low/remote virtual-band sectors. All these are constructed from **physical response operators**. If different k or external-field samples are interpreted in unrelated Wannier coordinate gauges without physical overlaps, an entirely artificial mixed Berry curvature may appear, even reversing its sign. This is a logically separate requirement from retaining remote-band **physical** virtual transitions and Feshbach contacts (Stage 7h).

The occupied-bundle mixed curvature is gauge-invariant as a physical observable *only when the Hilbert-space identification across parameter points is correct*. Numerical self-consistency of Wannier eigenvalues does NOT by itself certify this identification.

## 4. DFT/Wannier preflight and fail-closed status

**Preferred Route A: single fixed spinor-Wannier basis.** Supply same-frame \(H(\mathbf k)\), spin \(S_{x,y,z}(\mathbf k)\), position/Berry-connection or physical orbital overlap links, covariant velocity and full-BZ occupied gap. Insert the auxiliary \(\lambda S_i\) and spin-only Zeeman term in that fixed basis, and independently compare Wilson/projector against Kubo with window, mesh and source-step convergence.

**Route B: separately Wannierized source or magnetic field points.** Require physically meaningful \(k,\lambda\) cross-parameter basis overlap matrices and a common Hilbert-space operator/spin normalization. Without them, status **BLOCKED_CROSS_SOURCE_OVERLAP**: do not use bare coordinates to construct mixed Berry curvature. One can take a finite difference of individually gauge-invariant Wilson phases at different \(B\) if each phase itself is correctly computed; direct B-links are not needed for that specific operation, but consistent physical perturbations/gaps are still mandatory.

Both routes require:
- explicit fixed-stack vs domain-mapping symmetries, including C3 and reference layer-exchange constraints, *before* naming allowed tensor components;
- spinor/SOC model, occupation, full-BZ gaps and source-corner gap guards;
- Wannier remote-band window convergence, spin/velocity normalization and 2D sheet units;
- static response separated from driven-DC NMEE and actual switchable-ferroelectric energetics.

**Future physical gate:** select a genuinely existing gapped bilayer structural candidate and attempt a read-only Wannier operator/overlap preflight; only with that package complete proceed to VASP/QE response calculations. No DFT was run in Stage 7j.

## 5. Reproduce and archival rules

- Independent standalone code/data/provenance ZIP delivered as Stage7j research, not deposited in the PRL reproduction.
- Study module: \`bsf_nmee/wannier_gauge_links.py\`, tests: \`tests/test_stage7j_wannier_links.py\`, script: \`scripts/run_stage7j_wannier_gauge.py\`.
- Machine-readable gate: \`results/STAGE7J_WANNIER_GAUGE_GATE.json\`.
- Article-associated research pointer: \`papers/jia-2026-nonlin-magnetoelectric-edelstein/research/README.md\`. This is **navigation-only**; source reproduction, private correspondence, frozen tests and paper data are not modified.

**Scientific status:** mathematically exact gauge theorem plus numerical counterexamples in the synthetic static model. The true solid-state material response and a universal sign theorem remain **NOT ESTABLISHED**.
