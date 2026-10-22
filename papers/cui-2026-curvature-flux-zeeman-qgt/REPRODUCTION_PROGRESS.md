# Cui Flux — three-track progress checkpoint

**Date: 2026-10-22. Scientific full-paper status: PARTIAL.**
This is a gate-based ledger, not a percent-complete estimate.
“Passed tests” means the listed model protocol is numerically validated,
not that a real experimental observable or a disputed source claim is resolved.

| Workstream | Accepted gate | Still missing |
| --- | --- | --- |
| **1. Paper reproduction** | Eqs.(2)–(15) two-/multi-band Zeeman QGT, local Dirac flux and Wilson loop; Eq.(16) restricted clean interband causal response, fourfold sector expansion; local GHE/KME reciprocity; Appendix C conventional spin current, Hessian and Eq.(27) spin vorticity; 31-group local antiunitary Kubo covariance | Entire source Fig set and Table I six convention-sensitive entries not certified; full Eq.(19) material-normalized KME, genuine many-band degeneracies, contacts, collisions, SI response |
| **2. Physics capabilities** | Reusable nondegenerate Zeeman QGT; spectral Kubo plus independent Liouville; conventional spin-current derivative; finite-frequency antiunitary selection; static Jahn rank projector and all 31 local synthetic microscopic constraints | Conserving metallic electric drive in general Bloch/Wannier gauge, non-Abelian bands, full spin-torque and material-operator interfaces |
| **3. GitHub integration** | Merged PRs #47/#49/#52/#55/#57 plus physics-registry PRs #48/#50/#53/#56/#59; published source+tests+limited scientific-acceptance docs; Round16 PR #60 merged with 18 additional tests and acceptance documentation | Round16 PR #60 **merged**, double-version CI **118/118 PASS**; remaining offline Round9–11 research protocols deliberately not upstreamed as production |

## Verifiable scope counters

- 31/31 **MPG generator/real-spin projective representation + generic-k four-band local Kubo identities** (new).
- 62/62 **even/odd in-plane Jahn candidate tensor-space rank comparisons** (new).
- 7 representative groups independently tested via 16x16 orbital Liouville in Round16.
- 5 representative square-compatible groups with 8x8 BZ selection checks; C3/C6 hexagonal BZ **not validated**.
- Local tests: 100 previous + 18 new = 118; GitHub Actions run 39722106861: Python 3.11 and 3.12 jobs both success; actual logs each show Ran 118 tests.
- Full-paper acceptance still **PARTIAL**, irrespective of these counters.

## Real-material data gate

Historical AutoWann position-block candidates include HgI2, AlAgS2 and Hg3O3,
but no same-Wannier-gauge **physical spin** operator has been verified.
The Japan system search remains incompletely verified. Do not submit new VASP
runs until the nondegenerate physics interface is reviewed; a purposeful
simple-SOC control run is the next option if verified spin matrices remain absent.

## Next independent milestones

1. A crystallographically valid hexagonal reciprocal-space BZ for C3/C6
   with an explicitly periodic Hamiltonian and band symmetry checks.
2. Non-Abelian/Kramers-safe response, including TRIM and PT partners.
3. Integrate conserved finite-q longitudinal electric and gauge-covariant
   position-operator protocols from the earlier local Round9–11 research.
4. Reproduce a **material-specific** response using same-gauge HR, spin
   and position matrices; only then attempt macroscopic SI KME claims.

No global completion percentage is stated because test count, physics
generality, and true-paper observables do not share a defensible denominator.
