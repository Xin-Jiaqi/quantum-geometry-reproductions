# Research correction: distinguish abstraction level, not just shared phenomenon

**Date:** 2026-10-22  
**Purpose:** Correct the literature-kill mechanism and generate a plausible general-theory direction from specific stacking-response examples. This is a research hypothesis, not a novelty claim.

## 1. The user's conceptual correction

A paper on a concrete realization does not automatically pre-empt a later general theory. The examples raised by the user make the distinction clear:

- A material-specific sliding-ferroelectric result is evidence for a mechanism, not automatically a general theory of bilayer stacking ferroelectricity.
- A particular ferroic order or material realization does not automatically pre-empt a broader symmetry principle.
- Conversely, a broad principle may contain a concrete case as a corollary while leaving the material-specific prediction, magnitude, and experimental realization open.

Therefore **“direct competitor” must be defined relative to a proposition and its abstraction level**, not assigned from shared keywords. Do not infer bibliographic priority from the analogies alone; their role here is to correct the research logic.

## 2. Use two axes, not one linear novelty score

### Axis A — scope of the claimed result

| Scope | Example | What a matching prior paper would pre-empt |
|---|---|---|
| S0: realization | one material, stacking, observable, drive | the same material-specific prediction or measurement |
| S1: mechanism | a specified microscopic mechanism for a response | the same mechanism in the same relevant regime |
| S2: family classification | a class of layer groups / stackings / responses | the same classification and selection rules |
| S3: general principle | a theorem or covariance rule covering a broad family | the same theorem or a result that logically implies it |

### Axis B — type of contribution

- **Phenomenon / material prediction**
- **Microscopic mechanism**
- **Symmetry classification / theorem**
- **Computational method / database / algorithm**
- **Experimental protocol / device proposal**

The axes are independent: an algorithm may cover a broad family without introducing new physics; a single-material result may be scientifically important even if its mechanism is known.

### Claim-relative literature labels

1. **DIRECT-PREEMPTING:** a reliable source has already established the same proposition at the same abstraction level, or a more general result demonstrably implies it.
2. **FOUNDATIONAL / NEAREST PRIOR ART:** shares the framework or provides a lower-level case but does not itself establish the proposed general proposition.
3. **INTRODUCTION CONTEXT:** motivates the question, material, measurement, or technique but does not imply the exact claim.
4. **OPEN OVERLAP:** potentially relevant; compare equations, assumptions, scope and counterexamples before classifying.
5. **UNRELATED TO THIS CLAIM:** topical words overlap but the source, output, regime, and proposition do not.

A lower-level example may be crucial Introduction evidence and a test case for a higher-level theory without pre-empting that theory. A general theorem may pre-empt a purportedly new classification while leaving a new concrete prediction open. Keep these distinctions explicit.

## 3. What this means for our current LEE / stacking-response work

The supplied *Layer Edelstein Effect* notes document a substantial prior framework: ordinary Edelstein tensors, layer-resolved LEE, Type-I/Type-II classifications, and stacking routes for 80 layer groups, with MoSSe/MoTe2/WTe2 examples. Thus neither “layer Edelstein effect exists” nor “classify which stackings allow LEE” can be presented as new.

But this does **not** mean every new stacking-dependent response question is pre-empted. The right question is: what does the published theory classify, and what general statement is not already contained in it?

### Candidate generalization to investigate: a symmetry-covariant response field over stacking configuration space

Let \(\boldsymbol\delta\) denote the relative in-plane registry, and let \(K(\boldsymbol\delta,\omega)\) be a fully specified response tensor for one fixed drive, measured operator, and response protocol. For a unitary space operation \(g\) that maps one registry to another, the response should obey a covariance relation of the form
\[
K(g\!\cdot\!\boldsymbol\delta,\omega)
=D_O(g)\,K(\boldsymbol\delta,\omega)\,D_J(g)^{-1},
\]
where \(D_O\) includes the actual site/layer permutation as well as the operator's spin/orbital transformation. This equation is only a starting point; antiunitary operations require a causal-response/Onsager treatment and cannot be inserted as ordinary unitary matrices without qualification.

At a symmetry-fixed registry \(\boldsymbol\delta_*\), let \(H_*\) be its stabilizer. If \(K\) is differentiable there, the order-\(n\) Taylor coefficient is constrained by
\[
D^nK|_{\boldsymbol\delta_*}
\in\operatorname{Hom}_{H_*}
\left(\operatorname{Sym}^n U^*\otimes V_J,V_O\right),
\]
where \(U\) is the relative-registry representation, including layer exchange, and \(V_J,V_O\) are the source and output representations. The smallest \(n\) for which this space is nonzero is a candidate **symmetry-activation order**. Allowed does not mean microscopically nonzero.

**Important:** this covariance and tensor-product statement alone is standard symmetry mathematics. Renaming it a “generalized Neumann principle” would not be a contribution. The potential research value would have to come from a nontrivial result such as:
- a systematic layer-group / magnetic-layer-group classification of the *leading allowed registry derivative* for specified response channels;
- globally enforced zeros, parity relations, or sign constraints connecting symmetry-related registries across the configuration-space torus;
- a response channel for which layer permutations and antiunitary constraints produce a rule not already present in published LEE / stacking-response theory;
- a source-backed microscopic model showing the rule predicts a nontrivial response texture, beyond merely being symmetry-allowed.

These are hypotheses to test, not assertions that the literature lacks them.

## 4. First small test, not a full-paper commitment

Keep the existing idealized result as a seed:
\[
S_-^z(\boldsymbol\delta,\mathbf j)
=\gamma(u_xj_y-u_yj_x)+O(u^2j)
\]
near an ideal \(C_{3z}+m_x+m_z\) registry, where \(m_z\) exchanges the layers. The result only establishes one symmetry-allowed first derivative. The published minimal LEE Hamiltonian gives zero for this \(S_z\) channel, so it is not yet a viable microscopic test.

Before adding SOC or any free hopping, run two independent gates:

1. **Theory-overlap gate:** inspect the LEE PRL supplemental material for whether its 80-layer-group stacking construction already states the registry-dependent covariance or leading-derivative result. Record the exact equation/table and scope, not just the title.
2. **Hamiltonian gate:** continue the source-backed AA-MoS2 Slater–Koster matrix covariance test, including orbital mirror parity and rotations. Only after this passes should a source-supported SOC extension be considered.

If the general classification reduces to ordinary \(\operatorname{Hom}\)-space enumeration with no new consequence, retain it as a useful tool and move to the next concrete observable. If a nontrivial global constraint emerges, let that result dictate the paper's scope.

## 5. Change to our decision process

For every future idea, record separately:
- exact proposition / predicted observable;
- scope level S0–S3;
- contribution type;
- closest prior source and exact result that overlaps;
- what the prior work does **not** establish;
- smallest falsifying test;
- evidence status (symmetry allowed, model nonzero, material validated, or experimentally measurable).

Do not mark a direction as NO-GO because a related example exists. Mark NO_GO_CLAIM only when the registered proposition itself is pre-empted or disproved. A failed model test should be MODEL_NO_GO for that model, not for the whole physical direction.


## 6. A concrete global consequence worth testing: symmetry-enforced zeros on the registry torus

This is a more useful first target than merely restating the local tensor-product rule.

For two identical layers, the relative in-plane registry lives on the configuration-space torus
\[
\mathcal T=\mathbb R^2/\Lambda,
\]
where \(\Lambda\) is the in-plane Bravais lattice. If the horizontal mirror exchanges the two layers, its action is \(\boldsymbol\delta\mapsto-\boldsymbol\delta\) modulo \(\Lambda\).

For the layer-odd out-of-plane axial spin \(S_-^z=S_T^z-S_B^z\), the ordinary axial-spin factor under \(m_z\) is \(+1\), but the layer permutation contributes \(-1\). An in-plane polar current is \(m_z\)-even. Therefore the response must satisfy
\[
\chi_{S_-^z,\mathbf j_\parallel}(-\boldsymbol\delta)
=-\chi_{S_-^z,\mathbf j_\parallel}(\boldsymbol\delta).
\]

Any registry fixed by this involution, \(\boldsymbol\delta_*=-\boldsymbol\delta_*\pmod\Lambda\), obeys \(2\boldsymbol\delta_*\in\Lambda\). In a primitive 2D Bravais cell there are four such points: \(0\), \(\mathbf a_1/2\), \(\mathbf a_2/2\), and \((\mathbf a_1+\mathbf a_2)/2\). **If each point really preserves the layer-swapping mirror in the chosen structure and the response protocol respects the stated symmetry, the layer-odd \(S_z\) current response is symmetry-forced to vanish at each point.** Around each fixed point, the response is odd in local displacement, so every even Taylor order vanishes; the little group at that particular registry may forbid additional odd orders.

This is a candidate general-theory corollary: it predicts a *set* of response zeros across the registry torus, not just a zero at one reference stacking. It is still an elementary consequence of covariance, and it may already be known in related configuration-space work. Its value depends on whether the multi-registry result, correct layer permutation, and response-specific zero set are absent from the closest prior theory. First validate the symmetry of each registry, then compare the result with the LEE supplemental material and stacking-ferroelectricity configuration-space literature.
