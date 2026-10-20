# V31: Adiabatic finite-eta Laurent expansion in isolated energy clusters

**Paper:** Xiang, Jin, Wang, PRL 135, 146303 (2025). **Scientific acceptance: PARTIAL.**

## Precisely defined physics

For X = d_x f(H), Pi(Y) = sum_a P_a Y P_a, and Q_j(Y) = sum_{a != b} P_a Y P_b/[i(E_a-E_b)]^j, use the paper-inspired adiabatic recursion L1 = eta + i[H,.], L2 = 2 eta + i[H,.] and rho2 = L2^-1 d_x L1^-1 X.

At fixed k, rho2 = A/eta^2 + B/eta + C + O(eta), with

- A = 0.5 Pi d_x(Pi X)
- B = 0.5 Pi d_x(Q1 X) + Q1 d_x(Pi X)
- C = -0.5 Pi d_x(Q2 X) + Q1 d_x(Q1 X) - 2 Q2 d_x(Pi X)

Projectors are constructed onto **entire isolated degenerate energy clusters**, so no singular intra-doublet denominators or arbitrary individual-band U(2) gauges occur.

C is labeled **adiabatic Laurent constant**, never a universal **physical intrinsic tau^0 driven metallic DC response**. Identifying the two requires a distinct, validated collision/steady-state/limit-order protocol.

## Independent checks

- Existing PT 4-band, new global exact PT 6-band **three separate rank-2 doublets**, and nondegenerate 3/4/5 band models.
- Direct finite-eta orbital/band Liouville agrees with A/eta^2+B/eta+C with **O(eta)** matrix residual: at eta=0.08,0.04,0.02,0.01 the PT 6-band residuals are respectively **0.0020713, 0.0010377, 0.0005191, 0.0002596**.
- Additional population-gradient term independently derived from the Fermi derivatives and dH/dk agrees to **2.6e-14** on six bands.
- True intermediate third-energy-cluster matrix contribution is **nonzero**, maximum **0.002034903**. Intra-doublet degeneracy does not count as a third cluster.
- Against independently derived offline V19 model, PT-six-band Tr(O C) = **-0.00199453636063** (original **-0.00199453563048**), with O explicitly a **synthetic PT-even operator, not a real spin magnetization**.
- An explicitly labeled half-population candidate gives **-0.00123268434776**, difference **-0.00076185201287**. 48-point periodic one-dimensional BZ values: full C **-0.02921429403**, candidate **-0.02934600973**, difference **+0.00013171570**. These are not real-material numbers.
- Counterfactual outer rates 1 and 3 change the constant C and also match independent full finite-eta inversion, demonstrating sensitivity to adiabatic protocol. **Source paper outer rate is 2; neither counterfactual is attributed to it.**
- Random fixed orbital unitary covariance; fail-closed near crossings; numerical regression frozen.

## Boundaries

Not a material prediction, general non-Abelian covariant differential geometry interface, scattering theory, or a completed original SM Eq.(2) vs Eq.(6) four-source audit. Paper Fig.1(a) approximately 2x published amplitude, SM Eq.(55) angular reduction, SM Eq.(69)-(70), and metallic physical DC remain open. Do not force numerical agreement.

Usage (from this paper directory): python -m unittest discover -s tests -v; python scripts/reproduce_laurent_v31.py

Frozen outputs are in results/REFERENCE_V31.json. Previous paper benchmarks untouched.