# Spin QG × 堆叠对称性：A/B/C/D 并行 Brainstorm 第一轮

**日期：2026-10-21。性质：独立研究预研，不是原创性确认或论文结论。**

## 总体判断

不需要过早排除某条思路，但也不应把「有非零响应」「能写成几何张量」「有可切换符号」视为已经有原创物理。保留四条可证伪分支，每条只给一个首要问题和低成本判据：

| 分支 | 最重要的科学问题 | 最低成本证伪 | 当前状态 |
|---|---|---|---|
| A：隐藏非线性 layer-spin response | SOC 下有限空间波矢和边界是否将隐藏二阶层响应转变为有独立可读出的 spin density？ | 计算完整 `spin density + flux + torque` 连续性方程；先做可解的边界/有限 q 模型，而非再画 bulk spin-current 曲线 | 暂存可重启；与非线性 pure spin Hall / layer Hall 先例重叠 |
| B：真实滑移与 spin-related quantum geometry | 除传统 `k–δ` mixed Berry/metric 外，是否存在由层间**自旋不守恒与可测量力矩**定义的、不能由电荷理论平凡替换而来的量？ | 区分相位旋转与真实跃迁变化；关闭层间 hopping、检查 gauge/嵌入，并定义真正测量对象 | 本轮完成 toy-model null test |
| C：磁性双层 AC 层间转移 | 是否有已知 Néel SOT/普通广义力 Kubo **无法解释**的频率、对称性或边界谱关系？ | 重做完整 local+inter transfer；`m_zT` 相位表、Ward 恒等式；对照 Tang 2025 | C2 原命题 NO-GO；保留转向可观测量的假设 |
| D：BSF 连接算符与响应切换 | 同一双层两态之间多个独立 connection operators 对特定**磁性、层分辨、相位响应**是否导致可验证的零响应/互补响应？ | 先从 `D(g)` 对全部响应向量做群作用，找到互相矛盾的通道，再核对两种连接操作是否真的同时存在 | 数学检查可做；常规 Neumann 定理不能自称新颖 |

## B：本轮独立完成的两组计算

定义占据带投影 `P(k,δ)`，混合曲率采用固定约定

\[
F_{k\delta}=i\operatorname{Tr}\left[P\,[\partial_k P,\partial_\delta P]\right].
\]

### B-null：只有层轨道相位旋转

若 `H(k,δ)=U(δ) H(k,0) U†(δ)`，`U(δ)=exp(-i δG)`、`G` 与 `k` 无关，且占据带投影平滑周期，则

\[
\partial_\delta P=-i[G,P],\qquad
F_{k\delta}=\partial_k \operatorname{Tr}(GP),
\]

\[
\boxed{\displaystyle\oint_{BZ}\!dk\,F_{k\delta}=0.}
\]

这是一条对当前模型条件成立的**数学 null test**，不是有关任意物理滑移的普遍否定。若 `U` 只是固定局域轨道基的重相位，上述 `δ` 不是独立的物理结构变化。即使 `F(k,δ)` 局部不为零，整个 BZ 的积分依然为零。若存在非平凡 k-dependent embedding、扭转参数、非周期边界或者奇点，则必须另行审计。

在一个显式 2×2 绝缘玩具模型，取 `d=(0.7+0.2cos k, 0.25 sin k, 0.5+0.3cos k)`，对 `d_x,d_y` 作绕 z 的 `δ` 旋转；解析导数满足 `F=-1/2∂_k(d_z/|d|)`。数值：

- `max |F(k)| = 0.0777173575`，**局部混合曲率确实不为零**；
- `∫dk/(2π)F = 4.74×10^-18`，**整体积分为数值零**；
- 对解析公式的有限差分局部误差约 `1.3×10^-7`（网格 3000 点）。

### B-control：真实跃迁幅度随参数改变

考虑 Rice–Mele 类的 2×2 周期模型（`δ` 代表一个**人为设定的周期结构参数**，不等于某个材料真实的原子滑移路径）：

\[
\mathbf d(k,\delta)=\big[t_0+A\cos\delta+t_2\cos k,\ t_2\sin k,\ m\sin\delta\big]
\]

取 `t0=1.10, A=0.50, t2=0.90, m=0.35`，有真实 `∂δ` hopping/能量变化；在 δ=2.00，`∫dk/(2π) F_{kδ} = -0.232696767`，全周期 `∫dk dδ F/(2π) = -1.0000000000`（使用 3000×801 点）。

**解释边界：**后者是已知 Thouless 几何泵类型的玩具测试，不能据此宣传「新型滑移自旋泵」；它只证明我们用于下一步的区分方法能识别结构参数真正改变 Hamiltonian 的情况。

本轮可复验脚本：`b_geometry_nulltest.py`。

## B ↔ C：一个可检验的互易关系，暂不视为新物理

考虑带隙始终打开的三参数占据带投影 `P(k,δ,θ)`，其中 δ 是真实结构位移、θ 是作用在真实层间 hopping 的相对自旋 twist 源参数。标准 Berry 曲率满足 Bianchi 恒等式：

\[
\partial_\delta F_{k\theta}-\partial_\theta F_{k\delta}+\partial_k F_{\theta\delta}=0.
\]

对周期 BZ 积分后（不穿越带隙关闭或投影不光滑）：

\[
\boxed{\displaystyle\partial_\delta\int_{BZ} F_{k\theta}
=\partial_\theta\int_{BZ} F_{k\delta}.}
\]

含义：**滑移对电场诱导 spin-link generalized force 的调制**，与**相对自旋 twist 对电场诱导 sliding generalized force 的调制**来自同一个扩展参数空间几何。这是已知 Berry 连接微分几何的恒等式，并非新的普适定理。只有物理上明确的 θ 源、可测 sliding force / link transfer，以及对称性导致一侧有真正独有的信号，才值得继续。

在 C2 四带模型基础上对 interlayer hopping 做随 δ 的真实幅度调制，数值检查（29×29 k 网格，中心差分）：

- `∂δ ∫ F(k,θ) = +0.0005184472202`
- `∂θ ∫ F(k,δ) = +0.0005183388835`
- 二者差 `1.08×10^-7`，即约 `0.02%` 的有限差分误差；**说明数值实现与 Bianchi 关系一致，绝不证明新效应**。

脚本 `bc_three_parameter_bianchi.py`。后续可通过减小差分步长/加密 k 网格再做收敛性核对。

## A/C/D：下一道门槛与先例

- **A**：原四带模型中某些普通层差自旋流算符可被化简成 `I₄`，因而均匀绝缘体高阶响应为零。先尝试有限 q 或边界 spin accumulation，前提是 torque/source 连续性必须闭合。需重点排除 **PRB 112, 134421 (2025)** 中的 layer spin Hall，以及 **PRB 113, 125106 (2026)** 和 **PRB 114, L171108 (2026)** 中的非线性纯自旋 Hall。
- **C**：重新运行了 `c2_minimal.py`，`ω=0.2` 时 link/local response 分别为 `+0.00259065168/-0.00257773600`，两者之和与 `-iωχ_S` 一致到数值精度；标准 mixed curvature 又准确重现静态 link response。它是物理一致的现象，但无法凭此独立于 **Tang–Zhang–Cheng, Nat. Commun. 16, 7790 (2025)** 的绝缘体 NSOT。
- **D**：设两个连接操作 g₁,g₂ 都将状态 X 映到同一状态 Y，则对应响应向量必须同时满足 `R_Y=D(g₁)R_X=D(g₂)R_X`。故 `[D(g₁)-D(g₂)]R_X=0`；若某个标量通道对应 `+` 和 `−` 两种符号，它只能为零。这是 **Neumann principle / 稳定子约束的直接推论**；可用来验证磁性 BSF 与层分辨相位响应分类，但不应被称为新的群论定理。

## 本轮查看的直接研究先例（均需论文级别进一步阅读）

1. Zheng et al., *Interlayer Electric Multipoles Induced by In-Plane Field from Quantum Geometric Origins*, **Nano Lett. 24**, 8017 (2024). DOI: https://doi.org/10.1021/acs.nanolett.4c01657 — k–δ 量子几何、层间电多极矩已做。
2. Zhang, Gao, Xiao, *Topological charge pumping in twisted bilayer graphene*, **PRB 101**, 041410(R) (2020). DOI: https://doi.org/10.1103/PhysRevB.101.041410 — 动态滑移 / 几何泵已做。
3. Tang, Zhang, Cheng, *Néel spin-orbit torque in antiferromagnetic quantum spin and anomalous Hall insulators*, **Nat. Commun. 16**, 7790 (2025). DOI: https://doi.org/10.1038/s41467-025-63171-1 — 绝缘体交流电场 NSOT 与共振已做。
4. Wang et al., *Layer Hall and layer spin Hall effects in two-dimensional altermagnets induced by spin-layer coupling*, **PRB 112**, 134421 (2025). DOI: https://doi.org/10.1103/643d-wkc1 — layer spin Hall 对称性分类的直接先例。
5. Chen et al., *Ferroelectric Control of Layer-Resolved Spin and Orbital Transport in a Topological Heterobilayer*, **Nano Lett.** (2026-10-07). DOI: https://doi.org/10.1021/acs.nanolett.6c04169 — FE/层选择 SHE/OHE 已做。
6. *Mixed Weyl semimetals and low-dissipation magnetization control in insulators by spin–orbit torques*, **Nat. Commun.** (2017). DOI: https://doi.org/10.1038/s41467-017-01138-7 — k–magnetization mixed Berry curvature 与几何 SOT 先例。
7. *Sliding-tuned quantum geometry in moiré systems*, **PRB 113**, L121411 (2026). DOI: https://doi.org/10.1103/xhwb-5z6g — 滑移增强几何与非线性 Hall 已做。

## 现在不做的事

不改已有 Spin QG V32 冻结代码、不启用真实材料 DFT、不把单个非零 toy 数值作为创新证明、不写包装性论文摘要。下一次理论推进优先测试 **B/C 交叉互易关系的真实自旋/结构可观测意义**，以及 **D 的双连接反例能否推翻某个具体未经完整检查的候选响应**。如果都只是通用代数，就继续搜索其他物理问题。