# A 方向探索归档：隐藏非线性自旋流、层间几何与可观测性

> 记录日期：2026-10-21。性质：**scientific brainstorming archive / NOT a validated numerical release**。记录对话中的物理问题、可独立核对的代数与尚待复算的数值线索，避免未来新会话误将探索性结论宣传为发现。本轮没有修改论文复现仓库的冻结代码。

## 0. 三条独立研究方向

- **A**：中心对称双层的层反向非线性自旋流，之后尝试绝缘体有限频率、层间相干与可观测性。
- **B**：层间杂化/滑移坐标如何改变自旋相关量子几何；要求区分谱、波函数与表象效应。
- **C**：绝缘体中的交流 displacement spin Hall，以及磁/层/堆叠切换的严格对称性关系。

A 当前 **HOLD / 暂存**；B、C 重新评估，详见《2026-10-21-B-C量子几何选题深度研判.md》。

## 1. A1：整体反演对称不禁止所有层反向响应

选定层交换反演 \(\mathcal P=\tau_x\)（空间同时 \(\mathbf k\mapsto-\mathbf k\)），采用普通自旋流 \(J^i_a=\{v_a,\sigma_i\}/4\)，层反向流 \(J^i_{-,a}=\{\tau_z,J^i_a\}/2\)。反演下：电场 E 极矢量为奇、速度奇、自旋轴矢量偶、层差 \(\tau_z\) 奇；因此

- 总自旋流 \(J^i_a\) 的 \(n\) 阶电场响应系数在两个**严格反演相关的构型**间满足 \(\chi^{(n)}_\mathrm{tot}(-P)=(-1)^{n+1}\chi^{(n)}_\mathrm{tot}(+P)\)。
- 层反向自旋流 \(J^i_{-,a}\) 满足 \(\chi^{(n)}_-( -P)=(-1)^n\chi^{(n)}_-(+P)\)。
- 处于中心对称的**单一构型**时，总二阶自旋流必为零，而层反向二阶响应可以被允许；但剩余镜面、旋转、反幺正对称性还可能继续禁止具体分量。

这只是张量选择定则，不保证任何 \(\chi\) 非零，更不意味着可测净自旋积累。垂直门压 \(\Delta\tau_z\) 破坏反演并释放总二阶响应，属于合理调控思路，但与已发表的 LNHE/LEE 机制相近。

## 2. A2：最小连续双层 Rashba 与层间相干

\[
H(\mathbf k)=k^2/(2m)+t_\perp\tau_x+\alpha\tau_z(k_y\sigma_x-k_x\sigma_y),\quad
\varepsilon_\pm=k^2/(2m)\pm\sqrt{t_\perp^2+\alpha^2k^2}.
\]

该模型可拥有 P、T、PT 双重简并。在 t=0 时，层反向自旋流仍可能非零，表明 **nonzero hidden response ≠ interlayer-coherence mechanism**。打开 t 后的响应变化同样可能来自谱、费米面、算符期望值的改变；不能只用响应曲线判断 quantum metric 的因果作用。

对简并投影 \(P_a\)，可考察 \(iP_a[\partial_xP_a,\partial_yP_a]P_a\) 的非阿贝尔结构，但“存在非阿贝尔 Berry 曲率”既不自动说明它贡献了目标二阶自旋流，也不自动说明该几何对象新颖。需要完整有限 η/有限频率 Liouville、相同谱对照和自旋算符定义一致性。

## 3. A2 等能谱对照与关键的表象陷阱

对话中测试的家族是

\[
\begin{aligned}
H_\beta={}&\varepsilon_0I+t(\mathbf k)\tau_x+\alpha\sin k_y\tau_z\sigma_x\\
&-\alpha\sin k_x\cos\phi_\beta\tau_z\sigma_y+\alpha\sin k_x\sin\phi_\beta\tau_y,\\
\phi_\beta={}&\beta(1+\cos k_x).
\end{aligned}
\]

其本征值对 β 不变。然而模型间存在 \(k\)-依赖酉联系 \(H_\beta=U_\beta H_0U_\beta^\dagger\)，其中可取 \(U_\beta=\exp[-i\phi_\beta\tau_x\sigma_y/2]\)（符号依赖约定，重新使用前要核算）。

**两种解释不能混淆：**

1. 如果是同一个物理体系的 Bloch 表象改变，位置/Berry 连接、速度、自旋和层算符都应同时作正确变换，所有可观测量不变。
2. 如果是固定物理轨道基底中真实改变的 hoppings，则可以视为不同模型；但必须建立原子位移/跃迁微观来源，检查 \(v=\partial_kH\) 与位置算符的正确关系。

本轮之前的对话曾报告 β 扫描下 quantum metric 与二倍频数值发生差异，**但没有在可交付、版本固定的代码和完整微观算符中独立重跑**。这些数字的等级只能是 **UNVERIFIED EXPLORATORY OBSERVATION**，不得作为定量论文证据，更不能从“不同”推断“新自旋量子几何机制”。

## 4. A3：算符恒等式与反例（应优先独立验算）

令 \(v_y=\partial_{k_y}H_\beta\)，层反向普通自旋流

\[
J^x_{-,y}=\frac12\left\{\tau_z,\frac14\{v_y,\sigma_x\}\right\}.
\]

用 Pauli 反对易关系有候选代数恒等式

\[
J^x_{-,y}=\frac{\alpha}{2}\cos k_y I_4
+\frac12(\partial_{k_y}\varepsilon_0)\tau_z\sigma_x.
\]

当 \(\partial_{k_y}\varepsilon_0=0\) 时，算符正比单位矩阵。在保持动量处的单粒子迹恒定且完全绝缘填充的**特定均匀场 Liouville 方案**下，激发产生的各阶响应 \(\operatorname{Tr}(\rho^{(n)}J)\) 为零；因此一般性的“metric 非零 ⇒ 二倍频自旋流非零”是假命题。

这并不是普遍无响应定理：\(\varepsilon_0(k_y)\)、自旋流定义、Bloch embedding、有限空间 q、边界及 SOC torque 都会影响假设。

### 当前会话新增的独立最小矩阵核查（不等于响应函数复现）

取 Pauli tensor-product 4×4 矩阵及 \(k_x=0.73,k_y=0.48,t=0.45,\alpha=0.6,\phi=0.81(1+\cos k_x)\)：NumPy/SciPy 中直接求 \(U_\beta H_0U_\beta^\dagger\)，与所写 \(H_\beta\) 最大矩阵元差 \(7.85\times10^{-17}\)；\(\varepsilon_0=0\) 时，上述 \(J^x_{-,y}=\alpha\cos k_y I_4/2\) 的最大矩阵元差严格为 0（浮点输出）。**只验证这一处矩阵代数，不验证此前对话所有的 BZ 二倍频积分，也不构成 physical response PASS。**

## 5. 真实可观测性：连续性方程

对于 layer-odd spin density \(S^x_-\)，应从

\[
\partial_t S^x_-+\nabla\cdot \mathbf J^x_-=\mathcal T^x_-
\]

出发，其中 \(\mathcal T^x_-\) 包含 SOC 自旋力矩与上下层自旋转移/源汇项。若研究非均匀激发，则有限波矢下 \(S_-(2\omega,q)\) 可满足不同于均匀光场的空间选择定则；只有同时计算 torque、接触与边界条件才可讨论自旋积累或光学读出。

不能将 conventional \(\frac14\{v,\sigma\}\) 在 SOC 下直接当做严格守恒的输运流。

## 6. A 的科学否决与归档结论

| 子方向 | 否决/保留原因 | 状态 |
| --- | --- | --- |
| A1：隐藏二阶自旋流 + 栅压释放 | 机制与 LEE/LNHE/已发表非线性 pure spin Hall 工作重叠 | LOW PRIORITY |
| A2：层间耦合导致 QG 增强 | 几何随耦合变化不保证特定响应变化；等谱模型需排除表象伪效应 | HOLD |
| A3：绝缘双层二倍频 spin current | 简单算符可能严格零；finite-frequency spin current 的 gauge/torque/边界问题突出 | HOLD，未来可与 C 衔接 |

**重新启动 A 的严格触发条件**：存在被独立解析或数值验证的层间相干特有贡献；不存在任何两个解耦单层的重标记解释；提供 gauge-invariant formula、同谱控制、对称性许可、至少一种可读出观测量，且与相近已发表论文有明确不同。

## 7. 先例与边界（外部文献，仅作为 novelty check）

- Zhou et al., *Layer Edelstein Effect*, PRL 137, 076401 (2026), DOI `10.1103/qhvs-pfkq`；用户原始 Markdown 见同一项目文件。
- *Interlayer Electric Multipoles Induced by In-Plane Field from Quantum Geometric Origins*, Nano Letters (2024), DOI `10.1021/acs.nanolett.4c01657`。
- Sarkar, Das & Agarwal, *Symmetry-driven intrinsic nonlinear pure spin Hall effect*, PRB 113, 125106 (2026), DOI `10.1103/swzr-55sp`。
- Yue et al., *Ferroelectric switchable intrinsic nonlinear pure spin Hall current*, PRB 114, L171108 (2026), DOI `10.1103/8htr-q899`。
- Dong et al., *Band-Geometry-Driven Spin Photocurrent in Centrosymmetric Altermagnets*, PRL 136, 216702 (2026), DOI `10.1103/g85x-rgxm`。
- Xiang, Jin & Wang, *Spin Transport Revealed by Spin Quantum Geometry*, PRL 135, 146303 (2025), DOI `10.1103/14mp-263q`，原始 PDF/SM 为本项目上传文件。

## 8. 对下一位研究者的要求

不要在没有重新核算 H、J 和 covariant Liouville 的前提下继承本轮对话中的响应数值。先验证算符恒等式和一个完整 2D 周期格点模型，再讨论新物理。真正的新方向优先考虑 B、C。