# Stage 8o 独立方法审查：Stage 8n 根因尚未定位（2026-10-23）

> **证据：**直接阅读了真实推送到 GitHub 的 Stage8n commit `3505a3f2d6e9600e694ad5e6e6b3eef579404069`，尤其是 `stage8n_diag1.py`、`stage8n_diag4.py`、`stage8n_diag3.py`、`stage8n_kubo_core.py`、`stage8n_kubo_derivation.md` 及净化 JSON，并复查了 Stage8m `stage8m_plaq.py` / `stage8m_response.py`。本次**没有接触日本超算原始 MMN/SPN/WAVECAR**，零新增 PBS/VASP 作业，不改冻结的原 PRL 复现，不覆盖历史证据。
>
> **纠正后状态：`MATERIAL_RESPONSE_SYMMETRY_GATE_FAILED__STAGE8N_VASP_ROOT_CAUSE_NOT_ESTABLISHED__IDEA_OPEN`。** Stage8n 的 NO-GO 对继续*生产级数值预测*是对的，但 `ROOT_CAUSE_LOCATED_IN_DATA_SYMMETRY__ENGINE_AND_GRID_CLEAN` 并未被其代码支持。

## 1. 原始数值哪些仍然可信？

- WSe₂ 的原有 Stage8g/h 几何镜面和 Stage8i 已完成 SOC 采样能带、Stage8m 真实 VASP 后处理文件**不会因为这次审计而消失**。
- Stage8m 在 C₃v AB 得到的 Cartesian 2×2 总张量显著偏离理论所需的 \(M=\gamma I\)，且 glide 中点总响应应为零却显著非零。因此**真实材料响应尚未通过 symmetry/null-gate**。
- 现有 96-band 输出缺 band97，95/96 近简并；不能证明完整 96-band 子空间对称封闭。
- 这些事实都**不能唯一识别**误差产生于 VASP 的 ISYM=0、截断 SPN/MMN、倒空间几何或响应估计器。

## 2. Stage8n 的 0.32/1.03 mirror SPN 残差为何不构成根因证据？

**(a) 所用的 Cartesian mirror 不是实际六角晶格的对称操作。** 在真正给定的晶格 \(\mathbf a_1=(a,0),\ \mathbf a_2=(-a/2,\sqrt3a/2)\) 下，代码 `stage8n_diag4.py` 直接套入 Cartesian \(R=\begin{pmatrix}0&-1\\-1&0\end{pmatrix}\)，其 \(A^{-1} R A\) **不是整数 unimodular 变换**。同一矩阵如果属于**分数坐标**空间，则可能有效，但不能不经变基就当作 Cartesian 镜面。

**(b) k 配对错误。** `stage8n_diag1.py` 只读取 26 个 k，净化文件中其镜像匹配为 **0 对**。`stage8n_diag4.py` 虽读取 441 点，却用不考虑 reciprocal-periodic images 的 Cartesian 近邻查找，仅得到 **12 对**。真实晶体整数 \(R_{\rm frac}\) 应用 \(k'=R_{\rm frac}^{-T}k\bmod\mathbb Z^2\) 将完整 21×21 网格映射到自身，即 **441/441** 有定义（不能因此声称所有物理对称性已通过）。

**(c) 自旋方向也没有变换。** 自旋是**轴矢量**，真实 3D 空间操作 \(R_{\rm cart}\) 诱导 \(Q=\det(R_{\rm cart})R_{\rm cart}\)。应比对完整的
\[
 C_{ij}(k)=\operatorname{ReTr}[S_i^\dagger(k)S_j(k)],\qquad
 C(k')=Q C(k)Q^\mathsf T.
\]
若镜面交换 \(x,y\)，\(S_x^2(k)\) 应与 \(S_y^2(k')\) 对应而非同一分量。对截断的 96-band 子空间，甚至完整 Gram 也可能受窗口不封闭影响。

**(d) 真正的水平 glide \(\{m_z|\mathbf t\}\) 对所有 in-plane k **本来就不改变 k**，却在同一 k 使轴自旋 \(S_x,S_y\) 反号、\(S_z\) 不变。应直接检验**所有非简并 glide 本征态**的 \(S_x,S_y\) 对角期望是否为零；进一步的 sewing matrix 检验需要物理波函数。用另一个竖直 mirror 的 12 对 k 近似比较不能代替此关键物理零值。

**结论：**所报 0.32/1.03 由无效坐标映射、未旋转自旋等混合比较产生，不能叫“VASP 数据对称性破坏已证实”。

## 3. Stage8n 将已经除法的曲率误标为 raw Wilson phase

`stage8n_diag4.py` 的 `raw_plaq()` 调用 `eng.plaquette_phase()`，但 Stage8m 引擎返回的是
\[
{\arg\det(L_1L_2L_3L_4)\over \Delta k\,\Delta\lambda},
\]
即**已除以**步长的混合曲率估计，并非原始无量纲相位。因此报告的 \(\sim 0.004\) 不是“原始 rad”。

以其公开的 AB \(k=(0.44,0.12)\) 示例和 \(\Delta\lambda=10^{-4}\) eV、\(\Delta k\approx0.104\,\mathrm{\AA}^{-1}\) 为例，真正角度约 \(4.5\times10^{-8}\) rad，两个 B 值的真实角度差约 \(5\times10^{-12}\) rad。`|dP|/eps` 从报告中的 \(10^9\) 实际降到大约 \(2\times10^4\)（只是最理想的标量浮点 epsilon 比较）。36×36 极分解、行列式与 k 求和的**实际误差**仍没有估计，不能据此宣布引擎精度已过关。

## 4. Kubo 交叉验证确实没有通过

真实 `stage8n_AB_diag3.sanitized.json` 中，10 个 k/方向数据只有一个 rel≈0.037、另一个≈0.467，其余相对差 **9–275 倍**。当前 `stage8n_diag3.py` 只做 \(U^\dagger dM U\) 旋转，却忽略了随 k 变化的 \(U^\dagger\partial_k U\)；这个缺口可能由已有 MMN 链接中的 gauge-covariant 差分补齐，无需新 VASP。单凭其合成模型**对称零点双方都为零**，不能证明非零 Kubo 与 Wilson 的实现等价。`stage8n_kubo_derivation.md` 的投影交换子展开也应重新按 \(A^\dagger B-B^\dagger A\) 审阅，检查自旋矩阵元的共轭与指标、速度定义及整体符号。

该独立 Kubo 对照的结果**不能排除响应估计器或物理算符/坐标实现问题**；将差异唯一归咎于 eiegh eigenvector 相位也尚无证据。

## 5. 数值单位现存错误

以当前定义的*纯几何核*
\[
\partial_{B}\int d^2k/(2\pi)^2\,
i\operatorname{Tr}P[\partial_{\lambda}P,\partial_{k_j}P]
\]
而言，\(\partial_\lambda\) 的单位是 eV⁻¹，\(\partial_{k_j}\) 是 Å，\(d^2k\) 是 Å⁻²，\(\partial_B\) 是 T⁻¹，积分结果应为
\[
\boxed{\mathrm{\AA^{-1}\,eV^{-1}\,T^{-1}}},
\]
而**不是**此前 Stage8m/8n 标记的 \(1/(\mathrm{eV\,T\,\AA^2})\)。实际可观测响应仍需明确电荷、自旋与磁源转换系数及二维厚度，不能继续沿用先前 SI 换算或者与 toy 模型直接比较量级。

## 6. Stage8o 已在 GitHub 执行的工作

- 新增 `stage8o_symmetry_math_audit.py`：六角原胞真实坐标变换、完整 441 点整数对称 k 置换、轴自旋三分量 Gram covariance、真正 loop phase 缩放及二维单位检查。
- 新增 `stage8o_actual_spn_symmetry_screen.py`：读取现有 POSCAR/WIN/EIG/SPN，调用 spglib **实际原子 Seitz 操作**；对 441 点执行 reciprocal-fractional 对称配对及全部 \(3\times3\) spin Gram 变换，记录能谱残差，对含水平 glide 的结构检验**非简并态的面内自旋对角元零值**。其外部远端执行尚未发生，本报告不编造真实数据测试结果。
- 新增配套 synthetic/fail-closed 回归。**所有新增程序均不提交 VASP/PBS 作业，不读取许可赝势或提交大文件。**

## 7. 下一步只做 Stage8o-real，暂不新增 DFT

建议由现有 DeepSeek harness 在日本已保留的真实 Stage8m 目录**只读**执行新脚本，必须输出每个实际晶体 Seitz 操作的 441 点完整映射、真实 Hermitian spin Gram covariance、K/−K 和水平 glide 的 spin selection rule 结果，并归档净化 JSON 及源 SHA。并行修正原始 Wilson phase 实际角度预算和 Kubo/projector 的非零合成 benchmark。

如果正确的 real-data 检验真正表明完整 96-band 数据不符合对称性，再分离“96-band 上限/PAW spin vs MMN/wavefunction/VASP ISYM=0”的来源；**这时**才讨论是否需要重新后处理或新的 VASP 计算。

**目前最重要的科学判断：** Stage8n 停止生产级数值计算正确；但它对数据根因的断言未通过独立方法审查，**科学 idea 仍有开放空间**。
