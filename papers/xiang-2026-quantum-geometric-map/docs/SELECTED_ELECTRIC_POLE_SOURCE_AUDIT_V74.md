# v74 | 逐项定位原始 SM Eq. (15) 的电场绝热极点

**原论文**：Xiang et al., *Quantum Geometric Map of Magnetotransport*, PRB 113, L201406 (2026)。本审计的源公式依据此前已对原始论文和 16 页补充材料进行 SHA256 锁定的 v49、v50、v53 审计。**总体论文复现仍为 `PARTIAL`。**

## 一、为什么 v74 不重复 v50

v50 已从原始补充材料 SM Eqs. (4)–(10) 的双频率 Liouville 微扰提取与 SM Eq. (20) 相同的电场绝热极点系数，也已发现 SM Eq. (31) 的相反符号。v74 的独立增量是 **将电场带内极点、另一个电场带间部分、两种时序分开计算，证明为什么两条很大的时序响应必须相加，并检查调节量趋零的路径**。不能把代码版本数或测试数当成论文物理完成度。

## 二、严格的源模型和算符约定

使用 v49/v50 的连续倾斜 Dirac 两带模型、固定实际自旋算符 `sigma_x`，显式保留主文 Eq. (10) 和 SM Eq. (62) 两种相反旋向。`bar B=-g mu_B B`，因此输出电流按物理 B 的负号转换；模型设 `e=hbar=vF=g mu_B=1`，结果不称 SI 电导。原始磁微扰是 **spin Zeeman**，不是轨道 Peierls 位移。

写成原始源采用的 Liouvillian 形式：

\[
\mathcal L_z X=zX+i[H,X],\quad z_E=\eta_E-i\omega_E,\;z_B=\eta_B-i\omega_B,\;z_\Sigma=z_E+z_B.
\]

计算

\[
\rho_E=\mathcal L_{z_E}^{-1}(D_b\rho_0),\quad
\rho_B=\mathcal L_{z_B}^{-1}(-i[S,\rho_0]),
\]
\[
\rho_{B\to E}=\mathcal L_{z_\Sigma}^{-1}(D_b\rho_B),\qquad
\rho_{E\to B}=\mathcal L_{z_\Sigma}^{-1}(-i[S,\rho_E]).
\]

v74 **不使用 v50 的有限 k 差分**。它直接采用本征能带 Fermi 矩阵的 Fréchet 导数
\((D_b\rho_0)_{mn}=f^{[1]}(\varepsilon_m,\varepsilon_n)(v_b)_{mn}\)，其中带内是 \(f'_nv^b_{nn}\)，再对磁场响应使用 Liouvillian 的算符求导恒等式：

\[
\mathcal L_{z_B}(D_b\rho_B)=-i[S,D_b\rho_0]-i[v_b,\rho_B].
\]

这产生独立于 v50 数值微分的新计算方式，并保留全部两路时序。

电场 E→B 进一步分解为明确的 **带内输入极点**

\[
\rho_{E\to B}^{\mathrm{diag}}=
\frac{1}{z_E}\,\mathcal L_{z_\Sigma}^{-1}(-i[S,D_b\rho_0^{\mathrm{diag}}])
\]

及其余带间电场输入部分。**注意其余部分也可能包含输出频率之和的极点，因此不能笼统称其“解析正则项”。**

## 三、数值证据及关键负控制

参数为 \(\beta=1/8,\mu=.3,T=.03\)，64×88 Gauss/角度积分、UV cutoff 2.5，主文旋向，XY 通道，\(\omega_E=\omega_B=0,\eta_E=.003,\eta_B=.001\)。以物理 B 转换后的无量纲结果如下：

| \(\eta_E\) 乘后的独立贡献 | 数值 |
|---|---:|
| B→E | **+0.005619090327** |
| E→B | **−0.008125646806** |
| 两路完整合计 | **−0.002506556479** |
| 仅电场带内输入极点 | −0.002506387722 |
| 直接收缩源 SM Eq. (20) | **−0.002506260111** |
| 印刷的 SM Eq. (31) | **+0.002506260111** |

完整双场结果与独立 Eq. (20) 收缩的差为 \(2.96\times10^{-7}\)（当前有限调节量），相对约 0.0118%；v50 原始有限差分与 v74 Fréchet 解在完整 GitHub 仓库中还由独立回归检验。两路时序中的大项相互抵消，没有删除不喜欢的项或调参校正符号。

对 \(\eta_B/\eta_E=0.1,1,10\) 同时趋零，\(\eta_E\chi\) 的终点 \(\eta_E=.001\) 分别为 **−0.00250627558、−0.00250634886、−0.00250777426**；有限调节量造成明显路径误差，较小 \(\eta\) 时均接近 Eq. (20)，但**它们各自的 BE/EB 部分具有显著的路径依赖**。不能把有限的 4 个样本称作已严格证明所有路径的极限交换性。

若固定 \(\omega_E=.15,\omega_B=.04\) 再令 \(\eta\) 减小，乘以 \(\eta_E\) 的完整电流模逐渐减小（\(\eta_E=.001\) 约 \(1.18\times10^{-5}\)），符合与 DC 电场极点不同的有限真实频率限制。

另外，\(\mu\to-\mu\) 使有限温度两带源 SM Eq. (20) 的符号翻转；原文主/SM 两种旋向仍给出相反签名，不做未经作者证实的修复。

## 四、物理结论及剩余阻碍

- **VERIFIED**：给定原始 SM 的两带 spin-Zeeman 双频率形式，两路微扰与电场带内极点的独立算符分解；原始 Eq. (20) 的渐近系数；大项抵消；有限频率负控制；多个展宽路径及原来 Eq. (31) 不一致的符号记录。
- **CONDITIONAL**：连续 Dirac UV cutoff、指定 \(T,\mu,\beta\)、原始绝热无穷小而非真实碰撞积分；带内电场极点经现象学 \(1/\eta_E\to\tau\) 的替换只是原文的一项特定规定，v74 未推导微观耗散时间。
- **UNRESOLVED**：印刷 SM Eq. (31) 与 Eq. (20) 的整体符号、主文和 SM 相反的自旋旋向；完整源公式的直流输运物理与场预备态、顶角守恒、真实材料响应；原文 Fig. 2 SI 归一化；轨道 MNHE 的完整 Wannier 磁位移、SM Eq. (49)/(52) 因子以及三带遗漏。

严格结论：**v74 比 v50 更明确地定位并独立检验了已知的极点结构，但没有消除论文主线的源公式矛盾，也没有宣称原论文可被完整复现。论文状态仍为 `PARTIAL`。**

## 实现与验证

`src/qgmap/selected_electric_pole_v74.py`；`tests/test_selected_electric_pole_v74.py`（45 项）；`tests/test_selected_electric_pole_v74_frozen.py`（1 项）；`tests/test_source_v50_bridge_v74.py`（6 项，完整仓库 v49/v50 依赖）；`scripts/reproduce_v74.py` 和 `results/frozen/v74_selected_electric_pole.json`。本地独立运行 46/46 PASS；六项源 v50 交叉需在完整 GitHub 仓库中验收。v74 → v75 交接详见 `docs/HANDOFF_V74_TO_V75.md`。
