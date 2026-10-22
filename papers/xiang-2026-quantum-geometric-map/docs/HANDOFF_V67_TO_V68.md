# v67 → v68 独立科研交接｜2026-10-22

仓库：`Xin-Jiaqi/quantum-geometry-reproductions`，Xiang PRB 原论文科学状态仍为 **PARTIAL**。

## v67 的真正新内容

**不再复用 scalar gamma 模型当作第二个“微观验证”。** 新增两个明确 Hamiltonian 的半无限单通道紧束缚电极、对应的 retarded surface Green 函数与设备/电极的真实 Dyson 自能。固定能量 `mu=.13` 时，谱展宽直接来自可传播的电极通道，磁场依赖的电极耦合 `c(B)=U(B)c_0` 与 `d_B c` 一并加入接触自能和 Gamma 导数。器件重用 v65/v66 的同一有限 Peierls H(B)、轨道位置矩阵 Pa/Pb 和完整 d_B Pa/Pb。

测量的是**双端 reservoir-voltage Landauer 电导** `T(E_gate,B)`，以及独立静态 E_gate 和 B 控制下的 `∂B T`、`∂E_gate T`、`∂B∂E_gate T`。其中 `E_gate` 是器件门电势 Hamiltonian 参数 `H-E_gate Pb`，而驱动电流的是**电极间电化学势差**，两者不能混淆。这是与 v66 不同物理 bath 的 **conserving coherent elastic benchmark**，但不是原 Xiang 论文的 `sigma_xy^orbital(E,B)` 同观测量复现。

参考数值 `T=+0.900626086557784`，`d_B T=-0.180487543398844`，`d_E T=-0.557089765325712`，`d_B d_E T=+0.884057478695506`。三套磁/轨道基底完全一致。完整 Caroli 对应独立 scattering matrix 概率，`S†S=I`、磁场微分后 `∂B(S†S)=0`，spectral optical identity 及其磁导数，物理无电极连接 site 3 的局域连续性及其磁导数全部通过。独立 5 点 `dB/dE` 和 16 角四阶 `E/B` 差分验证完整磁接触，包括 Sigma_B/Gamma_B/d_B Pb。

**关键负控制**：错误地不变换电极耦合，基底改变就会得到不同的 T；忽略非交换轨道变换下的电极自能磁导数，`d_B d_E T` 从 `0.88405748` 错成 `0.61437901`；忽略 `d_B Pb` 错成 `0.78429578`。换电场原点时若不同时移动电极 onsite 和化学势，也会得到虚假的原点依赖。

## 科学验收分类

- **VERIFIED**：指定非相互作用弹性两电极模型中，完整 NEGF/Caroli、unitary on-shell scattering、磁场微分 Ward、电荷守恒、局域连续性、完整磁接触和轨道位置项，以及数值独立差分。v67 新增 45 项物理测试及 1 项冻结验收；独立 v60–v67 共有 **265 项**。
- **CONDITIONAL**：有限四轨道 Peierls 模型、外部提供的轨道间位置偶极矩、只施加在器件上的静态 E gate、理想非相互作用单通道电极。不是完整实际 Wannier 电磁嵌入或许多体散射。
- **UNRESOLVED**：原始 PRB 轨道 E/B 同观测量的微观守恒 Keldysh 与 bulk intrinsic orbital Hall 是否一致；完备的 Wannier 磁嵌入、真实电极/位移电流接触、有限频率与 DC 热力学极限、原文 SM 符号系数及三带遗漏。整篇 **PARTIAL**。

## v68 科学方向

优先选择原论文实际 orbital E/B 观测量的**同观测量** Keldysh 形式，把 v66 定义的 `J_a(E)=i[H-E Pb,Pa]` 接入 reservoir-charge 和 polarization/displacement continuity，并用 Ward 恒等式处理接触电流与守恒自能。不得从 v67 的 `∂E_gate ∂B T` 直接推出原论文 mixed E/B orbital conductivity。

若无可靠的相同观测量理论接口，应明确 `UNRESOLVED`，并继续研究真实 Wannier 嵌入与原文 SM Eq49/52、遗漏三带等问题。**禁止为了提高复现状态拟合符号、因子或冒充完整守恒微观碰撞理论。**

## 代码、CI 与交付

- `src/qgmap/reservoir_ward_v67.py`：半无限链微观自能、散射与微分 Ward、完整 Peierls/position B 接触；
- `tests/test_reservoir_ward_v67.py`：45 项独立物理验收；
- `tests/test_reservoir_ward_v67_frozen.py`：独立冻结数值重算；
- `scripts/reproduce_v67.py`、`results/frozen/v67_conserving_reservoir_ward.json`；
- `docs/CONSERVING_RESERVOIR_WARD_V67.md`：完整科学说明。

本轮发布必须以**最终 PR HEAD 的 GitHub Actions 全量测试**和**squash 合并后 main 上独立 GitHub Actions**均成功为准；每条记录保留提交 SHA、真实运行链接、测试数、`PARTIAL` 状态。独立 ZIP 包 SHA256 只能证明包内文件，不代表自动证明远端 GitHub 内容相同；新源码和测试文件另用 Git blob SHA 逐项比对。
