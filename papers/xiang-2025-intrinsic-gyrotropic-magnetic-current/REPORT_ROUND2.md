# IGMC 复现｜第二轮工作报告（2026-10-19）

目标：Xiang et al., *Intrinsic Gyrotropic Magnetic Current from Zeeman Quantum Geometry*, PRL 134, 116301 (2025)。继续第一轮 `MODEL_SUBSET_PASS / FULL_PAPER_PARTIAL`；**没有推送 GitHub，原 NMEE 基准和第一轮归档均未修改**。

## 本轮取得的实质进展

1. **反幺正张量规则独立验算。** 对 Dirac Hamiltonian 直接构造时间反演伙伴 `mass_sign=-1`，逐 k 检查 `Z` 为 T 偶，`Q` 为 T 奇；完整 Kubo 的实、虚响应对应地变换。独立张量 Neumann 求解器通过 `4mm1′`、`4m′m′`、反演三个典型群。正文 Eq.(14) 的共同 antiunitary `−1` 分支保留为待说明的原文协议，而不是自动修订原文。
2. **THz 约定从 SM §VI 找到明确依据。** 原文为 1 THz 写出的 `ℏω≈4×10^-3 eV` 对应 `ω=2π×10¹² rad/s`，程序显式区分 `cyclic_hz` / `angular_rad_s`。
3. **SI 单位已定量恢复且经过两重交叉验证。** 对 `g=2`，1 Å⁻¹ 自然响应对应 `281.7940 A/(m T)`；Rashba 的 Eq.(16) 平台预测 `3.88457 A/(m T)`，与原图接近 4 一致。Dirac C 在 `μ=1eV` 算得 `−33.49 A/(m T)`，与 Fig.2(b) 近 `−35` 一致。
4. **发现新的、可重复的定量不一致。** 在 `μ=0` 的绝缘体隙内，Dirac D 预测 `+0.017754 A/(m T)`（1 THz），而 Fig.2(d) 平台视觉读取约 `+0.04`。在 `μ=1eV`，预测 `−0.06419`，原图约 `−0.13` 至 `−0.14`。分别进行最高 `401²` 网格和完整 Kubo 校核，差异仍存在。**未引入任意×2比例因子。**
5. 本轮共 **23 项自动化测试通过**（首轮 11 项 + 新增 12 项）；旧的六个数值冻结目标仍不变。

## 关键结果状态

| 工作内容 | 状态 | 解释 |
|---|---|---|
| 反幺正 C/D channel parity | **PASS for clean-model symmetry** | Hamiltonian TR partner + Kubo + Neumann 三路线 |
| Rashba Fig.1(b) 绝对 SI | **PASS model/SI** | 3.88457，独立直接 SI 闭式相符 |
| Dirac Fig.2(b) 绝对 SI | **PARTIAL** | 选定 `μ=1eV` 峰值附近吻合；未全曲线数字化 |
| Dirac Fig.2(d) 绝对 SI | **OPEN / PARTIAL** | 峰值和隙内平台有约两倍差异，非网格/高阶频率造成 |
| SM 磁点群分类 | **PARTIAL** | 仅三个群测试；未枚举所有 MPG |
| 半经典及 Green 完整独立推导 | **NOT STARTED** | 本轮专注物理定义与量纲 |
| GitHub/CI 部署 | **NOT STARTED** | 仅本地验证，需后续审核后推送 |

## 下一轮唯一优先科学任务

检查 **Fig.2(d) 与正文 Eq.(5)/(8)/(12) 的幅度差异**，尤其作者的 `ω`、`g`、spin-rotation generator `σ/2` 及双带计数。在未找到独立证据前，禁止通过拟合倍数“修复”图像。完成后再做全化学势分段收敛和 SM 全 MPG 枚举。

详细物理推导、原始数据和复现路径见 `docs/ROUND2_TIME_REVERSAL_AND_UNITS.md`、`results/raw/round2_si_audit.json`、`docs/HANDOFF.md`。
