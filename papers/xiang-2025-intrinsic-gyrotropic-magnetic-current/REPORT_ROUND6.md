# 第六轮：用不依赖本征矢的解析方法审计 IGMC（2026-10-20）

## 核心增量

1. 对 Dirac 两带 `H=d(k)·σ` 用 Pauli 投影与矢量恒等式独立推导：`Z_-+=−∂_k n_β`，`2Q_-+/Δ=d·(∂_k d×e_β)/(2|d|³)`。完全不需要计算 Bloch 本征矢、Berry connection 或旧 Kubo 内核。全部 12 种（两磁态×2 种电流分量×3 种磁驱动方向）的随机 k 点测试最大误差低于 `2.7e−15`。
2. 从原文完整的正负频率分母出发获得精确的非共振有限频率复响应。两组时间反演相关的 Hamiltonian 满足 `K_TR(ω)=K_H(ω)*`，即 **conduction 的实象限 T-even，displacement 的虚象限 T-odd**。1/3 THz、12 个张量/频率组合与独立 Kubo 的最大差异为 `5.56e−17`（Å⁻¹）。
3. 不再依赖历史的均匀网格，改用 **Gauss–Legendre 整 BZ 积分**；128² 与 301² 均匀网格的 `D/ω` 差为 `2.01e−11 Å⁻¹ eV⁻¹`。

## 定量结论

| 项目 | 结果 |
|---|---:|
| Gauss 128²：绝缘体 `D/ω` | `0.01523407801563 Å⁻¹ eV⁻¹` |
| Fig. 2(d) `μ=0`, `f=1 THz`, `g=2`，正文公式 | `0.01775389319 A/(m·T)` |
| SM Eq. (34) 中间式按字面多乘 2，仅供诊断 | `0.03550778639 A/(m·T)` |
| 原文 Fig. 2(d) 肉眼估计 | 约 `0.040 A/(m·T)`，非原始数据 |

*结论：* 新的独立解析与不同数值积分再一次支持 **正文 Eq. (12) 的前因子**；没有科学依据把已经验证的公式乘以 2。Fig. 2(d) 的图源幅度问题 **仍未解决**，可能需要原始绘图数据或作者说明。Eq. (14) 的符号问题已用**完整复频率与两磁态**作出数值证据，但其作者具体约定仍有疑问。

## 交付与范围

- 理论审计：`docs/ROUND6_DVECTOR_FREQUENCY_AUDIT.md`
- 复现实验：`src/igmc_lab/d_vector_independent.py`, `scripts/audit_round6.py`, `tests/test_round6_dvector.py`
- 原始证据：`results/raw/round6_*`，新冻结基准：`results/frozen/round6_dvector_reference_v1.json`
- 可供用户决定是否联系作者的**未发送**问题：`docs/AUTHOR_DISCUSSION_BRIEF.md`

本轮不会假定已取得论文的原始图源数据；金属态不同化学势网格暂不作为严格收敛验收。没有发送邮件，也没有向 GitHub 推送。整篇复现状态仍为 `PARTIAL`。
