# IGMC 论文复现 · 第一轮研究报告（2026-10-19）

目标：Xiang et al., PRL 134, 116301 (2025)，DOI 10.1103/PhysRevLett.134.116301。工程参照：`Xin-Jiaqi/quantum-geometry-reproductions` 旧 NMEE 模块，远端快照 `78a34431bd4a165484ac3f714786b98120132f81`，未修改旧模块。

## 本轮实际完成

1. 逐图逐式验收地图：正文两张复合 Figure 共 8 panel，正文 Eq.(1)-(17)，SM I-VI 及 Tables I-III。SM 无新增数值图。
2. 两类独立模型计算：Rashba Eq.(15)/Eq.(16) 和 lattice Dirac Eq.(17)，包括 Zeeman quantum geometry、conduction IGMC、displacement IGMC、finite-frequency Kubo。
3. 生成 `results/figures/fig1_rashba_model.png` 与 `fig2_dirac_model.png`，可视化对应 8 panel 的主要机制，**未伪称绝对 SI 或像素一致**。
4. 11 个 unit tests 全通过，6 个模型数值冻结基准可重复验证；`scripts/verify_release.py` 返回 `PASS_MODEL_OFFLINE`。
5. 独立路线：Rashba 解析 Eq.(16) vs geometric k 积分；Dirac Fermi-sea Eq.(11) vs Fermi-surface Eq.(13)；Dirac low-frequency C+iD vs full-frequency Kubo；基变换、time-reversal 和网格收敛。

## 具体数值（单位为自定规范下的自然模型单位）

- Rashba λ_R=0.6 eVÅ, m=1.1m_e：μ≥0，去除 Eq.(16) e/ℏ prefactor 后 σC_xy=0.0137851385 Å⁻¹，直接数值积分差约 0.060%。μ=−0.015 eV 的直接积分误差约 0.115%。
- Dirac β=0.25 eV：Γ 能隙=0.5000 eV。μ=0、N=101²：σC_xy≈0，σD_xx/Ω≈0.015234078；Ω=0.001 eV 时 `Im Kubo=1.5234188588e−5` vs paper low-Ω `D=1.5234078036e−5`。
- Dirac μ=0.85 eV：N=201²，σC_xy=−0.0806328138，σD_xx/Ω=−0.0503994566。与 N=301² 的相对差分别 0.183%、0.042%。原论文未提供数字化曲线，无法据此称为与原图定量拟合完成。

## 重要科学问题

**正文 Eq.(14) 的反幺正符号不能未经审计直接用于两个响应。** Eq.(14) 把 `η_T=-1` 同时用于 C/D；若直接对有 T 对称性的 Rashba 使用，会禁止非零 C，但作者 Fig.1 和主文显然得到非零 C；且原文 Table I 的 Z 为 T-even、Q 为 T-odd，SM 两类 tensor 的 Jahn 记号不同。已从直接 `T H(k) T^-1=H(-k)` 和显式 Z/Q 积分验证这一冲突。需优先厘清 antiunitary 是否涉及不同 response channel 和 ω/complex conjugation。

## 仍未验收

- Dirac Fig. 2 的 **绝对 SI 数值**及原论文 Fig.2(d) 标称 ω=10¹² Hz 对应 2π 的问题；本轮 Fig.2(d) 明确画 `σD/Ω`。
- Rashba k=0 附近的自适应收敛、Dirac 费米面尖峰的局部积分误差；不能宣称所有 µ 点已收敛。
- SM §§III、IV 从头独立手推的完整代数验收、SM Tables I–III 全磁点群枚举。
- 作者图像数字化、实际 Wannier 材料测试、GitHub 新分支与 CI 验证（尚未推送）。

## 建议第二轮任务

先解开 **Eq.(14) 的 channel-specific T parity / frequency protocol**，并完善 Eq.(8)→Eq.(11)/(12) 的独立解析路线；随后转换 SI 单位和确认频率；最后做费米面峰值的自适应积分和原图对照。具体续跑命令、数据和未决问题已在 `docs/HANDOFF.md` 与 `docs/WORKLOG.md`。
