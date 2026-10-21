# BSF × NMEE | 三条并行研究主线（2026-10-21）

**导航性质，不改变任何冻结的 Jia 等 PRL 复现。** 所有物理结论仍以原始独立研究报告及数值证据为准。

| 主线 | 核心问题 | 现阶段已完成 | 下一真正科学门槛 |
|---|---|---|---|
| **① Idea / 新物理** | 何种 \(\Gamma^{\mathrm{stat}}_{iz;j}\) 在真实滑移结构中允许？层奇偶、磁场奇偶、静态/驱动条件为何？ | [Stage6–7 机制/群论](../RESEARCH_INDEX.md)、[Stage8f glide 理论](../material-screening/STAGE8F_GLIDE_HIDDEN_LAYER_PHYSICS.md)、[Stage9 spin conservation no-go](STAGE9_PHYSICAL_NOVELTY_AND_NO_GO_GATES.md) | **真实材料层差静态响应非零**；严格区分 driven-DC NMEE |
| **② Test / 可证伪计算** | 原子实际对称性、SOC 自旋混合和 \(\Gamma\) 的 full-BZ 可重复值分别如何？ | [Stage8h 日本 VASP 原始审计](../hpc-stage8g/STAGE8H_FINAL_SCIENTIFIC_ACCEPTANCE_20261021.md)：受限滑移中点 glide 保留、总张量禁止、无计算静态系数。 | [Stage8i 两份 SOC 单点计算](../hpc-stage8i/README.md) 由独立 DeepSeek harness 负责；不要因等作业而扩大跑量 |
| **③ Peer / 创新性校准** | 哪些“新”已被论文做过？我们必须多证明什么？ | **[Stage10 14 篇同行证据审阅](STAGE10_PEER_LANDSCAPE_20261021.md)**、[机器可读矩阵](STAGE10_PEER_MATRIX.json)；新加入 **PRB 2026 LNHE** 直接相邻竞争。 | 获得 P01–P05/P07 的原文 + SM 逐式比较，并在实际 \(\Gamma\) 数值出来后更新创新性定位 |

## 暂不作为主稿结论的内容

- **普通 LEE/层奇补偿、一般铁电 Rashba/EE 开关、spin QGT 的存在、mixed Berry curvature 在绝缘体可非零**：已有同领域先例。
- **半滑移中点可以存在隐藏层自旋数值信号**：目前只在**toy model**有存在性数值；真实 WSe₂ 层差响应尚未求出。
- **“约 24.1 meV/胞 的开关势垒”**：只是三点 non-SOC fixed-xy、z-only 受限 PES 中点能量差，不是 NEB。
- **静态 \(\Gamma\) 等于论文驱动 DC NMEE**：尚未证明。
- **SOC 电子结构已完成**：截至本文件，Stage8i GitHub 仅部署 agent-ready 输入与 PBS 程序；需外部 harness 的真实作业报告。

## 推荐沟通节奏

只有**真实材料自旋混合与带隙**、**层分辨静态响应数值**或**静态/驱动理论对照**中出现新的可核验结果，才考虑修改核心物理主张。同行文献调研随之更新即可，不要为了追求长文献列表而打断计算。
