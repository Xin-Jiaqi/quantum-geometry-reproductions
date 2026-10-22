# G3 继续 Brainstorm：真实 FE 切换路径而非静态端点张量

**2026-10-22。状态：SCIENCE SEARCH / NO PHYSICS DISCOVERY / ZERO NEW COMPUTE.** 本备忘录必须与[20 条淘汰主张和 8 个硬关口](BD_IDEA_ELIMINATION_GATES_20261022.md)一起阅读。原论文复现和 MnPS₃ R1 状态均不改变。

## 1. 从 G2 的负结果提炼出来的真正问题

普通 bulk tensor、Wyckoff 局域 induced representation、FE 端点连接和 Onsager 已有成熟理论。G2 p3 的 hidden local \(S_z\) 且全局零响应是正确的**对称性例子**，却不是新机制，也不是真实 ±P 铁电双态。

真实 FE 翻转往往不是两个均匀结构之间做一次刚体 Seitz 操作。2025 [Ke, Liu & Liu, PRL 135, 046201](https://doi.org/10.1103/jhlq-2dd7) 强调波状极化畴壁的迁移和有效电荷张量；2026 [Nano Lett 26, 1043](https://doi.org/10.1021/acs.nanolett.5c05482) 已比较不同 FE 壁型。因此**连接算符（端点对称关系）≠ 实际动态切换路径**。这是一个潜在的物理分岔，而不是新定理。

**候选问题：** 一个有真实滑移铁电畴壁及磁性子格序的双层，是否可能在两个端点的普通 E→总磁矩系数都禁戒的情况下，**移动壁**产生与不同物理磁模耦合的局域广义力，其可观测频率谱/壁方向选择不能被简单的磁性多极矩、壁局域点群选择定则或常规极化变化解释？

这比笼统声称“畴壁可产生磁响应”严格得多；**目前纯假说**。它可以独立于 MnPS₃ R1 进行公开结构及群论预检，不要求立即 DFT。

## 2. 可观测量必须定义清楚

对有结构 profile \(\delta(x,t)\) 和源认证磁性位点映射的壁，设每个 Mn/磁性 site 的交换力
\[
F_\nu(x,t)=-\partial_{\theta_\nu}H_{\rm ex}[\delta(x,t),\{m_i\}],
\]
并分别辨别：
- **fixed wall AC**：\(\chi^R_{F_\nu,E_b}(\omega,q)\)，壁结构固定时纯电子/磁力响应；
- **moving wall**：\(\delta F_\nu\) 来自 \(\partial_\delta F_\nu\,\delta\delta(x,\omega)\)，其中 \(\delta\delta\) 是经过验证的**壁动力学**；
- **FE current / bound charges**：极化变化产生普通位移/束缚电流，必须扣除；
- **spin-orbit / exchange / lattice torque**：三个实际电子—磁矩/晶格转移渠道分开，并检验连续性；
- **detector**：允许的近场磁光、空间分辨 SHG、局域 NV/磁探测或远场共振谱必须有可实现的耦合与分辨率，不能只画隐藏位点系数。

固定壁与运动壁不可混为同一 Kubo 响应；真实 \(q\) 和实验均匀场能否激励该模式也不可省略。非线性、迟滞与路径依赖进一步要求离开静态 equilibrium tensor 框架。

## 3. 直接竞争文献：先排除再计算

- [González-Hernández et al., Nat Commun 15, 7663 (2024)](https://doi.org/10.1038/s41467-024-51565-6)：子晶格局域 Edelstein/磁性力矩及无 SOC 效应。
- [*Local Magnetoelectric Effects as Predictors of Surface Magnetic Order*, PRX 15, 021094 (2025)](https://doi.org/10.1103/7brd-lynv)：bulk-site 磁性多极矩预测界面/表面磁序变化。
- [Shen et al., *Magnetoelectric Torque in Polar Magnetic Bilayers*, PRL 136, 016702 (2026)](https://doi.org/10.1103/l5fd-kpn5)：纯电场、无 SOC 也可驱动双层磁性力矩。
- [Skogvoll et al., *Local indirect magnetoelectric coupling at twin walls in CaMnO₃*, PRB 113, 224102 (2026)](https://doi.org/10.1103/6y99-279q)：中心对称母相中壁局域极化及磁化。
- [*Superlubric Motion of Wavelike Domain Walls in Sliding Ferroelectrics*, PRL 135, 046201 (2025)](https://doi.org/10.1103/jhlq-2dd7)，以及 [*Domain-Wall-Mediated Interfacial Ferroelectric Switching*, Nano Lett 26, 1043 (2026)](https://doi.org/10.1021/acs.nanolett.5c05482)：壁动力学本身已成熟。
- [*Octupole-driven spin torque switching of antiferromagnetic tunnel junctions*, Nature (2026)](https://doi.org/10.1038/s41586-026-11172-5)：反铁磁不同子格流与切换，近邻机制。

**当前可见文献已占据**“壁局域磁电”“纯 E 磁矩切换”“层/子格力矩”“移动铁电壁”的一般机制。尚未找到**已证实的独立新规律**，所以此处只作 OPEN_HYPOTHESIS。

## 4. 低成本 G3 决策（按先后顺序）

| Stage | 零/低成本检验 | 硬淘汰条件 |
| --- | --- | --- |
| 1 | 既有 BSF 结构数据或已发表真正 FE ±P 双态及**壁类型**，来源 DOI、位点对应、对称/极化与磁序 | 壁/磁序来源不明确 → BLOCKED |
| 2 | 各 FE 端点、壁截面和磁性结构 stabilizer/Seitz/antiunitary 检查，分别做输出局域磁性-mode 表示 | 只是普通降低点群使 tensor 允许 → NO_GO_CLAIM |
| 3 | 确定一项可读出 Fν 与实验驱动，分开固定壁/运动壁及磁-晶格耦合 | q/频率/局域探测不可能 → BLOCKED |
| 4 | 与 2024 local torque、2025 PRX surface ME、2026 PRL magnetoelectric torque 与壁动力学对**同一个量**比较 | 只是复现已有磁性多极矩/壁响应 → NO_GO_CLAIM |
| 5 | **只有经过 1–4** 才允许一个独立可重复、参数受约束的简化 Hamiltonian/Kubo 或少量材料计算 | 缺严格零对照、因果/时间协议 → 不启动任何生产作业 |

**并行科学选择：** 若 G3 domain-wall 联结也因已有 PRX/PRL 文献失去独立空间，保留更通用的「源认证 FE 端点 × site-mode observables × full SOC/antiunitary」软件审计为开源群论成果，不强行写成新物理文章。

**当前资源：** 零新 VASP、Slurm、WSL、日本超算任务；MnPS₃ R1 CLOSED、R2-P0 只读待真实数据回报。**此文不是计算授权。**
