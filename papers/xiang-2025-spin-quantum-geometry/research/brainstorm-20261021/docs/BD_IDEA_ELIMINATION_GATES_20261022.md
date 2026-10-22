# B+D 科研 Idea 淘汰机制 V1

**2026-10-22。** 科研证据与创新性审核规范；保留全部负结果、代码与最接近的同行来源。**淘汰具体论文主张，不删除技术资产或否定整个物理方向。** 独立研究不修改 Xiang 原论文复现代码、测试与验收，MnPS₃ R1 保持关闭，不授权任何计算重启。

## 证据状态（严格分离）

- **METHOD_PASS / MODEL_PASS**：群代数、Ward 恒等式、toy Kubo 或规范测试通过；绝不等同新发现。
- **NO_GO_CLAIM**：主张错误、只是不变量/重标记、或同机制已有论文；保存全部证据。
- **BLOCKED_EVIDENCE**：缺源结构、真实磁序、同规范观测算符、收敛计算或实验读出。
- **BLOCKED_NOVELTY**：同观测量最接近文献尚未逐公式比对；不以摘要未出现公式作为优先权证据。
- **PHYSICS_CANDIDATE**：经过明确的最小证伪、同行比较且有可测预测；仍不是被证明的新物理。
- **MATERIAL_ACCEPTED**：原生数值收敛且物理量独立验收；不自动给出论文新颖性。

## 强制八关，优先做零成本筛选

| Gate | 审查对象 | 必须立即停止的情况 |
| --- | --- | --- |
| K0 问题和算符 | 源扰动、测量算符、响应阶数、时间/frequency、空间层/位点/体、位置接触项 | 把局域自旋当输运自旋流；层反向当 Néel 模；遗漏定义 |
| K1 完整对称关系 | 区分结构连接、非相对论自旋群、SOC 完整幺正及反幺正连接 | 全幺正伙伴仅张量旋转却称新增独立通道；反幺正 AC 同序响应直接反号 |
| K2 Wyckoff/磁序 | 真 Seitz 站点交换、分数平移、磁子格、轨道诱导表示与置换基不变性 | 把一般轨道 3d 与三个特殊 1a/b/c 假定为严格母群对称 FE 伙伴 |
| K3 因果响应 | 完全谱求和/Onsager/接触/规范/Kramers-U2；核实响应定义 | 只用 Ward 恒等式、允许张量、群论维数作为新机制 |
| K4 反事实 | SOC-off、层间-off、磁序反转、k/eta/omega；区分 rank 与 Hom | 把 Hom 系数数目当可激发模式数、无物理实验驱动的有限 q |
| K5 同量同行 | DOI、原论文相同测量算符、相同源、机制、协议、近似、读出；比较最接近公式 | 同样的机制已在 2018–2026 文献出现，或还没核对就说“无人研究” |
| K6 可测性 | 激励频率/波矢，探测层/子格分辨、材料 FE 稳定性/磁序/模式 | 只预测一个观察不到的分量或不能生成目标波矢 |
| K7 算力 | 群代数 → 少参数独立 toy Kubo → 真实材料；用户独立批准 | SCF 未收敛却将 RUN_FINISHED 当 FE_VERIFIED；未经批准重新提交 |

**程序化检查：** [group-theory-operations-toolkit Draft PR #32 的 idea_elimination_gate.py](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/idea_elimination_gate.py) 可消费人工填报的候选 JSON，输出 NO_GO / BLOCKED / PHYSICS_CANDIDATE；程序不执行 HPC、不独立验证人工填的同行证据，也不授予创新或材料证书。该研究分支尚未通过保护仓库的合并审核。

## 14 条主张的累计决策登记

| ID | 具体主张 | 结论与理由 | 重新打开的必要新证据 |
| --- | --- | --- | --- |
| BD01 | 四带 PT 模型允许裸自旋相关层间 hopping | **NO_GO_CLAIM**：既有表示的 PT 禁戒；扩轨道是另一模型 | 物理反宇称轨道和 SOC downfolding |
| BD02 | toy 交换/层间力矩抵消就是新 Spin QG | **NO_GO_CLAIM**：连续性和交换力矩为已知；但 toy MODEL_PASS | 同可观测量区别于 NSOT 的量子几何核和实验 |
| BD03 | C3 稳定态线性面内电场产生宏观 Sz | **NO_GO_CLAIM**：C3 禁戒 | 真实低对称态或明确不同阶次/局域量 |
| BD04 | 完整幺正 FE 伙伴创造新的可激发秩 | **NO_GO_CLAIM**：等变性保持响应秩和非零奇异值 | 真正不等价 SOC 磁性态或固定外部环境 |
| BD05 | 反幺正伙伴 AC 同序响应简单反号 | **NO_GO_CLAIM**：Onsager 必须交换源/测量算符 | 证明额外同序互易，完备 Kubo 验证 |
| BD06 | p3 的 3d 与 1a+1b+1c 局域系数 4→0 是 FE 开关 | **NO_GO_CLAIM / GROUP PASS**：不同轨道并非严格 FE 伙伴 | 真实非共轭双相、能量/磁序与局域读出 |
| BD07 | hidden/local multipole spin 是全新效应 | **NO_GO_CLAIM**：2018 磁多极矩、2020 hidden effect 与 2026 HSP 已覆盖 | 独立算符与物理机制、同公式对照 |
| BD08 | moiré 极化涡旋导出位移光电流涡旋是新机制 | **NO_GO_CLAIM**：Agarwal 等 2025 npj 论文直接用非阿贝尔量子几何做了 | 新观测量/机制，非仅更换材料 |
| BD09 | 滑移参数空间 Berry 曲率导致拓扑电荷泵是新机制 | **NO_GO_CLAIM**：Su–Lin 等 2020 滑移合成维拓扑电荷泵已经完成 | 有新定义及保护规律的 spin/site proper pump |
| BD10 | 全 80 LG/528 MLG 普通允许张量/自旋层群枚举即创新 | **NO_GO_CLAIM**：2025 80 LG、2026 448 SLG、40,087 SLG 预印本 | 真正新的局域算符/可测选择定则 |
| BD11 | MnPS3 R1 有可用 FE/电子自旋物性 | **BLOCKED_EVIDENCE**：两项 2h 超时、首个 SCF 不收敛、零离子步 | 新授权、原生数值收敛和科学物性验收 |
| BD12 | 真 FE 连接 × magnetic-site mode × SOC/频率导致可测响应新类 | **OPEN_HYPOTHESIS**：G3 首选，但无独立新机制证据 | 源可靠 FE 伴侣、同观测量/互易/读出排除与模型 |
| BD13 | nonsymmorphic 区边磁模可以被均匀 E 直接激发 | **BLOCKED_EVIDENCE**：q 动量不匹配；相位代数本身正确 | 实际近场/声子/超胞折叠的 q 激励与读出 |
| BD14 | moiré 具有新的光电响应涡旋一般规律 | **BLOCKED_NOVELTY，高重叠**：2025 shift current vortex、2025–26 polar meron/node 已有 | 和已有非阿贝尔 Berry shift-current 公式的明确差分 |

**恢复规则：** 不得把旧 NO_GO 改名成更响亮的“新量子几何”。只有源算符、对称操作、可测量量或有新证据的物理机制发生实质变化，才能开新 idea ID，并链接旧结论。**负结果不删，PASS 不自动升为新颖性。**

## 2026-10-22 新增核对的最近同行

- [Agarwal et al., *Shift photocurrent vortices from topological polarization textures*, npj 2D Materials and Applications 9, 43 (2025)](https://doi.org/10.1038/s41699-025-00552-6)：用非阿贝尔 Berry 连接、TB/DFT 联结 polar merons 与光电流涡旋；直接阻止“首次 moiré shift-current vortices”。
- [Su & Lin, *Topological sliding moiré heterostructure*, PRB 101, 041113(R) (2020)](https://doi.org/10.1103/PhysRevB.101.041113)：滑移合成维的第一/第二 Chern、拓扑泵浦。2020 年还有相关 sliding pump 理论。
- [Park et al., *Topology-dependent node dynamics under mechanical manipulation in moiré ferroelectrics*, Nat Commun 17, 6834 (2026)](https://doi.org/10.1038/s41467-026-73621-z)：四/六/八支 moiré FE 节点动力学。
- [Lu et al., *Nonvolatile electrical control of spin splitting via altermagnets with sliding ferroelectricity*, PRB 114, 094423 (2026)](https://doi.org/10.1103/fvdw-j1mc)：Ca(CoN)2 中 FE/层极化与 AHE/MOKE。
- [Pan et al., *Multistate manipulation of charge-spin conversion in two-dimensional ferroelectric bilayers*, PRB 114, 115135 (2026)](https://doi.org/10.1103/np2g-86h3)：三态 EE/SHE 转换及 bilayer PtBi2。
- [Xu & Yao, *Thouless pumps and universal geometry-induced drift velocity in multisliding quasiperiodic lattices*, PRB 113, 205149 (2026)](https://doi.org/10.1103/dp19-2r4z)：滑移/合成维几何泵机制。

**核查深度**：以上反例适用于标题中所写的宽泛机制；未逐条穷尽论文 Supplemental Materials。对更窄的新方程应标记 BLOCKED_NOVELTY，而非因为没搜到就声称原创。

## G3 下一步：用真实 FE 连接区分等价与非等价

(1) 从既有 BSF/stacking catalog 取源可信 ±P 伴侣，保持**同一个物理位点、算符和磁序约定**；(2) 完整幺正/反幺正的等价对照必须分别过群变换与 Onsager；(3) 真正候选仅是结构相关但完整 SOC/magnetic 不相关的对照，或磁性稳定子/站点模不同的状态；(4) 先给可测预测与同算符前人区别，才准许新材料任务。

本文件是**故事与科学机制版本控制**，不是科研结论的美化工具。R1 CLOSED；R2-P0 仍为只读状态，无新计算授权。
