# B+D G3 | 暂缓畴壁 DFT：科学投入回报、低算力群论搜索、堆叠响应梯度

**2026-10-22 | DOI/仓库审阅 + 一个新的群论研究门槛。**  
**科学状态：** ROI_DECISION；REGISTRY_GRADIENT_METHOD_PROTOTYPE；NO_DISCOVERY；NO_NEW_COMPUTE_AUTHORIZED。R1 two-job Slurm TIMEOUT + no first-step SCF convergence；R1 CLOSED。之前 R2-P0 仅授权只读证据审计、尚未收到本地执行结果。绝不可据此提交日本超算/WSL/para 任务。

## 1. 直接结论：现在的 magnetic/FE domain-wall DFT 不划算

计算真实滑移铁电畴壁一般需要比端点原胞大得多的**双畴超胞**（有时多个不同壁型、距离、方向、表面/周期边界），同时控制层间不匹配、长程弹性与极化、电荷屏蔽、磁性位点排序、U/交换、SOC、费米占据和受限/非受限离子弛豫。如果要求 \(\chi_{F,E}(\omega)\)，还需要壁结构上的独立响应顶点、稠密 k/有限 q、局域投影/磁力以及因果/规范验收。**超胞尺寸与总 CPU 核时取决于晶格、实际壁宽、软硬件，不应凭经验编造定量作业成本。**

R1 **连两种原胞级磁性 SCF 都没有收敛**，此时放大到壁 DFT 几乎肯定会增加故障表面积、资源需求和科研解释难度。并且即使发现壁处局域自旋密度或 Edelstein，也不能自动宣称是新效应，因为[2026 PRB 已研究导电铁电壁 SOC/Edelstein](https://doi.org/10.1103/y87v-34d8)。该 PRB 模型为 **BaTiO3 的导电 90° 壁**，不是完全相同的磁性滑移双层；它否决的是“铁电壁有 Edelstein 就是新效应”的**宽泛首创主张**，不是所有其他材料的狭义机制。

**成本/信息价值等级（定性、并非排队时间估算）：**

| 层次 | 可做什么 | 计算规模 | 证据价值与风险 | 决定 |
| --- | --- | --- | --- | --- |
| **S0：精确结构/群论** | 从 BSF 已有连接、轨道/位点、磁序，分类端点与邻域的 \(\chi(\delta)\) 等变条件、\(d\chi/d\delta\) | 几乎零 DFT | 快速淘汰错误/已知选择定则；没有响应振幅 | **立即优先** |
| **S1：端点/几个代表位移的原胞取样** | 在已验证结构上，独立对称性的刚性层间位移；少量 SCF / 响应点，磁序固定且源可靠 | **低～中**，但先修复 R1 的 SCF | 拟合局域响应/能量，误差可检；不代表 FE 墙存在 | **仅 S0 有非平凡预测后** |
| **S2：对称约束的 registry 势能/响应 Fourier map + FK/弹性连续模型** | 用有限几何/DFT或文献采样拟合 \(V(\delta),P(\delta),K(\delta)\)，再最小化 continuum wall | **中/无额外 DFT（若能复用源数据）** | 壁宽/取向/路径可扫描，但 continuum、局域近似、梯度及电子非局域失效须审 | **可选** |
| **S3：大磁性畴壁 DFT/响应/动力学** | 位点 SOC+U、带/局域力矩全频响应、壁迁移 | **非常高** | 有真实材料能力，但壁/电子/磁晶格耦合难分；已知先例多 | **暂缓** |

理论上 \(E_{\mathrm{cont}}[\mathbf u(\mathbf r)]=\int d^2r[\frac12(\nabla\mathbf u):C:(\nabla\mathbf u)+V(\delta_0+\mathbf u)-\mathbf E\cdot\mathbf P(\delta_0+\mathbf u)]\) 是**已有的** stacking-energy/elastic/field 组合，不是新模型。早在 [Bennett & Remez, npj 2D Mater Appl 2022](https://doi.org/10.1038/s41699-021-00281-6) 已用 stacking energy、弹性和外场构建 moiré/FE 畴及 FK 路径计算；[Ke et al., PRL 2025](https://doi.org/10.1103/jhlq-2dd7) 和 [Nano Lett 2026](https://doi.org/10.1021/acs.nanolett.5c05482) 直接研究真实壁运动与不同壁型。因此不能把 Fourier/FK/continuum 本身称创新。

**Existing code reality:** \`extension-to-BSF\` 目前只生成 source-hashed POSCAR + grid / collision checks；**不计算 PES，不跑 VASP，且默认同层单层克隆**。 \`materials-structure-core\` 是坐标/哈希等几何合约；\`group-theory-operations-toolkit\` 已有 80 LG、实际 Wyckoff 和 Seitz/stacking FE，但**未认证局域结构动力学/电子响应**。不要因此声称 ready to run walls。

## 2. 进一步撞车的 2024–2026 直接近邻

1. [Nasir & Atkinson, *Spin-orbit coupling and the Edelstein effect at conducting ferroelectric domain walls*, PRB **113**, 035425 (15 Jan 2026)](https://doi.org/10.1103/y87v-34d8): FE 壁 SOC 与电流自旋响应；不是磁性滑移双层，但与“壁处 Edelstein”同观测。
2. [Dong et al., *Current-Induced Spin and Orbital Polarization in Magnetic Sliding Ferroelectrics*, phys. status solidi RRL (2024)](https://doi.org/10.1002/pssr.202400062): 双层 CrI3 中电流诱导 spin/orbital、磁电多极矩随 FE 变化；任何普通“磁性滑移 FE 中 spin/orbital Edelstein”不可当新机制。
3. [*Slow-phonon control of spin Edelstein effect in Rashba d-wave altermagnets*, PRB **113**, 184403 (May 2026)](https://doi.org/10.1103/4txy-gy8f): EPC/Holstein 导致自旋 Edelstein 变化；不能泛称首次 phonon-assisted Edelstein。
4. [*Coherent phonon control beyond amplitude saturation in a sliding ferroelectric*, Nature Physics (online **21 Oct 2026**)](https://doi.org/10.1038/s41567-026-03491-1): Td-WTe2 层间滑移模 **0.24 THz**、SHG 及复杂非线性晶格位移；特别指出滑移与严格意义上的剪切不同，Td↔1T' 实际路径涉及两个模，不能靠一条平移坐标解释所有相。
5. [*Linearly polarized light-induced magnetism in bilayer MoSe2*, Nature Communications (Oct 2026)](https://doi.org/10.1038/s41467-026-77987-y): 层间 breathing 与 Se stretching 声子/载流子耦合产生瞬态磁化；不是我们 hypothetical \(E\to\) spin torque，但**phonon→stacked transient magnetism**已先例。
6. [*Domain wall dynamics in trilayer sliding ferroelectric gamma-InSe*, npj Computational Materials (2026)](https://doi.org/10.1038/s41524-026-02296-3): DFT + machine-learning MD 的多层 FE 壁动力学；说明实际壁模拟可能升级为 ML 势与大规模 MD，而不是单独原胞 VASP。
7. [Shen et al., *Magnetoelectric coupling dependent on ferroelectric switching paths in 2D perovskite multiferroics*, PRB **103**, L220406 (2021)](https://doi.org/10.1103/PhysRevB.103.L220406): **不同 FE 切换路径给不同磁化结局**早有理论；路径依赖本身不是新的基本物理。

**结论：** 「磁性滑移双层 + 畴壁/声子 + Edelstein + 磁模」各个关键词均有强先例；缺的是**同一量、同一机制与同一可测协议层面的非平凡差分**。不能靠叠加关键词构成新故事。

## 3. S0 新的低投入切口：响应对堆叠位移的导数，先群论后电子

传统堆叠 FE 连接理论多比较静态端点 \(\delta_\pm\) 的响应 \(K(\delta_\pm)\)。先考虑**一个源认证高对称 registry** \(\delta_*\) 附近的微小、均匀相对位移 \(\mathbf u\)；如果该操作确实属于 \(H_*\) 的**规定结构表示** \(D_u(h)\)，那么局域线性电场响应
\[
\delta\langle O\rangle=K(\delta) E
\]
可展开为
\[
K(\delta_*+\mathbf u)=K_0+\sum_i u_i\Gamma_i+O(|\mathbf u|^2).
\]
完整变换关系 \(K(h\delta)=D_O(h)K(\delta)D_E(h)^{-1}\) 给出
\[
D_O(h)\Gamma_jD_E(h)^{-1}
=\sum_i\Gamma_i[D_u(h)]_{ij}.
\]
因此
\[
K_0\in\mathrm{Hom}_{H_*}(V_E,V_O),\qquad
\Gamma\in\mathrm{Hom}_{H_*}(V_u\otimes V_E,V_O).
\]

**最小严谨例子：** 假定 registry 邻域具有面内镜面 \(m_x\)，\(E=(E_x,E_y)\) 分量镜面特征 \((-1,+1)\)、轴自旋 \(S_z\) 为 \(-1\)，**且实际 U 也确实按面内极向量变换** \(u_x\mapsto-u_x,u_y\mapsto u_y\)。那么固定 registry 仅 \(E_x\to S_z\) 允许，但一阶 \(u_xE_y\to S_z\) 和 \(u_yE_x\to S_z\) 都允许。即：
- 高对称点：\(K_{zy}(\delta_*)=0\)；
- 一阶滑移邻域：\(K_{zy}(\delta_*+u_x\hat x)\sim \gamma u_x\)，**系数 \(\gamma\) 不由群论决定，可能为零**。
- 有限频率时可把 \(\gamma\) 视为 \(\omega\) 依赖的电子导数，但**声子动力学与非绝热电子-声子顶点不能直接用静态导数代替**。

此结果是标准 **Neumann/tensor-product selection**, **不是新机制**；[2026 PRB phonon-modulated Edelstein](https://doi.org/10.1103/4txy-gy8f) 和 [2026 Nature Physics coherent sliding phonons](https://doi.org/10.1038/s41567-026-03491-1) 对一般机制已有实质先例。它的价值在于：给现有 **80 LG + layer exchange + Wyckoff/site** 工具一个**可验收的第一导数筛选指标**，无需真实壁和磁性 SCF。

**Critical caveat:** 相对层滑移 \(\mathbf u=r_{\mathrm{top}}-r_{\mathrm{bottom}}\) 在**交换层次**的镜/反演操作下不一定按普通极矢量变换。必须从实际 site+layer permutation 生成 \(D_u(h)\)，不能统一套用普通 \(E\) 的表示。 [镜面显式 U 负对照脚本（Draft PR #32）](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/registry_gradient_gate.py) 仅测试上述**假设表示**的代数，不是一个 FE 材料计算，更不是新的 phonon-spectral 量。

## 4. 不扩建玩具 Hamiltonian：G3 的真正任务清单

**G3.1 零计算**：从 BSF 现有连接列表找两个具体、站点与层交换都可追踪的**实际堆叠路径或局部注册态**，将 \(\mathrm{Hom}_{H_*}(E,O)\) 与 \(\mathrm{Hom}_{H_*}(U\otimes E,O)\) 由真 Seitz 运算自动求出；一例无层交换负对照，一例含 layer exchange 正确处理 \(D_u\)。还必须说明读出 \(\hat O\)：宏观/层/磁子格、spin/orbital、场或电流驱动，不偷换。**如果只给出已知 phonon-assisted Edelstein 的 ordinary tensor activation，方法 PASS，但 NOVELTY NO-GO。**

**G3.2 低计算**：仅在 G3.1 对某一真正缺失的可观测预测有价值时，选择公开可靠单层/堆叠源，用 3–6 个**对称不等价原胞**对 \(V(\delta),K(\delta)\) 局部取样及误差外推；这只是候选性的极简采样范围，**不是规定的 DFT job 数或可靠 PES 拟合所需点数**。有了样本再决定需多少取点，且先解决既有磁性 SCF 问题。若用文献公开 PES 可做到零新增 DFT。

**G3.3 高成本，仅后续有新物理且用户重新批准**：结构连续模型/FK → 如有必要的磁性实空间有限壁超胞 → full causal site torque/电子+晶格。严守独立验证；不能从少数刚性堆叠原胞样本冒充真实弛豫壁、磁振子谱或开关电流。

## 5. 决策 / 新颖性停止线

- **已证 NO-GO：** 仅是域的坐标重定向、Onsager 互易、hidden local spin、C3 禁戒、滑移泵/位移光电流涡旋、普通 phonon-modulated Edelstein、铁电壁 Edelstein、普通 FE 改变磁矩/磁矩动力学。
- **保留的科研工具：** 同源 registry 响应梯度的 Seitz/site representation 计算；FE partner mapping 与频率 Onsager 交叉响应；诱导表示/Frobenius 局域响应；严格电子自旋/交换/SOC continuity gate。
- **仅存的论文候选问题：** 是否存在**源可信的真实局域结构状态+真实测量协议**，使响应梯度/模式耦合出现*已有磁点群张量、FE 滑移及 electron-phonon 理论都无法直接给出的强制约束*？目前没有给出这种构型/响应，因此应标 **OPEN_HYPOTHESIS / NOVELTY_UNRESOLVED** 而非 PHYSICS_CANDIDATE。

**Immediate decision:** ROI_PRIORITY = SYMMETRY_SITE_REGISTRY_GRADIENT. DOMAIN_WALL_DFT = DEFERRED. No material compute execution authorization.
