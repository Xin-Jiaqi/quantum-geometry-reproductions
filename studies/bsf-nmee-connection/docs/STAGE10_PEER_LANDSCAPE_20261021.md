# Stage 10 | BSF × NMEE 同行工作系统审阅与创新性定位（2026-10-21）

> **性质：原创研究的独立文献调研。** 不属于 Jia 等 PRL 论文复现结果，不代表已完成 WSe₂ 的 SOC/NMEE 计算。Stage8i 的两项 SOC 作业由外部 DeepSeek harness 负责，本轮不查询、不修改、不追加超算任务。证据主要来自 APS、Nature、arXiv 原始出版页/官方摘要，以及已有项目原始研究笔记；未获取全文的文章不把摘要以外的细节写作“确认”。文献覆盖截至 **2026-10-21**，并非声称穷尽所有预印本。

## 0. 首先锁定我们自己的科学问题

**明确区分三种响应**，避免同行对比时错位：

1. 普通/层 Edelstein 效应（EE/LEE）：电场或电流造成的**非平衡、通常费米面主导**的自旋积累（存在内禀带间部分时需进一步按作者具体定义区分），符号与弛豫参数依赖物理协议。
2. Jia 等人的 nonlinear magnetoelectric Edelstein（NMEE）：**电场 × 磁场**共同诱导自旋磁化，原文用 spin-space Berry curvature dipole / metric dipole 区分内禀/外禀，并明确声明**内禀项可在绝缘体存在**。项目中对原文最终公式的驱动 DC/静态极限仍有独立复现分歧。
3. **我们当前自己已在模型中验证的对象**：无耗散、带隙、\(T=0\)、均匀 E 与**spin-only** 纵向 Zeeman 场的**平衡/静态**混合曲率导数
\[
\Gamma^{\mathrm{stat},\pm}_{iz;j}=
\left.\partial_{B_z}\int_{\rm BZ}\frac{d^2k}{(2\pi)^2}
i\,\mathrm{Tr}P[\partial_{\lambda_i^\pm}P,\partial_{k_j}P]\right|_{B=0},
\quad i,j=x,y.
\]
其物理表达属于自旋源参数 / 电极化的混合静态几何，**不是经证实等于** Jia 的驱动 DC 内禀表达式，也不是有真实传导电流的 LEE。构造真实层算符时必须确保其在同一完整自旋电子 Hilbert 空间内的规范协变性。

在实际 3R-WSe₂ 来源结构上，目前**只**做过 3 个 non-SOC fixed-xy/z-only PBE+D3(BJ) 受限离子弛豫（Stage8g/h），已有远端原始输出的 DeepSeek harness 审计。**滑移中点真实弛豫原子坐标仍保留 \(\{m_z|\mathbf t\}\) 层交换滑移镜面**，所以均匀总响应对称禁止；层奇响应对称允许但**尚无实际非零数值**。SOC 电子结构和 \(\Gamma\) 还未从本项目获得。

## 1. 核心同行矩阵：哪里竞争最直接？

| 编号 | 论文与证据 | 作者实际完成了什么（不要移花接木） | 与我们的重合/差别 | 对论文立题的影响 |
|---|---|---|---|---|
| **P01** | **Jia et al., PRL 136, 226302 (2026)** [DOI](https://doi.org/10.1103/r3x2-kb93)，APS 原文摘要，A 级 | 提出 \(\mathbf E\mathbf B\) 协同诱导磁化的内禀/外禀 NMEE，spin-space Berry curvature/quantum metric dipoles；内禀可存在于绝缘体（Kane–Mele）。 | **物理刺激与 spin-quantum geometry 最近**；已有 *insulator + intrinsic* 主张。我们的新空间只可能在**滑移、层奇偶、真实材料、不同静态算符和机制**。 | **最高级直接竞争**；禁止“首次提出绝缘体内禀 NMEE”。必须写清静态/DC 关系。 |
| **P02** | **Zhou et al., PRL 137, 076401 (2026), Layer Edelstein Effect** [DOI](https://doi.org/10.1103/qhvs-pfkq)，APS，A | 普适非磁性双层**电流驱动**层反向自旋积累；用堆叠与对称性做判定，DFT 示例 MoSSe、MoTe₂、WTe₂。 | **层间抵消/可通过场揭示的概念直接重合**；但我们的 E×B 静态双轴量与它的**电流驱动 axial/polar 二阶**不同。 | 不能声称“首次提出 hidden layer Edelstein / 层反向自旋”或“首次系统堆叠判据”。 |
| **P03** | **Li et al., PRB 114, L111410 (2026), Layer Nonlinear Hall Effect** [DOI](https://doi.org/10.1103/2gyh-y52n)，APS，A | 2/3 阶**层奇非线性霍尔电流**，层交换导致宏观抵消、\(E_z\) 可揭示隐藏响应，**80 layer groups** 对称性分类，WTe₂/ReS₂ 计算。 | **隐藏层奇 + 总量为零 + 对称性分类**结构与我们的 glide 叙事高度同构；但是响应是非线性霍尔**电流**，不是层自旋对 \(B_zE_\parallel\) 的静态敏感度。 | 本轮**新增的最高风险邻近竞争**；“零总量但非零局域量”的一般构想不新。 |
| **P04** | **Pan/Jiang/Fabian, PRB 114, 115135 (2026)** [DOI](https://doi.org/10.1103/np2g-86h3)，APS，A | 铁电**金属 PtBi₂ 双层**的 EE/SHE 电荷自旋转换：极化翻转 EE 符号，反平行态将 EE 压为零，DFT 建立三态逻辑。 | **铁电翻转 + EE 开/关**近似直接重合；但我们不是金属费米面 EE，不应将其器件逻辑宣称为自己首次提出。 | “ferroelectric switching of Edelstein”的新意明显缩水。 |
| **P05** | **Lin et al., PRB 100, 155408 (2019)** [DOI](https://doi.org/10.1103/PhysRevB.100.155408)，APS，A | MX₂ (Mo/W, S/Se/Te) 双层电子型铁电耦合 Rashba 自旋分裂，铁电反转翻转 spin texture，部分应变下分裂显著。 | **材料家族与铁电-Rashba 直接交叉**；但它研究 spin texture，不是我们 \(\partial_{B_z}\mathcal F_{\lambda k}\) 的响应。 | “WSe₂ 滑移铁电产生 Rashba / 翻转自旋纹理”不新，必须提供真实的**不同响应算符**。 |
| **P06** | **Pan/Fabian, PRB 113, 165140 (2026)** [DOI](https://doi.org/10.1103/j5s5-m7j5)，APS，A | 金属铁电 PtBi₂ 的可调 Edelstein 效应，重点在无需额外掺杂的导电自旋转换。 | 与 P04 共同形成铁电自旋转换主线；我们研究静态绝缘体，应主动解释这一区分。 | “polar material with switchable EE”不可作为主稿主张。 |
| **P07** | **Xiang/Jin/Wang, PRL 135, 146303 (2025)** [DOI](https://doi.org/10.1103/14mp-263q)，APS，A | **spin quantum geometry / spin QGT / Zeeman QGT** 统一描述 spin current 与 spin magnetization，进一步提出非线性 Drude 自旋响应。 | **我们的量子几何算符语言并非首次创立**；它是 Jia NMEE 重要理论前驱。 | 不能声称“首创 spin-space QGT/mixed Berry 理论”；应说明利用并拓展在哪里。 |
| **P08** | **Huang et al., PRB 113, 035144 (2026)** [DOI](https://doi.org/10.1103/5v7q-7x7y)，APS，A | mixed quantum metric 产生 nonlinear Edelstein / intrinsic antidamping SOT，讨论注入、shift、anomalous 带间机制。 | 是 quantum-geometric spin conversion 的近邻，但其扰动与模型场景不同。 | “非线性自旋量子几何来源于远端带间相干”也是广泛概念，不能笼统算首创。 |
| **P09** | **Das/Agarwal, arXiv:2604.20249 (2026)** [arXiv](https://arxiv.org/abs/2604.20249)，**预印本、B 级** | 在叠层中电场/磁场驱动 layer–orbital quantum geometry，预测 intrinsic magnetoelectric **Hall** response，rhombohedral pentalayer graphene 例子，绝缘区仍可有限。 | \(\mathbf E\mathbf B\)+layer+mixed geometry **思想接近**，但测量对象是电荷**Hall conductivity**，磁耦合是**轨道**，未必需 SOC；不是我们的 \(\Gamma_{iz;j}^{\rm stat}\)。 | 是概念性直接邻居；区分 spin vs orbital source、\(q/\omega\) 极限及可观测量。**不要标记已发 Science/PRL**。 |
| **P10** | **Sliding Ferroelectric Control of Unconventional Magnetism in Stacked Bilayers**, PRL 135, 056801 (2025) [DOI](https://doi.org/10.1103/dmzg-ck2t)，APS，A | **connection operator + sliding FE + magnetic spin splitting/AHE**，系统 80 层群分类，AgF₂/Fe₂MoSe₄ DFT。 | 堆叠镜像/连接算符推导**方法论直接已被覆盖**；磁性 altermagnet/ferrimagnet、非相同自旋物理。 | 我们的群论程序是可靠工具，但不是独立物理创新的终点。 |
| **P11** | **Li et al., Nat. Commun. 17, 2457 (2026)** [DOI](https://doi.org/10.1038/s41467-026-68854-x)，Nature 全文，A | 3R-WSe₂ AB/BA **实验极化畴**、电场诱导 band switching、导带/价带边缘 \(\Lambda/K\)，非泛泛材料假设。 | 直接支持选材和真实域存在，**未直接给我们的静态 NMEE 数值**。 | WSe₂ 合适作为可实证材料，但不能把其 FE/开关声称为新发现；采样需全 2D BZ。 |
| **P12** | **Non-relativistic torque and Edelstein effect in non-collinear magnets**, Nat. Commun. (2024) [article](https://doi.org/10.1038/s41467-024-51565-6)，Nature，A | 线性 spin–E 响应拆分 intraband / interband，包括与 spin–momentum Berry 几何相关的清洁极限；文章讨论绝缘隙内依然存在带间响应。 | 提醒**“Edelstein 都只发生在金属/费米面”并不准确**。我们的独特性不只是“insulating response”。 | 必须比较作者实际的散射、带间磁响应定义，避免过度强化静态/输运二分。 |
| **P13** | **Mixed Weyl semimetals and low-dissipation magnetization control in insulators by spin–orbit torques**, Nat. Commun. (2017) [article](https://pmc.ncbi.nlm.nih.gov/articles/PMC5684220/)，全文，A | \(k\)–磁化方向的 **mixed Berry curvature** 与绝缘体内低耗散 SOT。 | **混合 Berry 曲率与绝缘体响应都有先例**。我们不能包装为首次发明参数空间 Berry curvature。 | 指向最接近的早期 Berry 自旋响应理论链，正文必须引用。 |
| **P14** | **Freimuth/Blügel/Mokrousov, Direct and inverse spin-orbit torques**, PRB 92, 064415 (2015) [DOI](https://doi.org/10.1103/PhysRevB.92.064415)，发表 PDF，A | 发展 direct/inverse SOT 的 Kubo 响应关系及微观 DFT 公式，包含磁性双层。 | 我们的静态 spin-source/MEE 是其方法后继，但响应协议、磁性背景和多场阶次需要具体区分。 | 普适 Kubo 带间几何/数值路径不是全新；只能凭新的层分辨响应、材料预测与算符验收形成独立贡献。 |

**证据等级：A = 出版社官方文章、摘要或可读全文直接支持此行所列事实；B = arXiv 预印本摘要与正文，未在本轮确认同行评审。** 对只读摘要的 A 级来源，**A 只表示对摘要明确陈述部分有高可信度**，不是对所有公式、方法和与我们细节的全文级审计。本报告对不存在的对照、未核实的公式号或精确能量不作推断。


## 1a. 用项目中已经精读的三篇原始 Markdown 再次校正定位

这一节补充前面按出版社摘要整理的矩阵，依据项目上传的现有完整阅读笔记，而不将笔记中的解释自动当成作者原文定理：

- **《Layer Edelstein Effect.md》**：正文明确以 \(\delta S_j=\chi_{ji}^{(j)}j_i\)（charge current）为输出；其 bilayer \(k\cdot p\) 以 layer pseudospin \(\tau\)、interlayer hopping 和外部 \(E_z\) detuning 为关键变量，SM 包含 80 layer groups 的 Type-I/II 分析。**不能把其 \(\chi(j)\) 直接视为我们 \(\Gamma^{\mathrm{stat}}(E,B)\)**，即便出现相似层奇偶性。
- **《Sliding Ferroelectric Control of Unconventional Magnetism in Stacked Bilayers.md》**：重点是 magnetic connection operator 对非相对论 spin splitting、SOC 下 AHE 的连接关系，以及 AM/CFiM 具体材料。它证明了堆叠连接算符已经是完整的同行研究路线；我们的新价值应集中在**不同物理响应的算符/实际量化**，而非重做 80 层群枚举。
- **《Intrinsic Magnetoelectric Hall Effect from Layer-Orbital Quantum Geometry.md》**：笔记清楚区分**三个**场：面内驱动 \(E_x\)、面外 gate \(\mathcal E_z\) 和面外轨道磁场 \(B_z\)，可观测量是
  \[
  j_y=\chi^{\mathrm{in}}_{yx;zz}\,E_x\mathcal E_z B_z.
  \]
  而我们是 \(\delta S_x=\Gamma^{\mathrm{stat}}_{xz;x}E_xB_z\)（在本身无净电流的绝缘体静态物理定义里）。因此 P09 不能被写成“已经研究同一张量”；它的 **layer dipole–orbital moment** mixed geometry 仍是概念上的重要参照。
  
这些上传笔记本身也可能包含教学习惯性的简化，因此**具体对比仍须返回对应文章与 SM**。原始上传笔记保留在 Project 中，Stage10 不复制整个笔记或改写旧文件。

## 2. 文献审阅中最值得纠正的四个“新意幻觉”

1. **“第一个层奇效应”** —— 错。LEE (P02) 和 LNHE (P03) 已大量覆盖。我们最多研究**新的具体 \((S_\parallel,B_z,E_\parallel)\) 静态三阶张量**如何分解成总/层差，并且必须给出非零实际材料系数及域奇偶关系。
2. **“第一个铁电翻转的 Edelstein”** —— 错。Lin (P05)、Pan (P04/P06) 已覆盖。我们的 AB/BA \(\gamma\) 变号若仅是对称性，也不够一个高水平主贡献。
3. **“首次绝缘体内禀 spin quantum geometry”** —— 错。Jia (P01) 明确宣称 insulating intrinsic NMEE，Xiang (P07) 已建立 spin QGT，2017 mixed Berry torkance (P13) 更早。准确的卖点必须是**新层/堆叠参数与自旋源耦合的限定机制**。
4. **“WSe₂ SOC 很强所以肯定能算出大 \(\Gamma_{xz;x}\)”** —— 错。确切的 no-go 是若 \([H,S_z]=0\)，横向自旋对 \(E_j\) 的交叉带间响应为 0；自旋非守恒、真实算符矩阵元/跨隙通道不建立前，不许将 huge Ising splitting 写作 spin mixing 证据。

## 3. 真正有价值且目前尚未被以上文献**直接给出**的组合

### A. **第一候选故事：静态磁电层隐藏自旋几何在真实 3R-WSe₂ 中的有限实现**

在合理定义的层自旋投影 \(\hat S^{(+/-)}\) 下，计算
\[
M^{(\pm)}_{ij}=\Gamma^{\mathrm{stat},\pm}_{i z;j},\quad i,j=x,y,
\]
并验证三个同时成立的条件：
- polar \(C_{3v}\) AB：\(M^{(+)}=\gamma I_2\)，BA：\(-\gamma I_2\)；数值 \(\gamma\neq0\)，且需要热力学静态定义及量纲一致。
- 在原始 DFT 实测保留的 half-slide glide：\(M^{(+)}=0\) 是严格空间选择定则；**\(M^{(-)}\neq0\) 必须真实计算证实**，不仅是对称性允许。
- 比较真实空间 top/bottom spin density 的分区方式，检验 \(M^-_i\) 不仅由 Wannier 轨道标号决定；给出可观测/可外部耦合的层敏感性。

**创新评价：条件性中高**。层交换/隐藏概念和 spin quantum geometry 都有前人，但它们在这一**响应定义、绝缘铁电翻转、真实材料、镜面中点**上的综合定量结果当前矩阵未发现已被直接报道。**“未找到”不等于证明首次**。

### B. **第二候选故事：微观自旋混合与 remote-band interference**

展示为何近 K/K′ 的 Ising splitting 不是可用预测因子；把真实 SOC spinor 的 \(\langle v|S_x|c\rangle\langle c|v_x|v\rangle\) 分解到 \(\Gamma,\Lambda,K\) 热区和不同能带窗口，说明何种层间杂化与 SOC 产生有效的静态响应。如果 full 8-band (或 DFT 全窗口) 与少带近似符号相反，且不来自任意 gauge/投影，就有更强的**具体材料微观新意**。

**创新评价：条件性中高**。量子几何/远端能带一般重要性并不新，值得发表的应是此响应的特定**自旋矩阵元与层堆叠干涉机制**，外加 gauge-covariant Wilson/Kubo 验证。

### C. **理论高风险高回报：精确定义静态 \(\Gamma\) 与 Jia driven-DC 的极限关系**

逐项检验 diagonal density-matrix / Fermi-sea / contact term、\(\omega\to0\) 与散射率 \(\eta\to0\) 的顺序、Zeeman spin-only versus physical orbital \(B\) 耦合。在同一八带模型内比较：
- 混合占据投影 \(\partial_{B_z}\mathcal F_{\lambda_x k_x}\)；
- clean driving Kubo / 密度矩阵完整表达；
- finite broadening + occupations 的极限。

如果两者不等价，**改名为滑移调控的“平衡态磁电自旋敏感度”**而不是强行与 NMEE 驱动 DC 混淆。若能得出明确条件下的恒等关系，可形成独立理论价值。但需小心 P01/P07 与既有 interband SOT 研究很接近。

**创新评价：潜在高但风险高；不该拖延已发出的 SOC 物理可行性检查。**

### D. **最容易被抢走新意的方向**

“普适层群分类 + switchable hidden response + 极化反向翻号”的宽泛论述最危险。P02、P03、P10 各自已给出类似结构的层群框架。我们的光滑路径群论更适合作为**材料-specific 可证伪的控制变量**，而不是另起一篇“普适隐藏层效应分类”的主张。

## 4. 建议的同行比较图及论文必须补充的数据

| 建议图 | 目的 | 当前证据 | 通过门槛 |
|---|---|---|---|
| **Fig. 1：响应/对称性物理地图** | 区分 EE、LEE、LNHE、NMEE、我们的静态 \(B_z E_\parallel\) spin response；明确**电流 vs 磁化、掺杂/带隙、odd/even layer、T/B** | 可立即根据模型和原文绘制 | 不给超出原文的静态/DC 等价宣称；指出 P01/P02/P03/P04 各自已知部分 |
| **Fig. 2：真实结构 + 受限 PES** | AB / s=.25 / glide midpoint 真实 DFT 对称性，能量相对 AB、域反演路径 | **Stage8h 已有三点及群论验证** | 标注 z-only、non-SOC、有限 k 网格，不写 CI-NEB barrier |
| **Fig. 3：全 2D BZ SOC gap 和 spin-mixing diagnostics** | 排除 pure Ising \([H,S_z]=0\) 的微观零响应风险 | **Stage8i 已安排/等待 harness** | 不仅画 \(\langle S_{x,y,z}\rangle\)，还须分析非简并 spinors 的**跨带**自旋/速度矩阵元 |
| **Fig. 4：真正的总 vs 层奇 \(\Gamma\)** | 证明 AB/BA 域奇偶，半滑移节点的 \(M^+=0\)、\(M^-\neq0\) | **尚无真实材料响应** | full BZ、物理 layer projector、绝缘占据、规范协变方法和结果单位 |
| **Fig. 5：微观来源 + 论文对比 benchmark** | P01 公式极限、P07 spin QGT、远端带抵消和矩阵元热区 | 模型已验证部分，真实 DFT 无 | k/能带分解、Wannier 窗口收敛、Wilson vs Kubo、真实 static/DC 边界 |

**建议尽早准备 Fig. 1（理论/同行定位）与 Fig. 2（真实 DFT 结构）**，但若 Fig. 3 的 SOC 自旋混合非常弱，应立刻重估 WSe₂ 的可行性，避免将 Fig. 4/5 的昂贵计算变成无效工作。

## 5. 审稿人最可能提出的五个实质问题（及回应所需证据）

1. **Jia (2026) 已经有绝缘体内禀 NMEE，你新增了什么？** 要用一张“作者的驱动 DC vs 我们平衡/静态 + 额外层/滑移定义”的**公式/极限表**回答，而不是泛称机制不同。
2. **Zhou LEE 和 Li LNHE 已建立 hidden layer 和 80 层群规则，为什么还算新？** 把非平衡电流/霍尔输运与 **\(E B\) 静态 layer spin** 的三阶响应张量、磁场奇偶和角向选择定则分开，**再给实际 \(\Gamma^-\neq0\) 数值**。
3. **WSe₂ Ising SOC 为什么产生面内 \(S_x\)？** 先回答 SOC spinor 的真正 off-diagonal spin/velocity 机制，并测量对更大 Wannier 窗口的收敛，而不是展示较大的 K 谷 spin splitting。
4. **层差 spin 是哪个可观测量？** 用空间 PAW/real-space layer windows、层间尾部敏感性和与物理场耦合/检测方式限定，不凭 Wannier 标号。
5. **24.117 meV/cell 就是切换势垒吗？** 不是；只是受限 z-only path 的 half-slide 相对能量。得有收敛 PES、完全离子弛豫、最小能量路径/CI-NEB 和切换动力学才可报告 barrier 或稳定第三态。

## 6. 研究流程：与用户明确提出的三条主线对齐

\`\`\`text
IDEA / 物理问题
    Stage6–7: 滑移 + 量子几何 + 总/层奇张量 + Sz-no-go
          ↓
TEST / 可证伪的材料验证
    Stage8e–8h: WSe2 真实数据、受限弛豫、半滑移 glide（已过）
    Stage8i: SOC gap + spin mixing（harness 执行中/待回报；本次不追问或新提交）
    StageNext: 仅在 SOC 可行时做物理层投影 + 全 BZ Γ
          ↕
PEERS / 始终校准创新与验证门槛
    本 Stage10: PRL NMEE / PRL LEE / PRB LNHE / FE Edelstein
    + Spin QGT / mixed Berry SOT / layer-orbital / experimental WSe2
\`\`\`

**研究建议：** 采用一条**更窄但可证伪**的中心命题：
> 真实 3R-WSe₂ 中，若自旋非守恒带间耦合足够大，滑移可能在同一材料内重分配平衡态 \(B_z E_\parallel\) 自旋敏感度的总量与层差量，产生一个真实结构保护的总响应节点、同时允许且实际可观测的层奇信号；量子几何与远端带间干涉决定它是否显著。

这是未来要检验的**假设**，不是已得到的 WSe₂ NMEE 结果。若后续 SOC/矩阵元检验否决假设，应优先调整材料/自旋分量，不应靠不断增加滑移 DFT 来保住故事。

## 7. 核心参考文献与来源策略

**优先仔细阅读原文（直接竞争）**：P01 [Jia](https://doi.org/10.1103/r3x2-kb93), P02 [Zhou LEE](https://doi.org/10.1103/qhvs-pfkq), P03 [Li LNHE](https://doi.org/10.1103/2gyh-y52n), P04 [Pan FE EE/SHE](https://doi.org/10.1103/np2g-86h3), P05 [Lin FE Rashba](https://doi.org/10.1103/PhysRevB.100.155408).

**方法原点**：P07 [Spin Quantum Geometry](https://doi.org/10.1103/14mp-263q), P08 [Metric NLEE/SOT](https://doi.org/10.1103/5v7q-7x7y), P12 [Nonrelativistic torque](https://doi.org/10.1038/s41467-024-51565-6), P13 [mixed Berry Weyl/SOT](https://pmc.ncbi.nlm.nih.gov/articles/PMC5684220/).

**滑移/材料**：P10 [Sliding FE / altermagnets](https://doi.org/10.1103/dmzg-ck2t), P11 [WSe2 experimental](https://doi.org/10.1038/s41467-026-68854-x).

**相邻预印本**：P09 [Das/Agarwal layer-orbital](https://arxiv.org/abs/2604.20249), 注意**截至本轮只确认 arXiv 预印本，不能冒称 Science 正式发表**。

> **可复现范围声明：** 本轮核验出版社原文标题/作者/日期/摘要和部分可读正文，并抽查了本 Project 中《Layer Edelstein Effect》《Sliding Ferroelectric Control》《Intrinsic Magnetoelectric Hall Effect》的完整中文精读笔记；并非逐式推导每篇完整 supplemental material。特别是 P01、P02、P03 与 P04 的有些公式、数值归一化、碰撞算符和 layer projector 实现，需要下一轮得到全文和 SM 后逐式细审。在没有这些证据前，任何“已经证明我们的 Eq. X 与同行 Eq. Y 不同”的断言都不成立。新增文献是否与我们有**完全相同的** \(\Gamma_{iz;j}\) 和 layer projection，仍属开放检索目标。
