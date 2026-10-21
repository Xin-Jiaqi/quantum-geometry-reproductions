# BSF × IGMC：同行工作复核、原创性边界与计算判别式

**日期：2026-10-21；性质：独立同行研究审查（source-audited research note），不是 DFT 结果；适用仓库：本目录 BSF×IGMC 独立项目。**

**更新目的：** 继承既有八带模型、AgF₂ 结构与磁连接预检，以及已提交一次的 G1 作业，不重启项目。将研究拆为 **I. Idea 提出；II. 物理/材料测试；III. 同行工作调研**，并以 III 反过来审查 I、决定 II 的最小下一步。

本文件对之前 2026-10-21 长篇审查报告的同行调查作二次核验，明确可核实出处、最相近竞争主张和必要修订。此处的数值分成（a）原论文作者结果，（b）仓库冻结模型结果，（c）尚无证据的预测，绝不相互替代。文末 DOI 链接经在线交叉核对；文献原文的具体量值仅在能从正文核查时使用。

## 0. 决策摘要

1. **值得继续，但以材料是否允许真实非零响应为核心科学判断。** 八带模型已有独立 Kubo / finite-Zeeman occupied Wilson 的受控符合，且无手工滑移奇 onsite 项；但其能量扫描不支持所选 ±0.19 为双稳态。八带模型是「可能性证据」，不是 ferroelectric 材料证据。
2. **原创性不能定位为「首次磁位移电流」「首次滑移调控磁性」「首次磁电极化 P∝B」或「首次连接对称性的奇偶规则」。** 这些基础分别已有 Xiang (2025)、Ghorai (2026)、Zhu (2025)、Peng (2026)、经典 magnetoelectricity 文献。
3. **强一些的潜在科学问题：** 真实有隙的滑移磁性双态中，层间杂化如何改变 mixed k–B 自旋/轨道量子几何，怎样在满足单态磁群的前提下产生可切换、可校准的低频 α^PB，并在极化态/零面外极化参照之间识别真正的机制。
4. **磁连接奇偶性是必要约束，不是新的普适响应律。** 自旋与轨道对同一物理磁场的静态平衡磁电贡献原则上具有相同的张量变换行为，不能把单纯的 m_z odd、m'_z even 宣传成只有自旋量子几何才有的新效应。
5. **AgF₂ 仍优先，只读核实 G1 后再分批请求计算批准。** CuF₂ 是有文献双稳态证据的明确备选，不因同行调研就擅自开启第二材料计算。

## 1. I｜Idea 的正确物理定位

### 1.1 最小可证伪命题

固定一个化学体系、同一个外参（包括应变、U/vdW、磁畴定义），存在由结构滑移联系的两个局部稳定有隙磁性构型 B1、B2：

\[
B_1\xrightarrow{\mathrm{slide}}B_2,\quad P_z(B_2)=-P_z(B_1),
\quad \alpha_{ij}^{PB}=\left.\partial P_i/\partial B_j\right|_{B=0}.
\]

在低频、静态、实数极限，因 P 为 T 偶极矢量、B 为 T 奇轴矢量，连接操作 (R,T^\theta) 给出

\[
\alpha^{PB}_{2}=(-1)^\theta\det(R)\;R\,\alpha^{PB}_{1}\,R^{\mathsf T}.
\]

该式是**普通静态磁电响应变换律**，不是从 Xiang 有争议的完整复频率 Eq.(14) 直接推广而来。对 \(R(m_z)=\mathrm{diag}(1,1,-1)\)，纵向 \(\alpha_{xx},\alpha_{yy}\) 在 unitary \(m_z\) 连接下 odd、antiunitary \(m'_z\) 连接下 even；AHE 的域间奇偶性相反。这里的奇偶性仅在对应映射、磁态和坐标约定成立时准确，不能推出响应数值非零。

低频 B(t)=B0 cos(ωt) 导致 \(j_i^D(t)=-\omega\alpha_{ij}^{PB}B_{0j}\sin(\omega t)+O(\omega^2)\)。**是 AC 极化电流，不是恒定 DC 传输；完整周期净电荷为零。** 需要与轨道 Zeeman 贡献、Faraday 感生电场和器件电容背景区分。

### 1.2 提升原创性的额外条件

仅证明磁连接奇偶性，属于将已知一般规则用于新系数的对称性例证。更重要的实质性增量至少包含：

- 两个离子弛豫的反向极化局部极小值及滑移路径，磁构型经自洽和 SOC 验证；
- 与 Wilson/finite-field 极化导数交叉验证的非零 α^PB（首先 spin-only），含绝缘隙、gauge、单位、k 收敛；
- 由轨道投影、层间跃迁和 k-resolved 对照支撑的**显微机制**，而不是只给两个符号；
- 非 FE 中间态或强制零 Pz 对照、其他磁畴/连接对照，检验 α 与静态 Pz 是否必然同步；
- 对 orbital contribution 和 induction background 的可信边界讨论。若暂不能计算 orbital，则标题/摘要必须明确 spin-Zeeman contribution。

**特别注意：** AgF₂ Type-I（0% 应变、Néel 反平行、幺正 m_z）与 Type-II（4% 双轴应变、Néel 平行、反幺正 m'_z）是两条独立的外参分支，不是无需应变调控即可原位切换的两个磁畴。

## 2. III｜最相关论文：事实、竞争程度、直接改变的计算

按与我们科学主张的直接关系分级，避免为了数量收集量子几何综述。

| 层级 | 文献（作者、年份、DOI） | 已完成的同领域物理 | 对本项目真正构成的约束 / 下一项行动 |
| --- | --- | --- | --- |
| **A：直接基础** | Xiang, Jia, Xu, Qiao, Wang, 2025，*Intrinsic Gyrotropic Magnetic Current from Zeeman Quantum Geometry*，PRL 134, 116301；[10.1103/PhysRevLett.134.116301](https://doi.org/10.1103/PhysRevLett.134.116301) | Zeeman 几何、传导 C / 位移 D 分解，允许绝缘体内禀 AC 电流 | 切勿声称 IGMC 或 Zeeman QGT 为我方发明；要验证 \(\partial P/\partial B\)，同时保留原复现 Fig.2(d) 幅值与一般反幺正复频率未决 |
| **A：直接材料/对称先例** | Zhu *et al.*, 2025，*Sliding Ferroelectric Control of Unconventional Magnetism in Stacked Bilayers*，PRL 135, 056801；[10.1103/dmzg-ck2t](https://doi.org/10.1103/dmzg-ck2t) | 80 层群 connection operators，AgF₂/Fe₂MoSe₄ 双层极化、磁性和 AHE 切换 | 既有堆叠与磁连接分类不能重新包装成新理论；用作者 AgF₂ 结果作复现目标，不冒充作者坐标已经公开 |
| **A：最直接的新响应竞争** | Ghorai, Sarkar, Agarwal, 2026，*Intrinsic gyrotropic magnetic current of orbital origin*，PRB **2026-09-02 已接收**；[10.1103/xr5v-gkv5](https://doi.org/10.1103/xr5v-gkv5)，[arXiv:2601.04787](https://arxiv.org/abs/2601.04787) | 密度矩阵导出 **orbital** IGMC，CuMnAs 中由 Néel reversal 翻转信号；作者在 CuMnAs 案例报告 orbital 可接近 spin 的十倍 | 当前 spin-only 不能等同物理总响应；真正创新需要结构可切换机制、实际材料信号和轨道渠道控制。CuMnAs 的大小比例不可机械移植到 AgF₂ |
| **A：最强相邻材料竞争** | Peng *et al.*, 2026，*Sliding Ferroelectricity Driven Spin-Layertronics in Altermagnetic Multilayers*，Advanced Science e76050；[10.1002/advs.76050](https://doi.org/10.1002/advs.76050) | 作者 DFT 给出 CuF₂ 双态 Pz=±1.23 pC/m，势垒 11.49 meV/f.u.，层锁定 altermagnetic splitting 和中间堆叠 | 已发表 FE＋磁性＋层输运控制；我们必须增加不同的磁场→极化几何观察量。**中间态只能据原文称 Pz=0，不能未经检验称总 P=0**；可设计 FE-I／IM／FE-II 控制 |
| **B：材料风险** | Derzsi, Tokár, Piekarz, Grochala, 2022，*Charge ordering mechanism in silver difluoride*，PRB 105, L081113；[10.1103/PhysRevB.105.L081113](https://doi.org/10.1103/PhysRevB.105.L081113) | bulk AgF₂ 中关联、磁性、晶格软模与 charge order 的竞争 | 对双层不是直接基态证明；但 G1 需审计 Ag moment、charge/PDOS、gap，不可把一次 DFT+U 收敛当唯一物理相 |
| **B：理论边界** | Malashevich, Souza, Coh, Vanderbilt, 2010，*Theory of orbital magnetoelectric response*，NJP 12, 053032；[10.1088/1367-2630/12/5/053032](https://doi.org/10.1088/1367-2630/12/5/053032) | 周期绝缘体轨道磁电极化率的 Kubo / Chern–Simons 项，有限场/解析理论检验 | α=∂P/∂B 并不是新型响应的命名发明；其 3D 公式不能直接当 2D spin-only Wannier 数值实现 |
| **B：金属/绝缘区分** | Zhong, Moore, Souza, 2016，*Gyrotropic Magnetic Effect and the Magnetic Moment on the Fermi Surface*，PRL 116, 077201；[10.1103/PhysRevLett.116.077201](https://doi.org/10.1103/PhysRevLett.116.077201) | 洁净金属费米面磁矩（spin+orbital）的传统 GME | 严格区分 metal conduction、insulating displacement 及 Faraday background |
| **B：数值工程** | Qiao, Zhou, Yuan, Zhao, 2018，*Calculation of intrinsic spin Hall conductivity by Wannier interpolation*，PRB 98, 214402；[10.1103/PhysRevB.98.214402](https://doi.org/10.1103/PhysRevB.98.214402) | MLWF 下 Hamiltonian、overlap、spin operator 联合插值 | 不是 IGMC 算法的直接证明；但提醒单有 hr.dat 和 centers 不足，必须审计同一 gauge 的 spin / Berry-connection / position 算符 |
| **C：层响应近邻** | Zhou, Zhou, Pan, Hu, Liu, Sun, 2026，*Layer Edelstein Effect*，PRL 137, 076401；[10.1103/qhvs-pfkq](https://doi.org/10.1103/qhvs-pfkq) | 电流诱发 layer-resolved spin accumulation；层交换与外电场调控 | LEE 的 Type-I/II 与 AgF₂ 磁连接 Type-I/II **不同定义**；借鉴先单态、再层关联、再总响应，不照搬 Boltzmann Fermi-surface 数值 |
| **C：层 Hall 近邻** | Li, He, Zhou, Pan, Tan, Sun, 2026，*Prediction of a layer nonlinear Hall effect in bilayer nonmagnetic or antiferromagnetic systems*，PRB 114, L111410；[10.1103/2gyh-y52n](https://doi.org/10.1103/2gyh-y52n) | layer-odd nonlinear Hall、层间抵消与门电场解锁，80 层群 | 计算 α 时可并行做 layer-projected 解读，但投影值可能有约定依赖；整体 gauge-invariant α 优先 |
| **C：重要的近邻对称论文** | Chen, Duan, Xiao, Zhao, Deng, Wang, 2026，*Nonlinear layer Hall effect and detection of the hidden Berry curvature dipole in PT-symmetric antiferromagnetic insulators*，PRB 114, L171110；[10.1103/k5f2-w2z7](https://doi.org/10.1103/k5f2-w2z7) | PT-AFM 的 hidden Berry curvature dipole / metric dipole，磁序 odd/even 与 Ez odd 的可辨响应 | 论文报告必须把「磁有序反转」「FE 电场/滑移反转」「测量张量变换」分开，不能将两种翻转混称 |
| **C：早期滑移—磁响应先例** | Jian Zhou, 2022，*Photo-magnetization in two-dimensional sliding ferroelectrics*，npj 2D Mater. Appl. 6, 15；[10.1038/s41699-022-00297-6](https://doi.org/10.1038/s41699-022-00297-6) | 滑移切换下光注入非平衡 spin+orbital 磁化反转 | 「滑移切换磁响应」不是首次；但光驱动磁化和本项目磁场驱动绝热电极化有不同因果方向与频率物理 |
| **C：最新结构—磁耦合近邻** | Luo, Tan, Chen, Yin, Zou, Zhong, 2026，*Switching of magnetic anisotropy via strong magnetoelectric coupling driven by sliding ferroelectrics*，PRB 114, 144423；[10.1103/c9vm-twc3](https://doi.org/10.1103/c9vm-twc3) | 双层 AlN 滑移极化调控 MnS₂ 磁易轴，界面电荷转移及轨道杂化 | 突出材料 SOC easy axis、结构畸变与 hybridization 需要协同检验；该文不是磁位移电流先例 |

**没有核实的近邻：** “Control of the Rashba-Edelstein Effect via Sliding Ferroelectricity in Bilayer Systems” 在先前用户研究笔记中为 submitted work。此次仍未得到可核实 DOI 或正式刊发页面，不能写为已发表或以不存在的量值比较。已检索的公开文献未见完全相同的“AgF₂ 稳定双态 + spin-only αxx/αyy 的 Wannier 验证 + magnetic connector parity”全链成果，但**检索未见不等于证明首创**。

### 2.1 最重要的 peer-paper 方法学启发

**从 Ghorai 得到：** 追问总信号中 orbital 是否主导，2D 实际有限厚度与 in-plane B 的轨道耦合不同于零厚度单层的理想极限。不能简单用某一原子局域 L 算符代表完整 orbital coupling；需要空间位置、速度、场作用和同规范处理。对所有 orbital 数值声明必须额外验证正确性。

**从 Zhu/Peng 得到：** 只报对称允许不够。适合强论文的图需要经过 relaxation 的 FE 两极小值、极化、相应磁序和磁群、能带，再连接实际响应。CuF₂ 的 IM 为 Pz=0 **不等于**其整个极化矢量必然为零（原文 IM 空间群 Pca2₁ 有可能允许其他方向的极化）；后续图注称“zero-out-of-plane-polarization intermediate”。

**从 Qiao/经典 ME 理论得到：** α 的 gauge-covariant 多带实现是独立研究瓶颈。应从一个 B1 态开始，先证明 direct spectral kernel 与 finite-field occupied Wilson polarization 在同一嵌入位置和 Zeeman 定义下相等，再算第二态，不允许用玩具模型数值类比 AgF₂ 电流。

## 3. II｜模型已有证据与下一个最省资源的数值判别

### 3.1 已有项目证据（与同行作者结果分开）

- 八带模型：两个真实空间不同轨道中心 × 双层 × spin；层间 \(u\) 仅通过显式 two-center hopping 的位移依赖进入，未加 \(\delta(u)\tau_z\)。在 \(u=\pm0.19\) 的 AF-like \(D_{xx}/\omega\approx\pm0.00308533\)，FM-like \(\approx+0.01937573\)，为 **reduced model units**。
- Direct interband kernel 与 finite-Zeeman occupied non-Abelian Wilson 极化导数显示网格收敛一致趋势，但有限 mesh 不完全相等。
- AF-like 模型中移除 spin-dependent interlayer hopping 后响应约下降 57 倍：这是既定模型的受控灵敏度，**非**材料普适增强规律；以后需分解 gap denominator 与矩阵元变化。
- 占据带 band-energy scan 的样本最小点在 \(u=0.5\)，\(\pm0.19\) 不是模型中已证明的 local minima；此模型**未证明 BSF**。
- \(u=0.5\) 的 FM-like \(\langle\tau_z\rangle\approx0\) 而 \(D_{xx}/\omega\approx0.00677\)；严格表述是「零 layer-occupation imbalance proxy 下存在 finite response」，**不是**完整晶体极化为零的结论。
- AgF₂：已有来自 C2DB 单层的构造 B1/B2，磁镜面映射经轴矢量规则核查；不是作者 relaxed POSCAR。真实 B1/B2 结构稳定性、磁基态与材料 α 均未独立建立。

### 3.2 短链 decision tree（不增加未经授权的计算）

| 阶段 | 要回答的二元科学问题 | 通过的最低证据 | 未通过时停止什么 |
| --- | --- | --- | --- |
| **P0：已提交 G1** | 作业真的运行且可信 SCF 吗？ | 真实 PBS history/OUTCAR/vasprun，收敛、Ag site moments、gap/occupancy、TOTEN、程序版本；解释 receipt timestamp 不一致 | **停止任何材料结果主张** |
| **P1：B1 结构磁态** | 有可重复的稳定、有隙磁性局部状态吗？ | 额外授权后极少量不同起始层距、磁种子 relaxation；力、能量、磁态 | **停止 B2/NEB/IGMC 投入** |
| **P2：双稳态** | B1/B2 都是局部极小值且 Pz 相反吗？ | 一致泛函/磁序下 relaxed pairs、2D 极化、少量滑移中间点的能量升高 | **停止「ferroelectric switching of IGMC」主张** |
| **P3：真实响应** | 非零 spin α、正确磁群映射与独立数值校验都过关吗？ | SOC 下 MPG、Wannier H+spin+position 同规范、Kubo/Wilson finite-B 交叉验证；首先 αxx/αyy | **停止增加应变/材料/频率扫描** |
| **P4：实验/投稿层级** | 量级与物理通道能支持可测 readout 吗？ | orbitals/induction/units、2D→device、相位、robustness 与 control | **不声称完整可测电流** |

每一步都先明确目标，不提前启动后续阶段。Type-II 4% 是另一个需独立批准、弛豫和验证的分支；不能由零应变 Type-I 的一个结果直接外推。

### 3.3 对论文定位的更新

- **低风险可说：** “a model-validated spin-Zeeman quantum-geometric polarization response with connection-dependent parity”; “a source-grounded AgF₂ candidate awaiting ab initio validation”。
- **不宜说：** “new IGMC mechanism from scratch”; “first ferroelectric magnetic displacement current”; “general spin-only odd/even law”; “measurable AgF₂ current”；“CuF₂ nonpolar in every component”。
- **有条件的投稿主线：** stable sliding FE magnetic bilayer → two symmetry-connected magnetic domains → same-gauge spin / orbital mixed geometry → quantitative low-frequency α → controlled switching and physically distinct readout, with explicit experimental backgrounds.
- **负结果有价值：** 若 AgF₂ 不满足磁绝缘双稳态，可将此作为有效筛选失败，保留八带机制，不人为提高 U 或使用不等价构型凑结论。CuF₂ 候选优先需要验证其**具体 magnetic point group 是否允许目标 α**，而非先进行高成本 Wannier。

## 4. 与上一版长报告相比的审查修订

1. **新增直接竞争**：Chen *et al.* 2026 PT-AFM nonlinear layer Hall 和 Luo *et al.* 2026 滑移磁各向异性，两者都已查到 APS 正式 DOI，但属于近邻而非同一 IGMC 通道。
2. **纠正 CuF₂ IM 的“nonpolar”措辞**：已核实作者只报告其面外 Pz 被 glide 对称性压制；Pca2₁ 不意味着全部矢量极化为零。只可称 “zero-Pz intermediate”，直至对称性及三分量极化另行证明。
3. **提高 Ghorai 竞争等级**：这不是仅供装饰的 orbital 综述，而是已经具有密度矩阵统一推导及反铁磁材料演示的 IGMC 直接先例；但它的 orbital/spin 数量级比仅针对 CuMnAs 的作者模型，不能直接认为 AgF₂ 必然如此。
4. **对 Novelty 降调且更可证伪**：磁连接 odd/even 是静态 ME 通用选律，新颖部分必须来自 switchable material + converged quantum geometry + mechanism/control。
5. **确认文献边界**：旧研究报告的 10 篇 DOI 大体成立；这里将“Methods”“Competing mechanism”“Materials/ferroic switching”“Layer-response analogies”明确分组，并添加未核查 submitted 工作状态。
6. **算力优先级不变**：G1 PBS 号 490620.cmspbs1 的旧 PR #25 不是实时状态；本次文献更新不能替代 read-only PBS/OUTCAR 审计，且**绝对不重复 qsub**。

## 5. 直接可用的阅读顺序（按决策收益）

**首先：** Ghorai 2026 (orbital IGMC) 的 Eq. 结构、spin/orbital 统一定义、CuMnAs 比值；Zhu 2025 (AgF₂) 的磁连接与结构；Peng 2026 (CuF₂) 的 FE-I/IM/FE-II、磁态、真实势垒。  
**第二层：** Malashevich 2010（静态 α 的已有理论边界）、Qiao 2018（Wannier 算符一致性）；Chen 2026（不同 switch parity 的实验可分辨分类）。  
**等待材料 gate 通过才深入：** 显微 orbital IGMC 实现、有限厚度 Faraday/background、多个滑移势阱的高精度 NEB 和 4% 应变 Type-II。

## 6. 证据和维护约束

- 外部 DOI 与期刊状态核查截至 2026-10-21，直接指向期刊/预印本。不对未收到的作者 relaxed POSCAR、自洽作业输出或尚未完成的材料 α 进行臆断。
- 关联现有仓库 [研究 README](../README.md)、[模型审查](SCIENTIFIC_ASSESSMENT.md)、[原论文未决事项](../../../docs/OPEN_PHYSICS_QUESTIONS.md)、[AgF₂ magnetic readout](../agf2_preflight/docs/STAGE5_MAGNETIC_POLARIZATION_READOUT.md) 和 [DFT gates](NEXT_DFT_GATES.md)。
- 此次仅新增同行审查文档及 README 文档索引。**没有**改动任何代码、模型冻结数值、VASP 参数、POTCAR、PBS 作业或 G1 PR #25。
