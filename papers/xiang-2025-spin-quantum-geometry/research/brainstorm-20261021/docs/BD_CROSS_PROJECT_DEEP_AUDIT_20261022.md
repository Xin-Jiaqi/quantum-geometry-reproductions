# 2026-10-22 科学总审计：复现、B+D Brainstorm、同行竞争与下一轮理论门槛

> **证据等级：** REPRO_PARTIAL / INDEPENDENT_METHOD_VERIFIED / TOY_MODEL_ONLY / REAL_INPUT_PREFLIGHT / REAL_NUMERIC_NOT_CONVERGED / OPEN_RESEARCH_HYPOTHESIS。  
> 这份交叉审计首先依据 2026-10-22 当前可访问的 GitHub 原始 README、\`RELEASE_STATUS.md\`、\`physics/README.md\`、各分支科学报告及 R1 私有验收记录；外部文献采用出版社正式页面或明确标注预印本。不能将测试脚本 PASS 推断为发表论文全篇复现，也不能把假定的系统与实验验证混淆。
>
> **执行边界：** 用户此前关闭 MnPS₃ R1；R2-P0 只授权阅读已有私有证据。没有新的 DFT/WSL/日本超算/Slurm 作业许可。项目原始 PRL 复现代码及科学验收不修改。

## I. 跨四篇文章的真实复现水平

| 复现文章/目录 | 已验证的独立能力 | 主要 unresolved | 科学状态 |
| --- | --- | --- | --- |
| Jia et al. PRL 136, 226302 (2026), \`jia-2026-nonlin-magnetoelectric-edelstein\` | Kane–Mele、PT 子空间/密度矩阵与若干静态几何数值测试 | intrinsic NMEE 的 source SM Eq.(9)/(11)/(13)/(14) 对角密度矩阵项与驱动 DC、CuMnAs 退简并、响应/接触项协议一致性 | **PARTIAL** |
| Xiang et al. PRL 134, 116301 (2025), \`xiang-2025-intrinsic-gyrotropic-magnetic-current\` | Rashba/Dirac、Zeeman 几何及部分频率响应 | published Fig.2(d) 归一化/磁场与电流符号约定 | **PARTIAL** |
| Xiang–Jin–Wang PRL 135, 146303 (2025), \`xiang-2025-spin-quantum-geometry\` | finite-η full-band Liouville、PT 项目子空间、U(2) 检查、SM Eq.(2) 四源项及 Laurent 模型门槛（截至所读 release 70 测试） | Fig.1(a) 绝对值近 2 倍、若干 SM 代数问题、η→0 /τ⁰ 金属 DC 极限依赖碰撞/制备协议；**C 与 η 调节相关，不能自动称 intrinsic driven-DC** | **PARTIAL** |
| Xiang et al. PRB 113, L201406 (2026), \`xiang-2026-quantum-geometric-map\` | 独立多带非简并 Berry 曲率/度规/轨道磁矩核、有限 torus Hall、Kubo/SCBA/Ward 有限子集、原文 16 页 SM 的源式审阅、一个 Ω×轨道磁矩子项 | SM Eq.(20)→(31) 符号、主文 Eq.(10) vs SM Eq.(62) 自旋螺旋度、Fig.2(b) 约 2π 量级及完整磁场/轨道位置 shift、非阿贝尔/真实 Wannier 接口 | **PARTIAL** |

**总体判断：** 方法库增长远快于整篇文献的物理验收。尤其是 PRB 的 v54 非简并 \`N\` 带几何 \`VERIFIED\`，不是全多带磁性输运 \`VERIFIED\`。复杂 QWZ 随机 Hall 的许多结果属于**额外方法开发**，不是原文自旋 PHE 的必要复现门槛。新的代码应优先补齐 **响应协议与观测算符的物理一致性**，而不是再增加随机测试或堆积未解释的项目。

**科研价值：** 原文符号/手性不一致的有源文件、公式级证据，是高质量独立审阅；但单篇论文中疑似公式问题不能替代新效应。若须继续复现，优先选 **NMEE 的 source-diagonal 项/有限 Zeeman 静态 oracle** 或 **PRB 两色/碰撞制备 Liouville 对照** 中一个，以精确封闭一个 physical gate。

## II. B+D 在物理上证明了什么、否定了什么

1. **Gate 1–3：** 四带 PT 表示无法合法容纳特定裸自旋层间 hopping；加入反宇称轨道与原子 \(L\cdot S\) 后可以通过投影得到有效自旋相关层间跃迁。**不是微扰方法或 SOC 机制的首创**。
2. **Gate 3 NO-GO：** 层反向 AC 自旋响应在 \(m=0\) 下仍有信号；磁性 odd 的总自旋静态响应是标准混合 \((k,B)\) Berry 曲率，不是新 spin QGT。
3. **Gate 4：** 含不等价 A/B 蜂窝基元的 toy 能产生等能、反号的层电子电荷差和局域双稳态；但 32 带模型中非极性 AA 能量比 toy AB/BA 更低，**没有真实铁电基态**。最重要的严格约束是 AB/BA 若保持 \(C_3\)，线性 \(E_{x,y}\to S_z\) 禁戒。
4. **Gate 5–6：** 全电子本征态 Kubo 的层自旋连续性恒等式成立。\(\omega=0.06,N_k=23^2\) 时层间 \(\sim+i1.021742\times10^{-4}\)、电子 exchange \(\sim-i1.008285\times10^{-4}\)、SOC \(\sim-i1.927658\times10^{-6}\)、净 \(\sim-i5.819214\times10^{-7}\)。这是**接近抵消的模型现象，不是新的角动量连续性定律或 Néel 磁矩动力学**。
5. **不可挪用的磁序差别：** 这些 toy 多为**层 A 型** \(m\tau_z\sigma_z\)，真实 MnPS₃ 的每个层内都有两个相反 Mn 子晶格。layer-odd spin 与真正的源验证 Néel 模不相同。
6. **R1 真实材料：** 文献派生未弛豫 \(\mathrm{MnPS_3}\) 两态、Cm/L2/Eq.13 的几何/磁性前检通过；para Slurm 两项分别在约 2h TIMEOUT、0 个完成的离子步、首个 SCF 未收敛。原始私有 \`result.json\` 为 \`RUN_FINISHED=true, NUMERIC_CONVERGED=false, FE_VERIFIED=false\`。**物性和原创假设仍未在真实材料上测试。**
7. **G0 全层群 toy：** 完整**幺正**磁性等价畴的静态实响应矩阵 \(K_B=D_Q(c)K_A D_E(c)^{-1}\)，保留秩、奇异值与对称性允许系数维数。固定实验室轴亮/暗切换可能只是坐标重新定向。**标准等变性负对照，不是新定理**。

## III. 本轮新增、独立执行的关键纠正：反幺正交流互易与“开关”

在 PR #32 追加 [**antiunitary_kubo_onsager_gate.py**](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/blob/research/fe-domain-mode-addressability-20261022/research/fe_domain_mode_addressability_20261022/antiunitary_kubo_onsager_gate.py)，用两个完备谱系独立计算以下 Kubo：
\[
\chi^R_{AB}(\omega)=\sum_{n,m}
\frac{(f_n-f_m)A_{nm}B_{mn}}{\omega+E_n-E_m+i\eta}.
\]

对于 **真正映射完整 Hamiltonian 的反幺正** \(\Theta=UK\)，有
\[
H_{\bar d}=\Theta H_d\Theta^{-1},\quad
\Theta A\Theta^{-1}=\eta_A A,\quad
\Theta B\Theta^{-1}=\eta_B B,
\]
在该算符/频率约定下满足经典 Onsager–Casimir **交叉**关系
\[
\boxed{\chi^{(\bar d)}_{AB}(\omega)=\eta_A\eta_B\,\chi^{(d)}_{BA}(\omega)}.
\]
不是天真地比较 **同序** \(\chi^{(\bar d)}_{AB}\overset{?}{=}\eta_A\eta_B\chi^{(d)}_{AB}\) 或简单复共轭。带真实 spin-\(\tfrac12\) 的 \(\Theta^2=-1\) 四带随机 Hermitian 例子（\(A=s_z\) 奇，\(B\)=轨道偶）在 \(\omega=.24,\eta=.03\)：
- \(\chi^{(d)}_{A B}=0.0324050660-0.0015220503i\)
- \(\chi^{(d)}_{B A}=0.0319998163+0.0016507039i\)
- \(\chi^{(\bar d)}_{A B}=-0.0319998163-0.0016507039i\)
- 互易残差 \(\sim10^{-17}\)；天真的同序负号失败。

这是**已知** Onsager 互易，原创性不在恒等式。**本项目此前的 FE-switching odd/even AC 声称必须按具体连接算符逐项核验**。尤其 MnPS₃ 来源 Eq.13 非相对论 \([E||M_zM_y]\) 不是自动的 SOC 反幺正操作。

**范围：** 测试的是普通两/四带完备 Lehmann (η>0)，不是实际磁层群的频率完整余表示、碰撞、热噪声、磁共振或 Wannier 偶联，也不是新增完整响应分类。

## IV. 新同行研究：直接收窄我们的三条主线

| 正式新文献 / 预印本 | 已做工作 | 对 B+D 的影响 |
| --- | --- | --- |
| Yuan et al., [PRL **137**, 106901 (2026)](https://doi.org/10.1103/12cl-b9jj) | 双圆偏振光探测 altermagnetic magnon 量子几何，非拓扑条件下仍有效 | **非平凡 magnon geometry + 光学探测本身不新** |
| Patra & Cocchi, [PRB **114**, 225412, 19 Oct 2026](https://doi.org/10.1103/ldws-7yvc) | 双层 MoS₂ 堆叠/扭转点群到 \(\chi^{(2)}\)，频率无关 SHG 角花瓣旋转 | **stacking/twist→二阶光学张量/旋转不新** |
| Schilberth et al., [npj Quantum Materials **11**, 26 (2026)](https://doi.org/10.1038/s41535-026-00857-9) | Co₂Mo₃O₈ 的 IR/Raman 声子选择定则实验，发现相对论磁群比非相对论 spin group 更符合磁序后观测 | **spin-group/SOC 选择定则差异已有实验先例** |
| Zhou et al., [PRL **137**, 076401 (2026)](https://doi.org/10.1103/qhvs-pfkq) | 完整 LEE、80 LG 堆叠筛选、layer-odd/current-induced spin、材料预测 | **任意层响应矩阵的通用 LEE 分类不新** |
| Zhang et al., [arXiv:2605.25484v2](https://arxiv.org/abs/2605.25484) PREPRINT | \(I_k\le12\) 统计 **40,087** 自旋层群及小表示 | **自旋层群枚举/基本 corep 不新** |
| Sun et al., [Nat. Commun. **17**, 3101 (2026)](https://doi.org/10.1038/s41467-026-69635-2) + Yang et al. [ACS Nano **20**, 23852–23860 (2026)](https://doi.org/10.1021/acsnano.6c06799) | FE 六态磁带/纯自旋光电流的连接关系 | **简单畴连接导致磁/光学响应反号或重定向不新** |
| Meguro et al., [PRB **114**, 194416 (2026)](https://doi.org/10.1103/gdnn-wrbc) | \((k,m)\) 混合空间量子几何介导非线性自旋驱动电流 | **混合动量-磁矩几何自身不新** |
| Li et al., [Nat. Commun. **16**, 5451 (2025)](https://doi.org/10.1038/s41467-025-60647-y) | twisted hBN/WSe₂/MoSe₂ 极化 meron/antimeron 纹理 | **moire FE 纹理和拓扑极化涡旋已有实验** |

**核查深度说明：** 上表针对出版商摘要/可访问正文与项目既有同行文件；没有对以上每篇 SM 逐条公式完成查重。因此不存在“经过证明的无人涉足的新效应”结论。其他 prior art（全层群 FE PRL、448 collinear spin LG 的 PRB、vdW SOT、MnPS₃ THz、NSOT、spin quantum geometry）已在 [主同行审阅](PEER_WORK_AUDIT_20261021.md) 中。

## V. 全层群理论的真实高收益问题与候选排序

**优先 A：site-/layer-resolved *physical observable* × FE switching graph × correct unitary/antiunitary response mapping。**  
不是新建 528 张矩阵表，而是输入 (a) 每个真实磁性 Wyckoff 轨道、(b) 有物理归一化的局域自旋/磁模力算符、(c) 磁层群实际诱导表示与连接算符、(d) 测量/驱动双算符，然后输出：
1. 每个态严格禁戒的观测通道；
2. 完整幺正等价带来的 invariant / fake lab-fixed on/off；
3. 反幺正互易带来的 **forward vs reciprocal** response 对偶关系，而不虚构全频域简单奇偶；
4. 结构连接但 SOC/磁序不连接时的可测模式/允许通道差异；
5. 哪些结论可由已有 SOT/LEE / IR 选择定则完全覆盖。

**科学问题：** 这样一个针对相互关联材料响应**而非能带表**的连接图，是否能导出已有 layer-magnetic/spin-LG 分类**没有直接提供的、可测的非平凡零通道/互易判别**？目前未证明。若只重新叙述 Neumann/Onsager/LEE，**停止**。

**备选 B：从 BSF 局域堆叠空间 \(T^2\) 到 moiré 真实空间的** *symmetry-equivariant local response field*。对每个局域 registry \(\boldsymbol\delta\) 定义明确的 tensor \(K(\boldsymbol\delta)\)，并在局域点群、畴连接和梯度校正下求响应纹理的严格受保护零点与实空间观测；与已经存在的 twisted FE polar merons、twist SHG rotation、moire 非线性 Hall 区分。先**群论驱动的零点/变换界限**，后结构弛豫/非局域响应；不能把局域 tensor 直接当作 moiré 真实测量值，不考虑梯度/电子非局域性。当前高风险但与原有 BSF、BPVE/SHG、材料数据库更贴合。

**暂缓 C：magnetic magnon FE-mode excitation 作为全层群主论文。** 最近 PRL magnon quantum geometry + Nature SSG magnon + 2026 FE magnon chirality + SOT six-class 分类，使它成为竞争最激烈、且真实磁性动力学/DFT 成本高的路线。可保留做 A 的代表性 observable，但不要反过来绑架 A 成“必须 MnPS₃ VASP 才能发表”。

**不投资 D：仅增加玩具八带/32 带数值、未经源结构验证的 MnPS₃ 频率谱、纯允许 tensor 表/完整预印本已给出的 spin-LG irreps、未处理 antiunitary 的奇偶图。**

## VI. 最小后续动作，保证有真结果

**G1 (完成)：** 普通幺正畴变换可证规范约束与 nonsymmorphic \(g^2=T_a\) 相位负对照，位于群论仓库 [Draft PR #32](https://github.com/Xin-Jiaqi/group-theory-operations-toolkit/pull/32)，仓库保护主分支，**未合并**。

**G1b (本轮完成)：** 增加独立完整 Lehmann 的 antiunitary **reciprocal** 交叉响应测试，实证证明 forward 同序奇偶推断不可靠。测试脚本与原理已加 PR #32；**不是群/频率全覆盖**。

**G2 (下一轮最小目标)：** 用群论仓库已有 Seitz/层交换数据，挑 2–3 个**具体且源可追踪的站点轨道/磁性群**：一个 unitary conjugate NO-GO，一个 full SOC mapping OBSTRUCTED，一个不同 magnetic stabilizer。实际计算 site-mode 诱导表示、Hom 的维数与 basis，再加入反幺正响应双算符表。绝不从群名推断真实磁模式；报告 \`PASS/NO-GO/UNKNOWN\`。**如找不到非平凡、可观测且不已知的区别，不扩展 528 群**。

**G3 (严格有条件)：** 如果 G2 真有可复验零通道/互易判别，用已有 v54 多带几何和 PT/U(2) 子空间方法作一个 4–8 带、全部接触项正确的**full Kubo oracle**，验证自然响应强度和带隙/频率。这是独立响应程序，不是再指定参数追数值。

**G4 (物理和公开性成熟后)：** 一个材料/moire 原型。必须先解决真实算符的 Bloch/Wannier 位置连接、磁性位点、材料稳定性和输出耦合；R1 不得自动重算。

## 结论

- **科研最强的资产：** BSF/层群空间几何＋磁群/层交换操作数据库＋逐论文自研量子几何/Kubo oracle＋独立反事实验收。
- **最大风险：** 构造许多物理正确但属于已知 response 的模型，并把软件 CI 或张量允许当作新发现。
- **投稿标准：** 首先出现一个**反驳/补足已发表分类的可测新选择规律或定量几何机制**，然后一个严谨模型验证，最后才需要材料。
- **目前只有候选 A、B，没有已证明的原创机制。**
