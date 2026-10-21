# Spin QG × Layer / Magnetism：A–C 科学 Brainstorming 交接

**日期：2026-10-21**  
**交接目的：**在一个全新的 ChatGPT 对话中，立即继续量子几何驱动的自旋输运原创选题筛选，主要推进 **C2：磁性双层中的交流层间自旋转移与几何响应**。本文同时回答“A–C 是源于 Spin QG 还是其他论文”这一未正式回答的问题。  
**工作性质：**科学讨论、分析性设想与部分模型代数；**不是新论文的已验证科学结果，也不是论文复现软件 release**。  
**用户当前决定：**A 暂存；B 作为备选；C 优先。先物理问题和独立推导，**不急于技术接口、代码重构、真实材料 DFT 或上传 GitHub**。强调 high-value/low-cost、严格可证伪、创新性查重、避免概念堆砌。中文直接、物理深度充分、不讲黑话、不重复解释基础。

---

## 0. 一段话交接给新 ChatGPT

> 已完成对 Xiang、Jin、Wang《Spin Transport Revealed by Spin Quantum Geometry》（PRL **135**, 146303, 2025）的核心计算方法分阶段复现，在 GitHub `Xin-Jiaqi/quantum-geometry-reproductions` 的独立论文目录中保存，科学状态仍为 PARTIAL。讨论过 A（隐藏层反向二阶自旋流）、B（层间杂化/滑移与自旋量子几何）、C（磁性绝缘体交流 displacement spin Hall）；A 发现概念拥挤、简单双层结果不需要真正的层间相干、算符恒等式可能强制响应为零以及 k 依赖 gauge 陷阱，已归档；B mixed (k,δ) geometry 存在先例，留作备选；C 对 Spin QG 原文线性 AC spin Berry curvature / spin metric 的相位、时间反演奇偶性，结合磁点群和层交换，得到一个待严验的 m_z T 交叉过滤选择定则。最新拟重点研究 **C2 层间自旋转移算符**，其在 t⊥→0 时严格为零，避免独立单层叠加的平凡机制；但这不自动意味着它有新几何或能产生可测 Néel torque。请先严格重推原文 AC Kubo、反幺正的复频选择定则、完整自旋连续性方程和最小绝缘模型。先给低成本 no-go 或明确 go/no-go，不要立即做 DFT/材料筛选/接口开发。

---

## 1. A–C 的论文来源：不是只来自 Spin QG

| 方向 | 核心科学出发点 | GitHub 其他论文能贡献什么 | 不能混为一谈的地方 |
| --- | --- | --- | --- |
| **A** 中心对称隐藏二阶自旋流、层相反贡献 | **Spin QG PRL**：非线性 Drude 自旋流与 spin quantum metric dipole（NDSC）；外部 **Layer Edelstein Effect、Layer Nonlinear Hall** 的层抵消/层选择逻辑 | **Jia 2026 NMEE** 的 PT 简并、自旋响应、printed vs static 协议；**Xiang 2026 PRB** 的有限 η Liouville、磁输运方法可以作独立数学工具 | NDSC、NMEE、非线性电荷 Hall 不是同一种响应；层差自旋流算符与自旋密度也不同 |
| **B** 层间耦合与滑移参数空间的几何 | Spin QG 的自旋权重量子几何和简并投影；外部 **interlayer electric multipoles / sliding quantum geometry** | **Xiang 2026 PRB** 的 spin/orbital 几何、Kubo/Liouville 模块可提供测试方法；**NMEE** 提供同规范自旋算符/层态参考 | 普通 quantum metric、spin quantum metric、mixed `(k,δ)` geometry、轨道响应彼此不能简单互换；必须证明可观测量，不是发明几何记号 |
| **C** 绝缘体交流自旋输运与磁性连接操作 | **最直接来源是 Spin QG PRL Eq.(7)、SM Eq.(59)–(62)**：线性 displacement spin Hall 分解成 spin Berry curvature 与 spin metric 两个相位通道 | **Jia NMEE**：磁序与 PT、层分辨 spin density / torque 的模型经验；**Xiang 2025 IGMC**：有限频率 Kubo、因果/接触项核查思路；**Xiang 2026 PRB**：Liouville、频率、磁输运/连续性方法，但不直接提供本文 AC spin Hall 公式；**Group Theory Operations Toolkit**：磁点群、反幺正与堆叠连接群操作 | **IGMC 是 B 驱动电流，不能直接拿来替代 E 驱动的 displacement spin Hall**；NMEE 是非线性 magnetoelectric，C 起点则是**线性、同频**交流自旋响应 |

**因此：**A–C 的物理主线始自 Spin QG，但我们的优势是**跨论文方法、对称性与层自由度的物理组合**；不能把 NMEE、IGMC、PRB 的响应公式直接粘到 Spin QG 公式上。真正新意要体现在可证明的新机制或可测关系，而不是引用更多热门概念。

### 公开代码与阅读入口

- 四论文总仓库：https://github.com/Xin-Jiaqi/quantum-geometry-reproductions
- Spin QG：https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/tree/main/papers/xiang-2025-spin-quantum-geometry
- NMEE：https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/tree/main/papers/jia-2026-nonlin-magnetoelectric-edelstein
- IGMC：https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/tree/main/papers/xiang-2025-intrinsic-gyrotropic-magnetic-current
- PRB：https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/tree/main/papers/xiang-2026-quantum-geometric-map
- 群论：https://github.com/Xin-Jiaqi/group-theory-operations-toolkit
- 现有仓库通用交接：https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/RESEARCH_HANDOFF.md
- 统一物理协议设计：https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/blob/main/docs/PHYSICS_INTERFACE_SPEC.md
- 本轮检查时复现仓库 main HEAD：`bc362053f0f643bcdda7eaef7176ccaa2ec32a57`；群论仓库 main：`551a0daf8d4ba5a5f24b2811694bfe438616de48`。**这是时间点快照，不代表新会话的最新状态，开始时重新查询。**
- 原文（合法 DOI）https://doi.org/10.1103/14mp-263q 。此项目文件中还提供正文 PDF `Xiang 等 - 2025 - Spin Transport Revealed by Spin Quantum Geometry(2).pdf` 与 SM `SuppleMater(3).pdf`。新对话若附件无法访问，应重新上传，不能假装已读取完整 PDF/SM。

---

## 2. Spin QG 原文中必须严格区分的两个物理量

**(i) 金属非线性 Drude 自旋流（NDSC）**：原文含 spin quantum metric dipole 的 `E^2` 响应，带有 Drude/τ² 性质。我们既有模型复现已可调用，但对物理金属 intrinsic DC 的解释不能随意推广。

**(ii) 绝缘体线性交流 displacement spin Hall**：原文 Eq.(7)、SM Eq.(59)–(62)，是 **E(ω) → J_s(ω)**，不是 E²→J(2ω)；包含与常规 DC intrinsic spin Hall 背景区分的**附加频率依赖部分**。在所用频率/时间约定下，写为

\[
\delta J_s(t)=E_0[\mathcal A_{\Omega}(\omega)\cos \omega t+\mathcal A_g(\omega)\sin\omega t].
\]

- `Ω` 来源于论文定义的 **spin Berry curvature**，在原文变换约定中对应 **T-even** 通道。
- `g` 来源于原文的 **spin quantum metric**，在其约定中对应 **T-odd** 通道。
- 低于绝缘体光学带间阈值、无其他奇异项时，原文位移响应分量有 `Ω ~ O(ω²)`、`g ~ O(ω)` 的预期标度。**这是原论文的内容，不是我们的新结果；实际可被积分/对称性消去，先以 SM 原式重推。**
- 对称性有 `T` 的非磁体系，整体 `g` 通道会被禁止；若存在其他如 `PT`、`C2z` 等，指定分量可能进一步全被禁止。**不是随便挑一套 AFM 就能有非零 spin metric**。
- 时间反演作用到复频响应时，一定同时考虑反幺正共轭、输入/输出算符奇偶、正负频率关系、场与磁序的翻转，不可以仅凭静态张量符号决定 `sin/cos` 部分。

**Spin QG V32 复现发布状态：**有限 η、多带、PT 投影、Laurent、原文 SM Eq.(2) 四源项均有受控方法代码与 70 项 scoped 测试；**整个论文科学状态仍 PARTIAL**（原 Fig.1(a) 大约 2 倍归一化、SM 若干方程代数、物理稳态 DC 极限未闭合）。这些主要是金属/非线性支线；**并不意味着我们已经拥有 C2 的完整 AC 层间自旋转移求解器**。

---

## 3. A 方向已经做过什么，为什么暂停

详细归档：`2026-10-21-A方向探索归档.md`（随本次会话可单独下载）。

### A1 隐藏层反向 E² 自旋流

在反演交换两层的中心对称双层，普通总二阶自旋流被 P 禁止，层差二阶自旋流 **可以** 被允许。垂直栅压破坏 P 可以释放总信号，或用对称外参调制隐藏信号。此类叙事与 LEE/LNHE、非线性 pure spin Hall 高度接近，**创新门槛不足**。

### A2 层间耦合与非阿贝尔自旋量子几何

最小双层 Rashba：

\[
H=k^2/(2m)+t_\perp\tau_x+\alpha\tau_z(k_y\sigma_x-k_x\sigma_y).
\]

` t⊥=0 ` 时隐藏自旋流已可能非零；所以 **隐藏 ≠ 层间相干**。一般量子度规/非阿贝尔 Berry 曲率随 `t⊥` 改变，与目标响应变化的相关性不是因果证明。

等谱模型 `Hβ(k)=Uβ(k)H0(k)Uβ†(k)` 存在严重表象陷阱：纯 Bloch 基变换必须同步变换位置、速度、自旋与层算符；否则变化可能是 gauge 伪效应。只有明确**固定原子轨道基中的真实 hopping 差异**，才能视为两个物理模型。

### A3 绝缘体二倍频层反向自旋流

对一个特定四带 Hamiltonian 有候选精确算符恒等式：

\[
J^x_{-,y}=\frac{\alpha}{2}\cos k_y I_4+
\frac12\partial_{k_y}\varepsilon_0\,\tau_z\sigma_x.
\]

若 `∂_{ky}ε0=0`，响应算符是单位阵倍数，在特定保持每 k 态粒子数的均匀场 Liouville/完整绝缘填充假设下，非平衡阶数响应为零。不能声称 quantum metric 非零就必有二倍频自旋流。SOC 下自旋流、力矩、有限 q 和边界读出尚未闭合。

**证据分级：**曾经口头讨论的绝缘双层二倍频数值/β 扫描，没有交付经过独立 CI/数值协议冻结的计算代码，**一律视为 UNVERIFIED EXPLORATORY**；不要在新对话援引为已验证发现。A 当前 HOLD。

---

## 4. B 方向仍保留什么

**目标：**不是“滑移改变普通量子度规”或定义新几何符号，而是能否证明真实滑移位置参数 `δ` 下，层间量子相干产生了特定可测自旋响应，并且无法由单独层、自旋纹理变化、谱重排或 Bloch embedding 伪效应解释。

可以用固定原子轨道规范的

\[
H(\mathbf k,\delta)=\begin{pmatrix}
 h_T(\mathbf k)&T(\mathbf k,\delta)\\
 T^\dagger(\mathbf k,\delta)&h_B(\mathbf k)
\end{pmatrix}
\]

分析混合 `(k,δ)` 投影几何，如 `Tr(∂_k P ∂_δ P)`，**但这类 mixed geometry 已有先例**（Nano Lett. 2024 interlayer multipoles 等）。任何 spin/layer weighted geometric tensor 必须从物理可测算符与 Kubo/力响应推导，不是插一个 τz、S 就算创新。若 `t⊥=0` 仍保留“层间几何信号”而无真实变化，优先排除 gauge。

**B 目前优先级中低**。如果 C2 的自旋转移几何推导确有希望，B 的 `(k,δ)` 框架可能作为**物理微观调控**重新加入，不要现在开发多层滑移数据库或 DFT。

---

## 5. C 已推进的脉络：名称与时间线不要混淆

最初 **C** = “绝缘体 AC displacement spin Hall，寻找层态/磁态的可切换关系”。

在 `2026-10-21-B-C量子几何选题深度研判.md` 中，曾用 **C2** 表示“P/T 相关 4 态连接选择表”；后续最新讨论**重新定义**了工作子题：

- **C1：`m_z T` 选择的 AC spin Berry curvature / spin metric 交叉过滤器。**
- **C2（最新主线）：真正依赖层间耦合的 AC layer-staggered spin-transfer / torque operator。**
- **C3：外加垂直栅压破坏磁性连接对称后，原本禁戒的相位通道被解锁。**

为了避免新对话混乱，后续一律按**最新 C1/C2/C3** 使用名称。旧报告的 4 态表可称“C-4state parity table”。

### C1：磁性反幺正对称的交叉过滤想法

研究对称操作 `m_z T`：`m_z` 将上下层交换，`T` 时间反演。对**面内驱动 E_x 与输出 J_y^x**，假定某单态磁群保留 `m_z T` 而不保留额外会杀掉这一分量的 `PT`、`C2z` 等，则按原论文的 `Ω`/`g` 时间奇偶约定，可得到**候选的**：

| 线性 AC 响应 | spin Berry `Ω` | spin metric `g` |
| --- | --- | --- |
| 总自旋流 `J_y^x / E_x` | 禁止 | 允许 |
| 层反向自旋流 `J_{-,y}^x / E_x` | 允许 | 禁止 |

解析直觉：水平镜面 m_z 对上述总流与场的比值给负号；层差的层交换增加一个负号；原文 spin Berry / metric 在时间反演下分别偶/奇。**但这不是已核验的完整反幺正复频 Kubo 选择定则**，必须对照原文复数响应的实虚部、ω→−ω 关系和正确响应定义独立推导。

此前讨论曾口头报告一个**有能隙的八带磁性层模型**：无量纲 `gap ≈4.06`，`ω≈0.36`、21×21 k 网格，总流的 `Ω≈0/g≈−2.07e−6`，层差流的 `Ω≈6.96e−7/g≈0`，翻磁序令 g 变号。**这些数值没有在当前 handoff 中配套完整 H(k)、算符、哈希、可执行脚本、独立 oracle 与频率收敛曲线，必须一律标记 UNVERIFIED / NOT REPRODUCED。** 不能作为模型已物理通过的证据。

**已有明确负面结果（待复验）：** `t⊥=0` 时上述允许通道仍可非零；因此 C1 更多是一种**磁对称性过滤规则**，不是必需层间相干的新机制。且 PRL 本文与邻近磁性光学/输运文章已很接近，选择表本身很难构成创新。

### C2：新主线，层间自旋转移算符

从层差自旋 `S_-^z = τ_z σ_z/2` 与物理的 interlayer hopping `H⊥ = t⊥ τ_x` 出发，**只定义其层间跃迁对时间导数的贡献**：

\[
\boxed{\mathcal T_{-,\mathrm{inter}}^z
=\frac{i}{\hbar}[H_\perp,S_-^z]
\overset{\hbar=1}{=} t_\perp\tau_y\sigma_z.}
\]

矩阵等式可单独重新验算。若 `t⊥=0`，**这个特定算符**严格为零，这解决了“两个单层之和造成所谓新响应”的第一道（但非充分）障碍。更一般的 hopping `T(k,δ)` 和自旋不守恒会让层间转移与 SOC spin torque 相互交织。

实际**层差自旋连续性**应为：

\[
\partial_t S_-^a+\nabla\cdot J_-^a
=\mathcal T_{-,\mathrm{inter}}^a+\mathcal T_{-,\mathrm{SOC}}^a
+\text{可能的其它局域源/汇}.
\]

概念挑战：**(a)** 层间转移不自动等于可测 Néel spin-orbit torque；**(b)** `z` 向转移不直接构成 `z` 向磁序的横向开关力矩；**(c)** `t⊥→0` 的零响应部分只来自算符被定义为含 `H⊥`，是“构造上的零”，并不意味着发现了新的量子几何定理；**(d)** 需检查自旋基、层投影、空间分辨连续性、磁性反幺正、物理电磁耦合与实际频率响应中的 diamagnetic/contact 项。

**真正需要证明的强假设**（目前未证明）：

> 某个 `m_z T` 允许的有限频率层间自旋转移（或横向 Néel 力矩）可以写成**规范不变、特定对称性筛选的自旋几何贡献**，它的谱、低频标度或实验可测相位与普通层分辨 spin Hall / spin-orbit torque 有严格区别，并且在真实 `H⊥→0` 的微观极限下消失。

这是我们接下来要**尝试证伪**的科学问题，不是已经发现的新响应。

### C3：电控解锁

额外 `Δ τ_z` 型垂直层势差一般破坏交换两层的对称性，可以使 C1 的“禁止”相位通道打开，对 `Δ→−Δ` 出现特定奇偶关系。但这大部分是普通磁/空间张量学预言。**只有物理读出明显优于已有工作或有独立几何谱关系时，才值得延伸。** 现阶段不作为主线，不写器件应用断言。

---

## 6. 最重要的立即执行任务（下一轮 ChatGPT 应从这里开始）

**建议先做解析与少量独立数值，而不是启动多日工程：**

### GATE 0：重查 Spin QG 原文（必需）

逐字核对 main Eq.(7)、Table I、SM Eq.(59)–(62)：

1. spin Berry curvature 与 spin quantum metric 的**自旋流算符定义**与对应频率核；`T-even / T-odd` 究竟作用于哪些张量分量；
2. time Fourier 规范、相位 `cos/sin`、频率符号以及现实物理读出；
3. `Ω∼ω², g∼ω` 所在的确切绝缘体范围，并显式分离 conventional intrinsic spin Hall 背景。

未能访问正文 PDF + SM 时应提出重新上传，或者仅按 GitHub/已有归档作**有限范围推论**。

### GATE 1：解析磁群与算符审计（关键）

1. 选定 `P, T, m_z, m_z T, PT, C2z` 在**轨道×层×自旋基底**中的矩阵表示；包括 layer permutation 与 k 映射、磁矩方向、外加场；
2. 从复频 Kubo / Liouville **重新推导** `m_zT` 对总 `J_y^x/E_x`、层差 `J_-^x/E_x` 的变换和实/虚分量禁戒，不仅借用静态奇偶号；列出最小允许磁点群和额外禁戒操作；
3. 对 `S_-^a` 与 `T_{-,inter}^a` 进行全部三方向自旋分量的镜面、时间反演变换。选择真正能驱动磁序或可读出的**横向**力矩（如相对于磁矩方向垂直的分量），避免只研究没有有效切换力矩的 `T^z`；
4. 由局域 Heisenberg 方程完整拆解自旋连续性（包括 SOC torque 与层间流），确立计算的物理观测量定义；
5. 在 `H⊥=0` 与 `M=0`、SOC=0 的极限给出严格 no-go，检查是否自动沦为已有普通 SOT。

**交付：**一张经过代数验证的 `observable × magnetic operation × Ω/g phase` 选择表；一个明确定义的 **可测或可观测转换通道**；若根本不能非零直接 STOP。

### GATE 2：最小绝缘模型（只在 GATE 1 有意义时）

1. 明确写出二维**周期**、有绝缘能隙的 4–8 带 Hamiltonian `H(k)`，每个 SOC、交换、层间跃迁项分别满足/破坏哪些操作；`m_zT` 允许时必须没有额外 `PT/C2z` 禁戒，**逐矩阵核对**；
2. 选真实算符 `S_layer`, `J_layer`, `T_inter`, `T_SOC` 和电耦合/velocity，保持统一 orbital embedding/gauge；
3. 从两个独立路线计算有限频率 AC response（标准因果 Kubo vs 独立 orbital Liouville），有已知的源接触项就按同一 scheme 包含；先低分辨率的小网格做对称性强制零，再做收敛；
4. 冻结 `t⊥→0`、`M→−M`、`Δ→−Δ`、`SOC→0`、k 网格、`ω→0`、`η→0` 和 gauge 变换，对比 Spin QG 原文 AC 线性参考；
5. 若只是普通已有 SOT Kubo 换了一个算符，且没有独立 geometrical term / sum rule / selectivity，则 **STOP C2**。

**交付：**一个清晰非零模型、至少一条 symmetry-forbidden 的严格零基线、一条可以证伪的几何关系，以及计算条件/单位/证据；优先 2–3 页物理笔记，再考虑写 code。

### GATE 3：独立文献查重与 Go / No-go

精准查：绝缘 AFM AC Néel SOT、interlayer spin-transfer torque、layer-selective spin Hall、spin quantum metric torque、`m_zT` magnetic antiunitary tensor constraints、spin displacement/torque dipole、geometry-based FMR readout。重点比较**相同的响应算符、同一时间/空间对称性、相同低频幂次和可读出方式**。

Go 条件：至少一个有可复现证据支持的非平凡新关系（受保护的谱特征/力矩通道/几何恒等式）、独立对照与明确 near-neighbor novelty distinction。No-go：仅由已知 Spin QG 线性响应 + 普通磁张量协变 + `tτyσ` 定义构成，或者无法产生可观测/区分信号。

---

## 7. 文献查重路线（**非穷尽**，请新会话核对原文）

- Xiang, Jin & Wang, *Spin Transport Revealed by Spin Quantum Geometry*, PRL 135, 146303 (2025), DOI `10.1103/14mp-263q`：本方向 **线性 AC 与二阶金属**的核心来源。
- Jia et al., *Nonlinear Magnetoelectric Edelstein Effect*, PRL 136, 226302 (2026), DOI `10.1103/r3x2-kb93`：反铁磁自旋密度 / finite-field / PT，但**不是相同的 AC spin current**。
- Xiang et al., *Intrinsic Gyrotropic Magnetic Current from Zeeman Quantum Geometry*, PRL 134, 116301 (2025)：因果有限频率磁响应工具参考，但**B-field 驱动 charge current**，不能套公式。
- Xiang et al., *Quantum Geometric Map of Magnetotransport*, PRB 113, L201406 (2026)：Liouville / spin-orbital / E-B mixing / Ward 工具参考，physical steady-state 仍 PARTIAL。
- *Layer Edelstein Effect*, PRL (2026)，及 *Layer Nonlinear Hall Effect*, PRB (2026)：layer hidden response 和层交换先例。
- *Interlayer Electric Multipoles Induced by In-Plane Field from Quantum Geometric Origins*, Nano Letters (2024), DOI `10.1021/acs.nanolett.4c01657`：B `(k,δ)` mixed QG 近邻。
- *Symmetry-driven intrinsic nonlinear pure spin Hall effect*, PRB (2026), DOI `10.1103/swzr-55sp`，以及 *Ferroelectric switchable intrinsic nonlinear pure spin Hall current*, PRB (2026), DOI `10.1103/8htr-q899`：A3、C 必须排除的非线性 spin Hall 先例。
- *Layer-dependent control of spin Hall effect in antiferromagnetic topological insulator*, PRB (2026), DOI `10.1103/cyqr-7sgt`：layer-specific spin Hall 先例。
- *Ferroelectric Control of Layer-Resolved Spin and Orbital Transport in a Topological Heterobilayer*, Nano Letters (2026), DOI `10.1021/acs.nanolett.6c04169`：FE layer readout 先例。
- *Néel spin-orbit torque in antiferromagnetic quantum spin and anomalous Hall insulators*, (2025)：绝缘体 AFM torque 的近邻，核对是否覆盖 C2。
- *Sliding Ferroelectric Control of Unconventional Magnetism in Stacked Bilayers*, PRL 135, 056801 (2025), DOI `10.1103/dmzg-ck2t`：磁性 connection operator / 堆叠 symmetry 方法。
- 更系统的磁性 layer groups 和 antiunitary selection rules，使用现有 Group Theory Operations Toolkit 数据，但所有反幺正**频率响应约定必须从物理定义重推**。

---

## 8. 科学与协作底线

- **不把未经独立复现的八带数值、先前二倍频系数及任何未冻结程序当作 PASS。** 已证实的只有列明的单次 Pauli 矩阵代数和原论文各自可验收的受控方法子集。
- **不把“可被对称性允许”说成“必定非零”，也不把“层间转移算符随 t⊥ 为零”说成“量子几何独特机制”。**
- 用户兴趣是**群论、对称性、偏物理的材料问题**；优先简洁清晰的定理/no-go、最小模型与真正可读出物理。现阶段不做材料库/DFT/MnBi2Te4 原文重复计算。
- 关于历史 PDF、SM、研究档案和 GitHub，若无法获得要直接说明；GitHub repo 可能由其他线程持续更新，不要以本文状态覆盖 live 状态。若需写 GitHub，先按独立 studies 目录、独立分支/PR 操作。
- 只有 user explicitly asks to publish / upload brainstorm notes 才上传 GitHub。**本 handoff 仅为本地 Markdown 文档，不代表已经推送远端。**

## 9. 建议新会话的第一条消息（直接复制）

> 请阅读随附的《SpinQG_ABC_to_C2_scientific_handoff_20261021.md》，以及 Spin QG 正文 PDF 与 SM（如能访问）。不要马上写接口、做 DFT 或宣传已有新发现。先回答当前方向 C2 的最关键物理问题：在保留 m_zT 且无额外 PT/C2z 禁戒的磁性绝缘双层中，是否存在被对称性允许、必须依赖层间跃迁、并能够被识别为不同于已有 Néel SOT 的**有限频率自旋转移/横向力矩**？请首先从复频 Kubo 及自旋连续性方程独立推导 observable 和选择定则，然后给出一个低成本最小模型与 NO-GO 检验。同步查重，不要将 C1 的探索性八带数字当作已复现数据。首先给一个真正关键的科学结论，而非推进代码架构。
