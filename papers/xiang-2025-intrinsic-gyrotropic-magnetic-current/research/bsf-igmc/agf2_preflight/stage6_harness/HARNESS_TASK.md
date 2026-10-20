# Execute with Harness: AgF₂ BSF × IGMC (Stage 6)

任务目标：获得首个**可审计的真实材料**结果，优先验证 AgF₂ 双层零应变 Type-I 两个铁电态的磁群、能隙和低频自旋 Zeeman 量子几何响应。模型几何响应与真实材料 DFT 结果必须严格分开。不要自动修改 `src/igmc_lab/`、历史冻结文件或其他论文目录。

## A. 先获取和锁定科学输入，不要盲算

1. 检索并保存 Zhu et al. PRL 135, 056801 (2025) 的**官方** SM 与公开结构材料，以及原作者提供的 AgF₂ 单层 LG17 和 B1/B2 双层原子坐标。记录 DOI/URL、版本、下载日、SHA256、单位、文件许可。`(0.5,-0.16)` / `(0.5,+0.16)` 只是相对滑移，**不能替代完整原子坐标**。禁止拿体相 AgF₂ POSCAR 冒充 2D LG17。
2. 从可核查的 SM/项目技术文档提取交换关联泛函、Ag 4d 的 DFT+U 约定（如采用）、vdW、PAW 版本、平面波截断、真空和偶极策略、磁序、k 网格。严格区分“原文给定”“Harness 假设”和“仍然缺失”；不可随意套用另一个 AgF₂ 计算的 U。
3. 如果没有双层结构但有**可靠的单层结构**，可以按作者的堆叠操作产生候选 B1/B2；必须报告生成方法、晶轴、层交换、原子映射、PBC 和滑移矢量的坐标约定，再经对称性/结构审查；只能称为“构造的待弛豫候选”，不是原作者的弛豫结构。不能用临时编造的原子坐标填补来源缺口。
4. 运行已有 `agf2_preflight.py` 输入门槛和原生 `group-theory-operations-toolkit`（固定明确版本），核对真实原子、空间群/层群及有自旋时的磁空间操作。注意 `m'_x` 与 `m_z`、`m'_z` 属于不同层级；不能用非磁空间群推定实际磁群。
5. 任一关键项缺失，**停止**并返回 `BLOCKED_SOURCE` 或 `BLOCKED_PROTOCOL` 及缺项清单。无需虚构 INCAR/CHGCAR/POTCAR，更不能声称运行了真实 VASP。

## B. 在明确允许且经审核的资源上做真实 DFT

第一批只计算 **0% 应变 Type-I** 的 B1/B2。独立于本仓库，在私有运行根目录中执行：

1. 对相同面内晶格和相同 DFT/PAW/U/vdW 参数，建立 B1/B2 结构；正确初始化各磁子晶格的 AF 自旋，检查总磁矩/层分辨磁矩和层群。放松离子时保持二维真空方向固定，不允许自动优化真空。
2. 非 SOC 下检查两种相对层磁序 Type-I/Type-II（同一结构、同一应变），从真实 OUTCAR/vasprun.xml 提取能量/力/自旋；通过 ENCUT、k、真空、SCF 和几何收敛检验。文献的 `0.26 meV/f.u.` 能差是比较目标，不是拟合指标；要提供误差条并检查 U/vdW 依赖。
3. 含 SOC 下分别比较 `[100]`、`[010]`、`[001]` 磁易轴，核对自旋方向、实际磁群与 partner mapping。对于同一磁态 B1/B2，统一坐标系、真空厚度和能量归一方式；检查 `E(B1)-E(B2)` 是否在收敛精度内相等。
4. 取得 SOC 能带：直接带隙、间接带隙、费米位置、是否绝缘。若金属态，不能以“绝缘磁电极化位移电流”来解释；须另立协议。建议先证明物理资格，再进入 Wannier。
5. Type-II 的 **4% 双轴应变**必须是独立 campaign。只有零应变 Type-I 经过初步验收后才启动，不得混合两个不同应变分支的总能或切换奇偶性。
6. 铁电切换势垒/极化：若有来源可信的 B1/B2 则可以计算二维 sheet P 与受控滑移路径（离子及电子贡献一致）；仅比较电子能带之和不等于 DFT 总能，也不等于独立证明双稳态。路径端点、晶胞、偶极与极化分支必须一致。

可复用 `wannier-workflow` 的**输入哈希锁定、运行状态机和 top-K**，但磁序、U、PAW 的设定要人为冻结。不要把 `nlo-response-workflow` 的光学核误当 IGMC 求和器。

## C. 真正的新物理：完整多带算符和响应

只有 DFT/SOC/磁群/能隙合格，才进入算符步骤：

1. 用明确声明的自旋量子化规范生成 SOC Wannier `H(k)`，并得到**与其处于同一 Wannier 规范下**的 `S_x(k),S_y(k),S_z(k)` 与完整位置 / Berry 联络矩阵。不得以 HR + Wannier centers 或对角 `S_z` 代替完整算符；不能把 Wannier 子空间不同的 B1/B2 排名由 IGMC 结果反向决定。
2. 在已知两带 benchmark 上先核对多带求和器的符号、磁场耦合和 `g μ_B` 归一化；再检查 gauge-covariant 随机幺正旋转、近简并处理、轨道投影完备性、矩阵 hermiticity、重建能带误差。
3. 计算同一实验室轴的 `alpha^{PB}_{xx}`, `alpha^{PB}_{yy}`，从其计算 AC `j_i=dot{P_i}`。记录响应的 **sheet 单位**、`g` 值、`ω` 是角频率还是 `f` 是周期频率、展宽、k 网格、温度、占据和频率低于直接带隙的范围。使用至少三组网格测试饱和性，报告 top-K 表示差异。
4. 检验 Type-I `alpha_xx(B2)≈-alpha_xx(B1)`、`alpha_yy(B2)≈-alpha_yy(B1)`。数值不显著时报告 `SYMMETRY_ALLOWED_BUT_UNRESOLVED`，不能把两侧几乎为零当成“可切换读出”。Type-II 4% 分支若通过前面门槛则核对**同号**关系。
5. 为可解释的微观机制，分解紧邻费米能级的带间矩阵元、SOC 变更、层间混合/反交叉与自旋/位置算符的变化；每个分解步骤需要规范敏感性检查。外场实际测量还要单独讨论轨道 Zeeman 与法拉第感生电场背景。

**Stage-6 的目标不是做得出漂亮的两条反号曲线，而是给出真实材料、同规范算符、可收敛的非零 `alpha^{PB}` 和物理来源。**

## D. 汇交内容（计算结束后交给 ChatGPT 分析）

请保留私有原始输入/OUTCAR/vasprun.xml/Wannier 文件不进入 GitHub；提供以下用户可审查产物：

- `stage6_receipt.json`：按照旁边 template 填写，`data_origin=real_dft` 只有在真实后端已运行时才填写；注明源码/输入哈希、算法及实际证据。
- `STRUCTURES_AND_METHODS.md`：已审查结构来源、原子映射、DFT 参数、软件版本、赝势 lock 的哈希；允许公开前应脱敏。
- `CONVERGENCE.csv`：每个状态的 ENCUT、k 网格、真空、能量 meV/f.u.、最大力、能隙、磁矩和相应误差条。
- `MAGNETIC_SYMMETRY.json`：真实磁矩及磁操作、`symprec`、layer/Cartesian setting、domain connector residual。
- `WANNIER_OPERATOR_AUDIT.json`：算符来源、h/r/spin 同一规范及 gauge covariance、数据哈希、带拟合质量。
- 如已通过响应门槛，再给 `IGMC_COMPONENTS.csv`、频率与 k 收敛和复响应相位图。

执行 `python stage6_harness/audit_receipt.py --receipt /private/receipts/stage6_receipt.json --root /private/campaign-root`。审计器会给出当前阶段性 `BLOCKED_*`/`REVIEW_REQUIRED`，**不会**仅凭 receipt 自动宣告物理 PASS。

如果第一轮 Harness 发现缺结构或参数，只返回证据和缺失项并停止重计算。用户随后把审查过的 receipt、表格和日志上传到对话，我会判断是否可以继续下一阶段。
