# DeepSeek Harness 研究执行合同：AgF₂ B1 的 CG 线搜索审计与下一计算决策

**版本：2026-10-22。角色：研究执行者 + 独立数值审计者。授权：只读现有 VASP 输出，运行轻量 Python 解析，向 GitHub 提交结果。明确禁止新 VASP/PBS 作业。**

## 0. 唯一目标与最重要的交付约定

直接在用户已配置的 WSL 环境中，对已经执行完毕的 AgF₂ B1_g3.0 单次 15 步非 SOC 结构弛豫做**原始 OUTCAR 取证、CG 优化器诊断、逐原子力与层间运动分解和下一最小计算决策**。不是重新提出规划，也不是重复运行模型或 VASP。

**研究者不必负责 AI 间转述。** 你应自行读 GitHub PR #38、连接既有 WSL、自动发现原始日志、解析与复核、制作审计文件，直接推送回 PR #38 原分支，并在 PR Conversation 留下最终结论。最后只给研究者 PR URL、审查状态和 2–3 句要点。ChatGPT 将直接从该 PR 继续研究。

- GitHub 仓库：https://github.com/Xin-Jiaqi/quantum-geometry-reproductions
- 本次唯一写入 PR：https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38
- 原分支：**research/agf2-b1-wsl-relax-20261021**，目标 main。先 fetch 远端并检查 HEAD，不依赖本文编写时的快照 SHA。禁止 force push、rebase 或覆盖别人的提交。
- 工作基础目录：**papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/**
- 历史 G1/PBS 提交：PR #25。旧日本超算 PBS job **490620.cmspbs1** 的旧 receipt 显示 Q，但不是实时状态。本任务无权 qsub、qdel、重试、自动轮询或改变该作业。
- 同行调研：PR #32。不得编辑 PR #25/#32 或冻结的 Xiang 复现、八带模型成果。
- 本 PR #38 已有完整评估、原始候选 POSCAR、15 步结果摘要、未收敛末态快照和独立几何审核，请先读，不要重做已有文献调研。

**权限边界：** 只能运行 grep/awk/Python 等低资源离线分析、校验与 Git 操作。不得调用任何 VASP 二进制、提交超算任务、启动后续结构弛豫/SCF、调整磁态或执行 B2/NEB/SOC/Wannier/IGMC。下一次真正计算须在本轮取证后经科研审核及另行授权。

## 1. 既有证据：必须保留的科学状态

- Ag4F8，12 原子，C2DB 衍生的 B1_g3.0 起始构型，不是 Zhu 等作者弛豫的结构。Type-I 层间 Néel 反平行种子。
- 本地 WSL VASP6.6.0，OpenMPI4.1.6，8 ranks；单次固定胞非 SOC 共线弛豫；输入 PBE、Ag-4d Ueff=4eV、D3 zero (IVDW=11)、ENCUT=500、EDIFF=1e-6、12×8×1 Gamma、ISMEAR=0、SIGMA=0.05；IBRION=2, ISIF=2, NSW=30, POTIM=0.25, EDIFFG=-0.02。
- 原研究者记录显示 15 步全部电子收敛，STOPCAR 软停，累计 9336.6 秒；开始/末步 TOTEN **-35.17580→-35.34932 eV**，最大力初始 0.1443、最低在 step 6 为 0.0510、step 14 为 0.696、step 15 为 **0.4973 eV/Å**。**力未收敛。**
- Ag 局域磁矩末态 [-0.566,+0.566,+0.566,-0.566] μB，净磁矩约0，Type-I 磁序在本次轨迹保留；无 SOC，不能决定易轴或磁基态。
- 最后占据阈值的近似 gap ~1.445eV，未做 k 收敛，不是定量材料能隙。
- Step 14 TOTEN -35.25755 eV，最近跨层 3D F–F **2.644 Å**；step 15 最低 TOTEN，但末态受力 0.4973。单独“能量最低”**绝不证明**局部极小。
- 公开原始 POSCAR 与 step15 快照的独立几何计算已确认：Ag 层间垂直分隔 4.969271→4.044265 Å；**内侧 F 平面垂直 gap 3.000000→2.211875 Å**；最近跨层 F–F **3.68056→3.05226 Å**；两 Ag 层相对 y 滑移 **约 0.00515 Å** 或 0.000876b。旧 receipt 的 F-plane gap=2.687 Å 与上述定义不一致；旧报告相对滑移=0.03 Å 存在换算错误。
- **未证明**：B1 结构局部极小值、Hessian 软本征模、正的滑移曲率、B1/B2 铁电双阱、竞争磁态基态、磁电极化率或材料 IGMC。

首先必读 PR #38 已有文件：
1. results/B1_RELAX_RUN_RECEIPT.json，results/B1_RELAX_SCIENTIFIC_ASSESSMENT.md；
2. results/G1_FORCE_AND_GEOMETRY_AUDIT.md，results/POSCAR.relaxed，results/INCAR.relax；
3. results/B1_POSTRUN_INDEPENDENT_REVIEW_AND_READONLY_TASK_20261022.md；
4. results/B1_INDEPENDENT_COORDINATE_AUDIT_20261022.md 和 results/audit_agf2_b1_positions.py。
以上 results 相对目录为 **stage2_real_dft/b1_wsl_relax/**。

## 2. 阶段 A：原始文件来源审计（先证明你读的是什么）

A1. 在用户已连接、授权的 WSL 研究目录中自动定位**这次 B1 15-step relaxation**，而不是 G1 单点或其他作业的 OUTCAR、vasprun.xml、OSZICAR、INCAR、KPOINTS、POSCAR、CONTCAR、stdout 和私有 POTCAR。核对晶胞、12 原子、计算时间、代码构建信息、离子步数与电子自洽结束标记。仅提交非敏感路径代号、文件大小、sha256，不向 GitHub 复制受限材料。

A2. 与锁定输入比对：original G1 POSCAR sha256 **2d323c3fe2353ef984330b3b99fc24a75a85337883bbaf36d08b20ed45e6d60a**；KPOINTS **4ab9ad06699176341f30e79f18b822be389d066f64e0295ce368e74b8b99420a**；POTCAR **8c0da262849433cad092946af66c1b584d457e46a70d4c1eb768ea68305f0fff**；本次 INCAR.relax **1d8ff5dcbdae142eef2bc06ace2a604be3dc73cc3ec821139ce2df3b795fe097**。POTCAR 只能验证 hash，禁止上传或显示内容。公开 step 15 POSCAR.relaxed 的 sha256 应为 **d49c82807fc2309575e5d3757904ab78caab348aef10e895ac2327a904ac2d68**。

A3. 本地原始输出若不存在、hash 不匹配、目录冲突无法消歧或 OUTCAR 不含所述 15 completed steps，立刻标记 **BLOCKED_SOURCE_INTEGRITY**；提交一个有证据的缺失报告，不猜测后续数值，不让用户复制日志到聊天窗口。

## 3. 阶段 B：CG line-search 完整诊断（核心）

B1. 编写简洁稳健的 Python 解析脚本，按 VASP 6.6.0 实际记录索引全部 ionic configurations；从 OUTCAR 提取每一步的 **12 原子 POSITION 和 TOTAL-FORCE**，相应的 TOTEN、energy without entropy、energy(sigma→0)、电子收敛标志，并尝试与 vasprun.xml 的 force/structure 和 OSZICAR 的 F/E0/mag 对照。**同一个数值必须对齐同一 completed ionic step**。原日志缺少的值填 null，不补造。

B2. 专门提取 CG line-search 标记 **trial / trialstep / g(F) / gam / ort / ZBRENT 等实际出现的字段**，保留行号和必要短上下文。VASP IBRION=2 有试探/修正过程，不要将 15 个 TOTEN 点视为等间距抽样势能曲线，也不要在没有 Hessian 时宣称“软模本征不稳定”。不同 VASP build 缺字段就报告 NOT_PRESENT。

B3. 输出 15-step 总表：**step、F(TOTEN)、E0、电子收敛、最大力矢量范数、最大 |Fz|、最大 |Fy|、最大力原子 ID、Ag 层 z 间距、内侧 F-plane 垂直间距、最近跨层 F–F 三维距离、Ag 层间 y registry、总磁矩、CG trial/corrector 状态**。全部单位明确，PBC x/y 使用最小像，z 真空方向不取跨盒虚假邻居；结构因 PBC wrapping 产生突跳时对真实原子位移 unwrap。

B4. 导出 **15×12 每原子力和坐标**，分解 top/bottom 全层、各层 Ag/F 子集的总力与平均力；投影两层相对 z 压缩、相对 y 滑移和层内 Ag–F buckling，给出清楚的向量定义/单位；不要把总力投影当作 Hessian 特征向量。突出解析 steps **6、13、14、15** 中哪些原子及 Fx/Fy/Fz 使最大力由 0.051 增至 0.696，为什么 step15 能量比 step14 低但仍有大力。无法根据现有记录证明的机制归为备选解释，而不是事实。

B5. 特别检查 STOPCAR 的真实停止语句与完整的最后 ionic step：提交的 POSCAR.relaxed 是否逐原子对应第15个 **completed** 的 OUTCAR force block？原 CONTCAR 是否确实未写出或全零？现有的“末态”是否有 predictor-only、trial-only 的疑点？将重启用结构与原始快照严格区分。本阶段不得写回原始私有输出文件。

## 4. 阶段 C：交叉核验与纠错

C1. 必须核验 OUTCAR 与 vasprun.xml 逐原子受力差异（最大绝对误差）、以及至少四步的能量、层距和最大力；若解析错误，对脚本和合成测试修复并重新解析**已有**原始日志，不能通过新跑 VASP 弥补。

C2. 必须逐项解释：
- 为什么旧 receipt 的 **vertical_F_plane_gap_A=2.687 Å** 与已经提交的第15步 POSCAR 计算的**内侧 F-plane gap=2.211875 Å** 不一致。可能是定义不同，也可能是笔误；无来源时标记原字段不可用于此定义，**保留历史记录、另加 dated erratum**。
- 旧报告的相对滑移 0.03 Å 不正确。由 -0.3400→-0.3391 乘以 b=5.87966575 Å 约为 **0.0053 Å**；用未四舍五入原子坐标则约 0.00515 Å。分别解释各层绝对移动与层间相对移动的定义。
- 从最低能量点、最后保存的 POSCAR 或 15/15 电子自洽，**都不能推出离子收敛**。第15步 max|F|=0.497 eV/Å，相对 EDIFFG=0.02 超标约25倍。

C3. 用审查级别区分：**RUNNER_REPORTED**（原始数值由 Harness 提交）、**INDEPENDENTLY_RECOMPUTED_FROM_COMMITTED_COORDINATES**、**CROSS_VERIFIED_FROM_RAW_WSL_OUTPUT**、**INFERRED_HYPOTHESIS**、**NOT_COMPUTED**。不要混用；磁基态、SOC easy axis、Berry Pz、滑移 FE 和 IGMC 目前都属 NOT_COMPUTED。

## 5. 阶段 D：下一次唯一最小 VASP 方案，只做纸面设计

先依据真实 trialstep 和原子力判断造成剧烈压缩的原因，以及第6/13/15步中哪一个**只适合做试探起点**。给出 **唯一首选** + **一个备选**：例如从可信的已完成 snapshot 进行更稳健的单次 B1 非 SOC 离子优化，或若 z 方向路径无法控制则小规模约束层间位移诊断。**不得自行运行任一方案。**

方案必须包含：为什么这样改，初始几何 SHA 和来源，建议 INCAR diff（不是偷偷改计算结果），力收敛与磁序检查，允许的 8 MPI ranks/WSL 约7.8GB RAM 资源边界，walltime 量级与安全 STOPCAR 规则、失败及异常接触停止条件、预期科学信息增量。历史 15-step 用时 2.6h 只能用于保守估算，不是新运行耗时保证。保留原 PBE、Ueff=4、IVDW=11、PAW 锁定、固定晶胞、无 SOC。优化器选择（IBRION/POTIM）必须由 line-search 证据支持，不可空泛地说“减小步长必然收敛”。

**下一条计算的目标只能是可证伪的 B1 力收敛局部极小值**，先于 B2、NEB、磁各向异性、Wannier、IGMC。不要因本次看到“能量降低 0.17 eV”和“磁序保留”便扩权。

## 6. 阶段 E：GitHub PR #38 自动交付，用户无需中继

在本 PR 的 **stage2_real_dft/b1_wsl_relax/results/** 内增量提交以下内容：
1. **B1_CG_LINESEARCH_AUDIT.md**：结构化科研审查，CG line-search 的源证据、步骤对齐、受力投影和机制判断，含结论边界。
2. **B1_CG_TRAJECTORY.csv**：15 步定量全表，可从原始文件离线重算；原记录缺数据留空/标 null。
3. **B1_CG_ION_FORCES.csv**：15×12 的 (atom index, species, xyz, Fx/Fy/Fz, step)，缺失不得插值。
4. **B1_CG_AUDIT.json**：status、输入与原输出 hashes、版本、逐项一致性检查、差异修正清单、科学结论证据等级。
5. **B1_NEXT_COMPUTATION_PROPOSAL.md**：唯一首选方案、一个备选及 no-go gate，明确 PROPOSED_NOT_EXECUTED。
6. 可复用 **scripts/parse_vasp_b1_cg.py**、最小合成单元测试、运行方法；只运行 Python，不需 VASP/POTCAR 参与。
7. 在**同一个 PR #38 Conversation** 留下最后摘要和上述新文件的相对链接，方便 ChatGPT 直接抓取审查。不要新开多个 PR 或把报告留在本地。

隐私限制：**不可上传 POTCAR 内容、WAVECAR、CHGCAR、完整 OUTCAR/vasprun.xml、私有主机名、IP/账户/路径、秘钥、收费账单等**。保留原始私有文件，只公开必要的派生数表、hash、物理输入记录和最小非敏感日志上下文。保留旧 B1 run receipt，**禁止静默覆写历史数值**。

## 7. 三种退出路径

**PASS_READONLY_AUDIT：** 真实日志追溯、trial/corrector/force 对齐、几何差异解释、单位与坐标复核、脚本测试、计算方案提出均完成，直接推送 PR #38。

**PARTIAL_READONLY_AUDIT：** 少量字段缺失或日志版本差异，但剩余数值有明确证据；完整标记缺失、提交有用成果和阻塞原因，不制造伪精度。

**BLOCKED_SOURCE_UNAVAILABLE：** WSL 原始 OUTCAR 不存在或无法合法读取/无法保证 provenance，只提交已有 GitHub 来源审查和需要恢复的最小文件清单；不提交假解析产物。

最后答复用户仅需：**“PR #38 已更新；审查状态 PASS/PARTIAL/BLOCKED；关键物理发现；推荐的下一单次计算方案（仅建议，未运行）。”**

**交接规则：GitHub PR #38 是 DeepSeek Harness 和 ChatGPT 的唯一可复查信息渠道，用户不承担中间复制粘贴工作。**
