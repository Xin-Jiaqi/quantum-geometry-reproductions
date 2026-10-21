# DeepSeek Harness 执行任务：AgF₂ B1_g3.0 首条真实结构弛豫（WSL 优先）

**日期：2026-10-21**  
**执行方式：单次授权的本地 WSL 任务，不是日本 PBS 追加 qsub**  
**目标研究：** BSF × IGMC，验证真实 AgF₂ Type-I 双层在当前 PBE+U+vdW 下能否保持磁有序、初步绝缘性并收敛到合理的局部结构。  
**执行结果直接回传 GitHub 分支：** \`research/agf2-b1-wsl-relax-20261021\`；如已有 PR，更新同一个 PR，**不再要求研究者在 DeepSeek 和 ChatGPT 之间转述/复制中间消息**。

> **对 Harness 的核心指令：这是一个执行型任务，不是计划书。** 直接使用研究者现有 WSL 连接、已成功的本地 VASP 6.6.0 环境、私有 POTCAR 与已完成 G1 单点的本地文件。自动发现并审计实际路径，先读取已有单点受力；若适合，再只运行**一次** B1_g3.0 非 SOC 固定晶胞离子弛豫。在原有结果不可获得或资格不满足时，用证据解释停止，不猜测结果。只给用户最后 GitHub PR 链接和一句科研结论。

## 0. 项目事实：不要重启、不要误读

- GitHub 仓库：<https://github.com/Xin-Jiaqi/quantum-geometry-reproductions>
- 旧任务 **PR #25**：<https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/25>，head 在本任务编写时为 \`ca5c5eefc8588a37fcb4c602f29fb195394d7d89\`。其中 \`G1_LOCAL_WSL_DIAGNOSTIC.json\` 记录本地 VASP 6.6.0 SCF 完成；**不能**把此事说成 PBS 成功。
- 同行研究 **PR #32**：<https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/32>，包含文献/原创性审查和 G1 科学状态更新。它与本次计算任务平行，**不要改动或重写**。
- 日本东北大学 IMR 服务器采用 **OpenPBS**。原 PBS Job \`490620.cmspbs1\` 已提交一次，其最后一次 GitHub receipt 是 **Q、未获真实 PBS OUTCAR**；该快照不是实时查询。**禁止新增 qsub、重投、取消或修改这个作业。** 若有现成访问可对该 Job 做一次只读 \`qstat -f/-x\`，不是本任务成功的前提。
- 本地 WSL 已运行相同输入：VASP 6.6.0 (06Mar2026 / 2026-09-02 build)、OpenMPI 4.1.6、8 MPI ranks、约 1195 秒；本地可用 RAM 约 7.8 GB，历史峰值约 3.2 GB。电子 EDIFF=1e-6 已到，27 Davidson，TOTEN −35.17579601 eV / 12-atom cell，Ag 顺序磁矩 [−0.567,+0.567,+0.567,−0.567] μB，净磁矩 0，粗略 occupation-gap ≈1.49 eV。**都是未弛豫的 non-SOC 单点，不是可发表材料结论。**
- 已构造的 Ag4F8 B1/B2 基于 C2DB 单层 \`2AgF2-1\` 和 Zhu et al. PRL 135, 056801 (2025)；B1 shift (0.5,−0.16)，B2 shift (0.5,+0.16)。3.0 Å 是人为 nominal interface F–F plane-gap **起始猜测**，不是实际最近 F–F 距离或平衡层距。Type-I 面内磁矩模式对应两层 Néel 反平行；不要再用极矢量镜面规则处理真实自旋。
- 原始冻结输入在 \`papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage2_real_dft/one_shot_g1/\` 及 \`../stage2_geometry/candidates/B1_g3.0/POSCAR\`。先实际读取路径，不可假定本地目录布局。G1 INPUT hashes：
  - POSCAR \`2d323c3fe2353ef984330b3b99fc24a75a85337883bbaf36d08b20ed45e6d60a\`
  - INCAR \`8ca091bbb4c59572b235fe201a71509c40d6ae84fba67b0f2bfbade34a025166\`
  - KPOINTS \`4ab9ad06699176341f30e79f18b822be389d066f64e0295ce368e74b8b99420a\`
  - POTCAR 合并 sha256 \`8c0da262849433cad092946af66c1b584d457e46a70d4c1eb768ea68305f0fff\`（**仅验证哈希，禁止上传内容/片段或无授权分发**）。
- 项目研究目标是后续 *sliding FE → magnetic connector → low-frequency spin-Zeeman α^PB → AC displacement response*；本次只验证 B1 结构和磁性的计算可行性，**不计算或声称任何材料 IGMC**。

## 1. 先做 G1 力和几何的**只读审计**（必做，不耗新 VASP）

1. 自动找到用户本地完成 G1 的 **原始** \`OUTCAR\`、\`vasprun.xml\`、\`OSZICAR\`、\`POSCAR\`、\`INCAR\`、\`KPOINTS\` 和 **私有** \`POTCAR\`。不可把别的材料文件认作此 G1；按 sha256 四项匹配。若源文件找不到、散列不符、VASP 未收敛或实测没有 12 原子，**STOP**，报告事实，无新计算。
2. 从最后有效电子步/对应 ionic step 提取 12 原子的 \`TOTAL-FORCE (eV/Angst)\`，输出每个原子 \((F_x,F_y,F_z)\)、\(\max_i |\mathbf F_i|\)、\(\max_i |F_{iy}|\)、作用于各层的总力/净转矩（若层映射可靠），检测最大的残余受力方向。对照 \`vasprun.xml\` 交叉检查力矩阵（若解析可得），有歧义就标记，不混用单位。
3. 用 POSCAR PBC 最小镜像计算**真实最近跨层 F–F 三维距离**，正确区分 nominal *plane gap* 3.0 Å 与最近距离；记录每层 z-centroid、层间偏移和构造时 \(t_y=-0.16\) 的对应。最近距离应基于真实笛卡尔坐标及 x/y 周期；不要穿越真空 z 边界找假近邻。
4. 从旧计算确认最后 Ag 局域自旋模式与 net M、occupations、带隙诊断和输入 lattice；整理 \`G1_FORCE_AND_GEOMETRY_AUDIT.md\`。该审计即便下一步被阻止，也必须作为独立有用结果提交。
5. 若 \(\max|\mathbf F|<0.02\) eV/Å，先报告“在本次筛选容差内近似静止”，**不再浪费资源重复同一几何优化**；真实 0.01 eV/Å 收敛和磁基态仍是下一关。如果存在显著近接触（如跨层 F–F 距离 <2.2 Å）或 \(\max|\mathbf F|>5\) eV/Å，先标记 \`STRUCTURE_SUSPECT\`，**停止计算**等待科学审查，而不是让 VASP 强行优化异常结构。以上阈值是本次执行安全阈值，不是论文物理阈值。

## 2. 条件满足后执行**一次** WSL B1 非 SOC relaxation

**这一条研究者指令授权的是：本地 WSL 上**仅 1 个** \`B1_g3.0\` 固定晶胞离子弛豫（不包含 PBS 新作业）。** 不需要反复再问“是否可以启动”；但如果实际需改用收费超算、扩大资源或超出以下边界，停止并说明需新授权。

从已校验的**原始 G1 POSCAR** 创建**新的隔离计算工作目录**（不要改原 G1 输入、原 OUTCAR/CONTCAR，也不要覆盖已冻结 results）：

- 固定面内晶格和包含真空的整体晶胞（\`ISIF=2\`），非 SOC、共线 \`ISPIN=2\`；先允许 12 个原子位置充分弛豫，**不预先固定滑移方向**，以检验 B1 是否会自然漂移。检查来源 POSCAR 是否已有 selective dynamics，避免无意冻结原子；保留绝对坐标系和结构来源。
- **只变更原 G1 INCAR 中与弛豫必需相关的设置**：\`SYSTEM=AgF2_B1_g3.0_TypeI_nosoc_WSL_relax1\`，\`IBRION=2\`，\`ISIF=2\`，\`NSW=30\`，\`POTIM=0.25\`，\`EDIFFG=-0.02\` (eV/Å，本次廉价试算目标，论文级将需 <0.01 eV/Å)。禁止 \`ISIF=3\`，不优化真空或晶格。
- 保留 \`GGA=PE\`、\`ENCUT=500\`、\`EDIFF=1E-6\`、\`ISPIN=2\`、\`MAGMOM=-1 1 1 -1 8*0\`、\`LDAUTYPE=2\` / Ag-4d \`U_eff=4 eV\`、\`IVDW=11\` (D3 **zero damping**) 、\`ISMEAR=0\`/\`SIGMA=0.05\`、\`ISYM=0\`、\`LORBIT=11\`、\`NCORE=4\`、\`LREAL=.FALSE.\`、\`ALGO=Normal\` 和所有原始势/元素顺序；首轮沿用 **12×8×1 Gamma screening**，明确 **不是**收敛生产网格。若系统版本要求细微 syntax 修正，留下精确 diff；绝不凭空改物理模型。
- 由于原 G1 \`LWAVE=.FALSE.\` / \`LCHARG=.FALSE.\`，默认从头电子自洽（\`ISTART=0\`/\`ICHARG=2\`）；**不要假装存在可复用的 WAVECAR**。
- 只用本地已验证合法 VASP6.6.0 二进制和 OpenMPI，最多 **8 MPI ranks**、约 7.8 GB WSL RAM，不在本地启动多个 VASP 并行任务。运行时采样内存与 CPU，若发生 OOM/失控内存、严重 SCF 发散、原子重叠、磁态坍缩或不合理结构，**停止且保留现场**；不要自动改变 U、MAGMOM、smearing 或重试。
- 使用 **tmux / systemd-run / 已验证可靠的持久会话**，避免 SSH 断线使 VASP 被杀；已发生过一次这种中断。可设置**最长约 8 小时**的安全观察边界；在达到该时间时优先通过该独立目录的 \`STOPCAR\` 让 VASP 完整结束当前 ionic step（必要时再人工受控终止），严格记录是否仅得到 partial relaxation。**禁止默默重启/连续尝试，禁止把达到 \`NSW=30\` 上限当作力收敛。** 可以中途向 PR 上传不含隐私的 RUNNING 状态，此时不要声明最终科学结论。
- 如果 WSL 实际内存/运行环境不安全或无法保障长时持久化，就只交付力审计和 \`BLOCKED_RUNTIME\`，**不要自动转去日本超算 qsub**。

## 3. 唯一的科学验收

本轮判据是：“初始 B1 是否在同一能量模型中走向结构合理的局部近似极小值，同时目标 Type-I Ag 磁性是否仍然保留？”

至少给出：

1. 初始/末态 TOTEN、\(E_0\)（若输出）、\(\Delta E\)（同设置且相同胞）、离子步数、所有电子步是否真的收敛、最大最终力及 force criterion 是否通过，是否发生明显层间距离/滑移 \(t_y\) 变化。
2. 初始/末态 Ag 逐原子磁矩、总磁矩、F 磁矩、初步绝缘 gap 与 occupations；磁序若变成另一种，**这是物理结果，不得重新初始化后偷偷替换**。非 SOC 结果不能决定 easy axis。
3. 相同晶胞的完整末态坐标、Ag/F 顺序、origin、两层 z-centroids、nearest interlayer F–F separation、是否出现反常近接触和横向滑移。
4. 结构对称性和磁群：在适当容差下可做初步空间映射；如果 VASP relaxation 打破原 \`pb11\`，明确给出数值精度/磁畴或结构改变的可能性，**不要强行 symmetrize 成预期 B1**。真实 spin symmetry 没有 SOC 的磁空间群需说明与后续有 SOC 的易轴选择不同。
5. 若优化失败、超时、撞到 \`NSW\`、电子/磁态失稳，必须给出 evidence-based **FAIL / PARTIAL / BLOCKED**，而不是继续拓展为三初猜、B2、4% 应变、NEB、SOC、Wannier 或 IGMC。

本次即使结构和局域磁性 PASS，也**不能宣称**磁基态、滑移铁电双势阱、极化数值或 IGMC 量级已经证明。最终 DFT production target 是 <0.01 eV/Å、多磁态、k 和数值收敛；不应在这个 cheap pilot 中一次做完。

## 4. GitHub 回传：让 ChatGPT 直接读取，用户零次中继

- 使用当前仓库本次独立分支 \`research/agf2-b1-wsl-relax-20261021\`。从远端 fetch 并确认该任务文件存在；不要 force push、rebase 或覆盖别人的提交。建议所有输出保留在：

  \`papers/xiang-2025-intrinsic-gyrotropic-magnetic-current/research/bsf-igmc/agf2_preflight/stage2_real_dft/b1_wsl_relax/results/\`

- **必须至少提交**：
  - \`G1_FORCE_AND_GEOMETRY_AUDIT.md\`：已有单点受力、几何真实检查。
  - \`B1_RELAX_RUN_RECEIPT.json\`：运行状态 (\`NOT_STARTED|RUNNING|CONVERGED|PARTIAL|FAILED|BLOCKED\`)、源输入 SHA、修改参数、实际 VASP/MPI 版本、MPI ranks、耗时、校验与出错标记；不能把任务提出时的计划当结果。
  - \`B1_RELAX_SCIENTIFIC_ASSESSMENT.md\`：起始/终态对比表、结构、磁序、gap、能量和符合/不符合的物理判断。
  - 新的**非敏感** \`INCAR.relax\`、\`KPOINTS.relax\` 和收敛时的 \`POSCAR.relaxed\`，可附简洁 force/magnetization 数据表；input hashes 必须实际计算。
- 不上传 **POTCAR（或任何其内容片段）、WAVECAR、CHGCAR、vaspout.h5、私有主机名/密钥/账户/日志敏感路径**，不上传大量原始 OUTCAR/vasprun 或未经检查的版权所有数据；原始完整输出只存本地，使用 sha256 + 已证实的受控路径形式留引用（公开 MD 中不要含私有绝对路径）。
- 若有权限：推送到同一个本次 branch，创建**一个** Draft PR（若已存在则更新），目标 \`main\`，描述本次只涉及 B1 单次 relax。更新结果即使 FAILED/PARTIAL 也应形成审计记录。不要动 PR #25 的独立 PBS receipt 或 #32 的同行文献，避免“结果同步”时篡改物理历史。
- 只需给用户**最终**一句话摘要 + GitHub PR URL；不要发送多段中途“请转告 ChatGPT”的消息。GitHub 上的报告就是交接。
- 若任务继续在后台运行，请先上传 \`RUNNING\` 收据并告知用户同一个 PR URL；**不得声称 Harness 会自行在聊天之外继续交互**。若已部署可靠本地后台进程/日志收集，明确目前实际存活的进程及可查询方式，不作未经验证的未来承诺。

## 5. 目标科学节奏（非本次执行授权）

G1-force → **单次 B1 relaxation** → 对局部极小值与 magnetic state 作 go/no-go → **另外批准之后**再比较 B1 多初猜及 B2 镜像双态和真实 Pz → SOC Wannier H/spin/position operator → spin-only αxx/αyy independent validation → orbital/Faraday background。不要跳关，也不要为“项目看起来完整”而一次性提交多作业。
