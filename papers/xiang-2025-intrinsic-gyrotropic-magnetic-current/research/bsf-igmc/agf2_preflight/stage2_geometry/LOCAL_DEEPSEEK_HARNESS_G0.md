# AgF₂ Stage 2 · DeepSeek Harness（本地 WSL → 日本 PBS）单次执行入口

**执行状态：G0 ONLY / NO JOB SUBMISSION。** 本文件供运行于研究者本机的 DeepSeek Harness 使用。它可以通过研究者**已有的 SSH 配置**连接日本超算进行只读预检，准备计算输入，但**不得提交或启动 DFT 作业**，除非研究者随后明确批准具体资源用量、队列与运行参数。请不要将仅一句“开始运行 Harness”解释为无限计算授权。

## 任务目标

把现有 Stage 2 的三个 B1 初始候选（界面 F–F 垂直间距 2.6、3.0、3.4 Å）转成**可审查但未提交**的日本 PBS + VASP 几何试算方案，并给 ChatGPT 一个可以通过 GitHub 独立复核的最小、非敏感报告。

首先阅读本目录：
- `JAPAN_PBS_WSL_RUNBOOK.md`（计算环境、物理边界与 G0）
- `CALCULATION_GATE.md`（G0/G1/G2/G3 关卡）
- `EXECUTOR_BRIEF.md`（最小物理任务）
- `README.md` 和 `candidates/GEOMETRY_MANIFEST.json`（构造结构及来源）
- `../stage6_harness/stage1_review/STRUCTURES_AND_METHODS.md`（Stage 1，磁连接：Type-I = 层间 Néel 反平行）

**只在本机已有的 WSL/SSH 环境执行；不得请求或回显 SSH 私钥、密码、登录令牌、PBS 账户信息，或上传 POTCAR。** 不允许修改已有复现代码、冻结数据及 Stage 1 的原始结构。

## Step A：本地环境及只读 SSH 探测（现在可执行）

1. 检查本机 WSL/SSH 工作环境，从现有 SSH 配置识别**用户已使用的日本 PBS 连接别名与双跳方式**。不要猜主机地址；若有多个可能目标，停止并向用户询问别名。SSH 失败后不得自行尝试提权或其他账号。
2. 用用户已有凭据通过 SSH 连接到**登录节点**。只读运行：`hostname`、`command -v qsub`、`command -v qstat`、`qstat -Q`（如果站点允许）、`type module`、`module avail vasp`、`command -v vasp_std`、`command -v vasp_ncl`、`command -v mpirun`。必要时读取 PBS 队列帮助、软件 module show；**不要**调用启动 VASP 计算的命令，不运行 `qsub`、`qdel`、`sbatch`。
3. 核实当前站点许可、VASP 实际可执行版本（历史记录 VASP 5.4.4/6.4.3 **不能直接当作当前事实**）、队列的节点/核数/内存/墙钟上限、PBS 资源语法。不要假设 Slurm。
4. 如能在**许可范围内合法访问** Ag、F PAW，**在私有端**检查 `TITEL`、`ENMAX`、元素顺序和 SHA-256；只保留版本/哈希/可用性用于对照，**永远不提交 POTCAR 全文或其原始物理路径**。不得将未知 PAW 版本“猜成”某一版。
5. 无法验证的事项保持 `null` 或 `BLOCKED_*`，不能为了让测试通过伪造数值。

## Step B：只准备可审查输入，不提交计算

- 仅以本目录 `candidates/B1_g2.6`、`B1_g3.0`、`B1_g3.4` 的未弛豫 POSCAR 为输入，验证每个文件的 SHA256、12 原子 Ag₄F₈、物种块顺序 `Ag F`、晶胞、20 Å 真空及结构来源。
- 物理目标限定在 **0% 应变、Type-I 层间 Néel 反平行、共线自旋、非 SOC**、固定面内晶格的离子弛豫；除非明确记载理由，原子滑移自由度**不得被人为固定**。自旋是轴矢量：对 `m_z` 镜面连接，面内磁矩反号。
- 文献约束仅冻结：PBE、Ag-4d Dudarev `U_eff=4 eV`、DFT-D3 方法族、ENCUT 500 eV、约 20 Å 真空、`EDIFF=1e-6 eV`、最大残余力 < 0.01 eV/Å。具体 PAW、`IVDW` 选项（D3 zero/BJ 差别）、`ISMEAR/SIGMA`、`MAGMOM` 大小、`ISYM`、`LDIPOL`、电子混合和队列脚本是**待核实的操作者选择**，请逐项列出依据或待确认状态。
- 使用 `12×8×1` 作为初步廉价 k 网格的**试算候选**，`14×9×1` 和 `16×10×1` 仅为后续收敛备选；不要宣称它们已经收敛。
- 在 WSL **独立的私有工作目录**中准备三份 `POSCAR/INCAR/KPOINTS` 和**未执行的** PBS 作业模板。若站点参数、POTCAR、实际 VASP 版本或队列约束不明确，模板保留 `TODO` 并标 `BLOCKED_PROTOCOL`；不要填造参数或拼接假的 PAW。
- 本阶段不运行 SOC、Type-II、Berry 相极化、NEB、Wannier90、IGMC，不产生任何声称真实的能量、带隙、磁矩或位移电流。

## Step C：交付与停止条件

生成并校验（不包含秘密）：
1. `SITE_G0_REPORT.json`：站点与 PBS 是否已验证、软件版本、Ag/F PAW 是否有来源/版本/哈希锁、可用资源边界、`actual_jobs_submitted=0`、是否得到用户的具体计算预算授权（本任务缺省 `false`）。
2. `INPUT_PROTOCOL_DRAFT.md`：三份输入文件摘要和 SHA-256，必须区分文献参数、选择参数与缺失参数；附 Type-I 初始磁矩排序及各原子映射依据。
3. `READY_FOR_QSUB_CHECKLIST.md`：首批三组作业的前置条件、预计核时、失败判据、建议先运行的最小单点/弛豫测试；**只建议、不运行**。
4. `G0_RECEIPT.json`：`status` 为 `PASS_INPUT_READY`（仅表示输入可供人工批准）或 `BLOCKED_*`，`data_origin=NOT_RUN`、`submitted_job_ids=[]`、`physics_approved_by_user=false`，不得声称 G1 计算已成功。

把报告（不含私人主机/IP、登录账号、路径、作业账户、许可文件等）放在本研究 `stage2_real_dft/preflight/` 的**新分支**并开 Draft PR，或者保存到本地待研究者手动上传；**不能直接合并 main**。私有 PBS/VASP/POTCAR/SSH 内容不得发布。未经明确批准不得运行 `qsub`。

**给用户的最终反馈应极简**：G0 是 PASS 还是 BLOCKED；缺哪项；哪个文件可以交给 ChatGPT 审查；是否真正提交了作业（应为“否”）。不要重新安排繁琐的 Harness↔ChatGPT 人工往返。
