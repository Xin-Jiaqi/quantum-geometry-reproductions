# DeepSeek Harness：AgF₂ B1 高性价比 idea 初筛（单次 WSL 计算）

**日期：2026-10-22；科研阶段：IDEA SCREENING，而非论文复现或生产计算。**

## 0. 直接执行，不重新规划整个项目

用户现阶段明确希望以低成本、较宽松的 DFT 筛选快速判断 BSF × IGMC 是否值得继续，不必投入完整高精度收敛。此任务由本对话的用户自行交给有 WSL 访问权限的 DeepSeek Harness 后生效：**授权仅一次本地 WSL、低成本、非 SOC 的 AgF₂ B1 离子优化初筛**，上限与安全门如下。它**不**授权任何日本超算/PBS 新任务，也不授权失败后重新提交或自动扩展到 B2、NEB、SOC 或 Wannier。

唯一协作入口：**https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38**，原分支 **research/agf2-b1-wsl-relax-20261021**。已有原始 SCF、15 离子步 CG 历史与全面 OUTCAR 取证。**不要重跑这些旧任务，不要用户复制粘贴中间输出**。最终将摘要和机器可读数据直接推回同一个 PR。

### 最直接的科研问题

在已有电子结构可自洽且 Type-I Néel 磁序保持的条件下，能不能以一次廉价且数值受控的 RMM-DIIS 尝试，找到**无明显异常接触、磁性保留、初步有隙、最大力降低到 0.08 eV/Å 以下**的 B1 筛选构型？

这个目标**不是**证明局部极小值、磁基态、滑移双势阱或 IGMC。若始终无法做到，明确判断继续优化 AgF₂ 的性价比是否偏低，而不是无限调参。

## 1. 不得混淆的科学基线

- [PR #38](https://github.com/Xin-Jiaqi/quantum-geometry-reproductions/pull/38) 先前的 **15-step B1 CG** 是 PARTIAL：电子 15/15 收敛，最终 max|F| 0.4973 eV/Å；历史最后 TOTEN −35.34931729 eV。Type-I Ag 磁矩 [-0.566,+0.566,+0.566,-0.566] μB，粗 gap ~1.445 eV（均仅既有记录）。
- 之前的 CG 遭遇 ZBRENT 大幅外推：第14步最近跨层 F–F ~2.644 Å，随后反弹，**不是 Hessian 本征软模证明**；完整 evidence 在 results/B1_CG_LINESEARCH_AUDIT.md，CSV 与脚本见同目录。
- 第15步快照为**完成了电子自洽与力计算但离子未收敛**的坐标；sha256 为 **d49c82807fc2309575e5d3757904ab78caab348aef10e895ac2327a904ac2d68**。第15步正确几何：Ag 面距 4.0443 Å；内 F plane gap **2.2119 Å**；最近跨层 F–F **3.0523 Å**；相对 y slip ~0.00515 Å。历史 2.687 Å 是 step12 旧数值，不能当末态。
- [此前独立审查](results/B1_NEXT_RUN_INDEPENDENT_REVIEW_20261022.md) 的严格纸面 P1 是 IBRION=1、POTIM=0.10、NELMIN=4、NSW=30、EDIFFG=-0.02。本文件明确将它**替换为更省资源的 brainstorm-screening 方案**；两者**不可分别执行两次**。
- 真实原始 VASP、POTCAR、CONTCAR、OUTCAR 等均在 DeepSeek 已成功使用的合法私有 WSL 环境中；**禁止向 GitHub 上传 POTCAR 内容、受限许可证、完整私有日志或私人账号/路径**。不访问/修改 PR #25 的已有 PBS 记录、PR #32 同行文献或 Xiang 原复现 baseline。

## 2. 本次批准的唯一 WSL 试算

### 2.1 起点与输入

1. 在 WSL 创建**一个新的隔离 run 目录**（绝不覆盖 first-relax 目录或 results），从已核验 SHA 的 step15 snapshot 创建 POSCAR；其 lattice、Ag/F 顺序及全部 12 个原子原样保留。原始构型并非作者 relaxed 坐标。
2. 与 first-relax 的 INCAR.relax（sha256 1d8ff5dcbdae142eef2bc06ace2a604be3dc73cc3ec821139ce2df3b795fe097）相比，只执行下列已审查的必要变更：

| INCAR 标签 | 历史值 | 本次探索值 | 理由 |
| --- | --- | --- | --- |
| SYSTEM | 原 first-relax 名称 | AgF2_B1_TypeI_idea_screen_1 | 明确避免混淆历史 |
| IBRION | 2 (CG) | **1 (RMM-DIIS)** | 避开先前 CG 大幅 ZBRENT 外推；可能失败，不保证 |
| POTIM | 0.25 | **0.10** | 更保守的迭代缩放，**非原子位移硬上限** |
| NELMIN | 未显式指定 | **4** | 支持更稳定的力历史 |
| NSW | 30 | **16** | 限制本轮计算量 |
| EDIFFG | −0.02 | **−0.08 eV/Å** | 仅用于低成本筛选，不算论文精度 |

3. **不改变** PBE、Ag-4d Dudarev Ueff=4 eV、IVDW=11 (D3 zero damping)、ENCUT=500、EDIFF=1e-6、ISPIN=2、原四 Ag 种子 MAGMOM、ISIF=2 固定晶胞、ISYM=0、12×8×1 Gamma KPOINTS、ISMEAR=0/SIGMA=0.05、LSORBIT=.FALSE.、LREAL=.FALSE.、原子顺序、合法私有 Ag/F POTCAR（concat sha256 8c0da262849433cad092946af66c1b584d457e46a70d4c1eb768ea68305f0fff）；保留原 ISTART=0/ICHARG=2，**不得假装能复用不存在的 WAVECAR**。禁止添加不适用于 IBRION=1 的 SMASS。
4. 提交前记录新 POSCAR、INCAR、KPOINTS、POTCAR **实际** SHA256（POTCAR 只提交 hash），审查是否旧 STOPCAR 残留及是否另有 VASP 运行占用资源。没有安全运行条件时标记 BLOCKED，不要强行启动。

### 2.2 资源与监测安全边界

- **仅一次 WSL VASP run**，最多 **8 MPI ranks**，一台本地用户机器，WSL 可见 RAM 历史约 7.8 GB；**最多 4 小时 wall-clock**（含优雅停止预留），约最多 32 rank-hours。本地持续会话必须可靠；若长时间无法监测/持久化，不启动。**日本 PBS Job 490620.cmspbs1 绝不重投/修改。**
- 监测每一个**完整离子步**的电子收敛、最大原子力范数、各层相对 z、最近跨层 F–F、内 F plane gap、局域 Ag 磁矩（如果输出提供）、总磁矩、TOTEN 和上一离子步最大位移。反复验证 site indexing 与 PBC。
- 如果发生**最近跨层 F–F <2.2 Å、max|F| >5 eV/Å、能量 NaN/单步非物理突变 >2 eV、电子不收敛、VASP OOM、Ag 磁序明显坍塌**或**连续多步剧烈力增长并伴随危险压缩**，仅在触发时向该 run 目录**创建** STOPCAR（LSTOP=.TRUE.）；**切忌在启动前预先创建有效 STOPCAR**。软停并记录实际完成的步骤；必要时当停机不可靠，应受控退出并明确最后完整数据块。
- 若接近 4 小时仍在运行，提前启动 STOPCAR 优雅停止；不能保证立即结束，严格记录 walltime 实况，不悄悄增加预算。不允许 second attempt、自动更换 IBRION、转日本超算，或继续直到正式 0.01 eV/Å。
- 如果同一初始构型发生单步反复能量暴涨/振荡，直接给 PARTIAL/NO-GO_FOR_THIS_NUMERICAL_ATTEMPT，**不能据此证明 B1 物理极小值不存在**。

### 2.3 通过标准：screening go/no-go

**SCREENING_PASS（不等于 stable minimum）** 必须同时：
- 电子在所有实际完成的离子步达到 EDIFF=1e-6；
- VASP 的离子停止行为与实际全原子 max norm force ≤0.08 eV/Å 对应，不把 NSW 耗尽当收敛；
- 层间无异常接触、结构未破坏；
- Ag 局域磁矩保持 Type-I sign pattern、净磁矩接近 0；
- 现有 k 网格下 occupations/带隙仍有初步有隙迹象（最后由已完成计算直接诊断）；若无法可靠确定则标 GAP_NOT_DETERMINED，不假称通过。

**SCREENING_PARTIAL**：虽未到 0.08 eV/Å，但力明显下降且结构磁序保持，可讨论进一步最小优化，而不默认继续。

**SCREENING_NO_GO_FOR_NOW**：首次 RMM-DIIS 再次出现失控、电子态/磁序不稳定、结构异常，或 16 步/4 小时仍无令人信服的趋稳信号。此时推荐优先改造物理候选筛选策略/备选材料，而不是反复修改优化器刷收敛；**no-go 仅针对投入优先级，不是材料不存在稳定相的证明**。

将结果放在已知原创性背景之下解释：Xiang 已有 IGMC，Ghorai 已有 orbital IGMC，Zhu/Peng 已有磁性滑移 FE。最终 idea 必须具备实质**可切换的材料磁电/量子几何效应**，不能仅靠静态 magnetic-connector tensor parity。

## 3. 结果直传 GitHub，帮助 ChatGPT 快速判断

把本次结果新增到同一 PR #38 的 **stage2_real_dft/b1_wsl_relax/idea_screen_1/** 下：

- **RUN_RECEIPT.json**：严格区分 NOT_STARTED、RUNNING、SCREENING_PASS、SCREENING_PARTIAL、SCREENING_NO_GO_FOR_NOW、BLOCKED；标出运行 VASP/MPI、8 ranks、实际 walltime、源/新输入 hash、实际结束原因、site magnetization、gap 来源与缺失值；未知字段 null，不杜撰。
- **TRAJECTORY.csv**：每完整离子步 TOTEN/max|F|/单点最大力的原子与分量/Ag z separation/inner-F vertical gap/nearest 3D F–F/relative y slip/net magnetization；若一项无法核验留空。
- **IDEA_SCREENING_VERDICT.md**：不超过 2 页，先用 3 句话归纳结果，再给决定性数据；评估这轮本地机时对于高性价比 idea 的信息增益。明确“尚未证明 B1 local minimum / FE double well / material IGMC”。
- **INCAR.screen + POSCAR.screen_input + POSCAR.screen_last_complete**（文本、非敏感），若 run 未完成也提供最后真实有效快照与对应 force；附 Python 解析/验证脚本或复用前一任务脚本的精确命令和 SHA。
- 若运行失败，同样提交源 hash、错误片段及原因，不因 FAIL 抹除不利证据。旧 15-step CG receipt / CG 审计 / 独立 reviewer 结论均保持不变。
- 只更新这个 PR 的原分支并在 PR Conversation 留一句摘要：**screening pass/partial/no-go、最终 max force、是否保持磁序/初步 gap、下一步唯一建议**。不要求用户把报告复制到 ChatGPT。用户下一次跟进时 ChatGPT 可以直接从 GitHub 检查这些文件。

## 4. 如何根据结果推进科研（供下一轮 ChatGPT 决策，不授权新计算）

- 如果 screening pass：**停止结构精度内卷**，先以该近似 B1 候选做**不耗 VASP**的群论/磁群/磁连接审查，明确是否允许目标 αxx/αyy、Pz 反转与 B2 连接，然后决定是否值得用一次轻量 B2/材料响应试算。
- 如果 partial：先判断数值导数、受力和几何是否在趋稳；如没有强 idea 证据，不投入长时间精修。
- 如果 no-go：将 AgF₂ 暂降优先级，对照已报道 CuF₂ 等源文献候选的实际磁群与磁电分量允许性做非 DFT 筛选；不直接开始新材料昂贵计算。

**本轮的成功不是论文结论，而是在明确的小时数预算内得到一个可信的继续/停止决策。**
