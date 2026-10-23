# AgF₂ B1：CG 原始审计的独立复核，以及下一条计算方案修正

**2026-10-22｜READONLY_AUDIT_ACCEPTED / NEXT_RUN_PROPOSAL_REVISIONS_REQUIRED / PROPOSED_NOT_EXECUTED**

审阅基线：PR #38，DeepSeek commit c4dec16dc5a536c3a39dd84d54a0bcc5c0e0453a。审阅了 B1_CG_LINESEARCH_AUDIT.md、B1_CG_TRAJECTORY.csv、B1_CG_ION_FORCES.csv、B1_CG_AUDIT.json 和 B1_NEXT_COMPUTATION_PROPOSAL.md。**本评审未运行 VASP，不授权新作业，未访问私有原始 OUTCAR。**

## 1. 本轮已独立确认的数字

提交的原子力 CSV 含 15 离子步 × 12 原子 = **180 行**。从已提交的 Fx/Fy/Fz 独立计算每个原子的欧氏受力范数，逐步取最大，结果与 15 行轨迹表一致，最大差值 **0.000041 eV/Å**（表格四位小数舍入）。15 步最大受力原子 ID 也全部一致。

| step | 独立重算 max norm (eV/Å) | 上层 ΣFz (eV/Å) | 下层 ΣFz (eV/Å) | max|Fy| (eV/Å) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 0.144234 | −0.21978 | +0.21979 | 0.04839 |
| 6 | 0.050990 | −0.22369 | +0.22369 | 0.00487 |
| 13 | 0.439463 | −0.32545 | +0.32545 | 0.24656 |
| 14 | 0.696259 | +0.53029 | −0.53029 | 0.48807 |
| 15 | 0.497281 | −0.31524 | +0.31524 | 0.29157 |

分层：Ag1/2 + F5–8 为上层，Ag3/4 + F9–12 为下层。上表证明已提交数据的内部数值自洽，以及第14步力向反转是**真实提交轨迹中的现象**。Harness 报告的 OUTCAR/vasprun 直接交叉符合约 7.67×10⁻⁷ eV/Å，本 ChatGPT 审查没有独立访问 WSL 原始文件，必须保留证据级别区别。

CG 外推过冲是最合适的数值机制假设；尚无 Hessian，不得声称本征软模或结构已达平衡。文献审查原文对次数有一处细节：ZBRENT increasing intervall **总共六次**（包含早期 3→4），其中 9→10、10→11、11→12、12→13、13→14 **连续五次**，不是“9–13 连续六次”。curvature 显示星号首先代表输出格式饱和，不足以直接诊断 Hessian 奇异。

## 2. 下一次计算提案的七项必要修订

1. **IBRION=1 可以尝试，但不保证稳定。** VASP RMM-DIIS 使用历史力构建近似逆 Hessian，接近极小值时高效、初态较远时可能失败。末态 max|F| 仍约 0.497 eV/Å，不能假定已接近 minimum。替换 CG 的 1D 外推只是合理的数值试验，不是成功保证。
2. **绝不能在启动前预先放入有效 STOPCAR。** 含 LSTOP=.TRUE. 的 STOPCAR 是停止指令，只应在实际观察到危险/超时并需要结束时创建。新隔离目录开跑前检查不存在 stale STOPCAR。
3. **NSW=80 与 5–6h 预算不相符。** 原有 15 步耗时约 2.6h，机械外推 80 步约 13.8h（不构成时间保证）。建议改为 NSW=30（实际上与原 INCAR.relax 一致），另设 8h 本地 walltime 监控边界，最多 8 MPI ranks，即本地最多 64 rank-hours，不自动改日本 PBS。
4. **SMASS default(0) 或“阻尼复位”不属于 IBRION=1 的合理调参。** VASP SMASS 控制离子动力学速度/热浴，RMM-DIIS 方案中保持不设置 SMASS；原文所谓“默认0”不正确。
5. **通过力阈值并不等于证明局部极小值。** EDIFFG=-0.02 仅表示所有离子力范数满足 0.02 eV/Å，仍可能是鞍点；需后续振动/位移稳定性和更精确验证。单次优化失败，也**不能证明材料不存在此局部极小值**。
6. **不能对未来作业要求固定“15/15 离子步电子收敛”。** 必须检查该次实际完成的全部离子步，无一处未收敛 SCF/不完整力块；事先不知道会有几步。
7. **电子力的准确度影响 RMM-DIIS。** VASP 建议合理的 NELMIN 下限，试运行可以注明 operator choice：新加 NELMIN=4，并保持 EDIFF=1E-6。MAGMOM 是电子态种子，不是对 Type-I 磁序的约束，必须审核每个 Ag 磁矩而不只看 net M≈0。

## 3. 修订版 P1：纸面协议，尚未获授权

**单一问题：** 给定同一非 SOC PBE+U+D3(0)、固定晶胞与 AgF₂ Type-I 种子，能否得到电子自洽且原子力范数全部不高于 0.02 eV/Å 的 B1 候选结构？

**起点：** PR #38 已核对的 ionic step-15 POSCAR.relaxed（**未收敛**快照），SHA256 d49c82807fc2309575e5d3757904ab78caab348aef10e895ac2327a904ac2d68；来自第15个完整 OUTCAR force block，不能称平衡结构。

**相对已提交 INCAR.relax 的唯一拟议变动：**

| 设置 | 原值 | 拟议新值 |
| --- | --- | --- |
| IBRION | 2 CG | 1 RMM-DIIS |
| POTIM | 0.25 | 0.10（缩小初始步幅尺度，**不是**硬位移上限） |
| NELMIN | 未显式设置 | 4 |
| NSW | 30 | **保持30，不设80** |
| SMASS | 未设置 | **仍不设置** |

原 PBE、Ag 4d Dudarev Ueff=4 eV、IVDW=11 即 D3 zero-damping、同一合法私有 POTCAR（锁定 SHA256 8c0da262849433cad092946af66c1b584d457e46a70d4c1eb768ea68305f0fff）、ENCUT=500、EDIFF=1e-6、EDIFFG=-0.02、ISIF=2、12×8×1 Gamma 筛选网格、ISPIN=2、原 MAGMOM 模式、无 SOC 等不变。不要从不存在的 WAVECAR/CHGCAR 默认接续。

**仅在用户明确授权之后执行：** 创建新的隔离 WSL 目录，检查源坐标和输入哈希、磁种子、算力、RAM、磁序、已有 STOPCAR 状态；VASP 6.6.0 + OpenMPI，最多8 ranks，一次性单作业，NSW≤30、壁钟监测不超过8小时（软停需留出余量）。所有 completed ionic steps 检查能量、力、层间 F–F、Ag 磁矩，安全异常发生时**才创建** LSTOP STOPCAR；不自动调整参数、不自动重投。原始文件留本地，向同一 PR #38 提交摘要/CSV/审计 receipt 和 hash；不改旧证据。

**仅当**全部有效离子步电子收敛且最终每个原子力范数 ≤0.02 eV/Å、无严重近接触、Type-I Ag 磁矩依然成立，才能称为 **force-converged B1 candidate**。随后才需要更严格 <0.01 eV/Å、k 网格收敛、几何正曲率和竞争磁态等；**绝非已经证明局部极小值、磁基态、BSF 双稳态或 IGMC**。

**失败路径：** 若 IBRION=1 发散/达到 NSW/软停，状态 PARTIAL 或 FAILED，仅说明此次数值方法未收敛；不能否定所有 B1 结构极小值。Damped IBRION=3 或 z 方向受控扫描是未来独立备选，不得在本次自动执行。

**算力估计仅供审批：** 最多 8 MPI ranks × 8 h，本地计算，不产生日本 PBS 任务；历史 15 步 2.6h 仅用于估计，真实 walltime 不可保证。**本文件不是计算授权。**

## 4. 依据

- VASP 官方 [Structure optimization](https://vasp.at/wiki/Structure_optimization) 与 [IBRION](https://vasp.at/wiki/IBRION)：RMM-DIIS 适用范围和 CG 机制。
- [EDIFFG](https://vasp.at/wiki/EDIFFG)：负阈值对全部力范数生效。
- [STOPCAR](https://vasp.at/wiki/STOPCAR)：LSTOP 是停止指令。
- [SMASS](https://vasp.at/wiki/SMASS)：离子 MD 的速度/热浴，而非该方案 RMM-DIIS 保证收敛的标签。

**结论：** 接受此次 DeepSeek 的只读 CG 取证为科研审计合格；新计算建议须按以上协议修正并单独授权。没有创建、启动、重试或取消任何计算。
