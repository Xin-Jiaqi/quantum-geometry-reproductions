# AgF₂ G1：单作业真实 VASP 验证，用户已批准

**授权日期：2026-10-21。** 用户在 ChatGPT 会话明确同意 `B1_g3.0` **32 MPI ranks、120 GB、最长 2 小时**的单点 VASP 试算，并要求尽快推进。该批准只覆盖**一次正式 PBS 提交**；不是 3 个几何弛豫、自动重试、Type-II、SOC、Wannier 或 IGMC 的授权。后者仍需单独科学与资源审查。

**任务：无需再次向用户请求同样的单作业资源批准。** 本地 DeepSeek Harness 可以在所有下述技术条件实际通过后，**直接提交这一个作业**，之后跟踪并交付证据。若有实质缺项，明确阻塞，不得自行越权扩大预算。

## 真正的执行对象

- 仅 `stage2_geometry/candidates/B1_g3.0/POSCAR`（GitHub main），SHA256 必须为 `2d323c3fe2353ef984330b3b99fc24a75a85337883bbaf36d08b20ed45e6d60a`；12 原子、Ag4F8、物种顺序 Ag/F；构造未弛豫，非论文作者 relaxed geometry。
- 日本东北大学 IMR PBS：`CP_001` 队列；`select=1:ncpus=32:mpiprocs=32:mem=120gb`，`walltime=02:00:00`。最大请求额 `32×2=64 requested core-hours`。**若站点按整节点计费导致额度超过用户批准范围，暂停提交并汇报收费口径**。
- PBE、Ag-4d Dudarev `U_eff=4 eV`、`ENCUT=500 eV`、D3、`EDIFF=1E-6`。采用 **`IVDW=11` D3 zero damping** 作为合理、明确记录的选择，区别于 `IVDW=12` D3(BJ)，不能把选择写成已复原作者私有 INCAR。
- **非 SOC、共线 `ISPIN=2`**；Type-I 层间 Néel 矢量相反，按 POSCAR Ag 顺序 `MAGMOM=-1,+1,+1,-1`（F 0）；该符号是相对自旋顺序，不代表非 SOC 磁矩具有 SOC 易轴。
- `ISMEAR=0, SIGMA=0.05 eV`（未知实际带隙时先 Gaussian），12×8×1 Gamma-centered k 点仅用于**便宜的 SCF 诊断**，不是收敛生产计算。
- 只做**单点**：`IBRION=-1, NSW=0`，不得误运行结构弛豫。原始建议 `ISMEAR=1, SIGMA=0.1` 是错误先例，不得沿用。
- 本目录的 `INCAR.smoke` / `KPOINTS.smoke` / `run.pbs.template` 是公开、非保密的**推荐最小输入**。正式 PBS 文件必须在私有工作目录中写入真实、经核实的站点 loader 和 CPU `vasp_std` 路径，不能让占位符进入 PBS 队列。

## 本地 Harness：一次执行完整闭环

1. 拉取 GitHub main 最新代码；复制冻结 `B1_g3.0` POSCAR 和本目录 `*.smoke` 文件到**新的私有单作业目录**，不复用旧版三组 relax INCAR。保留原文件 SHA256。仅保存可公开的 INCAR/KPOINTS/PBS 脚本摘要和 SHA，不提交受许可 POTCAR、私人主机/账号/作业日志原文。
2. 登录现有日本 PBS 站点，检查 `qsub/qstat` 与 `CP_001` 的实际约束，选择**一个**存在且可由当前账号合法使用的 CPU `vasp_std`（优先检查已发现的 6.6.1，若不兼容则使用另一站点支持版并记录实际版本）。在**不运行 VASP 生产计算的登录节点**查对执行权限、MPI launcher 版本、库依赖 `ldd`、站点现有加载脚本。PBS 脚本必须采用与所选 MPI 实现兼容的命令，不能不验证就沿用 `mpirun -ppn`。不要猜 MPI/软件路径。
3. 私下按 **Ag 然后 F** 拼接已获得合法授权的 PAW；在计算端核验 **完整拼接** SHA256=`8c0da262849433cad092946af66c1b584d457e46a70d4c1eb768ea68305f0fff`；另核对 TITEL：Ag 02Apr2005、F 08Apr2002，和 POSCAR 物种顺序。**不能上传 POTCAR 内容**。
4. 仅当 POSCAR/INCAR/KPOINTS、VASP bin、MPI、PBS 队列、PAW 均核对无误，且站点本次计费不超上限时，才 **`qsub` 一次**。记录 Job ID、作业脚本与四份输入摘要、耗时和状态；如队列超时或任务失败，**不要自动重新提交**。
5. 作业完成后私有保存 OUTCAR、vasprun.xml、stdout/stderr、完整 VASP 输入、PBS job history 及 sha256；公开只上传**经过去敏处理**的 `G1_RUN_RECEIPT.json`、`G1_SCF_AUDIT.md`、`G1_INPUT_FINGERPRINTS.json`，保留 SCF 是否达到 EDIFF、真实输出中 VASP 版本、最终局域 Ag 磁矩、总磁矩、总能量（eV/12-atom cell）、占据/是否有能隙的**初步**诊断、报错/中断状态。必要时保留 null。
6. 提交到此独立研究目录的新 Draft PR，目标 **main**；不要修改原 Xiang IGMC 冻结复现目录。允许发给 ChatGPT 的只有 Draft PR 链接/编号。若 Job 排队仍未开始，提交 `QUEUED` 的可审计 receipt 并停止主动反复催促；不要虚报完成。

## 审计必须区分的状态

`AUTHORIZED_ONE_SHOT` = 用户批准**一个**请求资源限额内的 G1 singlepoint，不等于 `physics_approved_by_user=true` 对 Stage6 全范围科研计算的批准。

`SCF_CONVERGED` = 真实 VASP 输出和收敛记录支持电子迭代结束；它**不是**能隙、磁基态或双稳态已经复现。即使 `SCF_CONVERGED`，仍需要磁矩、轨道占据、带隙及数值稳定性独立审查。

**禁止**自动扩展到 B1 其他间距、B2、三组弛豫、SOC、NEB、Wannier、IGMC。不得覆盖或修改任何尚未计算状态为已计算。保留所有失败/被阻塞证据以便修复。

## 参考

- 本仓库 `stage2_geometry/JAPAN_PBS_WSL_RUNBOOK.md`、`stage2_geometry/CALCULATION_GATE.md`、`stage2_real_dft/preflight/`。
- VASP 官方： https://vasp.at/wiki/ISMEAR 、 https://vasp.at/wiki/IVDW 、 https://vasp.at/wiki/ISPIN 。
