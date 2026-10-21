# AgF₂ Stage 2：日本 PBS 超算与 WSL 执行说明

**Status: SITE-ADAPTED-PLAN / NOT_SUBMITTED / NOT_RUN.** 本文只确定执行环境和准入步骤，不证明任何实际 VASP 作业已提交。原论文复现数据与研究 Stage-1 冻结文件保持不变。

## 1. 已知环境与职责

- **计算平台：** 日本东北大学金属材料研究所（Tohoku IMR）超算；历史作业采用 **PBS**，**不是 Slurm**。
- **本地端：** 研究者自己的 WSL/SSH 环境，承担 Git 同步、输入检查、远程会话、结果回传与分析。不要把 WSL 当作已核实有 VASP/POTCAR/高性能队列的计算节点。
- **历史记录（待登录时重新确认）：** 远程环境曾检查 VASP 5.4.4/6.4.3、Wannier90 接口可用；曾采用 1 node / 32 cores / 128 GB / 24 h 的 PBS 作业形状。这是**历史配置**，不是当前队列保证。
- **站点路径：** 工作根目录与 SSH 双跳别名由研究者的**私有配置**提供，公共仓库不收集 SSH 用户名、口令、私钥、PBS 账户编号或 POTCAR。
- **独立研究目录：** 本研究只写入 `research/bsf-igmc/agf2_preflight/stage2_geometry/` 和将来的 `stage2_real_dft/`；不改 Xiang IGMC 原论文复现目录。

## 2. 计算物理：唯一的第一批任务

第一批只使用现有、未弛豫的 `stage2_geometry/candidates/B1_g2.6`、`B1_g3.0`、`B1_g3.4` 三个初始结构。

- 0% 应变；**Type-I：层间 Néel 矢量反平行**，不是层间 FM；共线自旋、**非 SOC**。
- 第一轮固定面内晶胞矢量，让原子弛豫，**不人为冻结层间滑移自由度**，除非明确记录物理理由；监测是否滑入另一堆叠态。
- 文献约束：PBE、Ag 4d Dudarev U_eff=4 eV、500 eV cutoff、DFT-D3（具体零阻尼/BJ 变体应查证并显式记录）、约 20 Å 真空、`EDIFF=1e-6 eV`、最大力 <0.01 eV/Å。
- **待选择及验证的运行设置：** Ag/F PAW 数据集版本与 sha256，`IVDW` 的实际实现、`ISMEAR/SIGMA`、`MAGMOM`、`ISYM`、电子混合/优化器、偶极修正与约束；初始 k 网格 12×8×1 仅作廉价筛选，14×9×1 / 16×10×1 用于收敛审查。不得将这些 operator choices 写成原作者明确给出的输入。
- 仅 B1 三组初始间距，**不启动** Type-II、SOC、Berry 极化、NEB、Wannier、IGMC。B2 是后续经 `m_z` 轴矢量自旋规则连接并做独立验证的伙伴构型。

## 3. 日本站点 G0 检查（只读，不产生 DFT 作业）

从研究者已配置的 WSL SSH 别名，登录日本 PBS 环境（保留双跳，不在本仓库写密码）。在**真正的计算登录节点**检查：

```bash
hostname
pwd
command -v qsub
command -v qstat
qstat -Q
# 具体模块命令因站点而异；不要假设“module load vasp”一定存在
type module || true
module avail vasp 2>&1 || true
command -v vasp_std || true
command -v vasp_ncl || true
command -v mpirun || true
```

记录现场版本、可用 PBS 队列、PBS select/walltime 语法、预算余额，以及带软件许可的 VASP 可执行路径。私有 Ag/F POTCAR **只在站点安全目录下核对**数据集版本、元素排序、sha256；记录哈希和版本但**不上传文件本体**。确认不会向公共 GitHub 提交队列账户、SSH 信息或版权 PAW 文件。

**必须获得研究者对具体首批资源额度及真实提交的明确批准**，才可以由具备 SSH 权限的执行端运行 `qsub`。通过 GitHub CI 或本地几何检查 **不等于**获得算力审批。未满足任何条件保留 `BLOCKED_SOURCE_PROTOCOL` / `data_origin=NOT_RUN`。

## 4. PBS 执行原则

- 仅在 G0 成功后生成三组 B1 的 INCAR / KPOINTS / 作业脚本；确认 POSCAR 中物种顺序 **Ag → F** 与现场 POTCAR 顺序一致。
- 从站点查询到的 **PBS 资源格式**生成实际 `#PBS` 指令；**不得盲用 Slurm `sbatch`**。32 核/128 GB/24 h 只是历史先例，资源请求依当前队列与批准额度定。
- 作业目录独立，保持每组 `POSCAR, INCAR, KPOINTS, POTCAR-checksum-record, PBS script`、提交 Job ID、运行 stdout/stderr、`OUTCAR/vasprun.xml` 的私有存档与 SHA256；不要把原始大文件上传公共仓库。
- 对每组记录收敛性、初末结构、最短跨层原子距离、最大残余力、总能、Ag 局域磁矩、是否保持反平行 Néel 态、是否发生滑移以及失败原因。若磁序塌缩、SCF 不稳定或原子接触异常则 **STOP/REVIEW**。
- 三个候选收敛后，先做 14×9×1、16×10×1 复核；待确认结构/能量/磁态后，才讨论 B1/B2 两个局域极小值及滑移势垒。**镜像相关的等能本身不能证明铁电双稳态。**

## 5. 最小可提交的结果包

未来的 `stage2_real_dft/` 应给出：
- `SITE_G0_REPORT.json`：PBS/VASP/PAW 可用性、许可与批准情况、实际版本和非敏感哈希；
- `PILOT_CASES.json`：B1 三组运行条件、Job ID、初末文件 SHA256、单位、力和能量及磁序收敛结论；
- 独立的失败日志、收敛记录和经结构文件验证的几何摘要；
- 大文件只存超算私有路径；GitHub 仅上传可复算的小型报告与审计链。

**边界：ChatGPT 的 GitHub 连接可以维护本仓库代码与审查结果，但不提供 SSH/PBS 终端。实际登录、校验私有 POTCAR 和提交作业必须由授权的本地执行端完成；本文件本身不构成真实提交。**
