# BD36 Phase 0 | para 提交前预检报告（PREFLIGHT: FAIL, 2026-10-23）

状态：**未提交任何作业（无 Job ID）**。预检在真实 para 环境（bingxing-pod-0003,
ssh.paracloud.com:2222, USER=e6a000183）逐项执行，发现两项阻断性缺陷，
按授权规则"预检失败则停止、不擅自扩大授权"终止。

## 预检通过项
| 项 | 结果 |
|---|---|
| SSH/平台 | para 可达（Host para, Slurm 23.11.11, 分区 root.bingxing）✓ |
| 重复任务 | `squeue -u e6a000183` 为空，无运行中/排队同名作业 ✓ |
| MPI | `~/intel/oneapi/mpi/2021.6.0/bin/mpirun`（Intel MPI）✓ |
| wannier90.x | `~/zwy/wannier90.x` = 3.0.0（source oneapi setvars 后可运行）✓ |
| sbatch | /opt/slurm/slurm/bin/sbatch 正常 ✓ |
| 代码同步 | PR #99 研究分支已同步至 ac44ce0（含唯一运行目录/NCORE=1/严格 JSON 修正）✓ |

## 阻断性缺陷 1：para VASP 二进制不支持 .spn 写入（决定性）
候选构建全部核查（2026-10-23 实测）：
- `~/software-e6a000183/544-vtst+wannier/{vasp_ncl,vasp_std,vasp_gam}`（2024-12-19）
- `~/bin/vasp_std_vtst`（2025-09-20）

字符串级（大小写变体全查）：`LWRITE_SPN=0, WRITE_SPN=0, LSPN_FMTED=0, spn_formatted=0`。
接口模式（补丁结构判别，非仅字符串）：`LWANNIER90v2=0, VASP2WANNIER90v2=0, wannier_setup=0`；
`LWANNIER90` 上下文 = INCAR 关键字 + `LWANNIER90_RUN` + `LWRITE_WANPROJ` + `Writing WANPROJ`。
→ para 的 VASP 5.4.4 构建为 **VASP2WANNIER90 v1 接口（WANPROJ 写盘 + postw90 后处理）**，
**不含 VASP2WAN90_v2_fix 的 WRITE_SPN/.spn 写路径**（该功能 2020-06-25 补丁加入，字符串必然存在；
此处全零 + v1 接口结构双确认）。
→ INCAR_soc_w90_vasp544 的 `LWRITE_SPN=.TRUE.` 在 para 上无法产出 .spn，
Phase 0 验收（.spn 维数/Hermite/read_spn 解析）必然失败（EXIT-7/8），56 核时将白烧。

## 阻断性缺陷 2：本账户无 Pt/Se POTCAR 库
`~/` 下无 potpaw/potpaw_PBE/POTCAR_REC 目录；可读 POTCAR 均为其他用户项目文件
（仅 S/Cr 等物种），无 Pt、Se。`run_env.sh` 亦不存在（需 POTCAR_SOURCE 指向合法 Pt/Se
POTCAR）。→ 静态步（SOC 单点）同样无法启动。

## 结论与后续选项（需用户决策，未扩大授权）
1. para **不可行**（.spn 功能缺失 + POTCAR 缺失）→ Phase 0 在 para 无解。
2. 可选替代：IMR（VASP >6.4.2 原生 LWRITE_SPN，恒 formatted）——需用户确认版本≥6.4.3
   与工具链可达；或 WSL（VASP 6.6.0，Se 有/Pt 缺）——需补 Pt POTCAR。
3. 或在 para 重新编译带 VASP2WAN90_v2_fix 补丁的 vasp_ncl（需源码+编译器+授权，
   属环境工程，非 Phase 0 计算授权范围内）。

本报告不修改任何计算输入；研究分支仅新增本日志。
